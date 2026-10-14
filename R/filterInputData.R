#' Filter data quality, detection and prior associations
#'
#' Single entry point of the filtering step. In one call it (1) evaluates the
#' numeric validity (NA, NaN, Inf and zero variance) of the target features and
#' of every regulatory layer, (2) intersects the prior associations with the
#' features that survived, drops layers without usable priors and removes
#' features without priors, (3) counts the cells with a value of every target
#' and regulator in each group of the design and (4) prunes the prior pairs
#' according to those counts. Every excluded feature is recorded together with
#' the reason.
#'
#' @details
#' The quality metrics and the counts are computed with sparse-aware routines,
#' so sparse inputs (for example ATAC counts stored as \code{dgCMatrix}) are not
#' converted to dense matrices. Duplicated prior pairs are removed.
#'
#' Detection is evaluated on the original values (a cell "has a value" when it
#' is different from zero). The pruning is done in a single pass over the table
#' of prior pairs: a pair is kept when the target is usable in at least
#' \code{minGroups} groups and the regulator is estimable in at least one group
#' where the target is usable (see \code{\link{filterParams}}). Targets and
#' regulators left without pairs are excluded as a consequence; this needs no
#' further passes because each pair depends only on the counts of its target and
#' its regulator.
#'
#' With \code{interaction = FALSE} each regulator has one slope shared by all
#' groups, so the regulator is evaluated on the cells of the groups where each
#' target is usable: a pair is kept when the regulator has at least
#' \code{minCellsRegulatorGroup} cells with a value in total in those groups.
#'
#' @param targetData Matrix, sparse \code{Matrix} or data.frame of target
#'   features (features x cells).
#' @param regulatoryData Named list of matrices, sparse matrices or data.frames
#'   (regulators x cells), one per regulatory layer.
#' @param priorAssociations Named list of data.frames with the prior
#'   regulator-target pairs of each layer (first column = regulator, second
#'   column = target feature). Any additional columns are preserved.
#' @param metadata Data.frame with cell metadata. Row names must be cell IDs.
#'   If \code{NULL}, all cells form a single group.
#' @param fixedEffects Character vector of \code{metadata} columns used as
#'  fixed terms. They define the groups.
#' @param filters Detection filters. The output of \code{\link{filterParams}}, a
#'   named list with some of its elements or \code{NULL} to disable them.
#' @param interaction Logical. If \code{TRUE}, regulators are crossed
#'   with all fixed effects and their interactions.
#'
#' @return A list with the filtered \code{targetData}, \code{regulatoryData} and
#'   \code{priorAssociations}, the data.frame \code{excludedFeatures}
#'   (columns \code{feature}, \code{layer} and \code{exclusion_reason}) and
#'   \code{estimability}: a list with the group of each cell, the groups, the
#'   logical matrices \code{targetGroups} (targets x groups) and
#'   \code{regulatorGroups} (one per layer, regulators x groups) with the
#'   groups where each feature passes its threshold, the corresponding counts of
#'   cells with a value, and the filters used.
#'
#' @importFrom cli cli_abort cli_alert_info cli_alert_warning cli_alert_success
#' @keywords internal
filterInputData <- function(targetData,
                            regulatoryData,
                            priorAssociations,
                            metadata,
                            fixedEffects,
                            filters,
                            interaction) {

  filters <- resolveFilters(filters)

  excludedFeatures <- data.frame(
    feature = character(0),
    layer = character(0),
    exclusion_reason = character(0),
    stringsAsFactors = FALSE
  )

  # Excludes non-finite and zero-variance rows of one matrix
  qualityFilter <- function(mat, layer) {
    q <- rowQuality(mat)
    excludedFeatures <<- addExcluded(excludedFeatures, rownames(mat)[q$nonFinite],
                                     layer, "Contains NA, NaN, or Inf values")
    excludedFeatures <<- addExcluded(excludedFeatures, rownames(mat)[q$zeroVar],
                                     layer, "Zero variance across cells")
    if (any(q$nonFinite) || any(q$zeroVar)) {
      cli::cli_alert_info(
        "{.field {layer}}: excluded {sum(q$nonFinite)} feature(s) with NA/NaN/Inf and {sum(q$zeroVar)} with zero variance."
      )
    }
    mat[!q$nonFinite & !q$zeroVar, , drop = FALSE]
  }

  # ~~~~~~~~~~ 1. Zero variance and quality check on targets ~~~~~~~~~~ #
  targetData <- qualityFilter(targetData, "target")
  if (nrow(targetData) == 0) {
    cli::cli_abort("No valid target features remaining after filtering out NA values and zero-variance features.")
  }

  # ~~~~~~~~~~ 2. Zero variance and quality check on regulatory layers ~~~~~~~~~~ #
  cleaned_regulatoryData <- lapply(names(regulatoryData), function(m_name) {
    qualityFilter(regulatoryData[[m_name]], m_name)
  })
  names(cleaned_regulatoryData) <- names(regulatoryData)

  # ~~~~~~~~~~ 3. Validate layers and prior associations ~~~~~~~~~~ #
  valid_layers <- character(0)
  filtered_priors <- list()

  for (m_name in names(cleaned_regulatoryData)) {
    assoc <- priorAssociations[[m_name]]

    if (is.null(assoc) || nrow(assoc) == 0) {
      cli::cli_alert_warning("Layer {.field {m_name}}: {.arg priorAssociations} is NULL or empty. Excluding this regulatory omic.")
      next
    }

    if (ncol(assoc) < 2) {
      cli::cli_alert_warning("Layer {.field {m_name}}: {.arg priorAssociations} has fewer than 2 columns. Excluding layer.")
      next
    }

    # Intersect prior pairs with available valid features (col 1 = regulator, col 2 = target)
    reg_features <- rownames(cleaned_regulatoryData[[m_name]])
    target_features <- rownames(targetData)

    keep_idx <- (assoc[[1]] %in% reg_features) & (assoc[[2]] %in% target_features)
    valid_assoc <- unique(assoc[keep_idx, , drop = FALSE])
    rownames(valid_assoc) <- NULL

    if (nrow(valid_assoc) == 0) {
      cli::cli_alert_warning("Layer {.field {m_name}}: No matching target-regulator pairs found in expression matrices. Excluding layer.")
      next
    }

    filtered_priors[[m_name]] <- valid_assoc
    valid_layers <- c(valid_layers, m_name)
  }

  if (length(valid_layers) == 0) {
    cli::cli_abort("No regulatory layers remaining with valid prior associations and matching expression data.")
  }

  regulatoryData <- cleaned_regulatoryData[valid_layers]
  priorAssociations <- filtered_priors[valid_layers]

  # ~~~~~~~~~~ 4. Subset active targets and regulators ~~~~~~~~~~ #
  # Filter target features without any active prior association
  active_targets_in_priors <- unique(unlist(lapply(priorAssociations, function(df) df[[2]])))
  targets_without_regulators <- setdiff(rownames(targetData), active_targets_in_priors)

  if (length(targets_without_regulators) > 0) {
    cli::cli_alert_info("{length(targets_without_regulators)} target feature(s) have no prior regulatory associations and will be excluded.")

    excludedFeatures <- addExcluded(excludedFeatures, targets_without_regulators,
                                    "target", "No prior regulatory associations found")

    targetData <- targetData[rownames(targetData) %in% active_targets_in_priors, , drop = FALSE]
  }

  if (nrow(targetData) == 0) {
    cli::cli_abort("No target features remaining after filtering for prior associations.")
  }

  # Filter regulators without any active prior association in their respective layer
  for (m_name in names(regulatoryData)) {
    active_regs <- unique(priorAssociations[[m_name]][[1]])
    unused_regs <- setdiff(rownames(regulatoryData[[m_name]]), active_regs)

    if (length(unused_regs) > 0) {
      cli::cli_alert_info("Layer {.field {m_name}}: {length(unused_regs)} regulator(s) have no active target associations and will be excluded.")

      excludedFeatures <- addExcluded(excludedFeatures, unused_regs, m_name,
                                      "No active prior target associations found")

      regulatoryData[[m_name]] <- regulatoryData[[m_name]][active_regs, , drop = FALSE]
    }
  }

  # ~~~~~~~~~~ 5. Detection of targets and regulators in each group ~~~~~~~~~~ #
  cells  <- colnames(targetData)
  group  <- cellGroups(metadata, fixedEffects, cells)
  groups <- levels(group)
  nGroups <- length(groups)

  minGroups <- filters$minGroups
  if (minGroups > nGroups) {
    if (nGroups == 1) {
      cli::cli_alert_warning("The design defines a single group: {.arg minGroups} was set to 1.")
      minGroups <- 1L
    } else {
      cli::cli_abort(c(
        "{.arg minGroups} ({minGroups}) is larger than the number of groups of the design ({nGroups}).",
        "i" = "Groups: {.val {groups}}."
      ))
    }
  }

  cli::cli_alert_info("Groups defined by {.field {if (length(fixedEffects)) fixedEffects else 'none (single group)'}}: {nGroups} group(s).")
  if (!interaction) {
    cli::cli_alert_info("{.code interaction = FALSE}: regulators are evaluated on the cells of the groups where each target is detected, not group by group.")
  }

  # Cells with a value of each target / regulator in each group, and the groups that pass the thresholds
  tCounts <- countNonZeroByGroup(targetData, group)
  rCounts <- lapply(regulatoryData, countNonZeroByGroup, group = group)
  tOK     <- tCounts >= filters$minCellsTargetGroup
  rOK     <- lapply(rCounts, function(m) m >= filters$minCellsRegulatorGroup)
  tFeatureOK <- rowSums(tOK) >= minGroups

  # ~~~~~~~~~~ 6. Prune the prior pairs in a single pass ~~~~~~~~~~ #
  # A pair depends only on the counts of its target and of its regulator, so no iteration is needed
  pairKeep <- lapply(names(priorAssociations), function(m) {
    a  <- priorAssociations[[m]]
    r  <- as.character(a[[1]])
    t  <- as.character(a[[2]])
    tp <- tOK[t, , drop = FALSE]
    ok <- if (interaction) {
      # estimable in at least one group where the target is usable
      rowSums(tp & rOK[[m]][r, , drop = FALSE]) >= 1
    } else {
      # one slope for all groups: enough cells with a value in the groups where the target is usable
      rowSums(rCounts[[m]][r, , drop = FALSE] * tp) >= filters$minCellsRegulatorGroup
    }
    ok & tFeatureOK[t]
  })
  names(pairKeep) <- names(priorAssociations)

  # Targets without any kept pair
  keptTargets <- unique(unlist(lapply(names(priorAssociations), function(m) {
    as.character(priorAssociations[[m]][[2]])[pairKeep[[m]]]
  })))
  targetsOut <- setdiff(rownames(targetData), keptTargets)
  if (length(targetsOut) > 0) {
    reasonT <- ifelse(
      !tFeatureOK[targetsOut],
      sprintf("Detected in %d of %d groups (needs >= %d groups with >= %d cells with a value)",
              rowSums(tOK)[targetsOut], nGroups, minGroups, filters$minCellsTargetGroup),
      "No regulator estimable in the groups where the target is detected"
    )
    excludedFeatures <- addExcluded(excludedFeatures, targetsOut, "target", unname(reasonT))
    cli::cli_alert_info("{length(targetsOut)} target feature(s) excluded by the detection filters.")
  }

  # Regulators without any kept pair
  for (m in names(priorAssociations)) {
    keptRegs <- unique(as.character(priorAssociations[[m]][[1]])[pairKeep[[m]]])
    regsOut  <- setdiff(rownames(regulatoryData[[m]]), keptRegs)
    if (length(regsOut) > 0) {
      reasonR <- if (interaction) {
        ifelse(rowSums(rOK[[m]])[regsOut] == 0,
               sprintf("Fewer than %d cells with a value in every group", filters$minCellsRegulatorGroup),
               "No target in which it is estimable (groups where the target is detected)")
      } else {
        sprintf("Fewer than %d cells with a value in the groups where its targets are detected",
                filters$minCellsRegulatorGroup)
      }
      excludedFeatures <- addExcluded(excludedFeatures, regsOut, m, unname(reasonR))
      cli::cli_alert_info("Layer {.field {m}}: {length(regsOut)} regulator(s) excluded by the detection filters.")
    }
  }

  # Apply the pruning to priors and matrices
  nPairsIn  <- sum(vapply(priorAssociations, nrow, integer(1)))
  for (m in names(priorAssociations)) {
    priorAssociations[[m]] <- priorAssociations[[m]][pairKeep[[m]], , drop = FALSE]
    rownames(priorAssociations[[m]]) <- NULL
  }
  nPairsOut <- sum(vapply(priorAssociations, nrow, integer(1)))

  emptyLayers <- names(priorAssociations)[vapply(priorAssociations, nrow, integer(1)) == 0]
  if (length(emptyLayers) > 0) {
    cli::cli_alert_warning("Layer(s) {.field {emptyLayers}} have no regulator-target pair left after the detection filters and are excluded.")
    priorAssociations <- priorAssociations[setdiff(names(priorAssociations), emptyLayers)]
    regulatoryData    <- regulatoryData[names(priorAssociations)]
  }
  if (length(priorAssociations) == 0 || nrow(targetData) == length(targetsOut)) {
    cli::cli_abort(c(
      "No regulator-target pairs remaining after the detection filters.",
      "i" = "Consider relaxing {.arg filters}, for example {.code filterParams(minCellsTargetGroup = 3, minCellsRegulatorGroup = 3, minGroups = 1)}."
    ))
  }

  for (m in names(regulatoryData)) {
    active_regs <- unique(as.character(priorAssociations[[m]][[1]]))
    if (length(setdiff(rownames(regulatoryData[[m]]), active_regs)) > 0) {
      regulatoryData[[m]] <- regulatoryData[[m]][active_regs, , drop = FALSE]
    }
  }
  if (length(targetsOut) > 0) {
    targetData <- targetData[!rownames(targetData) %in% targetsOut, , drop = FALSE]
  }

  # ~~~~~~~~~~ 7. Estimability tables of what remains ~~~~~~~~~~ #
  estimability <- list(
    group            = group,
    groups           = groups,
    targetGroups     = tOK[rownames(targetData), , drop = FALSE],
    targetCounts     = tCounts[rownames(targetData), , drop = FALSE],
    regulatorGroups  = lapply(names(regulatoryData), function(m) rOK[[m]][rownames(regulatoryData[[m]]), , drop = FALSE]),
    regulatorCounts  = lapply(names(regulatoryData), function(m) rCounts[[m]][rownames(regulatoryData[[m]]), , drop = FALSE]),
    minGroups        = minGroups,
    filters          = filters,
    interaction      = interaction
  )
  names(estimability$regulatorGroups) <- names(estimability$regulatorCounts) <- names(regulatoryData)

  cli::cli_alert_info("{nPairsOut} of {nPairsIn} prior pair(s) kept; {sum(!estimability$targetGroups)} of {length(estimability$targetGroups)} target-group combinations are below {filters$minCellsTargetGroup} cells with a value (see {.code estimability}).")

  cli::cli_alert_success("Data quality filtering and prior matching completed successfully.")

  return(list(
    targetData        = targetData,
    regulatoryData    = regulatoryData,
    priorAssociations = priorAssociations,
    excludedFeatures  = excludedFeatures,
    estimability      = estimability
  ))
}
