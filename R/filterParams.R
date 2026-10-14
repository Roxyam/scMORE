#' Parameters of the detection filters
#'
#' Constructor of the \code{filters} argument of \code{\link{fitMORE}}. It groups
#' the thresholds that decide which target features, regulators and
#' regulator-target pairs are modelled, validates them and documents their
#' defaults in one place.
#'
#' @details
#' A group is a combination of the levels of the fixed effects (for example
#' \code{condition} x \code{cell_type}). A cell "has a value" when its value is
#' different from zero in the original data (before any residualization). The
#' filters work on the number of cells with a value in each group:
#' \itemize{
#'   \item A target is \emph{usable} in a group when it has at least
#'     \code{minCellsTargetGroup} cells with a value there.
#'   \item A regulator is \emph{estimable} in a group when it has at least
#'     \code{minCellsRegulatorGroup} cells with a value there.
#'   \item A target needs to be usable in at least \code{minGroups} groups.
#'     \code{minGroups = 2} keeps only targets that allow comparing groups; with
#'     \code{minGroups = 1} a target with a single usable group is modelled and
#'     tested, but contrasts between groups are not possible for it.
#'   \item A regulator-target pair is kept when the regulator is estimable in at
#'     least one group where the target is usable.
#' }
#' Targets and regulators without any kept pair are excluded and reported in
#' \code{excludedFeatures}.
#'
#' @param minCellsTargetGroup Numeric. Minimum number of cells with a value for a
#'   target to be usable in a group (whole number >= 0). By default, 3.
#' @param minCellsRegulatorGroup Numeric. Minimum number of cells with a value for
#'   a regulator to be estimable in a group (whole number >= 0). By default, 3.
#' @param minGroups Numeric. Minimum number of groups in which a target has to be
#'   usable (whole number >= 1). By default, 2.
#' @param rankGuard Logical. Reserved for future use. It will check, when the
#'   model of each target is built, that the columns of the model are linearly
#'   independent and will mark the dependent slopes as not estimable. It is stored
#'   but not applied yet. By default, \code{TRUE}.
#'
#' @return An object of class \code{scMoreFilters} (a named list).
#'
#' @examples
#' filterParams()
#' filterParams(minCellsTargetGroup = 10, minGroups = 1)
#'
#' @importFrom cli cli_abort
#' @export
filterParams <- function(minCellsTargetGroup    = 3,
                         minCellsRegulatorGroup = 3,
                         minGroups              = 2,
                         rankGuard              = TRUE) {

  checkCount <- function(x, nm, min = 0) {
    if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < min || x != round(x)) {
      cli::cli_abort("{.arg {nm}} must be a single whole number >= {min}.")
    }
  }
  checkCount(minCellsTargetGroup,    "minCellsTargetGroup")
  checkCount(minCellsRegulatorGroup, "minCellsRegulatorGroup")
  checkCount(minGroups,              "minGroups", min = 1)
  if (!is.logical(rankGuard) || length(rankGuard) != 1 || is.na(rankGuard)) {
    cli::cli_abort("{.arg rankGuard} must be {.code TRUE} or {.code FALSE}.")
  }

  structure(
    list(minCellsTargetGroup    = as.integer(minCellsTargetGroup),
         minCellsRegulatorGroup = as.integer(minCellsRegulatorGroup),
         minGroups              = as.integer(minGroups),
         rankGuard              = rankGuard),
    class = c("scMoreFilters", "list")
  )
}


#' Resolve the \code{filters} argument
#'
#' Accepts the output of \code{\link{filterParams}}, a (possibly partial) named
#' list, or \code{NULL}, and returns a validated \code{scMoreFilters} object.
#' Missing elements of a list take the default values. \code{NULL} disables the
#' detection filters (every threshold is 0 and \code{minGroups} is 1).
#'
#' @param filters \code{NULL}, a named list or a \code{scMoreFilters} object.
#'
#' @return A validated \code{scMoreFilters} object.
#'
#' @importFrom cli cli_abort
#' @keywords internal
resolveFilters <- function(filters) {
  if (is.null(filters)) {
    return(filterParams(minCellsTargetGroup = 0, minCellsRegulatorGroup = 0,
                        minGroups = 1, rankGuard = TRUE))
  }
  if (inherits(filters, "scMoreFilters")) return(filters)

  valid <- names(formals(filterParams))
  if (!is.list(filters) || (length(filters) > 0 && (is.null(names(filters)) || any(!nzchar(names(filters)))))) {
    cli::cli_abort("{.arg filters} must be {.code NULL}, the output of {.fn filterParams} or a named list.")
  }
  unknown <- setdiff(names(filters), valid)
  if (length(unknown) > 0) {
    cli::cli_abort(c("Unknown element(s) in {.arg filters}: {.val {unknown}}.",
                     "i" = "Valid elements: {.val {valid}}."))
  }
  do.call(filterParams, utils::modifyList(unclass(filterParams()), filters))
}
