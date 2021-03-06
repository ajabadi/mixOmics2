## -----------------------------------------------------------------------------------
## ---------------  custom stop to define specific error classes
## -----------------------------------------------------------------------------------
.stop <- function(.subclass, message, call = NULL, ...) {
  formals(stop)$call. <- FALSE
  err <- structure(
    list(
      message = message,
      call = call,
      ...
    ),
    class = c(.subclass, "error", "condition")
  )
  stop(err)
}

##TODO remove the first two
## ----------- for invalid signature
.inv_signature <- function(msg="incorrect input format") .stop("inv_signature", msg)
## ----------- for invalid data
.inv_mae <- function(data='data', msg=" is not a MultiAssayExperiment object") .stop("inv_mae", paste0(squote(data), msg))
## ----------- for invalid formula for single
.inv_sformula <- function(msg="'formula' must be a formula object of form Y~X where X and
                          Y are numeric matrices, or assay names from 'data'") .stop("inv_sformula", msg)
## ----------- for invalid formula for blocks
.inv_bformula <- function(msg="'formula' must be a formula object of form Y~X where Y is a
                          numeric matrix (or name of such an assay from 'data') and X is a
                          list of numeric matrices (or assay names)") .stop("inv_bformula", msg)
## ----------- for invalid X/Y
.inv_assay <- function(msg="invalid assay/colData name(s).") .stop("inv_assay", msg)

## -----------------------------------------------------------------------------------
## ---------------  custom warnings with specified class
## -----------------------------------------------------------------------------------
.warning <- function(.subclass, message, call=NULL, ...){
  formals(warning)$call. <- FALSE
  warn <- structure(
    list(
      message = message,
      call = call,
      ...
    ),
    class = c(.subclass, "warning", "condition")
  )
  warning(warn)
}

## -----------------------------------------------------------------------------------
## --------------- custom match.arg with call.=FALSE for stop()
## -----------------------------------------------------------------------------------
.matchArg <- function(arg, choices, several.ok = FALSE)
{
  if (missing(choices)) {
    formal.args <- formals(sys.function(sysP <- sys.parent()))
    choices <- eval(formal.args[[as.character(substitute(arg))]],
                    envir = sys.frame(sysP))
  }
  if (is.null(arg))
    return(choices[1L])
  else if (!is.character(arg))
    stop("'arg' must be NULL or a character vector")
  if (!several.ok) {
    if (identical(arg, choices))
      return(arg[1L])
    if (length(arg) > 1L)
      stop("'arg' must be of length 1")
  }
  else if (length(arg) == 0L)
    stop("'arg' must be of length >= 1")
  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
  if (all(i == 0L))
    stop(paste(match.call()$arg,gettextf("should be one of %s", paste(dQuote(choices),
                                                                      collapse = ", "))), call. = FALSE, domain = NA)
  i <- i[i > 0L]
  if (!several.ok && length(i) > 1)
    stop("there is more than one match in 'match.arg'")
  choices[i]
}
