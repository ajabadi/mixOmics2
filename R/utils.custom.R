## -----------------------------------------------------------------------------------
## ---------------  custom stop to define specific error classes
## -----------------------------------------------------------------------------------
.stop<- function(.subclass, message, call = NULL, ...) {
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
.inv_signature <- function() .stop("inv_signature", "incorrect input format")
## ----------- for invalid data
.inv_mae <- function(data) .stop("inv_mae", paste0(squote(data), " must be a MultiAssayExperiment object"))
## ----------- for invalid formula for single
.inv_sformula <- function() .stop("inv_sformula", "'formula' must be a formula object of form Y~X where X and Y are numeric matrices")
## ----------- for invalid formula for blocks
.inv_bformula <- function() .stop("inv_bformula", "'formula' must be a formula object of form Y~X where X is a list of numeric matrices")
## ----------- for invalid X/Y
.inv_assay <- function() .stop("inv_assay", "invalid assay/colData name(s).")

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
.matchArg <- function (arg, choices, several.ok = FALSE)
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
                                                                      collapse = ", "))), call.=FALSE, domain = NA)
  i <- i[i > 0L]
  if (!several.ok && length(i) > 1)
    stop("there is more than one match in 'match.arg'")
  choices[i]
}
