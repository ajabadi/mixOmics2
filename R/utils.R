## ---- custom stop to define error classes ----
stop_custom <- function(.subclass, message, call = NULL, ...) {
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


## ---- check that EXP can be evaluated and returned, if not throw error with error_message error  ----
internal_check_eval <- function(EXP, error_message = paste0(as.expression(substitute(EXP)), " cannot be evaluated")){
  formals(stop)$call. <-FALSE
  out <- tryCatch(EXP, error = function(e) e)
  if ("simpleError" %in% class(out)) stop(error_message)
  return(out)
}


## ---- get a MAE object, assay name/index, and the call list, and return the data matrix for MAE methods ----
## args.list contains X entry
#' @importFrom SummarizedExperiment assays
internal_mae2dm <- function(X, Assay){ ## MAE to data.matrix
  args.list <- match.call()[-1]
  ## ---------- get the assay name from either name  or index provided and check
  if (is.numeric(Assay)){
    if(Assay-floor(Assay)!=0) stop_custom(.subclass = "invalid_assay", message = paste0(Assay, " is not a valid assay index. Use an integer."))
    if(Assay<1 | Assay >length(assays(X))) stop_custom(.subclass = "invalid_assay", message = paste0("assay index must be positive integer smaller than or equal to the number of assays in ",
                                                                                                     args.list["X"]," (i.e. 1:",length(assays(X)),")"))
    Assay <- names(assays(X))[Assay]
  } else if(is.null(Assay)){
    stop_custom(.subclass = "invalid_assay", message = paste0("Please provide an assay from MultiAssayExperiment object ", args.list["X"]))
    } else if(is.na(Assay)){
      stop_custom(.subclass = "invalid_assay", message = "Assay cannot be NA ")
    } else if(!is.character(Assay)){
    ## if 'Assay' it none of acceptable forms
    stop_custom(.subclass = "invalid_assay", message = paste0("'Assay' must be either an assay name or index from ", args.list["X"]))
  }

  ## ---------- use assay name to get the data matrix
  if(! Assay %in% names(assays(X)))
    stop_custom(.subclass = "invalid_assay", message = paste0(Assay, " is not a valid assay from ","'",args.list["X"],"'"))
  ## transpose and create a data matrix
  X <- internal_check_eval(EXP = data.matrix(t(assay(X,Assay))),
                                        error_message = paste0("could not create a data matrix from assay '", Assay, "' in ", args.list["X"]) )
  if(!is.numeric(as.matrix(X))) stop_custom(.subclass = "invalid_assay", message = paste0("The ", Assay, " assay contains non-numeric values"))
  return(X)
}
