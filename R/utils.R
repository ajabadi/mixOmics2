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

##-----------------------------------------------------------------------------------------------------------------------
## ---- check that EXP can be evaluated and returned, if not throw error with error_message error  ----
internal_check_eval <- function(EXP, error_message = paste0(as.expression(substitute(EXP)), " cannot be evaluated")){
  formals(stop)$call. <-FALSE
  out <- tryCatch(EXP, error = function(e) e)
  if ("simpleError" %in% class(out)) stop(error_message)
  return(out)
}

##-----------------------------------------------------------------------------------------------------------------------
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



##-----------------------------------------------------------------------------------------------------------------------
##------------------------------------------
## ---- function to check whether a character Y_char matches to any and not both colData and assays of a MAE object, and return the numeric matrix Y if possible

Y_mae_check <- function(Y_char, data){
  args.list <- match.call()[-1]
  if(Y_char %in% names(colData(data))){
    if(Y_char %in% names(assays(data))) stop_custom("'Y' matches to both colData and assay name, change its name in one and continue.")
    Y= as.matrix(colData(data)[Y_char]) ## DataFrame to matrix
    ## if not numeric
    if(! typeof(Y)%in% c("numeric","integer")){
      if(typeof(Y)=="factor") {
        warning("The colData ", Y_char, " is factor, coercing to a numeric...")
      } else{
        stop_custom("invalid_coldata", paste0("The colData",Y_char," is not a numeric matrix (or a factor coercible to numeric)"))
      }
    }
    ## assay case
  } else if(Y_char %in% names(assays(data))){
    Y <- as.matrix(assay(data, Y_char))
  } else {stop_custom("invalid_assay", paste0(Y_char, " is not an assay or colData from ", args.list[["data"]] ))}
  return(Y)
}

## ---- function to check whether a character X_char matches to any assay of a MAE object, and return it if so

X_mae_check <- function(X_char, data){
  args.list <- match.call()[-1]
  if(!X_char %in% names(assays(data))) stop_custom("invalid_assay", paste0(X_char, " is not a valid assay from ", args.list[["data"]]))
  X <- as.matrix(assay(data, X_char))
  return(X)
}

##------------------------------------------
## ----- function to spls methods arguments plus the methods call mode and create internal level arguments
pls_methods_wrapper <- function(method.mode=c("xy", "formula", "formula_mae", "xy_mae"),...){
  formals(eval.parent)$n <- 2 ## within this function evluate formals 2 stacks higher in parent as it is a second level function
  ## eventually, mc must have all arguments required for internal - including the dots
  mc <- as.list(match.call()[-c(1,2)])
  match.arg(method.mode)
  ##------- if xy and y are provided no further work needed
  ##------- if only formula provided
  if(method.mode=="formula"){
    fterms <- as.list(eval.parent(mc$formula))[-1]
    if(length(fterms)!=2) stop_custom("invalid_formula", "The supplied formula must be of form Y~X")
    mc$X <- eval.parent(fterms[[2]])
    mc$Y <- eval.parent(fterms[[1]])
    mc$formula <- NULL
  }
  ##------- if formula = coldata/assay ~ assay and data=MAE provided
  else if(method.mode=="formula_mae"){
    mc$data <- eval.parent(mc$data)
    ## so that formula can also be stored in a variable - makes testing easier as well
    mc$formula <- eval.parent(mc$formula)
    fterms <- as.list(mc$formula)[-1]
    X_char <- as.character(fterms[[2]])
    Y_char <- as.character(fterms[[1]])
    mc$X <- X_mae_check(X_char = X_char, data=mc$data)
    mc$Y <- Y_mae_check(Y_char = Y_char, data = mc$data)
    mc$formula <- mc$data <-  NULL
  }
  ##------- if X=assay (symbol or char), Y=assay/colData (symbol or char) and data=MAE provided
  else if(method.mode=="xy_mae"){
    mc$data <- eval.parent(mc$data)
    ## so X and Y can be variables
    mc$X <- eval.parent(mc$X)
    mc$Y <- eval.parent(mc$Y)
    mc$X <- X_mae_check(X_char = mc$X, data=mc$data)
    mc$Y <- Y_mae_check(Y_char = mc$Y, data = mc$data)
    mc$data <-  NULL
  }
  do.call(.pls, mc)
}
