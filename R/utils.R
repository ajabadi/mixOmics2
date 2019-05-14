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
  } else {stop_custom("invalid_Y", paste0(Y_char, " is not an assay or colData from ", args.list[["data"]] ))}
  return(Y)
}

## ---- function to check whether a character X_char matches to any assay of a MAE object, and return it if so

X_mae_check <- function(X_char, data){
  args.list <- match.call()[-1]
  if(!X_char %in% names(assays(data))) stop_custom("invalid_X", paste0(X_char, " is not a valid assay from ", args.list[["data"]]))
  X <- as.matrix(assay(data, X_char))
  return(X)
}

##------------------------------------------
## ----- function to get X, Y, formula, and data, check for validity, and return matrix/numeric X and Y for internal
check_generic_args <- function(X,Y,formula,data, Expect=c("xy", "formula", "formula_mae", "xy_mae")){
  args.list <- match.call()[-1]
  ## expect 'Expect' to match given values
  Expect <- match.arg(Expect)
  ##------- if only matrix X and matrix/numeric Y expected
  if(Expect=="xy"){
    if(class(try(formula))!="NULL") stop_custom("args_conflict", message = "when formula is provided, X and Y must be NULL")
    if(class(try(data))!="NULL") stop_custom("args_conflict", message = "X and Y should be assay/colData names from data, or data should be NULL")
  }
  ##------- if only formula expected
  else if(Expect=="formula"){
    if(class(try(X))!="NULL" | class(try(Y))!="NULL") stop_custom("args_conflict", message = "when formula is provided, X and Y must be NULL")
    if(class(try(data))!="NULL") stop_custom("args_conflict", message = "when X and Y are matrices data should be NULL, or X and Y have to be assay name from data.")

    fterms <- as.list(formula)[-1]
    X=eval(fterms[[2]], envir = parent.frame())
    Y=eval(fterms[[1]], envir = parent.frame())
  }
  ##------- if formula = coldata/assay ~ assay and data=MAE expected
  else if(Expect=="formula_mae"){
    if(class(try(X))!="NULL" | class(try(Y))!="NULL") stop_custom("args_conflict", message = "when formula is provided, X and Y must be NULL")
    # if(class(try(data))!="MultiAssayExperiment") stop_custom("invalid_class", message = "''data' must be a MultiAssayExperiment object")
    fterms <- as.list(formula)[-1]
    X_char <- as.character(fterms[[2]])
    Y_char <- as.character(fterms[[1]])
    X <- X_mae_check(X_char = X_char, data=data)
    Y <- Y_mae_check(Y_char = Y_char, data = data)
  }
  ##------- if X=assay (symbol or char), Y=assay/colData (symbol or char) and data=MAE expected
  else if(Expect=="xy_mae"){
    if(class(try(formula))!="NULL") stop_custom("args_conflict", message = "when formula is provided, X and Y must be NULL")
    X_char <- as.character(substitute(X))
    Y_char <- as.character(substitute(Y))
    X <- X_mae_check(X_char = X_char, data=data)
    Y <- Y_mae_check(Y_char = Y_char, data = data)

  }
  return(list(X=X, Y=Y))
}

