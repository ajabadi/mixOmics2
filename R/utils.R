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
    if(Assay-floor(Assay)!=0) stop_custom(.subclass = "inv_xy", message = paste0(Assay, " is not a valid assay index. Use an integer."))
    if(Assay<1 | Assay >length(assays(X))) stop_custom(.subclass = "inv_xy", message = paste0("assay index must be positive integer smaller than or equal to the number of assays in ",
                                                                                                     args.list["X"]," (i.e. 1:",length(assays(X)),")"))
    Assay <- names(assays(X))[Assay]
  } else if(is.null(Assay)){
    stop_custom(.subclass = "inv_xy", message = paste0("Please provide an assay from MultiAssayExperiment object ", args.list["X"]))
    } else if(is.na(Assay)){
      stop_custom(.subclass = "inv_xy", message = "Assay cannot be NA ")
    } else if(!is.character(Assay)){
    ## if 'Assay' it none of acceptable forms
    stop_custom(.subclass = "inv_xy", message = paste0("'Assay' must be either an assay name or index from ", args.list["X"]))
  }

  ## ---------- use assay name to get the data matrix
  if(! Assay %in% names(assays(X)))
    stop_custom(.subclass = "inv_xy", message = paste0(Assay, " is not a valid assay from ","'",args.list["X"],"'"))
  ## transpose and create a data matrix
  X <- internal_check_eval(EXP = data.matrix(t(assay(X,Assay))),
                                        error_message = paste0("could not create a data matrix from assay '", Assay, "' in ", args.list["X"]) )
  if(!is.numeric(as.matrix(X))) stop_custom(.subclass = "inv_xy", message = paste0("The ", Assay, " assay contains non-numeric values"))
  return(X)
}


## ------------------------------------------
## ----- get MAE data and a call list containing character X and Y, check for validity of X (assay name) and Y (assay/coldata name)
## and return matrices of X and Y in mc$X and mc$Y
#' @importFrom MultiAssayExperiment complete.cases
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment assay

names2mat <- function(mc)({
  tryCatch({
    if(any(class(c(mc$X, mc$Y))!="character")) stop("mc must have character X and Y", call. = FALSE)
  }, error=function(e) message("class(mc$X/mc$Y) produced error"))
  if(!mc$X %in% names(assays(mc$data))) stop_custom("inv_xy", paste0(mc$X, " is not a valid assay from 'data'"))
  if(mc$Y %in% names(colData(mc$data))){
    if(mc$Y %in% names(assays(mc$data))) stop_custom("inv_xy", paste0(mc$Y, " matches to both colData and assay in 'data', change its name in one and continue."))
    ## ----- if Y is a colData column subset it using X samples
    Xcoldata <- suppressMessages( as.data.frame(colData(mc$data[,,mc$X]))) ## keep X assay coldata - DataFrame to data.frame
    mc$Y <- Xcoldata[,mc$Y] ## keep the coldata desired for Y
    mc$Y <- as.matrix(mc$Y) ## convert to matrix for spls
    ## if Y not numeric
    if(! typeof(mc$Y) %in% c("numeric","integer")){
      if(typeof(mc$Y)=="factor") {
        warning("'Y' is factor, coercing to a numeric...")
        mc$Y <- as.numeric(mc$Y)
      } else {
        stop_custom(.subclass = "inv_xy", message = paste0(" 'Y' is not a numeric/integer column (or a factor coercible to numeric)"))
      }

    }
    ## if all is well with Y
    mc$X <- as.matrix(assay(mc$data, mc$X))
    ## ----- If Y is assay name
  } else if(mc$Y %in% names(assays(mc$data))){
    mc$data <- mc$data[,complete.cases(mc$data[,,c(mc$X, mc$Y)])] ## keep complete data b/w two assays
    mc$X <- as.matrix(assay(mc$data, mc$X))
    mc$Y <- as.matrix(assay(mc$data, mc$Y))
  } else {stop_custom("inv_xy", paste0(mc$Y, " is not an assay or colData from 'data' " ))}
  mc$data <- NULL
  mc$X <- t(mc$X)
  if(!1 %in% dim(mc$Y)){ ## if Y is matrix transpose it
    mc$Y <- t(mc$Y)
  }
  return(mc)
})

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
    if(any(sapply(as.list(mc$formula), length)!=1)) stop_custom("inv_formula", "formula must be of form: Y~X")
    mc$X <- eval.parent(fterms[[2]])
    mc$Y <- as.matrix(eval.parent(fterms[[1]]))
    mc$formula <- NULL
  }
  ##------- if formula = coldata/assay ~ assay and data=MAE provided
  else if(method.mode=="formula_mae"){
    mc$data <- eval.parent(mc$data)
    mc$formula <- eval.parent(mc$formula)
    fterms <- as.list(mc$formula)[-1]
    ## if the length of each of the formula elements (LHS and RHS) is anything other than 1, stop
    if(any(sapply(fterms, length)!=1)) stop_custom("inv_formula", "formula must be of form: Y~X")
    ## so that formula can also be stored in a variable - makes testing easier as well
    mc[c("Y","X")] <- as.character(fterms)
    mc <- names2mat(mc = mc)
    mc$formula <- mc$data <-  NULL
  }
  ##------- if X=assay (symbol or char), Y=assay/colData (symbol or char) and data=MAE provided
  else if(method.mode=="xy_mae"){
    mc$data <- eval.parent(mc$data)
    ## if X and Y are variables
    mc$X <- eval.parent(mc$X)
    mc$Y <- eval.parent(mc$Y)
    mc <- names2mat(mc = mc)
    mc$data <-  NULL
  }
  do.call(.pls, mc)
}
