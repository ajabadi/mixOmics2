## ---- custom stop to define specific error classes ----
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
## ---- create custom warnings with specified class ----
warning_custom <- function(.subclass, message, call=NULL, ...){
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

## ------- put the character b/w two single quotation marks ----
squote <- function(char) paste0("'",char,"'")

## ------------------------------------------
## ----- get MAE data and a call list containing character X and Y, check for validity of X (assay name) and Y (assay/coldata name)
## and return matrices of X and Y in mc$X and mc$Y
#' @importFrom MultiAssayExperiment complete.cases
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment assay

names2mat <- function(mc)({ ## mc is list(data=MAE, X=X_name, Y=Y_name)
  mcc <- mc
  tryCatch({
    if(any(class(c(mc$X, mc$Y))!="character")) stop("mc must have character X and Y", call. = FALSE)
  }, error=function(e) message("class(mc$X/mc$Y) produced error"))
  if(!mc$X %in% names(assays(mc$data))) stop_custom("inv_xy", paste0(mcc$X, " is not a valid assay from 'data'"))
  ## --- if 'Y' is a colData
  if(mc$Y %in% names(colData(mc$data))){
    if(mc$Y %in% names(assays(mc$data))) stop_custom("inv_xy", paste0(mcc$Y, " matches to both colData and assay in 'data', change its name in one and continue."))
    ## ----- if Y is a colData column subset it using X samples
    Xcoldata <- suppressMessages( as.data.frame(colData(mc$data[,,mc$X]))) ## keep X assay coldata - DataFrame to data.frame
    mc$Y <- Xcoldata[,mc$Y] ## keep the coldata desired for Y
    ## if Y not numeric
    if(! typeof(mc$Y) %in% c("numeric","integer")){
      if(typeof(mc$Y)=="factor") {
        warning(paste0("The column data ", squote(mcc$Y)," is a factor, coercing to a numeric vector with names..."))
        mc$Y <- structure(as.numeric(mc$Y),
                          names=as.character(mc$Y), class="numeric")
      } else if(typeof(mc$Y)=="character"){
        ## if Y is a character colData and the number of unique terms are less than total,
        ## coerce it to factor and then numeric with a warning
        if(length(unique(mc$Y)) <  length(mc$Y) ){
          warning_custom("char_Y", message = paste0("The column data ",mcc$Y, " is character vector, coercing to factor and then named numeric for pls"))

          mc$Y <- structure(as.numeric(as.factor(mc$Y)),
                            names=mc$Y, class="numeric")
        } else {
          stop_custom(.subclass = "inv_xy", message = paste0(" 'Y' is not a numeric/integer column (or a factor coercible to numeric)"))
        }
      }

    }
    ## if all is well with Y
    mc$X <- assay(mc$data, mc$X)
    ## ----- If Y is assay name
  } else if(mc$Y %in% names(assays(mc$data))){
    mc$data <- mc$data[,complete.cases(mc$data[,,c(mc$X, mc$Y)])] ## keep complete data b/w two assays
    mc$X <- assay(mc$data, mc$X)
    mc$Y <- assay(mc$data, mc$Y)
  } else {stop_custom("inv_xy", paste0(squote(mcc$Y), " is not an assay or column data from the MAE object" ))}
  mc$data <- NULL
  mc$X <- t(as.matrix(mc$X))
  mc$Y <- as.matrix(mc$Y)

  if(!1 %in% dim(mc$Y)){ ## if Y is matrix transpose it
    mc$Y <- t(mc$Y)
  }
  return(mc)
})

##------------------------------------------
## ----- function to spls methods arguments plus the methods call mode and create internal level arguments
