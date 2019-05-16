############# build all documentation for the package and update namespace
file.remove(list.files('man/', pattern = '.Rd', full.names = TRUE))
roxygen2::roxygenise()
############# load the pks to access functions
devtools::load_all()
############# run ./tests tests which then redirects to testthat.R
library(mixOmics.data)
library(MultiAssayExperiment)
mae_data <- miniACC

X_index <- 1 ## RNASeq2GeneNorm
Y_index_assay <- 2 ## gistict
Y_index_coldata <- 2 ## "years_to_birth" - must be numeric
Y_index_invalid_coldata <- 1 ## "years_to_birth" - must be numeric
## ----------------------------------------------------- no need to change the following

X <- names(assays(mae_data))[X_index] ## X assay
Ya <- names(assays(mae_data))[Y_index_assay] ## Y assay
f_Ya<- as.formula(paste0(Ya, " ~ ", X)) ## formula with Y assay

Yc <- names(colData(mae_data))[Y_index_coldata] ## Y column data
f_Yc <- as.formula(paste0(Yc, " ~ ", X)) ## formula with Y column data

Xm_Yc <- t(as.matrix(assay(mae_data, X))) ## X matrix when Y column data
Xm_Ya <- t(as.matrix(assay(mae_data[,complete.cases(mae_data[,,c(X, Ya)])], X))) ## X matrix when Y is assay

Yam <- t(assay(mae_data[,complete.cases(mae_data[,,c(X, Ya)])], Ya)) ## Y assay matrix
Ycn <-  as.numeric(colData(mae_data[,,X])[,Yc]) ## Y column data numeric

Y_inv <- colnames(colData(mae_data))[Y_index_invalid_coldata] ## invalid coldata y
Y_inv.vec <- colData(mae_data)[,Y_index_invalid_coldata] ## vector


devtools::test()
############# test examples
devtools::run_examples()
## path to Rd files - will let you know if there are any errors/warnings
# testthat::test_example("man/parent_base_ext.Rd")
############# manual tests
## test files that cannot be included in standard testthat pipeline
# invisible(lapply(list.files("tests/manual", full.names = TRUE), source))
############# check
# devtools::check()
############# build it
## by default, it is made in parent directoy of the package, we manually set it.
## testhat has changed wd
# setwd('/Users/alabadi/Documents/_Projects/_Personal/someWrappers')

  # binary_dir <- "../__binary"
  # if(!dir.exists(binary_dir)){
  #   dir.create(binary_dir)
  # }
  # pkg <- devtools::build(path = binary_dir)
  # ############# install it in my library
  # install.packages(pkg, repos = NULL)

# CMD + SHift + B
