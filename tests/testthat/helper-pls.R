library(mixOmics.data)
library(MultiAssayExperiment)

mae_data <- miniACC

X_index <- 1 ## RNASeq2GeneNorm
Y_index_assay <- 2 ## gistict
Y_index_coldata <- 2 ## "years_to_birth" - must be numeric
Y_index_invalid_coldata <- 1 ## "patientID" - must be numeric
Y_index_char_coldata <- 7 ## pathologic_stage" - is in fact a valid class but not factor
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

Y_char <-colnames(colData(mae_data))[Y_index_char_coldata]
Y_char.vec <- colData(mae_data)[,Y_index_char_coldata]
