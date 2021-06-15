setwd<-'/local/workdir/mh865/GCA_SCA/CVData1920/output'
phenotype.file <- '../data/CrossMerge1920.csv'          # path to pehnotype file
mm.file        <- '../data/A.csv'          # path to covariates file
weight.file     <- NULL    # path to weight file if NULL no weights are used


ctr 		<- TRUE # TRUE/FALSE		# Valid only if mm.file is not NULL
std 		<- TRUE # TRUE/FALSE		# Valid only if mm.file is not NULL
weighting 	<- !is.null(weight.file) 	# Valid only if mm.file is not NULL



colIDy <- 3  # the P1 column. column in phenotype file that gives the IDs that link observations to covariates or grouping factor
# or NULL if ID is not used
# 4 #(basicallythe P2 column ---- 1126Notes)

