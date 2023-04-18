library(GEOquery)
library(affy) #to read all .cel file
library(limma)
library(Biobase)
library(gplots)

#3. Import GSE file to extract phenoData
bc <- getGEO("GSE137178")
bc #this is preprocess file

phenoData_bc <-phenoData(bc[[1]])
phenoData_bc

#4. Read .cel files into the environment.
bc_cel <-ReadAffy()
bc_cel #this is the raw file

#5. Normalise
eset_bc <-rma(bc_cel) #this is same like bc_cel
eset_bc

#6. Extract into expression dataset
exprs_bc <-exprs(eset_bc)
exprs_bc

#7. Variables in the phenoData, use limma
pheno_bc <-pData(phenoData_bc)

#8. Extract out column of variable needed
# in pheno_bc last column, CXCR3i is disease, DMSO is normal
#extract dubset from data frame, use square bracket
pData_bc <- pheno_bc[,35] #here use all rows and column 35
# here, it is string data type

#9. Describe model to be fitted
design_bc <- model.matrix(~pData_bc)
design_bc

#10. Fit each probeset to model
fit_bc <- lmFit(eset_bc, design_bc)
fit_bc 

#11. Empirical Bayes adjustment
efit_bc <- eBayes(fit_bc)
efit_bc #

#12. Select the top 50 genes, limma analysis
limma_bc <-topTable(efit_bc,number=50,coef=NULL,sort.by="P")
#logFC differentiate genes based on given disease (like teruk ke tak)
#downregulated, '-' sign in logFC.. upregulated, positive for logFC


#13. Extract out the gene names of differential expressed genes 
#from limma table formed in step 12.
gname <- rownames(limma_bc)
head(gname)

#14. Attach the expression values of each genes 
#for each sample to form another dataset.
bc_50 <- exprs_bc[c(gname),]
head(bc_50)

# 15. Conduct heatmap clustering, use gplots library
# (note: refer to the heatmap function 
# to play around with the scale, title, colour and etc.)
heatmap(bc_50)
# if mix, that two group are related to each other. if not, 
# to see expression of gene.. whether upregulated or downregulated. 
# darker colour is the genes are upregulated
# lighter colour, the genes regulated

# 16. questions
#explain about the heatmap.. 
#what colours represent up/down regulated