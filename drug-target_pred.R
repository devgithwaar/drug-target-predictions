# Purpose: Predicting drug-target interactions via R package "BioMedR"

# Author: Hui-Heng Lin, PhD . On 14th Aug.2020


# Original tutorial codes from PDF of "BioMedR: R/CRAN Package for generating various molecular representations for chemicals, proteins DNAs/RNAs and their interactions. ( BioMedR manual authored by Minfeng Zhu, Jie Dong, Dongsheng Cao, Package version: Release 1. 2019-07-03 "

# Annotations were added by Hui-Heng Lin.
require(BioMedR) # load library
gpcr=read.table(system.file('vignettedata/GRCR.csv', package= 'BioMedR') # load embedded dataset
protid=unique(gpcr[,1]); drugid=unique(gpcr[,2]) # protein and drug ID assignments to variable
protseq=BMgetProtSeqKEGG(protid, parallel=5) # retrieve protein seqeuence data from remote database, using IDs. Network connection is required

drugseq=c() # create a NULL variable for storging SMILES data of drugs
for (id in 1:length(drugid)){drugseq[id]=BMgetDrugSmiKEGG(drugid[id])} # remote retrieval of drug SMILES via drug IDs


"""
Alternatively, below codes could be used if your network connection is not good

protseq = readFASTA(system.file('vignettedata/GPCR_seq.fasta', package = 'BioMedR'))

drugseq = as.vector(read.table(system.file('vignettedata/GPCR_smi.txt', package = 'BioMedR'), col.names = 'SMILES'))

"""                


x0.prot=cbind(t(sapply(unlist(protseq), extrProtMoreauBroto)), t(sapply(unlist(protseq), extrProtCTDC))) # combine two featuresets of protein sequences

x0.drug=cbind(extrDrugGraphComplete(readMolFromSmi(textConnection(drugseq))), extrDrugPubChemComplete(readMolFromSmi(textConnection(drugseq)))) # executions failed. Functions had errors                
                
                

x.prot=matrix(NA, nrow = nrow(gpcr), ncol=ncol(x0.prot)) 
x.drug=matrix(NA, nrow = nrow(gpcr), ncol = ncol(x0.drug)) # creating matrices via specifying their row number and column numbers
for(i in 1:nrow(gpcr)) x.prot[i,] = x0.prot[which(gpcr[,1][i] == protid),]
for(i in 1:nrow(gpcr)) x.drug[i,]=x0.drug(gpcr[,2][i] == drugid)


y=as.factor(c(rep('1',nrow(gpcr)/2), rep('0', nrow(gpcr/2)))) # create label set

x=getCPI(x.prot, x.drug, type='combine') # generate drug-target interaction descriptors (combined descriptor matrix) using getCPI(). 
colnames(x) =paste('CCI', 1:dim(x)[2],sep='_')

require(caret)
x=x[,-nearZeroVar(x)]

# training set split
set.seed(20180808)
split_index=createDataPartition(y,p=0.75, list=F)

train_x=x[split_index,]
train_y=y[split_index]
test_x=x[-split_index,]
test_y=y[-split_index]


require(randomForest)
cv_result=rf.cv(train_x,train_y, cv.fold=5, type='classification', tree=500,mtry=30) # five-fold cross validation with random forest classifier

# train random forest classifier
rf.fit=randomForest(x=train_x, y=train_y, ntree=500, mtry=30, importance=TRUE)

# predict on the training set (for demonstration purpose only) and plot ROC curve
# predict on the test set (in fact the training set)
pre_res=predict(rf.fit,newdata=test_x, type='prob')[,2]

require(pROC) #  plot the Cross validation result and test result
require(RColorBrewer)
pal=brewer.pal(3, 'Set1')
opar<-par(no.readonly = TRUE)
par(mfrow=c(1,2))
plot.roc(train_y, cv_results$prob[,1], col=pal[2], grid=TRUE, print.auc=TRUE, main=' Cross Validation')
plotroc.(test_y, pre_res, col=pal[1], grid=T, print.auc=T, main='prediction')
par(opar)

                
                
                
                
                
                

# Debugging and modifications of above codes by Hui-Heng Lin


# References
""" > citation("caret") 

在出版物中使用程序包时引用‘caret’:

  Max Kuhn (2020). caret: Classification and
  Regression Training. R package version
  6.0-86.
  https://CRAN.R-project.org/package=caret

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {caret: Classification and Regression Training},
    author = {Max Kuhn},
    year = {2020},
    note = {R package version 6.0-86},
    url = {https://CRAN.R-project.org/package=caret},
  }

> citation("BioMedR")

在出版物中使用程序包时引用‘BioMedR’:

  Min-feng Zhu, Jie Dong and Dong-sheng Cao
  (2019). BioMedR: Generating Various
  Molecular Representations for Chemicals,
  Proteins, DNAs, RNAs and Their
  Interactions. R package version 1.2.1.
  https://CRAN.R-project.org/package=BioMedR

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {BioMedR: Generating Various Molecular Representations for Chemicals,
Proteins, DNAs, RNAs and Their Interactions},
    author = {Min-feng Zhu and Jie Dong and Dong-sheng Cao},
    year = {2019},
    note = {R package version 1.2.1},
    url = {https://CRAN.R-project.org/package=BioMedR},
  } """
