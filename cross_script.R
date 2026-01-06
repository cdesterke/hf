list.files()


load("voomGSE116250.rda")
gse116250<-edata
rm(edata)

load("voomGSE135055.rda")
gse135055<-edata
rm(edata)


cross2<-merge(gse116250,gse135055,by="row.names")

row.names(cross2)<-cross2$Row.names
cross2$Row.names<-NULL


load("voomGSE141910.rda")
gse141910<-edata
rm(edata)


cross3<-merge(cross2,gse141910,by="row.names")
row.names(cross3)<-cross3$Row.names
cross3$Row.names<-NULL

pheno<-read.table("phenotype_all.csv",h=T,sep="\t",row.names=2)

all(row.names(pheno)==colnames(cross3))

library(transpipe15)
pcatrans(cross3,pheno,group="gender",alpha=1,names=F)


library(sva)
batch = pheno$batch
mod = model.matrix(~1+as.factor(group), data=pheno)

edata = ComBat(dat=cross3,  batch=batch, mod=mod,par.prior=TRUE, prior.plots=TRUE)



pcatrans(edata,pheno,group="gender",alpha=1,names=F)

###ESET

matrix<-as.matrix(edata)
library(Biobase)
all(row.names(pheno)==colnames(matrix))
phenoData <- AnnotatedDataFrame(data=pheno)
eset <- ExpressionSet(matrix, phenoData=phenoData)

save(eset,file="HeartFailure_Voom_Combat_ESET.rda")


annot<-pData(eset)
data<-exprs(eset)





