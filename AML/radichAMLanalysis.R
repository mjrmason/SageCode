#amlSexVerification

library("data.table")
library("synapseClient")
library("parallel")
library("beadarray")
library("illuminaHumanv4.db")
library("beeswarm")
library(randomForest)
library(foreach)
library("doParallel")
library(gplots)
library(limma)
library(glmnet)
library(AUC)
synapseLogin()


# snp data
snpFiles = list.files("/home/apratap/projects/AML/external_data/dataset1-Raddich-lab/CARDINAL/SNP",full.names=T)
anno     = data.frame(fread("~/SageData/RadichAMLsnpData/HumanOmni5-4v1_C_Gene_Annotation.txt"))

xSnps = anno[anno$Chr == "X",1]
ySnps = anno[anno$Chr == "Y",1]

determineHomeoZ = function(x, xSnps, ySnps)
{
  tSnps = data.frame(fread(paste('zcat',x)))
  tXsnps = tSnps[as.character(tSnps[,1]) %in% xSnps,]
  tYsnps = tSnps[as.character(tSnps[,1]) %in% ySnps,]
  
  xhom = tXsnps[,"Allele1...AB"] == tXsnps[,"Allele2...AB"]
  yhom = tYsnps[,"Allele1...AB"] == tYsnps[,"Allele2...AB"] 
  yhom = yhom[tYsnps[,"Allele1...AB"] != "-" & tYsnps[,"Allele2...AB"] != "-"]
  
  return(list(xhom,yhom))
}


homs = mclapply(snpFiles, determineHomeoZ, xSnps, ySnps, mc.cores=10)


perchom = function(x,i){sum((x[[i]]))/length(x[[i]])}

percHx = sapply(homs, perchom,i = 1)
percHy = sapply(homs, perchom,i = 2)

##### expression data
load("~/SageData/BSDataWnWOweighting.Rda")
amlWeighted   = exprs(normaliseIllumina(BSDataWwts, transform = "none")); 
colnames(amlWeighted)   = gsub("_.*$","",colnames(amlWeighted))

amlClin = synTableQuery("SELECT * FROM syn5481721")@values
amlClin = amlClin[ !is.na(amlClin$PatientID) & amlClin$PatientID %in% colnames(amlWeighted),]
amlWeighted = amlWeighted[,match(as.character(amlClin$PatientID),colnames(amlWeighted))]

genes = as.data.frame(illuminaHumanv4SYMBOL)
chrs  = as.data.frame(illuminaHumanv4CHR); chrs = chrs[match(genes$probe_id, chrs$probe_id),]
genes =  genes[!is.na(chrs$chromosome),]
chrs  =  chrs[!is.na(chrs$chromosome),]


amlWeighted = amlWeighted[match(genes$probe_id, rownames(amlWeighted)), ]

xist = amlWeighted[genes$symbol == "XIST",]
yExp = amlWeighted[ chrs$chromosome == "Y",]
yPC  = prcomp(yExp)$rotation[,1]

plot(xist, yPC, pch = 16, col ="slategrey", cex = .5)

snpNames = gsub("^.*SNP/|_I.*zip$|I.*zip$","",snpFiles)
expNames = colnames(amlWeighted)[colnames(amlWeighted) %in% snpNames]
xistTemp = xist[colnames(amlWeighted) %in% snpNames]
yPCTemp  = yPC[colnames(amlWeighted) %in% snpNames]

percHxT = percHx[match(expNames,snpNames)]
percHyT = percHy[match(expNames,snpNames)]

boxplot(xist~amlClin$Gender, border="white")
boxplot(xist~amlClin$Gender, outline=F, border="grey60",add=T)
beeswarm(xist~amlClin$Gender, col="slategrey", pch=16, cex=.5, add=T)

boxplot(yPC~amlClin$Gender, border="white")
boxplot(yPC~amlClin$Gender, outline=F, border="grey60",add=T)
beeswarm(yPC~amlClin$Gender, col="slategrey", pch=16, cex=.5, add=T)


pdf("sexChromosomeHomozygosityNexpression.pdf")
par(mfrow=c(2,2))
tGend = amlClin$Gender[match(expNames,amlClin$PatientID)]

#boxplot(percHxT~tGend, border="white")
#boxplot(percHxT~tGend, outline=F, border="grey60",add=T)
#beeswarm(percHxT~tGend, col="slategrey", pch=16, cex=.5, add=T)

#boxplot(percHyT~tGend, border="white")
#boxplot(percHyT~tGend, outline=F, border="grey60",add=T)
#beeswarm(percHyT~tGend, col="slategrey", pch=16, cex=.5, add=T)

cols = rep("cyan",length(expNames)); cols[tGend == "Female"] = "orange"
plot(percHxT, percHyT, pch=16, col=cols,cex=.65,xlab = "% X Homoz",ylab = "% Y Homoz")
plot(percHxT,xistTemp, pch=16, col=cols,cex=.65,xlab = "% X Homoz", ylab="XIST expr")
plot(percHyT,yPCTemp, pch=16, col=cols,cex=.65,xlab = "% Y Homoz", ylab="PC Y genes expr")
# patient ID 6350 looks mail by snps and femail by expression
dev.off()


#### checking for batch/tissue effect

amlClin$cleanAge = sapply(amlClin$Age, function(x){if(grepl("-",x)){sum(as.numeric(unlist(strsplit(x,"-"))))/2;}else{as.numeric(x)}})

library(impute)

amlWeighted = impute.knn(amlWeighted)$data  

pcs = prcomp(amlWeighted,scale=T)


boxplot(pcs$rotation[,2]~amlClin$Sample_type, border="white")
boxplot(pcs$rotation[,2]~amlClin$Sample_type, outline=F, border="grey60",add=T)
beeswarm(pcs$rotation[,2]~amlClin$Sample_type, col="slategrey", pch=16, cex=.5, add=T)



cols = rep("slategrey", nrow(pcs$rotation)); cols[amlClin$Sample_type == "PB"] = "orange"
plot(pcs$rotation[,1], pcs$rotation[,2], col=cols, pch=16 )

cols = rep("slategrey", nrow(pcs$rotation)); 
cols[amlClin$Cooperative_Group == "SWOG"]  = "orange"
cols[amlClin$Cooperative_Group == "ECOG"]  = "blue"
cols[amlClin$Cooperative_Group == "CALGB"] = "cyan"

plot(pcs$rotation[,1], pcs$rotation[,2], col=cols, pch=16 )



###### actual modeling... via limma
trainD  = amlWeighted[, amlClin$response_class != "NA" & amlClin$PatientID != "6350"]
trainC  = amlClin[amlClin$response_class != "NA" & amlClin$PatientID != "6350",]
trainC$response_class = gsub("-","",trainC$response_class)
mm       = model.matrix(~response_class+Sample_type+mRNA_array_batch+Cooperative_Group+cleanAge+Gender,data=trainC)
fits     = lmFit(trainD,mm)
fits     = eBayes(fits)
ttResp   = topTable(fits, coef=2, adjust="BH", number=nrow(trainD), genelist=genes$symbol, sort.by="none")
ttTiss   = topTable(fits, coef=3, adjust="BH", number=nrow(trainD), genelist=genes$symbol, sort.by="none")
ttAge    = topTable(fits, coef=9, adjust="BH", number=nrow(trainD), genelist=genes$symbol, sort.by="none")
ttSex    = topTable(fits, coef=10, adjust="BH", number=nrow(trainD), genelist=genes$symbol, sort.by="none")
ttBatch1 = topTable(fits, coef=4, adjust="BH", number=nrow(trainD), genelist=genes$symbol, sort.by="none")
ttBatch2 = topTable(fits, coef=5, adjust="BH", number=nrow(trainD), genelist=genes$symbol, sort.by="none")
ttBatch3 = topTable(fits, coef=6, adjust="BH", number=nrow(trainD), genelist=genes$symbol, sort.by="none")
ttSite1  = topTable(fits, coef=7, adjust="BH", number=nrow(trainD), genelist=genes$symbol, sort.by="none")
ttSite2  = topTable(fits, coef=8, adjust="BH", number=nrow(trainD), genelist=genes$symbol, sort.by="none")

# checked if removing samples with low RIN or the 6350 sex anomolous sample had any improvement... it did not

# checked if weighting via voom helps... results are essentially uneffected

# roughly 70 genes have BH adj p < .01 for responder status... these genes show no gene ontology enrichment. (via a quick DAVID EASE analys)


###### use limma modeling to adjust for batch and site before creating classifier, 
trainD  = amlWeighted[, amlClin$response_class != "NA"]
trainC  = amlClin[amlClin$response_class != "NA",]
trainC$response_class = gsub("-","",trainC$response_class)
mm      = model.matrix(~mRNA_array_batch+Cooperative_Group,data=trainC)
tfits   = lmFit(trainD,mm)
tfits   = eBayes(tfits)
res     = residuals(tfits,trainD)

indParam = data.frame(t(res), tissue = trainC$Sample_type,age= trainC$cleanAge, sex = trainC$Gender)
y        = as.factor(trainC$response_class)

cl       = makeCluster(4); registerDoParallel(cl)
set.seed(13)
trainInds = sample(1:nrow(tDat), size=floor(.8*nrow(tDat)))
rf       = foreach(ntree=rep(2500, 4), .combine=combine, .packages='randomForest') %dopar% randomForest(x=indParam[trainInds,], y=y[trainInds], ntree=ntree)

plot( rf$importance[1:(nrow(rf$importance)-3),1], -log10(ttResp$adj.P.Val), pch = 16, col = "slategrey",cex =.5)


SD = apply(trainD,1,sd)
tInds = SD > .3
rf2    = foreach(ntree=rep(1000, 4), .combine=combine, .packages='randomForest') %dopar% randomForest(x=indParam[trainInds,c(tInds,rep(T,3))], y=y[trainInds], ntree=ntree)

plot( rf2$importance[1:(nrow(rf2$importance)-3),1], -log10(ttResp$adj.P.Val)[tInds], pch = 16, col = "slategrey",cex =.5)


out = data.frame(ttResp[tInds,],rfImp=rf2$importance[1:(nrow(rf2$importance)-3),1])

pred = predict(rf2, newdata = indParam[,c(tInds,rep(T,3))], type="prob")
boxplot(pred[,2]~y,border="grey60",outline=F)
beeswarm(pred[,2]~y,col="slategrey",add=T, pch=16, cex =.65)







#### check out lasso and elastic net for prediction
resp   = trainC$Respons
sex    = factor(trainC$Gender)
tissue = factor(trainC$Sample_type)
site   = factor(trainC$Cooperative_Group)
batch  = factor(trainC$mRNA_array_batch)
xFacs  = model.matrix(resp~sex+tissue+batch+site)[,-1]
xFacsSimp  = model.matrix(resp~sex+tissue)[,-1]
age    =  as.numeric(trainC$cleanAge)

#tDat  = as.matrix(data.frame(t(trainD), age, xFacs))
#tInds = ttResp$P.Value < .05 | ttTiss$P.Value < .05 | ttAge$P.Value < .05 | ttSex$P.Value < .05 
SD = apply(trainD,1,sd)
tInds = SD > .3
tDat  = as.matrix(data.frame(t(trainD[tInds,]), age, xFacs))
tDatSimp  = as.matrix(data.frame(t(trainD[tInds,]), age, xFacsSimp))
set.seed(13)
trainInds = sample(1:nrow(tDat), size=floor(.8*nrow(tDat)))


alphas = seq(0,1,by=.1)

glmNet = function(a,x, y, Seed)
{
  set.seed(Seed)
  tMod   = cv.glmnet(x, y, family="binomial", type.measure="auc", nfolds = 5,alpha=a)
  b      = tMod$glmnet.fit$beta; nparams = apply(b != 0, 2,sum)
  retMat = cbind(tMod$lambda,tMod$cvm, a, nparams[1:length(tMod$lambda)]); #retMat = retMat[-1,]
  return(retMat)
}

pdf("classifiersComparison.pdf")

modList1 = mclapply(alphas,glmNet, x = tDat[trainInds,], y = as.factor(resp[trainInds]),Seed = 5, mc.cores=10)
modList1Simp = mclapply(alphas,glmNet, x = tDatSimp[trainInds,], y = as.factor(resp[trainInds]),Seed = 5, mc.cores=10)
modList2 = mclapply(alphas,glmNet, x = tDat[trainInds,], y = as.factor(resp[trainInds]),Seed = 7, mc.cores=10)
modList3 = mclapply(alphas,glmNet, x = tDat[trainInds,], y = as.factor(resp[trainInds]),Seed = 13, mc.cores=10)
mL = modList3
plot(mL[[2]][-1,4],mL[[2]][-1,2], type = "l",ylim=c(.75,.92), xlab = "Number of parameters", ylab = "auc")
for(i in 3:length(mL)){lines(mL[[i]][-1,4],mL[[i]][-1,2], type = "l", col=i);print(alphas[i]);}

aucs1 = as.data.frame(lapply(modList1,function(x){return(x[,2])})); colnames(aucs1) = 1:ncol(aucs1)
aucs2 = as.data.frame(lapply(modList2,function(x){return(x[,2])})); colnames(aucs2) = 1:ncol(aucs2)
aucs3 = as.data.frame(lapply(modList3,function(x){return(x[,2])})); colnames(aucs3) = 1:ncol(aucs3)

heatmap.2(as.matrix(aucs1), Rowv=F, Colv=F, scale="none", trace="none", dendrogram="none", labCol=alphas, xlab = "alpha",ylab="lambda index")
heatmap.2(as.matrix(aucs2), Rowv=F, Colv=F, scale="none", trace="none", dendrogram="none", labCol=alphas, xlab = "alpha",ylab="lambda index")
heatmap.2(as.matrix(aucs3), Rowv=F, Colv=F, scale="none", trace="none", dendrogram="none", labCol=alphas, xlab = "alpha",ylab="lambda index")

x = tDat[trainInds,]; y = as.factor(resp[trainInds])
par(mfrow=c(2,2))
set.seed(7)
tMod   = cv.glmnet(x, y, family="binomial", type.measure="auc", nfolds = 8,alpha=0)
plot(tMod, main="alpha = 0")
set.seed(7)
tMod   = cv.glmnet(x, y, family="binomial", type.measure="auc", nfolds = 8,alpha=.3)
plot(tMod, main="alpha = .3")
set.seed(7)
tMod   = cv.glmnet(x, y, family="binomial", type.measure="auc", nfolds = 8,alpha=.6)
plot(tMod, main="alpha = .6")
set.seed(7)
tMod   = cv.glmnet(x, y, family="binomial", type.measure="auc", nfolds = 8,alpha=1)
plot(tMod, main="alpha = 1")
dev.off()

finModel  =  glmnet(x,y,alpha=.3,lambda=.27, family="binomial")
predicted = predict(finModel, newx=x, type="class")
table(predicted, y)
AUC  = auc(accuracy(predicted, y))

predicted = predict(finModel, newx=tDat[-trainInds,], type="class")
table(predicted, as.factor(resp[-trainInds]))
AUC  = auc(accuracy(predicted, as.factor(resp[-trainInds]))) # the auc is .81 for the small test set


##################################################
# checking if clinical info is skewed in test set
##################################################

amlClin$white = F; amlClin$white[amlClin$Race %in% c("White", "Caucasian")] = T
varOfInt = c("Gender", "Race","Cooperative_Group","Sample_type","RIN_number","Thaw_Batch","response_class","cleanAge","white")
is.na(amlClin$response_class)

amlClin$dataset = "train"; amlClin$dataset[amlClin$response_class == "NA"] = "test"
fisher.test(table(amlClin$dataset, amlClin$Sample_type))

tests = list()
for(var in varOfInt)
{
  if(class(amlClin[,var]) == "numeric")
  {
    tests[[var]] = wilcox.test(amlClin[amlClin$dataset == "train",var],amlClin[amlClin$dataset == "test",var])$p.value
  }
  else{tests[[var]] = fisher.test(table(amlClin$dataset, amlClin[,var]))$p.value}
}

# there are no strong misbalances in clinical parameters between the training and test sets

##################################################
# gene ontology enrichment
##################################################

library(package = "illuminaHumanv4.db", character.only = TRUE)
library(GOstats)

entrez = as.data.frame(illuminaHumanv4ENTREZID)
dat    = ttResp; dat$entrez = entrez[match(rownames(dat),entrez[,1]),2];

paramsBPup = new("GOHyperGParams",geneIds=dat$entrez[dat$adj.P.Val < .05 & dat$t >0], universeGeneIds=dat$entrez,annotation="illuminaHumanv4.db", ontology="BP", pvalueCutoff=0.01, conditional=T,testDirection="over")
paramsBPlo = new("GOHyperGParams",geneIds=dat$entrez[dat$adj.P.Val < .05 & dat$t <0], universeGeneIds=dat$entrez,annotation="illuminaHumanv4.db", ontology="BP", pvalueCutoff=0.01, conditional=T,testDirection="over")
paramsBP   = new("GOHyperGParams",geneIds=dat$entrez[dat$adj.P.Val < .05], universeGeneIds=dat$entrez,annotation="illuminaHumanv4.db", ontology="BP", pvalueCutoff=0.01, conditional=T,testDirection="over")
paramsCCup = new("GOHyperGParams",geneIds=dat$entrez[dat$adj.P.Val < .05 & dat$t >0], universeGeneIds=dat$entrez,annotation="illuminaHumanv4.db", ontology="CC", pvalueCutoff=0.01, conditional=T,testDirection="over")
paramsCClo = new("GOHyperGParams",geneIds=dat$entrez[dat$adj.P.Val < .05 & dat$t <0], universeGeneIds=dat$entrez,annotation="illuminaHumanv4.db", ontology="CC", pvalueCutoff=0.01, conditional=T,testDirection="over")
paramsCC   = new("GOHyperGParams",geneIds=dat$entrez[dat$adj.P.Val < .05], universeGeneIds=dat$entrez,annotation="illuminaHumanv4.db", ontology="CC", pvalueCutoff=0.01, conditional=T,testDirection="over")
paramsMFup = new("GOHyperGParams",geneIds=dat$entrez[dat$adj.P.Val < .05 & dat$t >0], universeGeneIds=dat$entrez,annotation="illuminaHumanv4.db", ontology="MF", pvalueCutoff=0.01, conditional=T,testDirection="over")
paramsMFlo = new("GOHyperGParams",geneIds=dat$entrez[dat$adj.P.Val < .05 & dat$t <0], universeGeneIds=dat$entrez,annotation="illuminaHumanv4.db", ontology="MF", pvalueCutoff=0.01, conditional=T,testDirection="over")
paramsMF   = new("GOHyperGParams",geneIds=dat$entrez[dat$adj.P.Val < .05], universeGeneIds=dat$entrez,annotation="illuminaHumanv4.db", ontology="MF", pvalueCutoff=0.01, conditional=T,testDirection="over")
goTestBPup = summary(hyperGTest(paramsBPup))
goTestBPlo = summary(hyperGTest(paramsBPlo))
goTestBP   = summary(hyperGTest(paramsBP))
goTestCCup = summary(hyperGTest(paramsCCup))
goTestCClo = summary(hyperGTest(paramsCClo))
goTestCC   = summary(hyperGTest(paramsCC))
goTestMFup = summary(hyperGTest(paramsMFup))
goTestMFlo = summary(hyperGTest(paramsMFlo))
goTestMF   = summary(hyperGTest(paramsMF))
#probeLis   = probeSetSummary(goTestBPup, ids = "ENTREZID", sigProbesets=rownames(dat)[dat$adj.P.Val < .05]) #example of getting the probes in each GO group

plotGO = function(x,N=10,bMar=5, yTop = 5,main="")
{
  x = x[1:(min(nrow(x),N)),]
  par(mar=c(c(bMar, 4, 4, 2) + 0.1))
  barplot(-log10(x$Pvalue), names.arg = x$Term, las=3,cex.names = .65,border =  "white", col="slategrey", ylim=c(0,yTop), ylab ="-log10(p)", main=main)
}

paramsKEGGup = new("KEGGHyperGParams",geneIds=dat$entrez[dat$adj.P.Val < .05 & dat$t >0], universeGeneIds=dat$entrez,annotation="illuminaHumanv4.db",  pvalueCutoff=0.05,testDirection="over")
paramsKEGGlo = new("KEGGHyperGParams",geneIds=dat$entrez[dat$adj.P.Val < .05 & dat$t <0], universeGeneIds=dat$entrez,annotation="illuminaHumanv4.db",  pvalueCutoff=0.05,testDirection="over")
paramsKEGG   = new("KEGGHyperGParams",geneIds=dat$entrez[dat$adj.P.Val < .05], universeGeneIds=dat$entrez,annotation="illuminaHumanv4.db",  pvalueCutoff=0.05,testDirection="over")
keggTestup   = summary(hyperGTest(paramsKEGGup))
keggTestlo   = summary(hyperGTest(paramsKEGGlo))
keggTest     = summary(hyperGTest(paramsKEGG))

pdf("~/SageOutput/AML/geneontologyBarPlots.pdf")
plotGO(goTestBPup, 20, 12, main="BP up")
plotGO(goTestBPlo, 20, 12, main="BP lo")
plotGO(goTestBP, 20, 12, main="BP both")
plotGO(goTestCCup, 20, 12, main="CC up")
plotGO(goTestCClo, 20, 12, main="CC lo")
plotGO(goTestCC, 20, 12, main="CC both")
plotGO(goTestMFup, 20, 12, main="MF up")
plotGO(goTestMFlo, 20, 12, main="MF lo")
plotGO(goTestMF, 20, 12, main="MF both")
plotGO(keggTestup, 20, 12, main="KEGG up")
plotGO(keggTestlo, 20, 12, main="KEGG lo")
plotGO(keggTest, 20, 12, main="KEGG both")
dev.off()

# gene ontology of model
coefs = as.matrix(coefficients(finModel))
coefs      = coefs[grep("^ILMN", rownames(coefs)),]
entrez = as.data.frame(illuminaHumanv4ENTREZID)
entrez = entrez[match(names(coefs),entrez[,1]),] 

paramsBPup = new("GOHyperGParams",geneIds=entrez[coefs > 0,2], universeGeneIds=entrez[,2],annotation="illuminaHumanv4.db", ontology="BP", pvalueCutoff=0.01, conditional=T,testDirection="over")
paramsBPlo = new("GOHyperGParams",geneIds=entrez[coefs > 0,2], universeGeneIds=entrez[,2],annotation="illuminaHumanv4.db", ontology="BP", pvalueCutoff=0.01, conditional=T,testDirection="over")
paramsBP   = new("GOHyperGParams",geneIds=entrez[coefs != 0,2], universeGeneIds=entrez[,2],annotation="illuminaHumanv4.db", ontology="BP", pvalueCutoff=0.01, conditional=T,testDirection="over")
paramsCCup = new("GOHyperGParams",geneIds=entrez[coefs < 0,2], universeGeneIds=entrez[,2],annotation="illuminaHumanv4.db", ontology="CC", pvalueCutoff=0.01, conditional=T,testDirection="over")
paramsCClo = new("GOHyperGParams",geneIds=entrez[coefs > 0,2], universeGeneIds=entrez[,2],annotation="illuminaHumanv4.db", ontology="CC", pvalueCutoff=0.01, conditional=T,testDirection="over")
paramsCC   = new("GOHyperGParams",geneIds=entrez[coefs != 0,2], universeGeneIds=entrez[,2],annotation="illuminaHumanv4.db", ontology="CC", pvalueCutoff=0.01, conditional=T,testDirection="over")
paramsMFup = new("GOHyperGParams",geneIds=entrez[coefs < 0,2], universeGeneIds=entrez[,2],annotation="illuminaHumanv4.db", ontology="MF", pvalueCutoff=0.01, conditional=T,testDirection="over")
paramsMFlo = new("GOHyperGParams",geneIds=entrez[coefs > 0,2], universeGeneIds=entrez[,2],annotation="illuminaHumanv4.db", ontology="MF", pvalueCutoff=0.01, conditional=T,testDirection="over")
paramsMF   = new("GOHyperGParams",geneIds=entrez[coefs != 0,2], universeGeneIds=entrez[,2],annotation="illuminaHumanv4.db", ontology="MF", pvalueCutoff=0.01, conditional=T,testDirection="over")
paramsKEGGup = new("KEGGHyperGParams",geneIds=entrez[coefs > 0,2], universeGeneIds=entrez[,2],annotation="illuminaHumanv4.db",  pvalueCutoff=0.05,testDirection="over")
paramsKEGGlo = new("KEGGHyperGParams",geneIds=entrez[coefs < 0,2], universeGeneIds=entrez[,2],annotation="illuminaHumanv4.db",  pvalueCutoff=0.05,testDirection="over")
paramsKEGG   = new("KEGGHyperGParams",geneIds=entrez[coefs != 0,2], universeGeneIds=entrez[,2],annotation="illuminaHumanv4.db",  pvalueCutoff=0.05,testDirection="over")
goTestBPup = summary(hyperGTest(paramsBPup))
goTestBPlo = summary(hyperGTest(paramsBPlo))
goTestBP   = summary(hyperGTest(paramsBP))
goTestCCup = summary(hyperGTest(paramsCCup))
goTestCClo = summary(hyperGTest(paramsCClo))
goTestCC   = summary(hyperGTest(paramsCC))
goTestMFup = summary(hyperGTest(paramsMFup))
goTestMFlo = summary(hyperGTest(paramsMFlo))
goTestMF = summary(hyperGTest(paramsMF))
keggTestup = summary(hyperGTest(paramsKEGGup))
keggTestlo = summary(hyperGTest(paramsKEGGlo))
keggTest = summary(hyperGTest(paramsKEGG))

pdf("~/SageOutput/AML/modelGeneontologyBarPlots.pdf")
plotGO(goTestBPup, 20, 12, main="BP up")
plotGO(goTestBPlo, 20, 12, main="BP lo")
plotGO(goTestBP, 20, 12, main="BP both")
plotGO(goTestCCup, 20, 12, main="CC up")
plotGO(goTestCClo, 20, 12, main="CC lo")
plotGO(goTestCC, 20, 12, main="CC both")
plotGO(goTestMFup, 20, 12, main="MF up")
plotGO(goTestMFlo, 20, 12, main="MF lo")
plotGO(goTestMF, 20, 12, main="MF both")
plotGO(keggTestup, 20, 12, main="KEGG up")
plotGO(keggTestlo, 20, 12, main="KEGG lo")
plotGO(keggTest, 20, 12, main="KEGG both")
dev.off()


### INPP4B looks super important in this classifier based on literature

cvfinModel  =  cv.glmnet(x,y,alpha=.3, family="binomial")
