#######################################################
## AZ drug combination challenge 2
## contains parallelized execution for Bayes Factor and
## empirical p-values
## Mike Mason (Sage),
#######################################################
rm(list=ls())
library(parallel)
library(data.table)
library(ppcor)
library(synapseClient)
library(plyr)
library("ROCR")
library("MatrixModels")

synapseLogin()

## load submission via file path (after downloading from synapse)
loadSubmission = function(x) # x is zipped file path
{
  cmd = paste('unzip -o \'',x,'\' -d \'', gsub('\\.zip$','',x),'\'',sep="")
  system(cmd)
  conf = read.csv(paste(gsub('\\.zip$','',x),'/confidence_matrix.csv',sep=""),row.names=1)
  pred = read.csv(paste(gsub('\\.zip$','',x),'/synergy_matrix.csv',sep=""),row.names=1)
  
  conf =  conf[order(rownames(conf)),]; conf =  conf[,order(colnames(conf))]; colnames(conf) = gsub("\\.","-",colnames(conf)); colnames(conf) = gsub("^X","",colnames(conf))
  pred =  pred[order(rownames(pred)),]; pred =  pred[,order(colnames(pred))]; colnames(pred) = gsub("\\.","-",colnames(pred)); colnames(pred) = gsub("^X","",colnames(pred))
  return(list("conf"=conf,"pred"=pred))
}

## get the subchallenge 2 predictions from each team submission (and confidence scores) from synapse
e           <- synGetEvaluation(4990389)
submissions <- synGetSubmissions(e@id, status='SCORED')
submissions <- synGetSubmissions(e@id, status='SCORED', limit = submissions@totalNumberOfResults)
subIds      <- sapply(submissions@results,function(x){x$id})
dat         <- lapply(subIds, function(x){print(paste("getting",x));synGetSubmission(x, downloadFile = T)})
fLocs2      <- sapply(dat,getFileLocation);
subData2    <- lapply(fLocs2, loadSubmission)


## get the true observations
obs          <- as.data.frame(fread("~/SageData/AZdrugComboDream/ch2_test.csv"))
obs          <- obs[, c("CELL_LINE", "COMBINATION_ID", "SYNERGY_SCORE")]
obs          <- obs[order(obs$COMBINATION_ID, obs$CELL_LINE),]
obs          <- obs[!obs$COMBINATION_ID == "AKT.AKT",] # filtering out self drug combo that results in NA values

## format the predictions
formatLikeObs <- function(x, y)
{
  calls  <- as.vector(unlist(x)); if(any(is.na(calls))){print(paste(sum(is.na(calls)), "predictions labled NA"))}
  calls[is.na(calls)] <- 0
  combos <- rep(rownames(x), ncol(x))
  cells  <- unlist(lapply(colnames(x), function(a){rep(a, nrow(x))}))
  df     <- data.frame(CELL_LINE = cells, COMBINATION_ID = combos, SYNERGY_SCORE = calls)
  df = df[match(paste(y$CELL_LINE,y$COMBINATION_ID,SEP="_"),paste(df$CELL_LINE,df$COMBINATION_ID,SEP="_")),]
  return(df[,3])
}

preds2        <- lapply(subData2, function(x){formatLikeObs(x[[2]],obs)}) 
confs2        <- lapply(subData2, function(x){formatLikeObs(x[[1]],obs)}) 
names(preds2) <- gsub("^.*/|\\.zip$","",fLocs2); names(confs2) <- names(preds2)
preds2        <- as.data.frame(preds2)
confs2        <- as.data.frame(confs2)

########################################################################################################
## with the data loaded and predictions in a single data frame, scoring metrics can be parallelized for 
## determining Bayes factor and empirical p-values
########################################################################################################

# --------------------------------------------------------------------------------------------
# main scoring function that enable parallelzation
# iteration number (iter) and indices (INDs) use as flags in addition to their main purpose
# --------------------------------------------------------------------------------------------

computeScores = function(iter=NA,INDs=NA,PREDS,OBS, bootStrap=F)
{
  if(!is.na(iter))
  {
    if(bootStrap == T){ PREDS = PREDS[INDs[,iter],]; }
    OBS  = OBS[INDs[,iter],]; 
    nr   = nrow(read.csv("./progressFile.csv"))/ncol(INDs); # for tracking progress
    print(paste(iter,": ", signif(nr,3),sep = ""))
  }
     
  resp       <- as.numeric(OBS$SYNERGY_SCORE)
  cl         <- as.factor(OBS$CELL_LINE)
  combo      <- as.factor(OBS$COMBINATION_ID)
  primary    <- apply(PREDS,2,function(X){X = as.integer(X); temp = coefficients(summary(lm(resp~0+cl+combo+X)));p = temp[rownames(temp) == "X",4];if(length(p)==0){return(0);}else{return(-log10(p))} }); # 3-way anova
  
  pos        <- resp > 20
  tieBreaker <- apply(PREDS,2,function(X){sens = sum(X==1 & pos)/sum(pos); spec = sum(X == 0 & !pos)/sum(!pos); return(mean(c(sens,spec)))})
  
  return(data.frame(primary=unlist(primary), tieBreaker))
}

# --------------------------------------------------------------------------------------------
# score challenge 2 submissions
# --------------------------------------------------------------------------------------------
date()
ch2Scores           = computeScores(NA, NA,PREDS=preds2, OBS=obs)
date()

# Bootstrapped scoring for computing Bayes factor
N        <- 1000 # should be 1000
cl       <- makeCluster(8,outfile='~/SageCode/progressFile.csv') # for tracking progress in shell, use "> watch tail -n 20 progressFile.csv"
clusterEvalQ(cl,library(plyr))
set.seed(13)

write.csv(matrix(,0,1), file = "./progressFile.csv", row.names=F,quote=F, col.names=F) # for tracking progress
bsPermInds <- matrix(1:nrow(obs), nrow(obs), N)
bsPermInds <- data.frame(apply(bsPermInds,2,sample, replace=T))
bsScores2   <- parLapply(cl,X= as.list(1:ncol(bsPermInds)), fun = computeScores, INDs = bsPermInds,PREDS=preds2, OBS = obs,bootStrap=T)

stopCluster(cl)
save(bsScores2, list=c('bsScores2'), file = "~/SageData/bsScores2.Rdat")

# compute Bayes Factor using the best primary score and boostrapped values above
bestI <- which(ch2Scores[,1] == max(ch2Scores[,1]))
M     <- as.data.frame(lapply(bsScores2,function(x){x[bestI,1] - x[,1]}))
BF2   <- apply(M,1,function(x){sum(x>0)/sum(x<=0)})

#################################################
# permutaion based empirical p-value computation
#################################################
N        <- 10000
cl       <- makeCluster(10,outfile='~/SageCode/progressFile.csv')
clusterEvalQ(cl,library(plyr))
set.seed(13)

write.csv(matrix(,0,1), file = "./progressFile.csv", row.names=F,quote=F, col.names=F) # for tracking progress
nullPermInds <- matrix(1:nrow(obs), nrow(obs), N)
nullPermInds <- data.frame(apply(nullPermInds,2,sample, replace=F)) # subtle difference from bootstrapping is the replace = F
bsScores2Null   <- parLapply(cl,X= as.list(1:ncol(nullPermInds)), fun = computeScores, INDs = nullPermInds,PREDS=preds2, OBS = obs, bootStrap=F)

stopCluster(cl)
save(bsScores2Null, list=c('bsScores2Null'), file = "~/SageData/bsScores2Null.Rdat")

# compute the empirical p-value using the permuted values above
nullPrimaryScoresMat <- as.matrix(as.data.frame(lapply(bsScores2Null, function(x){return(unlist(x[,1]))})))
nullPs <- apply((nullPrimaryScoresMat - ch2Scores[,1]) > 0,1, sum)/N
nullPs[nullPs == 0]  <- 1/N


####################################################
## some plotting
####################################################
pdf("~/SageOutput/AZdrugComboDream/subChallenge2FinalBayesFactor.pdf")

tBF = BF2; tBF[tBF == Inf] = 1.25*max(tBF[tBF != Inf])
plot( ch2Scores[,1], log10(tBF), pch =16, col = "slategrey", cex = 1, xlab = "Primary Scoring Metric", ylab = "log10 Bayes Factor", main="Subchallenge 2 R1")
abline(h= log10(5),col="cyan",lty="dashed",cex=3)
dev.off()

bsPrimaryScoresMat <- as.data.frame(bsPrimaryScoresMat)
colnames(bsPrimaryScoresMat) <- colnames(preds2)

pdf("~/SageOutput/AZdrugComboDream/subChallenge2FinalComboPlot.pdf", width=10,height=7 )
par(fig=c(0,.65, 0,1))
boxplot(bsPrimaryScoresMat[,Ord], horizontal=T, outline=F, border = "tan", col="tan", las = 2, cex.axis = .5, boxwex=.5,cex.main=.75, main="Primary Metric")
points(sort(ch2Scores[,1]), 1:ncol(bsPrimaryScoresMat), pch = 16, col = "slategrey", cex = .5)
par(fig=c(.55,.8, 0,1), new=T)
barplot(-log10(nullPs[Ord]),horiz=T, border="white", col = "tan", names.arg=NULL,width=.5, las = 2, cex.axis = .5, cex.main=.75,main="Permuted p-val")
par(fig=c(.75,1, 0,1), new=T)
barplot(ch2Scores[Ord,2],horiz=T, border="white", col = "tan", names.arg=NULL,width=.5, las = 2, cex.axis = .5, cex.main=.75,main="BAC")
dev.off()



