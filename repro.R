
#' Reproducibility/Repeatability
#'
#' Performs PCA and calculates the reproducibility of the analytical replicates based on the scores of the first two principal components
#'
#' @param data Data in the format of observations as rows and variables as columns
#' @param dsc This is the data starting column number, i.e. the column where the data starts after all the meta data
#' @param repgroup This is the column number with the groups of replicates. This must be in numeric form i.e. 1,2,3,4 etc...
#'
#' @return
#' @export
#'
#' @examples
#' repro(dataname, 6, 1)
#'
#' #Produces:
#'  #     PC1       PC2
#' #0.9224886 0.9176070
#' #[1] "Mean:  0.920047806692795"
repro<-function(data,dsc,repgroup){
data<-data[order(data[,repgroup], decreasing=FALSE),]
data[1:10,1:10]
mypca1=prcomp(data[,dsc:ncol(data)])
plot(mypca1$x[,1], mypca1$x[,2], col=data[,repgroup], xlab="PC1", ylab="PC2", pch=16)
pcascores=mypca1$x[,1:2]
#text(mypca1$x[,1], mypca1$x[,2], row.names(pcascores), cex=0.6, pos=4, col="black")


# Within Replicates
#Total Number of observations
N=nrow(data)


#Total number of Groups
k=max(as.numeric(data[,repgroup]))



#Scores by Group
scoresgrouplist=c()
scoresgroupM=c()
for (j in 1:k){
  x<-assign(paste("scoresgroup",j, sep=""), pcascores[as.numeric(data[,repgroup])==j,])
  scoresgrouplist=c(scoresgrouplist,list(x))
  scoresgroupM=rbind(scoresgroupM,scoresgrouplist[[j]])
}

groupmeanslist=c()
groupmeansM=c()
for (j in 1:k){
  mu=colMeans(scoresgrouplist[[j]])
  groupmeanslist=c(groupmeanslist, list(mu))
  groupmeansM=rbind(groupmeansM,groupmeanslist[[j]])
}




SSwlist=c()
SSwM=c()
for (j in 1:k){
  for (i in 1:nrow(scoresgrouplist[[j]])){
    x<-assign(paste("SSwgroup",j,sep=""), (scoresgrouplist[[j]][i,]-groupmeanslist[[j]])^2)
    SSwlist=c(SSwlist,list(x))
  }
}


SSwM=c()
for (i in 1:N){
  x=SSwlist[[i]]
  SSwM=rbind(SSwM, x)
}


numingroup=c()
for (j in 1:k){
  x<-nrow(scoresgrouplist[[j]])
  numingroup=c(numingroup,x)
}


groupcol=as.matrix(rep(1:k, numingroup))
colnames(groupcol)<-c("groupcol")
SSwM1=cbind(SSwM,groupcol)


SSwlist2=c()
for (j in 1:k){
  x<-assign(paste("SSwgroup2",j, sep=""), (subset(SSwM1,groupcol==j)))
  SSwlist2=c(SSwlist2,list(x))
}

SSw_first=c()
SSw_firstM=c()
for (j in 1:k){
  x=colSums(SSwlist2[[j]])
  SSw_first=c(SSw_first, list(x))
  SSw_firstM=rbind(SSw_firstM, SSw_first[[j]])
}

SSWrep=colSums(SSw_firstM)
SSWrep_PC1=SSWrep[1]
SSWrep_PC2=SSWrep[2]
SSWrep=SSWrep[-3]

Variation_SW2rep=c()
for (i in 1:length(SSWrep)){
  x=SSWrep[i]/N-k
  Variation_SW2rep=c(Variation_SW2rep,x)
}

########

# Between Groups
N=nrow(data)


#Total number of batches
k=max(as.numeric(data[,repgroup]))


#Scores by batch
scoresgrouplist=c()
scoresgroupM=c()
for (i in 1:k){
  x<-assign(paste("scoresgroup",i, sep=""), pcascores[as.numeric(data[,repgroup])==i,])
  scoresgrouplist=c(scoresgrouplist,list(x))
  scoresgroupM=rbind(scoresgroupM,scoresgrouplist[[i]])
}

groupmeanslist=c()
groupmeansM=c()
for (i in 1:k){
  mu=colMeans(scoresgrouplist[[i]])
  groupmeanslist=c(groupmeanslist, list(mu))
  groupmeansM=rbind(groupmeansM,groupmeanslist[[i]])
}

grandmean=colMeans(scoresgroupM)


SSBM=c()
SSBlist=c()
for (j in 1:k){
  x=(nrow(scoresgrouplist[[j]]))*((groupmeanslist[[j]]-grandmean)^2)
  SSBlist=c(SSBlist,list(x))
  SSBM=rbind(SSBM,SSBlist[[j]])
}

SSBgroup=colSums(SSBM)
SSBgroup_PC1=SSBgroup[1]
SSBgroup_PC2=SSBgroup[2]

Variation_SB2group=c()
for (i in 1:length(SSBgroup)){
  x=SSBgroup[i]/k-1
  Variation_SB2group=c(Variation_SB2group,x)
}


print(Variation_SB2group/(Variation_SB2group+Variation_SW2rep))


paste("Mean: ", mean(Variation_SB2group/(Variation_SB2group+Variation_SW2rep)))

}
