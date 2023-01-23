#######################################################################
# Purpose: functions needed for clustering categorical functional data
# Author: Xiaoxia Champon
# Date: 1/21/2023
# Update:
#####################################################################

##input categorical functional data n*t and output clustering results, latent curves, probability curves
catcluster=function(catfd,st,et,splines1D,M,knnum,pct,minPts,max.nc,min.nc){
  datapoints=dim(catfd)[2]
  tolcat=table(catfd)
  catorder=order(tolcat,decreasing = TRUE)
  numcat=length(catorder)
  refcat=catorder[numcat]
  refmat=catfd
  refmat[refmat!=refcat]=0
  refmat[refmat==refcat]=1
  nsub=dim(catfd)[1]
  ntime=dim(catfd)[2]
  subdata=array(data=NA,c(nsub,ntime,numcat))
  for (i in 1:numcat){
    datacopy=catfd
    datacopy[datacopy!=i]=0
    datacopy[datacopy==i]=1
    subdata[,,i]=datacopy
  }
  t=seq(st,et,length=datapoints)
  Zihat=array(data=NA,c(nsub,ntime,numcat))
  for (i in 1:numcat){
    datacopy=subdata
    Zihat[,,i]=Z_ihat(datacopy[,,i],t)
  }

  Zihatstar=array(data=NA,c(nsub,ntime,numcat-1))
  for (i in 1:(numcat-1)){
    datacopy=Zihat
    Zihatstar[,,i]=Zihat[,,i]+log(1+exp(Zihat[,,numcat]))-log(1+exp(Zihat[,,i]))-Zihat[,,numcat]
  }

  phatmat=phatf(Zihatstar)

  vecapply=matrix(1:(dim(Zihatstar)[3]),ncol=1)
  mfdataused=apply(vecapply,1,function(x) {mfundata(Zihatstar[,,x],t)})
  mvdata=multiFunData(mfdataused)

  uniexpan=list()
  # MFPCA based on univariate FPCA Z_ihat
  for (i in 1:(numcat-1)){
    uniexpan[[i]]=list(type = "splines1D", k = splines1D)
  }

  # MFPCA based on univariate FPCA Z_ihat
  uFPCA <- MFPCA(mvdata, M = M, uniExpansions = uniexpan)
  scores_z=uFPCA$scores

  dist=kNNdist(scores_z, k = knnum)
  distdata=data.frame(sort(dist))
  distdata$index=1:dim(distdata)[1]
  ninty5p=quantile(dist, probs = pct)
  dp <- ggplot(distdata,aes(index,sort.dist.)) + geom_line()+ggtitle(paste0(knnum,"-NN Distance Plot ",'\n',"(",dim(distdata)[1]," Subjects",")")) +
    xlab("Points sorted by Distance") + ylab("Distance")+ theme(plot.title = element_text(hjust = 0.5))+geom_hline(yintercept=ninty5p, color = "red")+
    geom_text(data=data.frame(round(ninty5p,2)),
              aes(x=dim(distdata)[1]/2,y=1.2*ninty5p,label=paste0("Distance at ",gsub("%$","",row.names(data.frame(round(ninty5p,2)))),"th percentile= ",round(ninty5p,2))))

  #7, 0.98 2 clusters  6, 88
  res <- dbscan(scores_z, eps =pct , minPts = minPts )

  clustertable=table(res$cluster)
  #tclustertable    #z score


  #########Kmeans
  reskmeans=NbClust(data = scores_z, diss = NULL, distance = "euclidean",
                    min.nc = min.nc, max.nc = max.nc, method = "kmeans")
  clustertablek=table(reskmeans$Best.partition)

  return(list("scores"=scores_z,"distfig"=dp,"dbcluster"=res$cluster,"dbtable"=clustertable,
              "kcluster"=reskmeans$Best.partition,"kmeantable"=clustertablek,
              "latentcurve"=Zihatstar,"meanfn"=uFPCA$meanFunction,"eigenvalue"=uFPCA$values,
              "eigenfn"=uFPCA$functions,"probcurve"=phatmat))
}
