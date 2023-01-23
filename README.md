# catfda
categorical functional data analysis
Package: catfda
Type: Package
Title: categorical functional data analysis using latent curves
Version: 0.1.0
Author: Xiaoxia Champon
Maintainer: The package maintainer <xzhao17@ncsu.edu>
Description: This package analyzes the categorical functional data. The catcluster function takes   	    
    categorical functional data as input and returns and latent curves, probability curves and 
    multivariate functional principal component scores, which are used to cluster the categorical 
    functional data. See the catfda document for more details.

				
License: What license is it under?
Encoding: UTF-8
LazyData: true
Imports:
    fda,
    refund,
    mgcv,
    funData,
    MFPCA,
    dbscan,
    fossil,
    NbClust,
    ggplot2
    
    
    To use the functions in the package:
  1) library(devtools)
     install_github("XiaoxiaChampon/catfda")
  2) library("catfda")
     #to cluster the categorical functional data
     #catfd is the n*t matrix, where n is the number of subjects, t is observatioal time points
     #st: first observational time, et: last observational time
     #splines1D is the number of basis expansions, it is usually set to 25
     #M: the number of principal component chosen, usually it's set to 2
     #knnum,pct,minPts: parameters from DBSCAN, usuallly set to 3, 0.9, 3
     #max.nc,min.nc: parameters from Kmeans clustering, usually set to maximum number of cluster 5, minimum 2.
     catcluster(catfd,st,et,splines1D,M,knnum,pct,minPts,max.nc,min.nc)
   3) sample code 
    #######################################################################
# Purpose: simulate data and clustering categorical functional data sample code
# Author: Xiaoxia Champon
# Date: 1/21/2023
# Update:
#####################################################################




#######################################################################################################
#test the code
#catsim=function(seed=123,datapoints,n,sparse1,sparse2,scorevar1,scorevar2,ps,st,et,k,propc1,propc2,propnoise)

dpoints=1000. #1000 time points for 1000 observations per subject
nsub=100. #the number of the subject
testcatdata=catsim(seed=123,dpoints,nsub,4,1,1,1,1,0.01,0.99,3,0.7,0.3,0.04)

#catcluster=function(catfd,st,et,splines1D,M,knnum,pct,minPts,max.nc,min.nc)
st=Sys.time()
testcluster=catcluster(testcatdata,0,1,25,3,3,0.99,3,5,2)
et=Sys.time()
telapse=et-st
telapse
testcluster$distfig

#randfn=function(dbcluster,kcluster,nsub)
randfn(testcluster$dbcluster,testcluster$kcluster,nsub)

#plot_cluster=function(scores_z,dbcluster,kcluster,st,et,datapoints)
plot_cluster(testcluster$scores,testcluster$dbcluster,testcluster$kcluster,0,1,dpoints)

#plot_latent=function(zlatent,phat,st,et,cluster,labelnum,argval)
plot_latent(testcluster$latentcurve,testcluster$probcurve,0,1,
            testcluster$dbcluster,0,seq(0,1,length=dpoints))  #testcluster$meanfn

plot_latent(testcluster$latentcurve,testcluster$probcurve,0,1,
            testcluster$dbcluster,1,seq(0,1,length=dpoints))
plot_latent(testcluster$latentcurve,testcluster$probcurve,0,1,
            testcluster$dbcluster,2,seq(0,1,length=dpoints))

plot_latent(testcluster$latentcurve,testcluster$probcurve,0,1,
            testcluster$kcluster,3,seq(0,1,length=dpoints))

plot_latent(testcluster$latentcurve,testcluster$probcurve,0,1,
            testcluster$kcluster,1,seq(0,1,length=dpoints))

plot_latent(testcluster$latentcurve,testcluster$probcurve,0,1,
            testcluster$kcluster,2,seq(0,1,length=dpoints))
#####################################################################################################
