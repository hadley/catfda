# catfda
categorical functional data analysis
Package: catfda
Type: Package
Title: categorical functional data analysis using latent curves
Version: 0.1.0
Author: Xiaoxia Champon
Maintainer: The package maintainer <xiachampon@gmail.com>
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
   3) sample code is in the folder: R
