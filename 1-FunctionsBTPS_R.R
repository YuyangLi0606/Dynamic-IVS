#####Functions for Bivariate Tensor Product B-splines

##Required Libraries

require(mgcv)
require(msm)
require(stats)
# require(rgl)
require(MASS)
require(reshape2)
require(plyr)
require(ggplot2)
########################Function Names and purpose
#fitBTPS-Fits a penalized bivariate tensor product B-spline using According to Section 3
#         Then estimates the first and second partial derivative with respect to the x values
#
#yieldDataSim-creates simulated data in accordance to section 4
#
#compute_bandwidth-Calculate bandwidth for kernel density estimation
#                   Helper function for VarianceEstimator
#
#predictDensity-calculate a bivariate kernel density, helper function for VarianceEstimator, uses compute_bandwidth
#
#VarianceEstimator-Calculates the variance for a fitBTPS object according to equation 22
#
#kernelDensity
#
#fitKernel

#######################################



##inputs: y- the dependent variable for the penalized BTPS
#x-the independent variable that will have it's derivative taken
#z-the second independent variable
#knots-the number of knots used for each axis
#degree-the degree of the spline used
#penalty the order of the penalty used
#tol1st- the tolerance when taking the first derivative
#tol2nd-the tolerance when taking the second derivative
#degree.sd-degree of the spline used when estimating sigma
#penalty.sd=c(1,1)-the penalty used when estimating sigma
#newXLim-the range of x values used to predict on
#newZLim-the range of z values used to predict on
#nonDeg-applied the adjustment to keep  positive on second derivatives along the x-axis
##(note that nonDeg is for a very specific case working only with quadratic terms)


##Returns a list with the following:
##outputs-origY a vector of the original y values used,
#origX-a vector of the original x values used
#origZ-a vector of the original z values used
#fitOrig-fitted values for the original data points
#newX-new x values to predict
#newZ-new z values to predict
#Fit-the predicted y of the new x and z values
#firstDerv- the estimated first derivatives at the new x and z values
#secondDerv- the estimated second derivatives at the new x and z values
#int.knots-the number of interior knots used for each axis
#smoothedSigma-the estimated sigma value, used for estimating the variance
#model-the output from the gam function



# fitBTPS_zc<-function(y,x,z,knots=c(0,0),degree=c(5,3),penalty=c(1,1),tol1st=.0000001,tol2nd=.1,degree.sd=c(3,3),penalty.sd=c(1,1),newXLim=c(min(x),max(x)),newZLim=c(min(z),max(z)),
#                      X_norm_stats = X_norm_stats, new_X_Z_grid = c(100, 50),  nonDeg=F, penalty_bool = TRUE){
#   
#   
#   ##Get a value for the number of B-spline Basis needed for the set up based on degree and knots
# #  myK<-c(knots+degree+2) # ??? should be knots+degree+1 
#   myK<-c(knots+degree+1) 
#   ##create the degree and penalty list for the function
#   m<-list(c(degree[1],penalty[1]),c(degree[2],penalty[2]))
#   
#   ##fit the BTPS using the gam function in r
#   if(penalty_bool == TRUE){
#     theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F) # ??? should add  fx = TRUE in te() to turn off the penalty ?
#   }else{
#     theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK,fx = TRUE),drop.intercept=F) # ??? should add  fx = TRUE in te() to turn off the penalty ?
#   }
#   
#   
#   
#   
#   ###create new data and fit the new data
#   # xNew<-seq(newXLim[1],newXLim[2],length.out = 50)##values of length.out might not work on all computers due to memory limitations
#   # zNew<-seq(newZLim[1],newZLim[2],length.out = 50)
#   xNew_unnorm<-seq(newXLim[1],newXLim[2],length.out = new_X_Z_grid[1])##values of length.out might not work on all computers due to memory limitations
#   zNew<-seq(newZLim[1],newZLim[2],length.out = new_X_Z_grid[2])
#   
#      
#   
#   hold.dataframe<-expand.grid(xNew_unnorm,zNew)
#   
#   xNew_unnorm<-hold.dataframe$Var1
#   # normalization 
#   xNew <- (xNew_unnorm - X_norm_stats[1]) / X_norm_stats[2]
#   zNew<-hold.dataframe$Var2
#   Fit<-as.numeric(predict(theFit,data.frame(x=xNew,z=zNew)))
#   
#   
#   yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
#   
#   ##Estimate the First derivative fit
#   Derv1Fit<--(predict(theFit,data.frame(x=xNew,z=zNew))-predict(theFit,data.frame(x=xNew+tol1st,z=zNew)))/tol1st
#   
#   ##Estimate the Second derivative fit
#   Derv2Fit<-(predict(theFit,data.frame(x=xNew-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew,z=zNew))+predict(theFit,data.frame(x=xNew+tol2nd,z=zNew)))/tol2nd^2
#   # browser()
#   Derv2FitNonNeg<-NA
# 
#   ##Ensure the fit is positive
#   if(nonDeg){
#     number<-sum(((predict(theFit,data.frame(x=xNew-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew,z=zNew))+predict(theFit,data.frame(x=xNew+tol2nd,z=zNew)))/tol2nd^2)<0)
#     if(number>0){
#       theFitInitial<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F )
#       coefMatrix<-matrix(coefficients(theFit),ncol  = myK[1])
#       
#       Dmatrix<-matrix(0,nrow  = myK[1],ncol=myK[2])
#       for(i in 1:(nrow(coefMatrix)-2)){
#         Dmatrix[i+1,c((i):(i+2))]<-c(1,-2,1)
#         
#       }
#       Dmatrix<-t(Dmatrix)
#       
# 
#       value<-diag(as.numeric(colSums( Dmatrix* matrix(coefficients(theFit),nrow  = myK[1]))<0))
#       theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F,paraPen =list(x=list(sp=value)) ) # penalize on x only ? sp: a vector of smoothing parameter values
#       
#       
#     }
#     
#     
#     yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
#     
#     Fit<-as.numeric(predict(theFit,data.frame(x=xNew,z=zNew)))
#     
#     ##Estimate the First derivative fit
#     Derv1Fit<--(predict(theFit,data.frame(x=xNew,z=zNew))-predict(theFit,data.frame(x=xNew+tol1st,z=zNew)))/tol1st
#     
#     ##Estimate the Second derivative fit
#     Derv2FitNonNeg<-(predict(theFit,data.frame(x=xNew-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew,z=zNew))+predict(theFit,data.frame(x=xNew+tol2nd,z=zNew)))/tol2nd^2
#     #
#     
#   }
#   
#   ###Now for the error estimation
#   
#   Residuals.V<-as.numeric(theFit$residuals)^2
#   m.sd<-list(penalty.sd,degree.sd)
#   
#   ##Fit the residuals as in equation 20
#   fittedRes<-gam(Residuals.V~te(x,z,bs=c("ps","ps"),m=m.sd,k=myK),drop.intercept=F)
#   residualsEst<-abs(as.numeric(predict(fittedRes,data.frame(x=xNew,z=zNew))))
#   lambdaValues<-theFit$sp
#   
#   
#   
#   
#   ##data returned
#   return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,newX_unnorm = xNew_unnorm, newX=xNew,newZ=zNew,Fit=Fit,firstDerv=as.numeric(Derv1Fit),secondDerv=as.numeric(Derv2Fit),int.knots=knots,smoothedSigma=residualsEst,model=theFit,NonNeg_Derv2Fit=Derv2FitNonNeg[Derv2FitNonNeg<0]<-0))
#   
# }


# concise version of fitBTPS_zc(): provide the option for if_pred_new option, ignore the derivative related part and residual related part, input x is normalized moneyness
fitBTPS_zc<-function(y,x,z,knots=c(0,0),degree=c(5,3),penalty=c(1,1),tol1st=.0000001,tol2nd=.1,degree.sd=c(3,3),penalty.sd=c(1,1),newXLim=c(min(x),max(x)),newZLim=c(min(z),max(z)),
                    X_norm_stats = X_norm_stats, new_X_Z_grid = c(100, 50),  nonDeg=F, penalty_bool = TRUE, if_pred_new = TRUE){

  ##Get a value for the number of B-spline Basis needed for the set up based on degree and knots
   myK<-c(knots+degree+2) # should be knots+degree+1: according to the ref book: degree here seems to be the (degree - 1) in our understanding, so m =2 represents a cubic spline
  # myK<-c(knots+degree+1) 
  ##create the degree and penalty list for the function
  m<-list(c(degree[1],penalty[1]),c(degree[2],penalty[2]))
  
  ##fit the BTPS using the gam function in r
  if(penalty_bool == TRUE){
    # browser()
    theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ? 
  }else{
    theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK,fx = TRUE, np = FALSE),drop.intercept=F) # should add  fx = TRUE in te() to turn off the penalty ?
  }
  
  
  
  if(if_pred_new == TRUE){
    ###create new data and fit the new data
    # xNew<-seq(newXLim[1],newXLim[2],length.out = 50)##values of length.out might not work on all computers due to memory limitations
    # zNew<-seq(newZLim[1],newZLim[2],length.out = 50)
    xNew_unnorm<-seq(newXLim[1],newXLim[2],length.out = new_X_Z_grid[1])##values of length.out might not work on all computers due to memory limitations
    zNew<-seq(newZLim[1],newZLim[2],length.out = new_X_Z_grid[2])
    
    
    
    hold.dataframe<-expand.grid(xNew_unnorm,zNew)
    
    xNew_unnorm<-hold.dataframe$Var1
    # normalization 
    xNew <- (xNew_unnorm - X_norm_stats[1]) / X_norm_stats[2]
    zNew<-hold.dataframe$Var2
    Fit<-as.numeric(predict(theFit,data.frame(x=xNew,z=zNew)))
    
  }
  
  
  yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
  
  # ##Estimate the First derivative fit
  # Derv1Fit<--(predict(theFit,data.frame(x=xNew,z=zNew))-predict(theFit,data.frame(x=xNew+tol1st,z=zNew)))/tol1st
  # 
  # ##Estimate the Second derivative fit
  # Derv2Fit<-(predict(theFit,data.frame(x=xNew-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew,z=zNew))+predict(theFit,data.frame(x=xNew+tol2nd,z=zNew)))/tol2nd^2
  # # browser()
  # Derv2FitNonNeg<-NA
  
  ##Ensure the fit is positive
  if(nonDeg){
    number<-sum(((predict(theFit,data.frame(x=xNew-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew,z=zNew))+predict(theFit,data.frame(x=xNew+tol2nd,z=zNew)))/tol2nd^2)<0)
    if(number>0){
      theFitInitial<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F )
      coefMatrix<-matrix(coefficients(theFit),ncol  = myK[1])
      
      Dmatrix<-matrix(0,nrow  = myK[1],ncol=myK[2])
      for(i in 1:(nrow(coefMatrix)-2)){
        Dmatrix[i+1,c((i):(i+2))]<-c(1,-2,1)
        
      }
      Dmatrix<-t(Dmatrix)
      
      
      value<-diag(as.numeric(colSums( Dmatrix* matrix(coefficients(theFit),nrow  = myK[1]))<0))
      theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F,paraPen =list(x=list(sp=value)) ) # penalize on x only ? sp: a vector of smoothing parameter values
      
      
    }
    
    
    yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
    
    Fit<-as.numeric(predict(theFit,data.frame(x=xNew,z=zNew)))
    
    ##Estimate the First derivative fit
    Derv1Fit<--(predict(theFit,data.frame(x=xNew,z=zNew))-predict(theFit,data.frame(x=xNew+tol1st,z=zNew)))/tol1st
    
    ##Estimate the Second derivative fit
    Derv2FitNonNeg<-(predict(theFit,data.frame(x=xNew-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew,z=zNew))+predict(theFit,data.frame(x=xNew+tol2nd,z=zNew)))/tol2nd^2
    #
    
  }
  
  ###Now for the error estimation
  
  # Residuals.V<-as.numeric(theFit$residuals)^2
  # m.sd<-list(penalty.sd,degree.sd)
  # 
  # ##Fit the residuals as in equation 20
  # fittedRes<-gam(Residuals.V~te(x,z,bs=c("ps","ps"),m=m.sd,k=myK),drop.intercept=F)
  # if(if_pred_new == TRUE){
  #   residualsEst<-abs(as.numeric(predict(fittedRes,data.frame(x=xNew,z=zNew))))
  # }
  
  
  lambdaValues<-theFit$sp
  
  
  
  
  ##data returned (AT: return origx = x where z is the input(normalized moneyness))
  if(if_pred_new == TRUE){
  return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,newX_unnorm = xNew_unnorm, 
              newX=xNew,newZ=zNew,Fit=Fit,int.knots=knots,model=theFit))
  }else{
    return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,
                int.knots=knots,model=theFit))
  }
  
}


# # simply change bs=c("ps","ps") to bs=c("bs","bs") in fitBTPS_zc()
# # But the model matrix actually NOT FULL RANK (rank deficient worse than 'ps')
# fitBTPS_zc_unnorm_bs<-function(y,x,z,knots=c(0,0),degree=c(5,3),penalty=c(1,1),tol1st=.0000001,tol2nd=.1,degree.sd=c(3,3),penalty.sd=c(1,1),newXLim=c(min(x),max(x)),newZLim=c(min(z),max(z)),
#                       new_X_Z_grid = c(100, 50),  nonDeg=F, penalty_bool = TRUE, if_pred_new = TRUE){
#   
#   ##Get a value for the number of B-spline Basis needed for the set up based on degree and knots
#   myK<-c(knots+degree+2) # should be knots+degree+1: according to the ref book: degree here seems to be the (degree - 1) in our understanding, so m =2 represents a cubic spline
#   # myK<-c(knots+degree+1) 
#   ##create the degree and penalty list for the function
#   m<-list(c(degree[1],penalty[1]),c(degree[2],penalty[2]))
#   
#   ##fit the BTPS using the gam function in r
#   if(penalty_bool == TRUE){
#     # browser()
#     # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
#     theFit<-gam(y~ti(x,z,bs=c("bs","bs"),m=m,k=myK, np = FALSE, ),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
#     # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
#     # theFit<-gam(y~te(x,z,bs=c("bs","bs"),m=m,k=myK),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
#     
#       }else{
#     theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK, np = FALSE, fx = TRUE),drop.intercept=F) # should add  fx = TRUE in te() to turn off the penalty ?
#   }
#   
#   
#   
#   if(if_pred_new == TRUE){
#     ###create new data and fit the new data
#     # xNew<-seq(newXLim[1],newXLim[2],length.out = 50)##values of length.out might not work on all computers due to memory limitations
#     # zNew<-seq(newZLim[1],newZLim[2],length.out = 50)
#     xNew_unnorm<-seq(newXLim[1],newXLim[2],length.out = new_X_Z_grid[1])##values of length.out might not work on all computers due to memory limitations
#     # zNew<-seq(newZLim[1],newZLim[2],length.out = new_X_Z_grid[2])
#     zNew<-seq(newZLim[1],newZLim[2],by = 1/new_X_Z_grid[2]) # this specification is used so that we can plot 2D plot on exact ptmday
#     
#     
#     
#     hold.dataframe<-expand.grid(xNew_unnorm,zNew)
#     
#     xNew_unnorm<-hold.dataframe$Var1
#     # normalization 
#     # xNew <- (xNew_unnorm - X_norm_stats[1]) / X_norm_stats[2]
#     zNew<-hold.dataframe$Var2
#     # Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew)))
#     Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew * 365))) # modified on 05.10: to plot the 2D plot surface
#     
#   }
#   
#   
#   yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
#   
#   # ##Estimate the First derivative fit
#   # Derv1Fit<--(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))-predict(theFit,data.frame(x=xNew_unnorm+tol1st,z=zNew)))/tol1st
#   # 
#   # ##Estimate the Second derivative fit
#   # Derv2Fit<-(predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2
#   # # browser()
#   # Derv2FitNonNeg<-NA
#   
#   ##Ensure the fit is positive
#   if(nonDeg){
#     number<-sum(((predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2)<0)
#     if(number>0){
#       theFitInitial<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F )
#       coefMatrix<-matrix(coefficients(theFit),ncol  = myK[1])
#       
#       Dmatrix<-matrix(0,nrow  = myK[1],ncol=myK[2])
#       for(i in 1:(nrow(coefMatrix)-2)){
#         Dmatrix[i+1,c((i):(i+2))]<-c(1,-2,1)
#         
#       }
#       Dmatrix<-t(Dmatrix)
#       
#       
#       value<-diag(as.numeric(colSums( Dmatrix* matrix(coefficients(theFit),nrow  = myK[1]))<0))
#       theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F,paraPen =list(x=list(sp=value)) ) # penalize on x only ? sp: a vector of smoothing parameter values
#       
#       
#     }
#     
#     
#     yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
#     
#     Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew)))
#     
#     ##Estimate the First derivative fit
#     Derv1Fit<--(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))-predict(theFit,data.frame(x=xNew_unnorm+tol1st,z=zNew)))/tol1st
#     
#     ##Estimate the Second derivative fit
#     Derv2FitNonNeg<-(predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2
#     #
#     
#   }
#   
#   ###Now for the error estimation
#   
#   # Residuals.V<-as.numeric(theFit$residuals)^2
#   # m.sd<-list(penalty.sd,degree.sd)
#   # 
#   # ##Fit the residuals as in equation 20
#   # fittedRes<-gam(Residuals.V~te(x,z,bs=c("ps","ps"),m=m.sd,k=myK),drop.intercept=F)
#   # if(if_pred_new == TRUE){
#   #   residualsEst<-abs(as.numeric(predict(fittedRes,data.frame(x=xNew_unnorm,z=zNew))))
#   # }
#   
#   
#   lambdaValues<-theFit$sp
#   
#   
#   
#   
#   ##data returned (AT: return origx = x where z is the input(normalized moneyness))
#   if(if_pred_new == TRUE){ # AT: newX and newX_unnorm both set to xNew_unnorm
#     return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,newX = xNew_unnorm, newX_unnorm = xNew_unnorm, newZ=zNew,Fit=Fit,int.knots=knots,model=theFit))
#   }else{
#     return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,
#                 int.knots=knots,model=theFit))
#   }
#   
# }

fitBTPS_zc_unnorm_bs_own<-function(y,x,z,knots=c(0,0),degree=c(5,3),penalty=c(1,1),tol1st=.0000001,tol2nd=.1,degree.sd=c(3,3),penalty.sd=c(1,1),newXLim=c(min(x),max(x)),newZLim=c(min(z),max(z)),
                               new_X_Z_grid = c(100, 50),  nonDeg=F, penalty_bool = TRUE, if_pred_new = TRUE){
  
  ##Get a value for the number of B-spline Basis needed for the set up based on degree and knots
  myK<-c(knots+degree+2) # should be knots+degree+1: according to the ref book: degree here seems to be the (degree - 1) in our understanding, so m =2 represents a cubic spline
  # myK<-c(knots+degree+1) 
  ##create the degree and penalty list for the function
  m<-list(c(degree[1],penalty[1]),c(degree[2],penalty[2]))
  
  ##fit the BTPS using the gam function in r
  if(penalty_bool == TRUE){
    # browser()
    # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    theFit<-gam(y~ti(x,z,bs=c("bs","bs"),m=m,k=myK, np = FALSE, ),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    # theFit<-gam(y~te(x,z,bs=c("bs","bs"),m=m,k=myK),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    
  }else{
    theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK, np = FALSE, fx = TRUE),drop.intercept=F) # should add  fx = TRUE in te() to turn off the penalty ?
  }
  
  
  
  if(if_pred_new == TRUE){
    ###create new data and fit the new data
    # xNew<-seq(newXLim[1],newXLim[2],length.out = 50)##values of length.out might not work on all computers due to memory limitations
    # zNew<-seq(newZLim[1],newZLim[2],length.out = 50)
    xNew_unnorm<-seq(newXLim[1],newXLim[2],length.out = new_X_Z_grid[1])##values of length.out might not work on all computers due to memory limitations
    # zNew<-seq(newZLim[1],newZLim[2],length.out = new_X_Z_grid[2])
    zNew<-seq(newZLim[1],newZLim[2],by = 1/new_X_Z_grid[2]) # this specification is used so that we can plot 2D plot on exact ptmday
    
    
    
    hold.dataframe<-expand.grid(xNew_unnorm,zNew)
    
    xNew_unnorm<-hold.dataframe$Var1
    # normalization 
    # xNew <- (xNew_unnorm - X_norm_stats[1]) / X_norm_stats[2]
    zNew<-hold.dataframe$Var2
    # Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew)))
    Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew * 365))) # modified on 05.10: to plot the 2D plot surface
    
  }
  
  
  yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
  
  # ##Estimate the First derivative fit
  # Derv1Fit<--(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))-predict(theFit,data.frame(x=xNew_unnorm+tol1st,z=zNew)))/tol1st
  # 
  # ##Estimate the Second derivative fit
  # Derv2Fit<-(predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2
  # # browser()
  # Derv2FitNonNeg<-NA
  
  ##Ensure the fit is positive
  if(nonDeg){
    number<-sum(((predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2)<0)
    if(number>0){
      theFitInitial<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F )
      coefMatrix<-matrix(coefficients(theFit),ncol  = myK[1])
      
      Dmatrix<-matrix(0,nrow  = myK[1],ncol=myK[2])
      for(i in 1:(nrow(coefMatrix)-2)){
        Dmatrix[i+1,c((i):(i+2))]<-c(1,-2,1)
        
      }
      Dmatrix<-t(Dmatrix)
      
      
      value<-diag(as.numeric(colSums( Dmatrix* matrix(coefficients(theFit),nrow  = myK[1]))<0))
      theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F,paraPen =list(x=list(sp=value)) ) # penalize on x only ? sp: a vector of smoothing parameter values
      
      
    }
    
    
    yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
    
    Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew)))
    
    ##Estimate the First derivative fit
    Derv1Fit<--(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))-predict(theFit,data.frame(x=xNew_unnorm+tol1st,z=zNew)))/tol1st
    
    ##Estimate the Second derivative fit
    Derv2FitNonNeg<-(predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2
    #
    
  }
  
  ###Now for the error estimation
  
  # Residuals.V<-as.numeric(theFit$residuals)^2
  # m.sd<-list(penalty.sd,degree.sd)
  # 
  # ##Fit the residuals as in equation 20
  # fittedRes<-gam(Residuals.V~te(x,z,bs=c("ps","ps"),m=m.sd,k=myK),drop.intercept=F)
  # if(if_pred_new == TRUE){
  #   residualsEst<-abs(as.numeric(predict(fittedRes,data.frame(x=xNew_unnorm,z=zNew))))
  # }
  
  
  lambdaValues<-theFit$sp
  
  
  
  
  ##data returned (AT: return origx = x where z is the input(normalized moneyness))
  if(if_pred_new == TRUE){ # AT: newX and newX_unnorm both set to xNew_unnorm
    return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,newX = xNew_unnorm, newX_unnorm = xNew_unnorm, newZ=zNew,Fit=Fit,int.knots=knots,model=theFit))
  }else{
    return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,
                int.knots=knots,model=theFit))
  }
  
}


# version of fitBTPS_zc() with unnormalized moneyness, input x is unnormalized moneyness
fitBTPS_zc_unnorm<-function(y,x,z,knots=c(0,0),degree=c(5,3),penalty=c(1,1),tol1st=.0000001,tol2nd=.1,degree.sd=c(3,3),penalty.sd=c(1,1),newXLim=c(min(x),max(x)),newZLim=c(min(z),max(z)),
                            new_X_Z_grid = c(100, 50),  nonDeg=F, penalty_bool = TRUE, if_pred_new = TRUE){
  
  ##Get a value for the number of B-spline Basis needed for the set up based on degree and knots
  myK<-c(knots+degree+2) # should be knots+degree+1: according to the ref book: degree here seems to be the (degree - 1) in our understanding, so m =2 represents a cubic spline
  # myK<-c(knots+degree+1) 
  ##create the degree and penalty list for the function
  m<-list(c(degree[1],penalty[1]),c(degree[2],penalty[2]))
  
  ##fit the BTPS using the gam function in r
  if(penalty_bool == TRUE){
    # browser()
    theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    # theFit<-gam(y~te(x,z,bs=c("bs","bs"),m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    
    # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    
  }else{
    theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK, np = FALSE, fx = TRUE),drop.intercept=F) # should add  fx = TRUE in te() to turn off the penalty ?
  }
  
  
  
  if(if_pred_new == TRUE){
    ###create new data and fit the new data
    # xNew<-seq(newXLim[1],newXLim[2],length.out = 50)##values of length.out might not work on all computers due to memory limitations
    # zNew<-seq(newZLim[1],newZLim[2],length.out = 50)
    xNew_unnorm<-seq(newXLim[1],newXLim[2],length.out = new_X_Z_grid[1])##values of length.out might not work on all computers due to memory limitations
    # zNew<-seq(newZLim[1],newZLim[2],length.out = new_X_Z_grid[2])
    zNew<-seq(newZLim[1],newZLim[2],by = 1/new_X_Z_grid[2]) # this specification is used so that we can plot 2D plot on exact ptmday
    
    
    
    hold.dataframe<-expand.grid(xNew_unnorm,zNew)
    
    xNew_unnorm<-hold.dataframe$Var1
    # normalization 
    # xNew <- (xNew_unnorm - X_norm_stats[1]) / X_norm_stats[2]
    zNew<-hold.dataframe$Var2
    print('Note: if back to simu, need to change code here')
    # Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))) # for simulation study 
    Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew * 365))) # modified on 05.10: to plot the 2D plot surface for real data
    
  }
  
  
  yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
  
  # ##Estimate the First derivative fit
  # Derv1Fit<--(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))-predict(theFit,data.frame(x=xNew_unnorm+tol1st,z=zNew)))/tol1st
  # 
  # ##Estimate the Second derivative fit
  # Derv2Fit<-(predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2
  # # browser()
  # Derv2FitNonNeg<-NA
  
  ##Ensure the fit is positive
  if(nonDeg){
    number<-sum(((predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2)<0)
    if(number>0){
      theFitInitial<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F )
      coefMatrix<-matrix(coefficients(theFit),ncol  = myK[1])
      
      Dmatrix<-matrix(0,nrow  = myK[1],ncol=myK[2])
      for(i in 1:(nrow(coefMatrix)-2)){
        Dmatrix[i+1,c((i):(i+2))]<-c(1,-2,1)
        
      }
      Dmatrix<-t(Dmatrix)
      
      
      value<-diag(as.numeric(colSums( Dmatrix* matrix(coefficients(theFit),nrow  = myK[1]))<0))
      theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F,paraPen =list(x=list(sp=value)) ) # penalize on x only ? sp: a vector of smoothing parameter values
      
      
    }
    
    
    yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
    
    Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew)))
    
    ##Estimate the First derivative fit
    Derv1Fit<--(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))-predict(theFit,data.frame(x=xNew_unnorm+tol1st,z=zNew)))/tol1st
    
    ##Estimate the Second derivative fit
    Derv2FitNonNeg<-(predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2
    #
    
  }
  
  ###Now for the error estimation
  
  # Residuals.V<-as.numeric(theFit$residuals)^2
  # m.sd<-list(penalty.sd,degree.sd)
  # 
  # ##Fit the residuals as in equation 20
  # fittedRes<-gam(Residuals.V~te(x,z,bs=c("ps","ps"),m=m.sd,k=myK),drop.intercept=F)
  # if(if_pred_new == TRUE){
  #   residualsEst<-abs(as.numeric(predict(fittedRes,data.frame(x=xNew_unnorm,z=zNew))))
  # }
  
  
  lambdaValues<-theFit$sp
  
  
  
  
  ##data returned (AT: return origx = x where z is the input(normalized moneyness))
  if(if_pred_new == TRUE){ # AT: newX and newX_unnorm both set to xNew_unnorm
    return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,newX = xNew_unnorm, newX_unnorm = xNew_unnorm, newZ=zNew,Fit=Fit,int.knots=knots,model=theFit))
  }else{
    return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,
                int.knots=knots,model=theFit))
  }
  
}

# version of fitBTPS_zc() with unnormalized moneyness, input x is unnormalized moneyness
# For simulation study IV_simu_step1_nongrid_data_unbinned.Rmd 
# the only diff from fitBTPS_zc_unnorm_simu:     Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))) # for simulation study
fitBTPS_zc_unnorm_simu<-function(y,x,z,knots=c(0,0),degree=c(5,3),penalty=c(1,1),tol1st=.0000001,tol2nd=.1,degree.sd=c(3,3),penalty.sd=c(1,1),newXLim=c(min(x),max(x)),newZLim=c(min(z),max(z)),
                            new_X_Z_grid = c(100, 50),  nonDeg=F, penalty_bool = TRUE, if_pred_new = TRUE){
  
  ##Get a value for the number of B-spline Basis needed for the set up based on degree and knots
  myK<-c(knots+degree+2) # should be knots+degree+1: according to the ref book: degree here seems to be the (degree - 1) in our understanding, so m =2 represents a cubic spline
  # myK<-c(knots+degree+1) 
  ##create the degree and penalty list for the function
  m<-list(c(degree[1],penalty[1]),c(degree[2],penalty[2]))
  
  ##fit the BTPS using the gam function in r
  if(penalty_bool == TRUE){
    # browser()
    theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    # theFit<-gam(y~te(x,z,bs=c("bs","bs"),m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    
    # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    
  }else{
    theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK, np = FALSE, fx = TRUE),drop.intercept=F) # should add  fx = TRUE in te() to turn off the penalty ?
  }
  
  
  
  if(if_pred_new == TRUE){
    ###create new data and fit the new data
    # xNew<-seq(newXLim[1],newXLim[2],length.out = 50)##values of length.out might not work on all computers due to memory limitations
    # zNew<-seq(newZLim[1],newZLim[2],length.out = 50)
    xNew_unnorm<-seq(newXLim[1],newXLim[2],length.out = new_X_Z_grid[1])##values of length.out might not work on all computers due to memory limitations
    # zNew<-seq(newZLim[1],newZLim[2],length.out = new_X_Z_grid[2])
    zNew<-seq(newZLim[1],newZLim[2],by = 1/new_X_Z_grid[2]) # this specification is used so that we can plot 2D plot on exact ptmday
    
    
    
    hold.dataframe<-expand.grid(xNew_unnorm,zNew)
    
    xNew_unnorm<-hold.dataframe$Var1
    # normalization 
    # xNew <- (xNew_unnorm - X_norm_stats[1]) / X_norm_stats[2]
    zNew<-hold.dataframe$Var2
    print('Note: if back to simu, need to change code here')
    Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))) # for simulation study
    # Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew * 365))) # modified on 05.10: to plot the 2D plot surface for real data
    
  }
  
  
  yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
  
  # ##Estimate the First derivative fit
  # Derv1Fit<--(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))-predict(theFit,data.frame(x=xNew_unnorm+tol1st,z=zNew)))/tol1st
  # 
  # ##Estimate the Second derivative fit
  # Derv2Fit<-(predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2
  # # browser()
  # Derv2FitNonNeg<-NA
  
  ##Ensure the fit is positive
  if(nonDeg){
    number<-sum(((predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2)<0)
    if(number>0){
      theFitInitial<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F )
      coefMatrix<-matrix(coefficients(theFit),ncol  = myK[1])
      
      Dmatrix<-matrix(0,nrow  = myK[1],ncol=myK[2])
      for(i in 1:(nrow(coefMatrix)-2)){
        Dmatrix[i+1,c((i):(i+2))]<-c(1,-2,1)
        
      }
      Dmatrix<-t(Dmatrix)
      
      
      value<-diag(as.numeric(colSums( Dmatrix* matrix(coefficients(theFit),nrow  = myK[1]))<0))
      theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F,paraPen =list(x=list(sp=value)) ) # penalize on x only ? sp: a vector of smoothing parameter values
      
      
    }
    
    
    yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
    
    Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew)))
    
    ##Estimate the First derivative fit
    Derv1Fit<--(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))-predict(theFit,data.frame(x=xNew_unnorm+tol1st,z=zNew)))/tol1st
    
    ##Estimate the Second derivative fit
    Derv2FitNonNeg<-(predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2
    #
    
  }
  
  ###Now for the error estimation
  
  # Residuals.V<-as.numeric(theFit$residuals)^2
  # m.sd<-list(penalty.sd,degree.sd)
  # 
  # ##Fit the residuals as in equation 20
  # fittedRes<-gam(Residuals.V~te(x,z,bs=c("ps","ps"),m=m.sd,k=myK),drop.intercept=F)
  # if(if_pred_new == TRUE){
  #   residualsEst<-abs(as.numeric(predict(fittedRes,data.frame(x=xNew_unnorm,z=zNew))))
  # }
  
  
  lambdaValues<-theFit$sp
  
  
  
  
  ##data returned (AT: return origx = x where z is the input(normalized moneyness))
  if(if_pred_new == TRUE){ # AT: newX and newX_unnorm both set to xNew_unnorm
    return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,newX = xNew_unnorm, newX_unnorm = xNew_unnorm, newZ=zNew,Fit=Fit,int.knots=knots,model=theFit))
  }else{
    return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,
                int.knots=knots,model=theFit))
  }
  
}


# version of fitBTPS_zc() with unnormalized moneyness, input x is unnormalized moneyness, plus for binned data in 'IV_real_binned_unnormalized_surface_final_step1.r'
# Note that difference with fitBTPS_zc_unnorm: now all use ptmday_frac (but actually scale-invariant in gam(te())): it normalize the data by default
fitBTPS_zc_unnorm_binned<-function(y,x,z,knots=c(0,0),degree=c(5,3),penalty=c(1,1),tol1st=.0000001,tol2nd=.1,degree.sd=c(3,3),penalty.sd=c(1,1),newXLim=c(min(x),max(x)),newZLim=c(min(z),max(z)),
                            new_X_Z_grid = c(100, 50),  nonDeg=F, penalty_bool = TRUE, if_pred_new = TRUE){
  
  ##Get a value for the number of B-spline Basis needed for the set up based on degree and knots
  myK<-c(knots+degree+2) # should be knots+degree+1: according to the ref book: degree here seems to be the (degree - 1) in our understanding, so m =2 represents a cubic spline
  # myK<-c(knots+degree+1) 
  ##create the degree and penalty list for the function
  m<-list(c(degree[1],penalty[1]),c(degree[2],penalty[2]))
  
  ##fit the BTPS using the gam function in r
  if(penalty_bool == TRUE){
    # browser()
    theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    
  }else{
    theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK, np = FALSE, fx = TRUE),drop.intercept=F) # should add  fx = TRUE in te() to turn off the penalty ?
  }
  
  
  
  if(if_pred_new == TRUE){
    ###create new data and fit the new data
    # xNew<-seq(newXLim[1],newXLim[2],length.out = 50)##values of length.out might not work on all computers due to memory limitations
    # zNew<-seq(newZLim[1],newZLim[2],length.out = 50)
    xNew_unnorm<-seq(newXLim[1],newXLim[2],length.out = new_X_Z_grid[1])##values of length.out might not work on all computers due to memory limitations
    # zNew<-seq(newZLim[1],newZLim[2],length.out = new_X_Z_grid[2])
    zNew<-seq(newZLim[1],newZLim[2],by = 1/new_X_Z_grid[2]) # this specification is used so that we can plot 2D plot on exact ptmday
    
    
    
    hold.dataframe<-expand.grid(xNew_unnorm,zNew)
    
    xNew_unnorm<-hold.dataframe$Var1
    # normalization 
    # xNew <- (xNew_unnorm - X_norm_stats[1]) / X_norm_stats[2]
    zNew<-hold.dataframe$Var2
    # Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew)))
    # Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew * 365))) # modified on 05.10: to plot the 2D plot surface
    Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))) # no need for * 365: b/c input is ptmday for x now
    
  }
  
  
  yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
  
  # ##Estimate the First derivative fit
  # Derv1Fit<--(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))-predict(theFit,data.frame(x=xNew_unnorm+tol1st,z=zNew)))/tol1st
  # 
  # ##Estimate the Second derivative fit
  # Derv2Fit<-(predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2
  # # browser()
  # Derv2FitNonNeg<-NA
  
  ##Ensure the fit is positive
  if(nonDeg){
    number<-sum(((predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2)<0)
    if(number>0){
      theFitInitial<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F )
      coefMatrix<-matrix(coefficients(theFit),ncol  = myK[1])
      
      Dmatrix<-matrix(0,nrow  = myK[1],ncol=myK[2])
      for(i in 1:(nrow(coefMatrix)-2)){
        Dmatrix[i+1,c((i):(i+2))]<-c(1,-2,1)
        
      }
      Dmatrix<-t(Dmatrix)
      
      
      value<-diag(as.numeric(colSums( Dmatrix* matrix(coefficients(theFit),nrow  = myK[1]))<0))
      theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F,paraPen =list(x=list(sp=value)) ) # penalize on x only ? sp: a vector of smoothing parameter values
      
      
    }
    
    
    yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
    
    Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew)))
    
    ##Estimate the First derivative fit
    Derv1Fit<--(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))-predict(theFit,data.frame(x=xNew_unnorm+tol1st,z=zNew)))/tol1st
    
    ##Estimate the Second derivative fit
    Derv2FitNonNeg<-(predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2
    #
    
  }
  
  ###Now for the error estimation
  
  # Residuals.V<-as.numeric(theFit$residuals)^2
  # m.sd<-list(penalty.sd,degree.sd)
  # 
  # ##Fit the residuals as in equation 20
  # fittedRes<-gam(Residuals.V~te(x,z,bs=c("ps","ps"),m=m.sd,k=myK),drop.intercept=F)
  # if(if_pred_new == TRUE){
  #   residualsEst<-abs(as.numeric(predict(fittedRes,data.frame(x=xNew_unnorm,z=zNew))))
  # }
  
  
  lambdaValues<-theFit$sp
  
  
  
  
  ##data returned (AT: return origx = x where z is the input(normalized moneyness))
  if(if_pred_new == TRUE){ # AT: newX and newX_unnorm both set to xNew_unnorm
    return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,newX = xNew_unnorm, newX_unnorm = xNew_unnorm, newZ=zNew,Fit=Fit,int.knots=knots,model=theFit))
  }else{
    return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,
                int.knots=knots,model=theFit))
  }
  
}


# one dimension version of fitBTPS_zc() with unnormalized moneyness, input x is unnormalized moneyness
fitBTPS_zc_unnorm_m<-function(y,x,knots=c(0),degree=c(5),penalty=c(1),tol1st=.0000001,tol2nd=.1,degree.sd=c(3),penalty.sd=c(1),newXLim=c(min(x),max(x)),
                              new_X_grid = c(100),  nonDeg=F, penalty_bool = TRUE, if_pred_new = TRUE){
  
  ##Get a value for the number of B-spline Basis needed for the set up based on degree and knots
  myK<-c(knots+degree+2) # should be knots+degree+1: according to the ref book: degree here seems to be the (degree - 1) in our understanding, so m =2 represents a cubic spline
  # myK<-c(knots+degree+1) 
  ##create the degree and penalty list for the function
  m<-list(c(degree[1],penalty[1]))
  
  ##fit the BTPS using the gam function in r
  if(penalty_bool == TRUE){
    # browser()
    theFit<-gam(y~te(x,bs=c("ps"),m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ? 
  }else{
    theFit<-gam(y~te(x,bs=c("ps"),m=m,k=myK,np = FALSE, fx = TRUE),drop.intercept=F) # should add  fx = TRUE in te() to turn off the penalty ?
  }
  
  
  
  if(if_pred_new == TRUE){
    ###create new data and fit the new data
    # xNew_unnorm<-seq(newXLim[1],newXLim[2],length.out = new_X_Z_grid[1])##values of length.out might not work on all computers due to memory limitations
    xNew_unnorm<-seq(newXLim[1],newXLim[2],length.out = new_X_grid[1])##values of length.out might not work on all computers due to memory limitations
    
    hold.dataframe<-expand.grid(xNew_unnorm)
    
    xNew_unnorm<-hold.dataframe$Var1
    # normalization 
    # xNew <- (xNew_unnorm - X_norm_stats[1]) / X_norm_stats[2]
    Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm)))
    
  }
  
  
  yPred<-as.numeric(predict(theFit,data.frame(x=x)))
  
  # ##Estimate the First derivative fit
  # Derv1Fit<--(predict(theFit,data.frame(x=xNew_unnorm))-predict(theFit,data.frame(x=xNew_unnorm+tol1st)))/tol1st
  # 
  # ##Estimate the Second derivative fit
  # Derv2Fit<-(predict(theFit,data.frame(x=xNew_unnorm-tol2nd))-2*predict(theFit,data.frame(x=xNew_unnorm))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd)))/tol2nd^2
  # # browser()
  # Derv2FitNonNeg<-NA
  
  
  ###Now for the error estimation
  
  # Residuals.V<-as.numeric(theFit$residuals)^2
  # m.sd<-list(penalty.sd,degree.sd)
  # 
  # ##Fit the residuals as in equation 20
  # fittedRes<-gam(Residuals.V~te(x,z,bs=c("ps","ps"),m=m.sd,k=myK),drop.intercept=F)
  # if(if_pred_new == TRUE){
  #   residualsEst<-abs(as.numeric(predict(fittedRes,data.frame(x=xNew_unnorm,z=zNew))))
  # }
  
  
  lambdaValues<-theFit$sp
  
  
  
  
  ##data returned (AT: return origx = x where z is the input(normalized moneyness))
  if(if_pred_new == TRUE){ # AT: newX and newX_unnorm both set to xNew_unnorm
    return(list(origY=y,origX=x,fitOrig=yPred,newX = xNew_unnorm, newX_unnorm = xNew_unnorm, Fit=Fit,int.knots=knots,model=theFit))
  }else{
    return(list(origY=y,origX=x,fitOrig=yPred,
                int.knots=knots,model=theFit))
  }
  
}




# version of fitBTPS_zc() with unnormalized moneyness, input x is unnormalized moneyness
fitBTPS_zc_unnorm_test<-function(y,x,z,knots=c(0,0),degree=c(5,3),penalty=c(1,1),tol1st=.0000001,tol2nd=.1,degree.sd=c(3,3),penalty.sd=c(1,1),newXLim=c(min(x),max(x)),newZLim=c(min(z),max(z)),
                            new_X_Z_grid = c(100, 50),  nonDeg=F, penalty_bool = TRUE, if_pred_new = TRUE){
  
  ##Get a value for the number of B-spline Basis needed for the set up based on degree and knots
  myK<-c(knots+degree+2) # should be knots+degree+1: according to the ref book: degree here seems to be the (degree - 1) in our understanding, so m =2 represents a cubic spline
  # myK<-c(knots+degree+1) 
  ##create the degree and penalty list for the function
  m<-list(c(degree[1],penalty[1]),c(degree[2],penalty[2]))
  
  ##fit the BTPS using the gam function in r
  if(penalty_bool == TRUE){
    # browser()
    # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    # theFit<-gam(y~te(x,z,bs=c("cc","cc"), m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    # theFit<-gam(y~te(x,z,bs=c("cc","cc"), m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    # if remove the m=m (actually no use -> b/c default cubic spline d_m = d_tau = 3)
    # theFit<-gam(y~te(x,z,bs=c("cc","cc"), m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    theFit<-gam(y~te(x,z,bs=c("cc","cc"), k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    
    # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    
  }else{
    # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK, np = FALSE, fx = TRUE),drop.intercept=F) # should add  fx = TRUE in te() to turn off the penalty ?
  }
  
  
  
  if(if_pred_new == TRUE){
    ###create new data and fit the new data
    # xNew<-seq(newXLim[1],newXLim[2],length.out = 50)##values of length.out might not work on all computers due to memory limitations
    # zNew<-seq(newZLim[1],newZLim[2],length.out = 50)
    xNew_unnorm<-seq(newXLim[1],newXLim[2],length.out = new_X_Z_grid[1])##values of length.out might not work on all computers due to memory limitations
    # zNew<-seq(newZLim[1],newZLim[2],length.out = new_X_Z_grid[2])
    zNew<-seq(newZLim[1],newZLim[2],by = 1/new_X_Z_grid[2]) # this specification is used so that we can plot 2D plot on exact ptmday
    
    
    
    hold.dataframe<-expand.grid(xNew_unnorm,zNew)
    
    xNew_unnorm<-hold.dataframe$Var1
    # normalization 
    # xNew <- (xNew_unnorm - X_norm_stats[1]) / X_norm_stats[2]
    zNew<-hold.dataframe$Var2
    # Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew)))
    Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew * 365))) # modified on 05.10: to plot the 2D plot surface
    
    # For the purpose of getting model matrix for the fine grids
    z_find_grid = zNew * 365
    theFit_new<-gam(Fit~te(xNew_unnorm, z_find_grid ,bs=c("cc","cc"), k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    
    
  }
  
  
  yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
  
  # ##Estimate the First derivative fit
  # Derv1Fit<--(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))-predict(theFit,data.frame(x=xNew_unnorm+tol1st,z=zNew)))/tol1st
  # 
  # ##Estimate the Second derivative fit
  # Derv2Fit<-(predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2
  # # browser()
  # Derv2FitNonNeg<-NA
  
  ##Ensure the fit is positive
  if(nonDeg){
    number<-sum(((predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2)<0)
    if(number>0){
      theFitInitial<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F )
      coefMatrix<-matrix(coefficients(theFit),ncol  = myK[1])
      
      Dmatrix<-matrix(0,nrow  = myK[1],ncol=myK[2])
      for(i in 1:(nrow(coefMatrix)-2)){
        Dmatrix[i+1,c((i):(i+2))]<-c(1,-2,1)
        
      }
      Dmatrix<-t(Dmatrix)
      
      
      value<-diag(as.numeric(colSums( Dmatrix* matrix(coefficients(theFit),nrow  = myK[1]))<0))
      theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F,paraPen =list(x=list(sp=value)) ) # penalize on x only ? sp: a vector of smoothing parameter values
      
      
    }
    
    
    yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
    
    Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew)))
    
    ##Estimate the First derivative fit
    Derv1Fit<--(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))-predict(theFit,data.frame(x=xNew_unnorm+tol1st,z=zNew)))/tol1st
    
    ##Estimate the Second derivative fit
    Derv2FitNonNeg<-(predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2
    #
    
  }
  
  ###Now for the error estimation
  
  # Residuals.V<-as.numeric(theFit$residuals)^2
  # m.sd<-list(penalty.sd,degree.sd)
  # 
  # ##Fit the residuals as in equation 20
  # fittedRes<-gam(Residuals.V~te(x,z,bs=c("ps","ps"),m=m.sd,k=myK),drop.intercept=F)
  # if(if_pred_new == TRUE){
  #   residualsEst<-abs(as.numeric(predict(fittedRes,data.frame(x=xNew_unnorm,z=zNew))))
  # }
  
  
  lambdaValues<-theFit$sp
  
  
  
  
  ##data returned (AT: return origx = x where z is the input(normalized moneyness))
  if(if_pred_new == TRUE){ # AT: newX and newX_unnorm both set to xNew_unnorm
    return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,newX = xNew_unnorm, newX_unnorm = xNew_unnorm, newZ=zNew,Fit=Fit,int.knots=knots,model_new=theFit_new))
  }else{
    return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,
                int.knots=knots,model=theFit))
  }
  
}


# version of fitBTPS_zc() with unnormalized moneyness, input x is unnormalized moneyness
fitBTPS_zc_unnorm_test2<-function(y,x,z,xz, knots=c(0,0),degree=c(5,3),penalty=c(1,1),tol1st=.0000001,tol2nd=.1,degree.sd=c(3,3),penalty.sd=c(1,1),newXLim=c(min(x),max(x)),newZLim=c(min(z),max(z)),
                                 new_X_Z_grid = c(100, 50),  nonDeg=F, penalty_bool = TRUE, if_pred_new = TRUE){
  
  ##Get a value for the number of B-spline Basis needed for the set up based on degree and knots
  # myK<-c(knots+degree+2) # should be knots+degree+1: according to the ref book: degree here seems to be the (degree - 1) in our understanding, so m =2 represents a cubic spline
  myK<-c(knots+degree+1)  # for s(bs = "bs") setup: https://stat.ethz.ch/R-manual/R-patched/library/mgcv/html/smooth.construct.bs.smooth.spec.html 
  ##create the degree and penalty list for the function
  # m<-list(c(degree[1],penalty[1]),c(degree[2],penalty[2]))
  m<-degree # for s(bs = 'bs')
  
  ##fit the BTPS using the gam function in r
  if(penalty_bool == TRUE){
    # browser()
    # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    # theFit<-gam(y~te(x,z,bs=c("cc","cc"), m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    # theFit<-gam(y~te(x,z,bs=c("cc","cc"), m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    theFit <- gam(y~s(x, k = myK, m = m, bs = "bs") + s(z, k = myK, m = m, bs="bs") + s(xz, k = myK, m = m, bs= "bs"),
               drop.intercept=F)
    # theFit <- gam(y~s(x, k = myK, m = m, bs = "bs") + s(z, k = myK, m = m, bs="bs"), drop.intercept=F) 
    
    # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    
  }else{
    # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK, np = FALSE, fx = TRUE),drop.intercept=F) # should add  fx = TRUE in te() to turn off the penalty ?
  }
  
  
  # browser()
  if(if_pred_new == TRUE){
    ###create new data and fit the new data
    # xNew<-seq(newXLim[1],newXLim[2],length.out = 50)##values of length.out might not work on all computers due to memory limitations
    # zNew<-seq(newZLim[1],newZLim[2],length.out = 50)
    xNew_unnorm<-seq(newXLim[1],newXLim[2],length.out = new_X_Z_grid[1])##values of length.out might not work on all computers due to memory limitations
    # zNew<-seq(newZLim[1],newZLim[2],length.out = new_X_Z_grid[2])
    zNew<-seq(newZLim[1],newZLim[2],by = 1/new_X_Z_grid[2]) # this specification is used so that we can plot 2D plot on exact ptmday
    
    
    
    hold.dataframe<-expand.grid(xNew_unnorm,zNew)
    
    xNew_unnorm<-hold.dataframe$Var1
    # normalization 
    # xNew <- (xNew_unnorm - X_norm_stats[1]) / X_norm_stats[2]
    zNew<-hold.dataframe$Var2
    
    xzNew = hold.dataframe$Var1 * hold.dataframe$Var2
    Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew, xz = xzNew))) # we are using ptmday_frac now 
    # Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew * 365))) # modified on 05.10: to plot the 2D plot surface
    
  }
  
  # prediction of the original data (insample)
  yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z, xz =xz)))
  
  # ##Estimate the First derivative fit
  # Derv1Fit<--(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))-predict(theFit,data.frame(x=xNew_unnorm+tol1st,z=zNew)))/tol1st
  # 
  # ##Estimate the Second derivative fit
  # Derv2Fit<-(predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2
  # # browser()
  # Derv2FitNonNeg<-NA
  
  ##Ensure the fit is positive
  if(nonDeg){
    number<-sum(((predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2)<0)
    if(number>0){
      theFitInitial<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F )
      coefMatrix<-matrix(coefficients(theFit),ncol  = myK[1])
      
      Dmatrix<-matrix(0,nrow  = myK[1],ncol=myK[2])
      for(i in 1:(nrow(coefMatrix)-2)){
        Dmatrix[i+1,c((i):(i+2))]<-c(1,-2,1)
        
      }
      Dmatrix<-t(Dmatrix)
      
      
      value<-diag(as.numeric(colSums( Dmatrix* matrix(coefficients(theFit),nrow  = myK[1]))<0))
      theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F,paraPen =list(x=list(sp=value)) ) # penalize on x only ? sp: a vector of smoothing parameter values
      
      
    }
    
    
    yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
    
    Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew)))
    
    ##Estimate the First derivative fit
    Derv1Fit<--(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))-predict(theFit,data.frame(x=xNew_unnorm+tol1st,z=zNew)))/tol1st
    
    ##Estimate the Second derivative fit
    Derv2FitNonNeg<-(predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2
    #
    
  }
  
  ###Now for the error estimation
  
  # Residuals.V<-as.numeric(theFit$residuals)^2
  # m.sd<-list(penalty.sd,degree.sd)
  # 
  # ##Fit the residuals as in equation 20
  # fittedRes<-gam(Residuals.V~te(x,z,bs=c("ps","ps"),m=m.sd,k=myK),drop.intercept=F)
  # if(if_pred_new == TRUE){
  #   residualsEst<-abs(as.numeric(predict(fittedRes,data.frame(x=xNew_unnorm,z=zNew))))
  # }
  
  
  lambdaValues<-theFit$sp
  
  
  
  
  ##data returned (AT: return origx = x where z is the input(normalized moneyness))
  if(if_pred_new == TRUE){ # AT: newX and newX_unnorm both set to xNew_unnorm
    return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,newX = xNew_unnorm, newX_unnorm = xNew_unnorm, newZ=zNew,Fit=Fit,int.knots=knots,model=theFit))
  }else{
    return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,
                int.knots=knots,model=theFit))
  }
  
}

# This function is for real data, when if_pred_new = TRUE, it is predicting the out-of-sample surface
# version of fitBTPS_zc() with unnormalized moneyness, input x is unnormalized moneyness
# xz, penalty_bool: no use
fitBTPS_zc_unnorm_test2_v2<-function(y,x,z,xz, knots=c(0,0),degree=c(5,3),penalty=c(1,1),tol1st=.0000001,tol2nd=.1,degree.sd=c(3,3),penalty.sd=c(1,1),newXLim=c(min(x),max(x)),newZLim=c(min(z),max(z)),
                                  new_X_Z_grid = c(100, 50),  nonDeg=F, penalty_bool = TRUE, if_pred_new = TRUE){
  
  ##Get a value for the number of B-spline Basis needed for the set up based on degree and knots
  # myK<-c(knots+degree+2) # should be knots+degree+1: according to the ref book: degree here seems to be the (degree - 1) in our understanding, so m =2 represents a cubic spline
  myK<-c(knots+degree+1)  # for s(bs = "bs") setup: https://stat.ethz.ch/R-manual/R-patched/library/mgcv/html/smooth.construct.bs.smooth.spec.html 
  ##create the degree and penalty list for the function
  # m<-list(c(degree[1],penalty[1]),c(degree[2],penalty[2]))
  m<-degree # for s(bs = 'bs')
  
  ##fit the BTPS using the gam function in r
  if(penalty_bool == TRUE){
    # browser()
    # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    # theFit<-gam(y~te(x,z,bs=c("cc","cc"), m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    # theFit<-gam(y~te(x,z,bs=c("cc","cc"), m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    # theFit <- gam(y~s(x, k = myK, m = m, bs = "bs") + s(z, k = myK, m = m, bs="bs") + s(xz, k = myK, m = m, bs= "bs"),
                  # drop.intercept=F)
    theFit <- gam(y~s(x, k = myK, m = m, bs = "bs") + s(z, k = myK, m = m, bs="bs"), drop.intercept=F)
    
    # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    
  }else{
    # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK, np = FALSE, fx = TRUE),drop.intercept=F) # should add  fx = TRUE in te() to turn off the penalty ?
  }
  
  
  
  if(if_pred_new == TRUE){
    ###create new data and fit the new data
    # xNew<-seq(newXLim[1],newXLim[2],length.out = 50)##values of length.out might not work on all computers due to memory limitations
    # zNew<-seq(newZLim[1],newZLim[2],length.out = 50)
    xNew_unnorm<-seq(newXLim[1],newXLim[2],length.out = new_X_Z_grid[1])##values of length.out might not work on all computers due to memory limitations
    # zNew<-seq(newZLim[1],newZLim[2],length.out = new_X_Z_grid[2])
    zNew<-seq(newZLim[1],newZLim[2],by = 1/new_X_Z_grid[2]) # this specification is used so that we can plot 2D plot on exact ptmday
    
    
    
    hold.dataframe<-expand.grid(xNew_unnorm,zNew)
    
    xNew_unnorm<-hold.dataframe$Var1
    # normalization 
    # xNew <- (xNew_unnorm - X_norm_stats[1]) / X_norm_stats[2]
    zNew<-hold.dataframe$Var2
    
    # xzNew = hold.dataframe$Var1 * hold.dataframe$Var2
    Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))) # use ptmday_frac now for both simu and real data
    # Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew * 365))) # No longer use like fitBTPS_zc_unnorm, use ptmday_frac now for both simu and real data

    # For the purpose of getting model matrix for the fine grids
    theFit_new <- gam(Fit~s(xNew_unnorm, k = myK, m = m, bs = "bs") + s(zNew, k = myK, m = m, bs="bs"), drop.intercept=F)
    
  }
  
  # prediction of the original data (insample)
  yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
  
  # ##Estimate the First derivative fit
  # Derv1Fit<--(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))-predict(theFit,data.frame(x=xNew_unnorm+tol1st,z=zNew)))/tol1st
  # 
  # ##Estimate the Second derivative fit
  # Derv2Fit<-(predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2
  # # browser()
  # Derv2FitNonNeg<-NA
  
  ##Ensure the fit is positive
  if(nonDeg){
    number<-sum(((predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2)<0)
    if(number>0){
      theFitInitial<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F )
      coefMatrix<-matrix(coefficients(theFit),ncol  = myK[1])
      
      Dmatrix<-matrix(0,nrow  = myK[1],ncol=myK[2])
      for(i in 1:(nrow(coefMatrix)-2)){
        Dmatrix[i+1,c((i):(i+2))]<-c(1,-2,1)
        
      }
      Dmatrix<-t(Dmatrix)
      
      
      value<-diag(as.numeric(colSums( Dmatrix* matrix(coefficients(theFit),nrow  = myK[1]))<0))
      theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F,paraPen =list(x=list(sp=value)) ) # penalize on x only ? sp: a vector of smoothing parameter values
      
      
    }
    
    
    yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
    
    Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew)))
    
    ##Estimate the First derivative fit
    Derv1Fit<--(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))-predict(theFit,data.frame(x=xNew_unnorm+tol1st,z=zNew)))/tol1st
    
    ##Estimate the Second derivative fit
    Derv2FitNonNeg<-(predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2
    #
    
  }
  
  ###Now for the error estimation
  
  # Residuals.V<-as.numeric(theFit$residuals)^2
  # m.sd<-list(penalty.sd,degree.sd)
  # 
  # ##Fit the residuals as in equation 20
  # fittedRes<-gam(Residuals.V~te(x,z,bs=c("ps","ps"),m=m.sd,k=myK),drop.intercept=F)
  # if(if_pred_new == TRUE){
  #   residualsEst<-abs(as.numeric(predict(fittedRes,data.frame(x=xNew_unnorm,z=zNew))))
  # }
  
  
  lambdaValues<-theFit$sp
  
  
  
  
  ##data returned (AT: return origx = x where z is the input(normalized moneyness))
  if(if_pred_new == TRUE){ # AT: newX and newX_unnorm both set to xNew_unnorm
    return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,newX = xNew_unnorm, newX_unnorm = xNew_unnorm, newZ=zNew,Fit=Fit,int.knots=knots,model_new=theFit_new)) # model=theFit is the insample fit model
  }else{
    return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,
                int.knots=knots,model=theFit))
  }
  
}



# This function is for simulation data, when if_pred_new = TRUE, it is predicting the in-sample surface
# Diff fitBTPS_zc_unnorm_test2_v2: 
# + Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew * 365))) # modified on 05.10: to plot the 2D plot surface
# + if_pred_new = TRUE, return (..., model=theFit) # for the in-sample surface
# xz, penalty_bool: no use
fitBTPS_zc_unnorm_test2_v2_simu<-function(y,x,z,xz, knots=c(0,0),degree=c(5,3),penalty=c(1,1),tol1st=.0000001,tol2nd=.1,degree.sd=c(3,3),penalty.sd=c(1,1),newXLim=c(min(x),max(x)),newZLim=c(min(z),max(z)),
                                     new_X_Z_grid = c(100, 50),  nonDeg=F, penalty_bool = TRUE, if_pred_new = TRUE){
  
  ##Get a value for the number of B-spline Basis needed for the set up based on degree and knots
  # myK<-c(knots+degree+2) # should be knots+degree+1: according to the ref book: degree here seems to be the (degree - 1) in our understanding, so m =2 represents a cubic spline
  myK<-c(knots+degree+1)  # for s(bs = "bs") setup: https://stat.ethz.ch/R-manual/R-patched/library/mgcv/html/smooth.construct.bs.smooth.spec.html 
  ##create the degree and penalty list for the function
  # m<-list(c(degree[1],penalty[1]),c(degree[2],penalty[2]))
  m<-degree # for s(bs = 'bs')
  
  ##fit the BTPS using the gam function in r
  if(penalty_bool == TRUE){
    # browser()
    # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    # theFit<-gam(y~te(x,z,bs=c("cc","cc"), m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    # theFit<-gam(y~te(x,z,bs=c("cc","cc"), m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    # theFit <- gam(y~s(x, k = myK, m = m, bs = "bs") + s(z, k = myK, m = m, bs="bs") + s(xz, k = myK, m = m, bs= "bs"),
    # drop.intercept=F)
    theFit <- gam(y~s(x, k = myK, m = m, bs = "bs") + s(z, k = myK, m = m, bs="bs"), drop.intercept=F)
    
    # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    
  }else{
    # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK, np = FALSE, fx = TRUE),drop.intercept=F) # should add  fx = TRUE in te() to turn off the penalty ?
  }
  
  
  
  if(if_pred_new == TRUE){
    ###create new data and fit the new data
    # xNew<-seq(newXLim[1],newXLim[2],length.out = 50)##values of length.out might not work on all computers due to memory limitations
    # zNew<-seq(newZLim[1],newZLim[2],length.out = 50)
    xNew_unnorm<-seq(newXLim[1],newXLim[2],length.out = new_X_Z_grid[1])##values of length.out might not work on all computers due to memory limitations
    # zNew<-seq(newZLim[1],newZLim[2],length.out = new_X_Z_grid[2])
    zNew<-seq(newZLim[1],newZLim[2],by = 1/new_X_Z_grid[2]) # this specification is used so that we can plot 2D plot on exact ptmday
    
    
    
    hold.dataframe<-expand.grid(xNew_unnorm,zNew)
    
    xNew_unnorm<-hold.dataframe$Var1
    # normalization 
    # xNew <- (xNew_unnorm - X_norm_stats[1]) / X_norm_stats[2]
    zNew<-hold.dataframe$Var2
    
    # xzNew = hold.dataframe$Var1 * hold.dataframe$Var2
    Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))) # use ptmday_frac now for both simu and real data
    # Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew * 365))) # No longer use like fitBTPS_zc_unnorm, use ptmday_frac now for both simu and real data
    
    # # For the purpose of getting model matrix for the fine grids
    # theFit_new <- gam(Fit~s(xNew_unnorm, k = myK, m = m, bs = "bs") + s(zNew, k = myK, m = m, bs="bs"), drop.intercept=F)
    
  }
  
  # prediction of the original data (insample)
  yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
  
  # ##Estimate the First derivative fit
  # Derv1Fit<--(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))-predict(theFit,data.frame(x=xNew_unnorm+tol1st,z=zNew)))/tol1st
  # 
  # ##Estimate the Second derivative fit
  # Derv2Fit<-(predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2
  # # browser()
  # Derv2FitNonNeg<-NA
  
  ##Ensure the fit is positive
  if(nonDeg){
    number<-sum(((predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2)<0)
    if(number>0){
      theFitInitial<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F )
      coefMatrix<-matrix(coefficients(theFit),ncol  = myK[1])
      
      Dmatrix<-matrix(0,nrow  = myK[1],ncol=myK[2])
      for(i in 1:(nrow(coefMatrix)-2)){
        Dmatrix[i+1,c((i):(i+2))]<-c(1,-2,1)
        
      }
      Dmatrix<-t(Dmatrix)
      
      
      value<-diag(as.numeric(colSums( Dmatrix* matrix(coefficients(theFit),nrow  = myK[1]))<0))
      theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F,paraPen =list(x=list(sp=value)) ) # penalize on x only ? sp: a vector of smoothing parameter values
      
      
    }
    
    
    yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
    
    Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew)))
    
    ##Estimate the First derivative fit
    Derv1Fit<--(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))-predict(theFit,data.frame(x=xNew_unnorm+tol1st,z=zNew)))/tol1st
    
    ##Estimate the Second derivative fit
    Derv2FitNonNeg<-(predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2
    #
    
  }
  
  ###Now for the error estimation
  
  # Residuals.V<-as.numeric(theFit$residuals)^2
  # m.sd<-list(penalty.sd,degree.sd)
  # 
  # ##Fit the residuals as in equation 20
  # fittedRes<-gam(Residuals.V~te(x,z,bs=c("ps","ps"),m=m.sd,k=myK),drop.intercept=F)
  # if(if_pred_new == TRUE){
  #   residualsEst<-abs(as.numeric(predict(fittedRes,data.frame(x=xNew_unnorm,z=zNew))))
  # }
  
  
  lambdaValues<-theFit$sp
  
  
  
  
  ##data returned (AT: return origx = x where z is the input(normalized moneyness))
  if(if_pred_new == TRUE){ # AT: newX and newX_unnorm both set to xNew_unnorm
    return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,newX = xNew_unnorm, newX_unnorm = xNew_unnorm, newZ=zNew,Fit=Fit,int.knots=knots,model=theFit)) # model=theFit is the insample fit model
  }else{
    return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,
                int.knots=knots,model=theFit))
  }
  
}


# version of fitBTPS_zc() with unnormalized moneyness, input x is unnormalized moneyness
fitBTPS_zc_unnorm_test3<-function(y,x,z,xz, knots=c(0,0),degree=c(5,3),penalty=c(1,1),tol1st=.0000001,tol2nd=.1,degree.sd=c(3,3),penalty.sd=c(1,1),newXLim=c(min(x),max(x)),newZLim=c(min(z),max(z)),
                                  new_X_Z_grid = c(100, 50),  nonDeg=F, penalty_bool = TRUE, if_pred_new = TRUE){
  
  ##Get a value for the number of B-spline Basis needed for the set up based on degree and knots
  # myK<-c(knots+degree+2) # should be knots+degree+1: according to the ref book: degree here seems to be the (degree - 1) in our understanding, so m =2 represents a cubic spline
  myK<-c(knots+degree+1)  # for s(bs = "bs") setup: https://stat.ethz.ch/R-manual/R-patched/library/mgcv/html/smooth.construct.bs.smooth.spec.html 
  ##create the degree and penalty list for the function
  # m<-list(c(degree[1],penalty[1]),c(degree[2],penalty[2]))
  m<-degree # for s(bs = 'bs')
  
  ##fit the BTPS using the gam function in r
  if(penalty_bool == TRUE){
    # browser()
    # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    # theFit<-gam(y~te(x,z,bs=c("cc","cc"), m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    # theFit<-gam(y~te(x,z,bs=c("cc","cc"), m=m,k=myK, np = FALSE),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    # theFit <- gam(y~s(x, k = myK, m = m, bs = "bs") + s(z, k = myK, m = m, bs="bs") + s(xz, k = myK, m = m, bs= "bs"),
    #               drop.intercept=F) 
    # theFit <- gam(y~ti(x, k = myK, m = m, bs = "ps") + ti(z, k = myK, m = m, bs="ps") + ti(x, z, k = c(myK,myK), m = c(m,m), bs= "ps"), drop.intercept=F)
    theFit <- gam(y~ti(x) + ti(z) + ti(x, z), drop.intercept=F)
    
    # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F) # supply knots knots=list(x=seq(0,2,knots),z=seq(0,1,length=knots)) ?
    
  }else{
    # theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK, np = FALSE, fx = TRUE),drop.intercept=F) # should add  fx = TRUE in te() to turn off the penalty ?
  }
  
  
  
  if(if_pred_new == TRUE){
    ###create new data and fit the new data
    # xNew<-seq(newXLim[1],newXLim[2],length.out = 50)##values of length.out might not work on all computers due to memory limitations
    # zNew<-seq(newZLim[1],newZLim[2],length.out = 50)
    xNew_unnorm<-seq(newXLim[1],newXLim[2],length.out = new_X_Z_grid[1])##values of length.out might not work on all computers due to memory limitations
    # zNew<-seq(newZLim[1],newZLim[2],length.out = new_X_Z_grid[2])
    zNew<-seq(newZLim[1],newZLim[2],by = 1/new_X_Z_grid[2]) # this specification is used so that we can plot 2D plot on exact ptmday
    
    
    
    hold.dataframe<-expand.grid(xNew_unnorm,zNew)
    
    xNew_unnorm<-hold.dataframe$Var1
    # normalization 
    # xNew <- (xNew_unnorm - X_norm_stats[1]) / X_norm_stats[2]
    zNew<-hold.dataframe$Var2
    
    xzNew = hold.dataframe$Var1 * hold.dataframe$Var2
    Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew, xz = xzNew))) # we are using ptmday_frac now 
    # Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew * 365))) # modified on 05.10: to plot the 2D plot surface
    
  }
  
  # prediction of the original data (insample)
  yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z, xz =xz)))
  
  # ##Estimate the First derivative fit
  # Derv1Fit<--(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))-predict(theFit,data.frame(x=xNew_unnorm+tol1st,z=zNew)))/tol1st
  # 
  # ##Estimate the Second derivative fit
  # Derv2Fit<-(predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2
  # # browser()
  # Derv2FitNonNeg<-NA
  
  ##Ensure the fit is positive
  if(nonDeg){
    number<-sum(((predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2)<0)
    if(number>0){
      theFitInitial<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F )
      coefMatrix<-matrix(coefficients(theFit),ncol  = myK[1])
      
      Dmatrix<-matrix(0,nrow  = myK[1],ncol=myK[2])
      for(i in 1:(nrow(coefMatrix)-2)){
        Dmatrix[i+1,c((i):(i+2))]<-c(1,-2,1)
        
      }
      Dmatrix<-t(Dmatrix)
      
      
      value<-diag(as.numeric(colSums( Dmatrix* matrix(coefficients(theFit),nrow  = myK[1]))<0))
      theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F,paraPen =list(x=list(sp=value)) ) # penalize on x only ? sp: a vector of smoothing parameter values
      
      
    }
    
    
    yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
    
    Fit<-as.numeric(predict(theFit,data.frame(x=xNew_unnorm,z=zNew)))
    
    ##Estimate the First derivative fit
    Derv1Fit<--(predict(theFit,data.frame(x=xNew_unnorm,z=zNew))-predict(theFit,data.frame(x=xNew_unnorm+tol1st,z=zNew)))/tol1st
    
    ##Estimate the Second derivative fit
    Derv2FitNonNeg<-(predict(theFit,data.frame(x=xNew_unnorm-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew_unnorm,z=zNew))+predict(theFit,data.frame(x=xNew_unnorm+tol2nd,z=zNew)))/tol2nd^2
    #
    
  }
  
  ###Now for the error estimation
  
  # Residuals.V<-as.numeric(theFit$residuals)^2
  # m.sd<-list(penalty.sd,degree.sd)
  # 
  # ##Fit the residuals as in equation 20
  # fittedRes<-gam(Residuals.V~te(x,z,bs=c("ps","ps"),m=m.sd,k=myK),drop.intercept=F)
  # if(if_pred_new == TRUE){
  #   residualsEst<-abs(as.numeric(predict(fittedRes,data.frame(x=xNew_unnorm,z=zNew))))
  # }
  
  
  lambdaValues<-theFit$sp
  
  
  
  
  ##data returned (AT: return origx = x where z is the input(normalized moneyness))
  if(if_pred_new == TRUE){ # AT: newX and newX_unnorm both set to xNew_unnorm
    return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,newX = xNew_unnorm, newX_unnorm = xNew_unnorm, newZ=zNew,Fit=Fit,int.knots=knots,model=theFit))
  }else{
    return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,
                int.knots=knots,model=theFit))
  }
  
}

##inputs: n- the number of data points to simulate
#YieldMeanFunc-the mean structure to use for the simulation
#"Linear": (-25+1.3*z) "Quad":  1.2*(z-50)+(z-150)^2/200
#YieldVarFunc-the variance of the mean structure distribution
#"Const":1 "NonConst":|z|^.2
#YieldError-the distribution of the yield structure
#"Normal":10N(0,1) or "Beta":50[Beta(alpha,Beta)-alpha/(alpha+beta)]
#alphaE-alpha value if the beta yield error is used
#betaE-beta value if the beta yield error is used
#p-value used for p in equation 2
#meas_error-the measurement error added to the estimated premium

##Returns a list with the following:
##Yield_Curr-simulated yield values
##Yield_Hist-simulated historical yield values (APH)
##Cov_Rate-simulated coverage rate values
##Premium-simulated premium values
##Obs_Prem-premium values with the simulated measurement error
##OverallSigma-sigma values for the simulation


yieldDataSim<-function(n,YieldMeanFunc="Linear",YieldVarFunc="Const",YieldError="Normal",alphaE=5,betaE=3,p=1,meas_error=.001){
  
  ####simulate x and z values
  xSim<-sample(seq(.55,.95,by=.05),n,replace=T)
  zSim<-rtnorm(n,200,25,100,300)
  
  ##create the yield mean based on linear or quadratic
  if(YieldMeanFunc=="Linear"){
    yieldMean<--25+1.3*zSim
  }
  
  if(YieldMeanFunc=="Quad"){
    
    yieldMean<-1.2*(zSim-50)+(zSim-150)^2/200
  }
  
  ##create the sigma based on non-constant and constant
  if(YieldVarFunc=="NonConst"){
    yieldSigma<-abs(yieldMean)^.2
    
  }
  
  if(YieldVarFunc=="Const"){
    
    yieldSigma<-rep(1,n)
  }
  
  ###input the YieldError (Normal or Beta) and estimate the premium value 
  premium<-rep(NA,n)
  if(YieldError=="Normal"){
    error<-rnorm(n,0,1)
    
    
    if(YieldMeanFunc=="Linear"){s.factor<-10}
    if(YieldMeanFunc=="Quad"){s.factor<-10}
    
    myIntFunc<-function(x,mu,sigma){
      
      
      return(pnorm((x-mu)/sigma,0,1))
      
    }
    
    
    for(i in 1:n){
      premium[i]<-p*integrate(myIntFunc,lower=zSim[i]+yieldSigma[i]*s.factor*-5,upper=xSim[i]*zSim[i],mu=yieldMean[i],sigma=yieldSigma[i]*s.factor)$value
    }
    
  }
  
  if(YieldError=="Beta"){
    nonCentral<-alphaE/(alphaE+betaE)
    error<-rbeta(n,alphaE,betaE)-nonCentral
    
    if(YieldMeanFunc=="Quad"){s.factor<-50}
    if(YieldMeanFunc=="Linear"){s.factor<-50}
    #print(s.factor)
    
    myIntFunc<-function(x,alpha,beta,mu,sigma){
      
      
      return(pbeta(alpha/(alpha+beta)+(x-mu)/sigma,alpha,beta))
      
    }
    bottomInt<--alphaE/(alphaE+betaE)
    
    for(i in 1:n){
      premium[i]<-p*integrate(myIntFunc,lower=zSim[i]+yieldSigma[i]*s.factor*bottomInt,upper=xSim[i]*zSim[i],mu=yieldMean[i],sigma=yieldSigma[i]*s.factor,alpha=alphaE,beta=betaE)$value
      
      
      
    }
  }
  
  ##Calculate the current yield
  
  currentYield<-yieldMean+yieldSigma*s.factor*error
  
  
  ###Return List
  return(list(Yield_Curr=currentYield,Yield_Hist=zSim,Cov_Rate=xSim,Premium=premium,Obs_Prem=premium+rnorm(length(premium),0,meas_error),OverallSigma=yieldSigma*s.factor))
}

#####

##inputs: bw- a charter string of the bandwidth selector to use for Gaussian kernels
#x-the values of which to estimate the bandwidth for


## outputs:returns vector of bandwidths evaluated at x
compute_bandwidth <- function(bw, x) {
  if (is.numeric(bw)) return(bw)
  if (is.character(bw)) bw <- match.fun(bw)
  
  if (is.function(bw)) {
    bw <- bw(x)
    message("Using bandwidth ", format(bw, digits = 3))
    bw
  } else {
    stop("Invalid bandwidth")
  }
}

##inputs: xOrig-x values used for fitting a bivariate kernel density
#yOrig-y values used for fitting a bivariate kernel density
#(xOrig, yOrig) are a paired observation
#xNew-x values on which to evaluate the fitted a bivariate kernel density
#yNew-y values on which to evaluate the fitted a bivariate kernel density
#(xNew, yNew) are a paired observation

## outputs: returns a list with a vector of density values for the (xNew, yNew) observations
predictDensity <- function(xOrig, yOrig,xNew,yNew) {
  
  n_out <- length(xNew)
  n_in<-length(xOrig)
  xh <- suppressMessages(compute_bandwidth("bw.nrd", xOrig))
  yh <- suppressMessages(compute_bandwidth("bw.nrd", yOrig))
  ax <- outer(xNew, xOrig, "-") / xh
  ay <- outer(yNew, yOrig, "-") / yh
  myDensity<- rowSums(matrix(stats::dnorm(ax) * stats::dnorm(ay), n_out, n_in) * 
                        1 / (length(xOrig) * xh * yh))
  
  return(list(myDensity))
}

##inputs: BTPS.obj-x values used for fitting a bivariate kernel density
#Derv: in vector form, the partial derivatives the variance should be calculated for 

## outputs:a list with the following values
#xValues: a vector of the corresponding x values of the variance values
#zValues: a vector of the corresponding z values of the variance values
#variance: a vector of the estimated variance for the points (xValues,zValues)

VarianceEstimator<-function(BTPS.obj,Derv=c(2,0),CropData=T){
  
  #calculate n_1 and n_1 as well as knots
  n_ind<-c(length(unique(BTPS.obj$origX)),length(BTPS.obj$origY)/length(unique(BTPS.obj$origX)))
  knots<-BTPS.obj$int.knots
  intHvalues<-c(NA,NA)
  hValues<-c(NA,NA)
  
  #get the needed info from the gam object
  penalties<-c( BTPS.obj$model$smooth[1][[1]]$margin[[1]]$null.space.dim, BTPS.obj$model$smooth[1][[1]]$margin[[2]]$null.space.dim)
  smoothingPenalties<-c(max.(as.numeric(BTPS.obj$model$sp[1]),1),max(as.numeric(BTPS.obj$model$sp[2]),1))
  for(i in 1:2){
    
    
    # function for calculating H as seen in appendix A
    Hfunction<-function(x,dervs=2,penalty=2){
      if(penalty==1){equ<-expression((1/(2)*exp(-x)), "x")}
      if(penalty==2){equ<-expression(1/(2*sqrt(2))*exp(-1/(sqrt(2))*x)*(sin(x/sqrt(2))+cos(x/sqrt(2))), "x")}
      if(penalty==3){equ<-expression(1/(6)*exp(-x)+exp(-1/(2)*x)*(cos(x*sqrt(3)/2)/6+sin(x*sqrt(3)/2)*sqrt(3)/6), "x")}
      if(penalty==4){equ<-expression(exp(-.9239*x)*(.231*cos(.3827*x)+.0957*sin(.3827*x))+exp(-.3827*x)*(.0957*cos(.9239*x)+.231*sin(.9239*x)), "x")}
      
      if(dervs==0){equ<-equ[1]}
      if(dervs>0){
        for(i in 1:dervs){
          equ<-D(equ,"x")
          
        }
      }
      
      
      return(eval(equ)^2)
      
    }
    
    ##integration of the H function
    intHvalues[i]<-integrate(Hfunction,lower=.00001,upper=Inf,dervs=Derv[i],penalty=penalties[i])$value
    
    ##calculation of the h values as given in lemma 1
    hValues[i]<-((((smoothingPenalties[i]*knots[i]/n_ind[i]))^(1/(2*penalties[i])))/knots[i])^(2*Derv[i]+1)
    n<-1
    
  }
  finalH<-1
  
  #combine the h_1 and h_2 values into h
  if(CropData){
    finalH<-hValues[1]*hValues[2]
    
    
  }
  
  #estimate f hat by a bivariate kernel density estimator
  kernelDensity<-predictDensity(BTPS.obj$origX,BTPS.obj$origZ,BTPS.obj$newX,BTPS.obj$newZ)
  
  #combine into the final variance
  EstVar<-4*intHvalues[1]*intHvalues[2]*BTPS.obj$smoothedSigma/unlist(kernelDensity)/(finalH*length(BTPS.obj$origX))
  return(list(xValues=BTPS.obj$newX,zValues=BTPS.obj$newZ,variance=EstVar))
}




##inputs: w: the dependent variable on which to fit the spline
#z: the independent variable on which to fit the spline
#m: a vector of c(degree,penalty) used to fit the p-spline
#knots: number of interior knots to use 
#fixedZ:the z values at which to predict the fitted kernel estimator with
#predictW: the w values at which to predict the fitted kernel estimator with

## outputs: fHat a dataframe  of the fitted kernel estimator values at points  with columns fixedZ and predictW rows)

kernelEstimator<-function(w,z,m,knots,fixedZ,predictW){
  
  ##fit the p-spline to the data
  theFit<-gam(w~s(z,bs="ps",m=m,k=knots),drop.intercept=F)
  
  
  #fit the squared residuals to get sigma^2
  sigmaSquare<-gam(residuals(theFit)^2~s(z,bs="ps",m=m,k=knots),drop.intercept=F)
  
  
  #Get the estimated mean and sigma^2
  Predicted_Mu<-as.numeric(predict(theFit,data.frame(z=fixedZ)))
  
  Predicted_SigmaSquare<-abs(as.numeric(predict(sigmaSquare,data.frame(z=fixedZ))))
  
  
  
  
  
  
  epsilon<-residuals(theFit)/sqrt(predict(sigmaSquare))
  
  #Iterate through for each value of w and z
  fHat<-data.frame(predictW)
  names(fHat)<-"w"
  for(i in 1:length(fixedZ)){
    
    minVal<-(min(fHat$w,na.rm=T)-Predicted_Mu[i])/sqrt(Predicted_SigmaSquare[i])
    maxVal<-(max(fHat$w,na.rm=T)-Predicted_Mu[i])/sqrt(Predicted_SigmaSquare[i])
    ifBad<-tryCatch(density(epsilon,n=length(predictW),from=minVal, to=maxVal,na.rm=T),error=function(e){"ERROR"})
    if(ifBad[1]=="ERROR"){browser()}
    
    mydensity<-density(epsilon,n=length(predictW),from=minVal, to=maxVal,na.rm=T)
    mydensity$y
    
    fHat<-cbind(fHat,mydensity$y/sqrt(Predicted_SigmaSquare[i]))
    names(fHat)[ncol(fHat)]<-fixedZ[i]
    
  }
  
  return(fHat)
}

##inputs: w: the dependent variable on which to fit the spline
#z: the independent variable on which to fit the spline
#degree: the degree of penalty used in the p-spline
#knots: number of interior knots to use 
#penalty: The penalty of the p-spline used
#fixedZ: the z values at which to predict the fitted kernel estimator with
#NumW: the number of w values to estimate the kernel density estimate at
#wRange: the range of w values to estimate the kernel density estimate at


## outputs: list with the following values
#estimate: The kernel density estimate
#SD: the standard error of the kernel density estimate
#fixedZvalue: the z values that the kernel density was estimated at
fitKernel<-function(w,z,knots=1,degree=2,penalty=2,fixedZ=NA,NumW=1000,wRange=NA){
  
  #put knots kernel and degree and penalty into form to use gam function
  myK<-c(knots+degree+1)
  
  m<-c(degree,penalty)
  
  
  
  ####Kernel density estimation
  #set up the values at which to estimate at (default 5 z values)
  if(is.na(wRange[1])){
    
    new_w<-seq(min(w),max(w),length.out=NumW)
  }else{   new_w<-seq(min(wRange),max(wRange),length.out=NumW)}
  
  if(is.na(fixedZ[1])){
    
    predictionValues<-seq(min(z),max(z),length.out = 5)
    
    
  }else{predictionValues<-fixedZ}
  
  ##get the kernel density estimates
  finalFhat<-as.matrix(kernelEstimator(w,z,m,myK,predictionValues,new_w))
  
  
  ###Now for the error estimation using the jackknife method
  standardErrorMat<-array(,dim=c(NumW,length(predictionValues)+1,length(w)))
  for(i in 1:length(w)){
    
    standardErrorMat[,,i]<-(as.matrix(kernelEstimator(w[-i],z[-i],m,myK,predictionValues,new_w))-finalFhat)^2
    
    
  }
  
  standardErrorF<-apply(standardErrorMat, c(1,2), sum,na.rm=T)*(length(w)-1)/length(w)
  
  
  ##return the values
  return(list(estimate=finalFhat,SD=sqrt(standardErrorF),fixedZvalue=predictionValues))
  
}



##inputs: w: the dependent variable on which to fit the spline
#z: the independent variable on which to fit the spline
#t: the year of the historical yield value
#m: a vector of c(degree,penalty) used to fit the p-spline
#knots: number of interior knots to use 
#fixedZ:the z values at which to predict the fitted kernel estimator with
#predictW: the w values at which to predict the fitted kernel estimator with

## outputs: fHat a dataframe  of the fitted kernel estimator values at points  with columns fixedZ and predictW rows)
kernelEstimatorEmp<-function(w,z,t,m,knots,fixedZ,predictW){
  
  t2<-t^2
  t3<-t^3
  theFit<-gam(w~t+t2+t3+s(z,bs="ps",m=m,k=knots),drop.intercept=F)
  
  
  
  sigmaSquare<-gam(residuals(theFit)^2~s(z,bs="ps",m=m,k=knots),drop.intercept=F)
  
  Predicted_Mu<-  as.numeric(predict(theFit,data.frame(t=rep(0,length(fixedZ)),t2=rep(0,length(fixedZ)),t3=rep(0,length(fixedZ)),z=fixedZ))) 
  
  Predicted_SigmaSquare<-abs(as.numeric(predict(sigmaSquare,data.frame(t=rep(0,length(fixedZ)),t2=rep(0,length(fixedZ)),t3=rep(0,length(fixedZ)),z=fixedZ))))
  
  
  
  epsilon<-residuals(theFit)/sqrt(predict(sigmaSquare))
  
  
  fHat<-data.frame(predictW)
  
  names(fHat)<-"w"
  
  for(i in 1:length(fixedZ)){
    
    
    
    minVal<-(min(fHat$w,na.rm=T)-Predicted_Mu[i])/sqrt(Predicted_SigmaSquare[i])
    
    maxVal<-(max(fHat$w,na.rm=T)-Predicted_Mu[i])/sqrt(Predicted_SigmaSquare[i])
    
    ifBad<-tryCatch(density(epsilon,n=length(predictW),from=minVal, to=maxVal,na.rm=T),error=function(e){"ERROR"})
    
    if(ifBad[1]=="ERROR"){browser()}
    
    
    
    mydensity<-density(epsilon,n=length(predictW),from=minVal, to=maxVal,na.rm=T)
    
    
    # fHat<-cbind(fHat,(predictW-Predicted_Mu[i])/sqrt(Predicted_SigmaSquare[i])/var(mydensity$x),mydensity$y*var(mydensity$x))
    
    fHat<-cbind(fHat,mydensity$y/sqrt(Predicted_SigmaSquare[i]))
    names(fHat)[ncol(fHat)]<-fixedZ[i]
    
  }
  
  return(fHat)
  
}


##inputs: w: the dependent variable on which to fit the spline
#z: the independent variable on which to fit the spline
#t: the year of the historical yield value
#degree: the degree of penalty used in the p-spline
#knots: number of interior knots to use 
#penalty: The penalty of the p-spline used
#fixedZ: the z values at which to predict the fitted kernel estimator with
#NumW: the number of w values to estimate the kernel density estimate at
#wRange: the range of w values to estimate the kernel density estimate at
#baseYear: the year that the yield density is being estimated for

## outputs: list with the following values
#estimate: The kernel density estimate
#SD: the standard error of the kernel density estimate
#fixedZvalue: the z values that the kernel density was estimated at
fitKernelEmp<-function(w,z,t,knots=1,degree=2,penalty=2,fixedZ=NA,NumW=1000,wRange=NA,sd=F,baseYear=2009){
  t<-baseYear-t
  myK<-c(knots+degree+1)
  
  m<-c(degree,penalty)
  ####Kernel density estimation
  
  if(is.na(wRange[1])){
    
    new_w<-seq(min(w),max(w),length.out=NumW)
    
  }else{   new_w<-seq(wRange[1],wRange[2],length.out=NumW)}
  
  if(is.na(fixedZ[1])){
    predictionValues<-seq(min(z),max(z),length.out = 5)
  }else{predictionValues<-fixedZ}
  finalFhat<-as.matrix(kernelEstimatorEmp(w,z,t,m,myK,predictionValues,new_w))
  ###Now for the error estimation
  
  standardErrorF<-1
  standardErrorMat<-array(,dim=c(NumW,length(predictionValues)+1,length(w)))
  
  if(sd){
    for(i in 1:length(w)){
      
      standardErrorMat[,,i]<-(as.matrix(kernelEstimatorEmp(w[-i],z[-i],t[-i],m,myK,predictionValues,new_w))-finalFhat)^2
      
      if(i%%10==0){print(i)}
    }
    
    standardErrorF<-apply(standardErrorMat, c(1,2), sum,na.rm=T)*(length(w)-1)/length(w)
    
  }
  
  
  
  
  ##fit the residuals
  
  return(list(estimate=finalFhat,SD=sqrt(standardErrorF),fixedZvalue=predictionValues))
  
  
  
}
