# clustreg (Clusterwise Regression Function)
# Joshua Peng
# The Data Incubator Challenge Q3

# Inputs:
# d = data set with numerical variables
#     first column in d must be the response variable
# k = number of clusters
# t = number of procedure iterations
# q = set seed for random sampling: set.seed(q)

clustreg = function(d,k,t,q) {
  
  rsq = c()       # array of R-Squared values
  d0 = d[,-1]     # data set without response variable
  n = dim(d)[1]   # number of observations
  p = dim(d)[2]   # number of columns
  ybar = mean(d[,1])
  xnames = names(d) 
  
  # Setting up multiple regression formula
  # formula = paste(xnames[1],'~')
  # for(i in 2:p) {
  #   if(i==2) {
  #     formula = paste(formula,xnames[i])
  #   } else {
  #     formula = paste(formula,'+',xnames[i])
  #   }
  # }
  formula = paste(xnames[1],'~.') # this single line is equivalent to the above code
  
  # Initializing variables
  s = floor(n/k)  # initial size of each cluster
  r = n%%k        # remainder, included in last cluster
  set.seed(q)     # set seed for reproducible results
  ss = sample(1:n,n,replace=FALSE)    # randomized indices
  
  
  # Setting up initial cluster indices 
  # Initial clusters each have n/k observations from the random sample array
  # The last cluster starts with n/k+r observations if n/k is not an integer
  c = list()
  for(i in 1:k) {
    if(i==k) {
      c[[i]] = ss[((i-1)*s+1):(n)]
    } else {
      c[[i]] = ss[((i-1)*s+1):(i*s)]
    }
  }
  
  # Setting up initial clusters
  cc = list()         # data in k clusters
  ccc = list()        # data in k clusters without response variable
  for(i in 1:k) {
    cc[[i]] = d[c[[i]],]
    ccc[[i]] = d0[c[[i]],]
  }
  
  for(ii in 1:t) {    # runs loop t times
    
    # Storing multiple regression models in a list
    fit = list()
    for(i in 1:k) {
      fit[[i]]=lm(formula,data=cc[[i]])
    }
    
    # Calculating Squared Residual Error (SRE) for all observations
    pred = matrix(data = NA, n, k)
    res = matrix(data = NA, n, k)
    actual = matrix(data = NA, n, k)
    sre = matrix(data = NA, n, k)
    for(i in 1:k) {
      pred[,i] = predict(fit[[i]], d0)
      actual[,i] = d[,1]
      res[,i] = actual[,i]-pred[,i]
      sre[,i] = res[,i]^2
    }
    
    # Reclassifying observations to the cluster with the minimum SRE
    # Computing overall R-Squared value for all cluster regression models combined
    a = c()
    rsqnum = c()
    rsqden = c()
    for(i in 1:n) {
      for(j in 1:k) {
        if(min(sre[i,])==sre[i,j]) {
          a[i] = j
          rsqnum[i] = (pred[i,j]-ybar)^2
          rsqden[i] = (actual[i,j]-ybar)^2
        }
      }
    }
    rsq[ii] = sum(rsqnum)/sum(rsqden)
    
    # Reforming new clusters for next iteration
    # Refitting k linear regression models
    c = list()
    fit = list()
    crsq = c()
    for(i in 1:k) {
      c[[i]] = which(a==i)
      cc[[i]] = d[c[[i]],]
      ccc[[i]] = d0[c[[i]],]
      fit[[i]] = lm(formula,data=cc[[i]])
      crsq[i] = summary(fit[[i]])$r.squared
    }
    
    # If the overall R-Squared value is the same as in last iteration,
    # the procedure has converged and we will break out of the loop
    if(ii>1) {
      if(rsq[ii]==rsq[ii-1]) break
    }
    
    # End of iteration loop
  }
  
  #output = list()
  #output[[1]] = fit  # return final regression models for each cluster
  #output[[2]] = rsq  # return overall R-Squared value for each iteration (combined from all clusters)
  #output[[3]] = crsq # return R-Squared value for the final regression model for each cluster
  #output[[4]] = c    # return cluster list of observation assignments
  #output[[5]] = a    # return observation list of cluster classifications
  #output[[6]] = ii   # total number of iterations
  output = list("fit" = fit, 
                "rsq" = rsq, 
                "crsq" = crsq, 
                "clusters" = c, 
                "obs" = a, 
                "i" = ii)
  return(output)
  # End of clusterwise regression function
}


# Example use of clusterwise regression function

d = read.csv('wages.csv', header=T) # assign data set to d
d = d[,c(11,3,4,6,7,9)] # remove columns with non-numerical data, move response variable to 1st column
d = na.omit(d) # remove missing data
colnames(d) <- c("ratio", "numfemale", "avgfemalelongevity", "nummale", "avgmalelongevity", "total")
d$ratio = as.numeric(gsub('%', '', d$ratio))
k = 3           # number of clusters
t = 50          # number of procedure iterations
q = 11          # set seed for reproducible results

m1 = clustreg(d,k,t,q)
m1$fit          # return final regression models for each cluster
m1$rsq          # 0.07030112 0.33305159 0.41986002 0.45735178 0.47366001 0.55663713 0.57654005 0.60526745 0.60143458 0.61910452
                # 0.72812594 0.72827728 0.71870780 0.74677250 0.75051372 0.75095160 0.75115813 0.75139433 0.75214741 0.75399931
                # 0.75568833 0.75406967 0.75859343 0.76985319 0.75896561 0.75589443 0.75603665 0.75609102 0.75609102
m1$crsq         # average cluster R-squared: 0.07530181 0.72338881 0.76866784
m1$clusters     # return cluster list of observation assignments
m1$obs          # return observation list of cluster classifications
m1$i            # converged on 29 iterations

plot((1:m1$i), m1$rsq, ylim=c(0,1), type="l", col=4, 
     main = "Variance Accounted For with 3 Clusters", 
     ylab = "Variance Accounted For", xlab="Number of Iterations")

k = 5           # number of clusters
m2 = clustreg(d,k,t,q)
m2$fit          # return final regression models for each cluster
m2$rsq          # 0.1016157 0.3562966 0.5164390 0.6384988 0.7706268 0.8240548 0.8511634 0.8560073 0.8548829 0.8556782 0.8583321
                # 0.8583232 0.8594951 0.8587382 0.8587581 0.8587581
m2$crsq         # average cluster R-squared: 0.6910960 0.7727334 0.8075314 0.9055918 0.3631383
m2$clusters     # return cluster list of observation assignments
m2$obs          # return observation list of cluster classifications
m2$i            # converged on 16 iterations

plot((1:m2$i), m2$rsq, ylim=c(0,1), type="l", col=4, 
     main = "Variance Accounted For with 5 Clusters", 
     ylab = "Variance Accounted For", xlab="Number of Iterations")