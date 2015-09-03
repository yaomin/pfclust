pfcluster.estVar <- function(d,dat,method=c("recursive","cut"),n.min=10,h,debug=F)
  ## estimate variance of the data 'dat' based on dendrogram 'd'.
  ## recursive method refers to iterative estimating for nodes with size > n.min
  ## cut method refers to pooled variance estimate using the clusters formed by cutting the d at h
  {

    ##
    method <- match.arg(method)
    if(method=="recursive") {
      if(is.leaf(d)) attr(d,'var') <- NA
      else if(attr(d,'members') <= n.min) #attr(d,'var') <- var(dat[labels(d),])
        ret <- var(dat[labels(d),])
      else if(is.leaf(d[[1]]) & !is.leaf(d[[2]])) {
        if(debug) cat('R\n')
        ret <- pfcluster.estVar(d[[2]],
                                        dat[labels(d[[2]]),],
                                        n.min=n.min)
        #attr(d,'var') <- attr(d.2,'var')
      }
      else if(is.leaf(d[[2]]) & !is.leaf(d[[1]])) {
        if(debug) cat('L\n')
        ret <- pfcluster.estVar(d[[1]],
                                        dat[labels(d[[1]]),],
                                        n.min=n.min)
        #attr(d,'var') <- attr(d.1,'var')
      }
      else {
        if(debug) cat('RL\n')
        var1 <- pfcluster.estVar(d[[1]],
                                        dat[labels(d[[1]]),],
                                        n.min=n.min)
        var2 <- pfcluster.estVar(d[[2]],
                                        dat[labels(d[[2]]),],
                                        n.min=n.min)
        #var1 <- attr(d.1,'var')
        #var2 <- attr(d.2,'var')
        n1 <- attr(d[[1]],'members')
        n2 <- attr(d[[2]],'members')
        ##print(paste(var1,var2,n1,n2))
        ret <- ((n1-1)*var1+(n2-1)*var2)/(n1+n2-2)
        #ret <- var0
      }
      ##ret <- attr(d,'var')
    }
    else if(method=="cut") {
      d.cut <- cut(d,h=h)
      var.denom <- 0
      var.num <- 0
      for(i in seq(length(d.cut$lower))) {
        n.members <- attr(d.cut$lower[[i]],'members')
        if(n.members > n.min) {
          var.num <- var.num+(n.members-1)*var(dat[labels(d.cut$lower[[i]]),])
          var.denom <- var.denom+(n.members-1)
        } else {
          var.num <- var.num
          var.denom <- var.denom
        }
        ret <- var.num/var.denom
      }
    }
    ret 
  }
