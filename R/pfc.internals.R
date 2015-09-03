common.expandIdx <- function(n,m,out.as.vector=F,self.included=T)
  # Function to calculate the unique combination of 1:n & 1:m.
  # so it will list all (1,1), (1,2), (1,3), ... but each one is a unique combination ( (1,2) and (2,1) are duplicate combination)
{
    comb <- expand.grid(1:n,1:m)
    if(self.included) comb.unique1 <- comb[comb[,1]>=comb[,2],]
    else comb.unique1 <- comb[comb[,1]>comb[,2],]
    comb.unique2 <- comb[comb[,2]>n,]
    return(as.matrix(rbind(comb.unique1,comb.unique2)))
}

heatColors <-colorRampPalette(c("#00007F", "blue",  "#007FFF",
                                "cyan", "#7FFF7F", "yellow",
                                "#FF7F00", "red", "#7F0000"))

## filter out the genes that have less interesting proflies
## filtering is based on 'fun.summary(abs(x-fun.reference(unlist(x))))'
## genes with high above values will be kept, e.g. top 20% keep
##
pfc.filter <- function(geneDt,
                       pct.keep,
                       fun.summary=max,
                       fun.reference=mean,
                       fun.transform=log,
                       return.summary=F)
  # geneDt: gene data with the col's the samples and the rows the genes
  # fun.summary: the function to calculate the summary measure of the magnitude of the profile,
  #              the options are 'max' or 'median'
  # fun.reference: the funciton to caculate the average reference line of the profile which is used
  #                to calculate the magnitude above, the options are 'mean' or 'median', the magnitude
  #                is calculated as the absolute difference between individual observations and this reference.
  {
    geneDt <- as.data.frame(geneDt)
    summary.measure <- apply(geneDt,1,
                             function(x,fun.ref) fun.summary(abs(x-fun.ref(unlist(x)))),
                             fun.ref=fun.reference)
    if (!is.null(fun.transform)) summary.measure <- fun.transform(summary.measure)
    toKeep <- names(summary.measure[summary.measure > quantile(summary.measure,probs=1-pct.keep)])
    
    filteredGene <- geneDt[toKeep,]
    if(return.summary) return(list(filteredGene=filteredGene,summary.measure=summary.measure))
    else return(filteredGene)
  }

pfc.getCoherentGroups <- function(d.coh)
  ## get the groups members
  {
    grps <- vector(mode="character")
    grpPos <- list()
    getgid <- function(X) {
      if(!is.leaf(X)) {
        if(!is.null(gid <- attr(X,"grpid"))) {
          grps <<- c(grps,gid)
          grpPos <<- c(grpPos,list(attr(X,"pos")))
          
        }
        else {
          for (j in seq(length=length(X))) {
            X[[j]] <- getgid(X[[j]]) 
          }
        }
      }
      X
    }
    getgid(d.coh)
    
    grpMembers <- vector(mode="list",length=length(grps))
    grpCohereIndex <- vector(mode="numeric",length=length(grps))
    i <- 1
    getmembers <- function(X) {
      if(!is.leaf(X)) {
        if(!is.null(attr(X,"grpid"))) {
          grpMembers[[i]] <<- labels(X)
          grpCohereIndex[i] <<- attr(X,"coherenceIndex")
          i <<- i+1
          
        }
        else {
          for (jj in seq(length=length(X))) {
            ## length(X) will always be 2
            X[[jj]] <- getmembers(X[[jj]]) 
          }
        }
      }
      X
    }
    getmembers(d.coh)
    list(grpid=grps,grpPos = grpPos,members=grpMembers,cohereIndex=grpCohereIndex)
  }
    

pfc.grpRepresentatives <- function(PfCluster,
                                   method=c('average','lowess','spline'),
                                   sp.spar=NULL,
                                   lw.f = 2/3
                                   )
  ## get the representatives for each profile groups.
  {
    method = match.arg(method)
    grp.n <- length(PfCluster$cohGroups$grpid)
    grp.rep <- list()
    
    for ( i in 1:grp.n) {
     
      grp.dt <- PfCluster$dt[PfCluster$cohGroups$members[[i]],]
      if (method == 'average') {
        grp.rep[[i]] <- unlist(colMeans(grp.dt,na.rm=TRUE),use.names=F)
      }
      else if (method == 'lowess')  {
        y <- unlist(grp.dt,use.names=F)
        x <- rep(1:dim(grp.dt)[2],each=dim(grp.dt)[1])
        grp.rep[[i]] <- lowess(x,y,f=lw.f)
      }
      else if (method == 'spline') {
        y <- unlist(grp.dt,use.names=F)
        x <- rep(1:dim(grp.dt)[2],each=dim(grp.dt)[1])
        grp.rep[[i]] <- smooth.spline(x,y,spar=sp.spar)
      }
      else {
        stop("Wrong method")
      }
    } 

    grp.rep
  }

pfc.pairwise.FUN <- function(X,MARGIN=1,self.included=F,FUN='-',...)
  ## by row: by =1
  ## by column: by =2
  {
    this.func <- function(this.index,x,fun,which.margin,...) {
      to.eval <- switch(which.margin,
                        '1'=do.call(fun,
                          list(x1=x[this.index[1],],
                               x2=x[this.index[2],],
                               ...)),
                        '2'=do.call(fun,
                          list(x1=x[,this.index[1]],
                               x2=x[,this.index[2]],
                               ...))
                        )
      eval(to.eval)
    }

    n <- dim(X)[MARGIN]
    idx <- common.expandIdx(n,n,self.included=self.included)

    ret <- apply(idx,
                 MARGIN,
                 this.func,
                 x=X,
                 fun=FUN,
                 which.margin=MARGIN,
                 ...)
    ret
  }

pfc.plotProfile <- function(x,
                            order=NULL,
                            toScale=F,
                            add = F,
                            plotProfile = T,
                            ylimit=NULL,
                            v=NULL,vtype="l",vlty="solid",vcol="blue",vlwd=1,
                            pcol="black",plty=1,plwd=1,p.label.cex=1,p.label.col="black",
                            ...)
  ## Function to plot the profile
{
    # assume the row index is gene and col index is samples, the profile is over the samples
    if(toScale) x <- t(scale(t(x)))
    ncase <- ncol(x)
    maxx <- max(x)
    if(is.null(ylimit)) ylimit <- c(min(x),max(x))
    if (!add) {
      plot(1,type='n',xlim=c(1,ncase*(1+.1)),ylim=ylimit,xaxt="n",...)
      axis(1,at=1:ncase,labels=1:ncase)
    }
    if (plotProfile) {
      to.labels <- dimnames(x)[[1]]
      text(rep(ncase+0.1,nrow(x)),x[,ncase],labels=to.labels,cex=p.label.cex,pos=4,col=p.label.col)
      if (is.null(order)) order <- 1:ncase
      else order <- match(order,dimnames(x)[[2]])
      if (dim(x)[1] > 1) apply(x[,order],1,lines,x=1:ncase,col=pcol,lty=plty,lwd=plwd)
      else lines(1:ncase,x[order],col=pcol,lty=plty,lwd=plwd)
      if(!is.null(v)& !add) {
        for (i in 1:length(v)) {
          abline(v=v[i],type=vtype,lty=vlty,col=vcol,lwd=vlwd)
        }
      }
    }
} 

pfc.profileContrast <- function(p.dim,type=1)
  # Function to generate the contrast matrix
  # This function generate the contrast needed later.
  # p.dim: dimetion for the contrast
  # type: type of contrast, 1 = parallel contrast; 2 = anti-parallel..
{

    if (p.dim <2) stop("p.dim has to be greater than 1")
    
    if (type == 1) {
        output <-   cbind(rep.int(0,p.dim-1),-diag(p.dim-1)) + cbind(diag(p.dim-1),rep.int(0,p.dim-1))
        return(output)
    }
    else if(type == 2) {
      output <-   cbind(rep.int(0,p.dim-1),diag(p.dim-1)) + cbind(diag(p.dim-1),rep.int(0,p.dim-1))
      return(output)
    }
    else stop("Sorry! The contrast you requested has not been implemented yet")
}

pfc.setBranchEdgePar <- function(d,bl=NULL,edgePar=list(col="red",lty="solid",lwd=1))
  {
    #edgePar <- c(attributes(d)$edgePar,edgePar)
    if (is.null(bl))
      {
       d <- dendrapply(d,function(d,edgePar) {attributes(d)$edgePar <- edgePar;d},edgePar=edgePar)
      }
    else
      {
        if (! is.list(bl)) bl <- list(bl)
        branch.n <- length(bl)
        for ( i in 1:branch.n) {
          if(attr(d[[bl[[i]]]],'members') > 1) {
            d[[c(bl[[i]],1)]] <- dendrapply(d[[c(bl[[i]],1)]],
                                            function(d,edgePar) {attributes(d)$edgePar <- edgePar;d},
                                            edgePar=edgePar)                        
            d[[c(bl[[i]],2)]] <- dendrapply(d[[c(bl[[i]],2)]],
                                            function(d,edgePar) {attributes(d)$edgePar <- edgePar;d},
                                            edgePar=edgePar)
          }
        }
      }
    d
  }
    

pfc.setCoherentGroups <- function(d,
                                  d.D,
                                  p.cutoff,
                                  coh.method=c('chisq','f'),
                                  edge.coherent=list(col="black",lty="solid",lwd=1),
                                  edge.noncoherent=list(col="gray50",lty="solid",lwd=.5),
                                  cohIndex.df
                                  )
  ## d: an object of dendrogram
  ## p.cutoff: cutoff value for coherence.
  ### P values that are greater than cutoff will be defined as coherent
  ## edge.cohrent: edgepar for coherent nodes
  ## edge.noncoherent: edgepar for noncoherent nodes
  {
    coh.method <- match.arg(coh.method)
    d <- dendrapply(d,pfc.setEdgePar.helper,edge.noncoherent)
    d <- dendrapply(d,pfc.setMaxD.helper,d.D)
    if(coh.method=='chisq')     d <- dendrapply(d,pfc.setCoherenceIndex.Chisq,cohIndex.df)
    else if (coh.method=='f')   d <- dendrapply(d,pfc.setCoherenceIndex.F,cohIndex.df)
    else stop(paste('coherence method:',coh.method,'is not supported!'))

    if (attr(d,'coherenceIndex') > p.cutoff) {
      cat("-->Note: Root is coherent\n")
      attr(d,'pos') <- 1
      attr(d,"grpid") <- 1
      attr(d,"edgetext") <- 1
      return(d)
    }
    
    i <- 1
    colEdge <- function(X) {
       if(!is.leaf(X)) {
        p.cutoff <- p.cutoff
        if(!is.na(attr(X,"coherenceIndex")) & attr(X,"coherenceIndex") > p.cutoff) {
          attr(X,"grpid") <- i
          attr(X,"edgetext") <- as.character(i)
          i <<- i+1
          for (j in seq(length=length(X))) {
            X[[j]] <- dendrapply(X[[j]],pfc.setEdgePar.helper,edge.coherent)
          }
        }
        else {
          for (jj in seq(length=length(X))) {
            attr(X[[jj]],'pos') <- c(attr(X,'pos'),jj)
            X[[jj]] <- colEdge(X[[jj]])
          }
        }
      }
      X
     }
    d <- colEdge(d)
    d
  }
                            

pfc.setCoherenceIndex.Chisq <- function(d,
                                               cohIndex.df
                                     )
  {
    ##    d.pvalue <- pfc.clustPvalues(d.pvalues,labels(d))
    
    maxDii <- attr(d,'maxD')
    p <- attr(d,'members')
    coherenceIndex <- (1-pchisq(maxDii,cohIndex.df))*p*(p-1)/2
    attr(d,"coherenceIndex") <- if(coherenceIndex > 1) 1 else coherenceIndex
    d
}

pfc.setCoherenceIndex.F <- function(d,
                                          cohIndex.df
                                          )
  {
    ##
    maxDii <- attr(d,'maxD')
    p <- attr(d,'members')
    if(p>cohIndex.df) {
      coherenceIndex <- (1-pf(maxDii*(p-cohIndex.df)/(cohIndex.df*(p-1)),
                              cohIndex.df,p-cohIndex.df))*p*(p-1)/2
      attr(d,"coherenceIndex") <- if(coherenceIndex > 1) 1 else coherenceIndex
    } else {
      attr(d,"coherenceIndex") <- NA
    }
    d
}

pfc.setMaxD.helper <- function(d,
                                     d.D)
  {
    idx <- match(labels(d),row.names(d.D))
    attr(d,"maxD") <- max(d.D[idx,idx])
    d
  }

pfc.setEdgePar.helper <- function(d,
                                     edgePar=list(col="gray",lty="solid",lwd=1)
                                     )
  ## Warning: edgePar will replace the attr(d,"edgepar")
  {
    attributes(d)$edgePar <- edgePar
    d
  }

pfc.getDendrAttr <- function(d,attributes=c("height"))
  ## Recursively get the values of a attribute.
  {
    val <- NULL
    if(is.leaf(d)) {
      return(val)
    }
    else {
      for (j in seq(length=length(d))) {
        val <- rbind(val,
                     sapply(seq(along=attributes),
                            function(x,d,attributes) attr(d[[j]],attributes[x]),
                            d=d,attributes=attributes),
                     pfc.getDendrAttr(d[[j]],attributes=attributes)
                     )
      }
    }
    val
  }

pfc.estM <- function(d,dat,h,debug=F,tol=1e-40,which=c('M','diag'))
  ## estimate the best var that satisfies common variance in the dendrogram
  {
    which <- match.arg(which)
    d.cut <- cut(d,h=h)
    n <- attr(d,'members')
    k <- length(d.cut$lower)
    W <- 0
    Si <- rep(0,k)
    ni <- rep(0,k)
    k.na <- 0
    for(i in seq(k)) {
      assign(paste('S',i,""),var(dat[labels(d.cut$lower[[i]]),]))
      assign(paste('n',i,''),attr(d.cut$lower[[i]],'members'))
      if (get(paste('n',i,""))<=1) k.na <- k.na+1
      else W <- W + (get(paste('n',i,""))-1)*get(paste('S',i,""))
    }
    S <- if(n>k) W/(n-k-k.na) else 0
    avg.diag.S <- mean(diag(S))
    M.1 <- 0
    if(debug) print(k)
    for(i in seq(k)) {
      if(get(paste('n',i,""))>1)
        M.1 <- M.1+(get(paste('n',i,""))-1)*log(diag(get(paste('S',i,""))))
    }
    M <- sum(((n-k-k.na)*log(diag(S))-M.1))
    M
    switch(which,M=M,diag=avg.diag.S)
  }

