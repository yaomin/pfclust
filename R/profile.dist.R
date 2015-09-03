
profile.dist <- function(x,
                         contr,
                         contr.type=1,
                         FUN='pfc.pairDiss',
                         ginv=F,
                         diss.type=1,
                         a=1,
                         b=0,
                         x.vcov=NULL
                         )
  ## Function to calculate the distance for each pair of column of a matrix
  ## assume nrow(x) is the genes, and ncol(x) is the samples.Each row define a profile
  ## We cluster the genes using the profiles defined by the rows
{
    n.gene <- nrow(x)
    n.cond <- ncol(x)
    if(missing(contr)) contr <- pfc.profileContrast(n.cond,
                                                          type=contr.type)
    b.ii <- contr%*%pfc.pairwise.FUN(as.matrix(x),
                                           self.included=F,
                                           FUN=FUN,
                                           a=a,
                                           b=b
                                           )
    b.bar <- rowMeans(b.ii)
    n.b <- ncol(b.ii)
    

    if(is.null(x.vcov))      x.vcov <- var(x)

    Sigma.b <- 2*contr%*%x.vcov%*%t(contr)

    
    if(ginv) {
      Sigma.b.inv <- ginv(Sigma.b)
    }
    else {
      Sigma.b.inv <- solve(Sigma.b)
    }
    
    bb.ii <- apply(b.ii,2,function(x) t(x)%*%x)
    D.ii <- apply(b.ii,2,function(x,Sigma.inv) t(x)%*%Sigma.inv%*%x,Sigma.inv=Sigma.b.inv)

    dist <- switch(diss.type,'1'=bb.ii,'2'=D.ii)
    diss2 <- matrix(0,nrow=n.gene,ncol=n.gene)
    diss2[lower.tri(diss2)] <- dist
    diss2[upper.tri(diss2)] <- t(diss2)[upper.tri(diss2)]

    diss <- sqrt(diss2)
    row.names(diss) <- row.names(x)

    diss2 <- matrix(0,nrow=n.gene,ncol=n.gene)
    diss2[lower.tri(diss2)] <- D.ii
    diss2[upper.tri(diss2)] <- t(diss2)[upper.tri(diss2)]
    D.ii <- diss2
    row.names(D.ii) <- row.names(x)
    
    
    list(diss=diss,bvcov=Sigma.b,b.ii=b.ii,bb.ii=bb.ii,D.ii=D.ii)
 
}

### auxiliary functions
pfc.pairDiss <- function(x1,
                         x2,
                         a=1,
                         b=0
                         )
  ## Calculate the distance between two profiles.
  ## Definition of distance: d(t)=y(t)-a*z*(t+b) ( currently implemented with b=0, no shift.)
  ## Parameters: b = 0 by default, which means no shift between two profiles,
  ##             the shift is defined on the final profile for comparasion
  ##             a = 1 by default, which means parallel profile
  ##             a = -1 means anti-reflective profile
  ##             arbitrary a between -1 and 1 and
  ##               b = 0 are for fan-shaped profile relation between two profiles
  ##             diss.type = 1 or 2, defines how to calculate the final meassure
  ##               based on the profile distance vector
  ##               1 for sum(b^2)(or t(b)%*%b) and 2 for  t(b)%*%solve(b.vcovt)*b
{
  d <- x1-a*x2
  return(d)    
}

    
