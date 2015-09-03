pfclust <- function(dt,
                      dist.type=2,
                      a = 1,
                      b=0,
                      linkage.method='ave',
                      cohIndex.df = dim(dt)[2]-1,
                      cohIndex.cutoff = 0.05,
                      coh.method="chisq",
                      scale = F,
                      scale.center = T,
                      scale.scale = T,
                      smooth.method = 'spline',
                      smooth.lw.f = 1/3,
                      smooth.sp.spar = NULL,
                      x.vcov=NULL,
                      ...
                      )
  ## Hierarchical PfCluster, which uses our dist measure
  {
    dt.org <-  dt
    if (scale) {
      dt = t(scale(t(dt.org),center = scale.center, scale = scale.scale))
    }

    PfCluster.dist <- profile.dist(dt,diss.type=dist.type,a=a,b=b,x.vcov=x.vcov)
    PfCluster.hclust <- hclust(as.dist(PfCluster.dist$diss),method=linkage.method,...)
    
    PfCluster.dendrogram <- as.dendrogram(PfCluster.hclust)

    PfCluster.dendrogram <- pfc.setCoherentGroups(PfCluster.dendrogram,
                                                  PfCluster.dist$D.ii,
                                                  p.cutoff = cohIndex.cutoff,
                                                  coh.method = coh.method,
                                                  cohIndex.df=cohIndex.df)
    PfCluster.coherentGrps <- pfc.getCoherentGroups(PfCluster.dendrogram)
    PfCluster.hcluster <- list(dt.org = dt.org,
                               dt = dt,
                               dist = PfCluster.dist,
                               hclust = PfCluster.hclust,
                               dendrogram = PfCluster.dendrogram,
                               cohGroups = PfCluster.coherentGrps
                               )
    PfCluster.representatives <- pfc.grpRepresentatives(PfCluster.hcluster,
                                                              method=smooth.method,
                                                              sp.spar=smooth.sp.spar,
                                                              lw.f=smooth.lw.f)
    PfCluster.hcluster <- c(PfCluster.hcluster,representatives=list(PfCluster.representatives))
    class(PfCluster.hcluster) <- 'pfclust'
    PfCluster.hcluster
  }
                                                        
