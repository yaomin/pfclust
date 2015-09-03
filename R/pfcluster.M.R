pfcluster.M <- function(d,dat,debug=F,tol=1e-40,h=NULL,n.node=25)
  ## d: dendrogram
  ## dat: data used in d
  ## h: threshold to use for computing M's above it
  ## n.node: minimum node size for calculating M's, an alternative to threshold h.
  {
    d.heights <- pfc.getDendrAttr(d,attributes='height')
    d.heights.order <- order(d.heights,decreasing=T)
    if(is.null(h)) d.heights <- d.heights[d.heights.order[seq(n.node)]]
    else d.heights <- sort(d.heights[d.heights>h])
    d.Ms <- sapply(d.heights,pfc.estM,d=d,dat=dat,debug=debug,tol=tol,which='M')
    avg.diag <- sapply(d.heights,pfc.estM,d=d,dat=dat,debug=debug,tol=tol,which='diag')
    d.inf <- c(d.inf=d.heights[which(d.Ms==min(d.Ms))], M.min=min(d.Ms))
    ret <- list(h=d.heights,M=d.Ms,S.diag=avg.diag,d.inf=d.inf)
    class(ret) <- 'Mvh'
    ret
  }
