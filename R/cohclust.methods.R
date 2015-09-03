members.cohclust <- function(cohclust)
  ## members method for class cohclust
  {
    labels(cohclust$dendrogram)
  }

nmembers.cohclust <- function(cohclust)
  ## nmembers method for class cohclust
  {
    length(labels(cohclust$dendrogram))
  }

plot.cohclust <- function(cohclust,
                          leaflab ='none',
                          dendrogram.title = 'Dendrogram',
                           edgePar.cohGrps.base =list(col = 'blue4',lty='solid',lwd=1),
                          ...
                          )
  ## plot method for class cohclust for dendrogram
  {
    d <- cohclust$dendrogram
    d <- pfc.setBranchEdgePar(d,edgePar=edgePar.cohGrps.base)
    plot(d,leaflab=leaflab,main=dendrogram.title,...)
  }

plotprofile.cohclust <- function(cohclust,
                                 plotRepresentatives = T,
                                 plotIndivProfile = T,
                                 rep.line.col='cyan3',
                                 rep.line.lty=1,
                                 rep.line.lwd=2,
                                 ...
                                 )
  ## plotprofile method for class cohclust
{
  if (plotRepresentatives | plotIndivProfile) {
    if (plotIndivProfile){
      pfc.plotProfile(cohclust$dat,
                      main = paste("Group:",attr(cohclust,'grpid')),
                      sub = paste("Coherence Index =",round(attr(cohclust,'cohIndex'),digits=3)),
                      ...
                      )
    }
    if (plotRepresentatives) {
      lines(cohclust$representative,col=rep.line.col,lty=rep.line.lty,lwd=rep.line.lwd)
    }
  }
}

summary.cohclust <- function(cohclust)
  ## summary method for class cohclust
  {
    cat("A coherent cluster object\n")

    sumry <- as.data.frame(cbind(attr(cohclust,'grpid'),
                                 nmembers(cohclust),
                                 attr(cohclust,'cohIndex')))
    cat("Summary:\n")
    cat(paste(" ","grpid:",sumry[1],sep="\t"))
    cat("\n")
    cat(paste(" ","nmembers:",sumry[2],sep="\t"))
    cat("\n")
    cat(paste(" ","cohIndex:",sumry[3],sep="\t"))
    cat("\n")

    names(sumry) <- c('grpid','nmembers','cohIndex')
    invisible(sumry)
  }
    


