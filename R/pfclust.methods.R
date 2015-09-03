plot.pfclust <- function(PfCluster,
#                           pdfFilename = NULL,
#                           pdf.pointsize = 6,
                           leaflab ='none',
#                           plotProfile = T,
                           scaleProfile = T,
                           scale.center = T,
                           scale.scale = T,
                           dendrogram.title = 'Dendrogram',
#                           profile.mfrow=c(3,2),
#                           plotRepresentatives = T,
#                           plotIndivProfile = T,
                           edgePar.base = list(col = 'gray30',lty='solid',lwd=.5),
                           edgePar.cohGrps.base =list(col = 'blue4',lty='solid',lwd=1),
                           edgePar.cohGrps.highlight = list(col='cyan3',lty='solid',lwd=1),
                           highlight.p.col = 'black',
                           highlight.t.col = 'white',
#                           rep.line.col='cyan3',
#                           rep.line.lty=1,
#                           rep.line.lwd=2,
                           grpid = NULL,
                           ...
                           )
  ### plot PfCluster.hclust
  {
    highlight <- F
    if(scaleProfile) {
      dt2plot <- t(scale(t(PfCluster$dt.org),center = scale.center, scale = scale.scale))
    } else {
      dt2plot <- PfCluster$dt.org
      plotRepresentatives = F
    }
    dt.max <- max(unlist(dt2plot))
    dt.min <- min(unlist(dt2plot))

    if (is.null(grpid)) gidx <- 1:length(PfCluster$cohGroups$grpid)
    else {
      gidx <- match(grpid,PfCluster$cohGroups$grpid)
      highlight <- T
    }
    d <- PfCluster$dendrogram
    d <- pfc.setBranchEdgePar(d,edgePar=edgePar.base)
    d <- pfc.setBranchEdgePar(d,PfCluster$cohGroups$grpPos,edgePar=edgePar.cohGrps.base)
    if (highlight) {
      for (i in gidx ) {
        d <- pfc.setBranchEdgePar(d,PfCluster$cohGroups$grpPos[[i]],edgePar=edgePar.cohGrps.highlight)
        attributes(d[[PfCluster$cohGroups$grpPos[[i]]]])$edgePar$p.col <- highlight.p.col
        attributes(d[[PfCluster$cohGroups$grpPos[[i]]]])$edgePar$t.col <- highlight.t.col
      }
    }
#    pdf(file=pdfFilename,pointsize=pdf.pointsize,height = 8.5, width = 11)
    plot(d,leaflab=leaflab,main=dendrogram.title)
    
##     if ( plotProfile & (plotRepresentatives | plotIndivProfile)) {
##       par(mfrow = profile.mfrow)
##       for (i in gidx ) {
##         if (plotIndivProfile){
##           pfc.plotProfile(dt2plot[PfCluster$cohGroups$members[[i]],],
##                                 main = paste("Group:",PfCluster$cohGroups$grpid[i]),
##                                 ylimit = c(dt.min,dt.max),
##                                 sub = paste("Coherence Index =",round(PfCluster$cohGroups$cohereIndex[i],digits=3))
##                                 )
##         }
##         if (plotRepresentatives) {
##           lines(PfCluster$representatives[[i]],col=rep.line.col,lty=rep.line.lty,lwd=rep.line.lwd)
##         }
##       }
##     }
##     dev.off()
   
  }

summary.pfclust <- function(pfcluster)
  ## summary of pfcluster
{
  grpdt <- pfcluster$cohGroups
  sumry <- as.data.frame(cbind(sapply(grpdt$members,length),
                               grpdt$cohereIndex),
                         row.names=grpdt$grpid)
  names(sumry) <- c("n_members","CohIndex")
  sumry.nr <- nrow(sumry)
  cat(paste("Number of coherent clusters:", length(grpdt$grpid),"\n"))
                                        #cat(paste(" ","ID","n_members","CohIndex",sep="\t"))
  cat("Summary:\n")
  if(sumry.nr>0) {
    ## for(i in seq(sumry.nr)) {
    ##       cat(paste(" ",sumry[i,1],sumry[i,2],sumry[i,3],sep="\t"))
    ##       cat("\n")
    ##     }
    print(sumry)
  }
  else cat("No coherent group found\n")
  invisible(sumry)
}


'[.pfclust' <- function(pfcluster,which)
 ## subsetting pfcluster
 {
   cohGroups <- pfcluster$cohGroups
   which.idx <- match(as.character(which),cohGroups$grpid)
   if(is.na(which.idx)) return(NA)
   dendr <- pfcluster$dendrogram[[cohGroups$grpPos[[which.idx]]]]
   membs <- labels(dendr)
   n.membs <- length(membs)
   cohIndex <- cohGroups$cohereIndex[which.idx]
   represent <- pfcluster$representatives[[which.idx]]
   dat <- pfcluster$dt[membs,]

   cohclust <- list(dat=dat,
                    dendrogram=dendr,
                    representative=represent)
   
   attr(cohclust,'grpid') <- which
   attr(cohclust,'cohIndex') <- cohIndex
   class(cohclust) <- 'cohclust'
   cohclust
 }
     
