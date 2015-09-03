pfc.plotHclust <- function(PfCluster,
                                 pdfFilename = NULL,
                                 pdf.pointsize = 6,
                                 leaflab ='none',
                                 plotProfile = T,
                                 scaleProfile = T,
                                 scale.center = T,
                                 scale.scale = T,
                                 dendrogram.title = 'Dendrogram',
                                 profile.mfrow=c(3,2),
                                 plotRepresentatives = T,
                                 plotIndivProfile = T,
                                 edgePar.base = list(col = 'gray30',lty='solid',lwd=.5),
                                 edgePar.cohGrps.base =list(col = 'blue4',lty='solid',lwd=1),
                                 edgePar.cohGrps.highlight = list(col='cyan3',lty='solid',lwd=1),
                                 highlight.p.col = 'black',
                                 highlight.t.col = 'white',
                                 rep.line.col='cyan3',
                                 rep.line.lty=1,
                                 rep.line.lwd=2,
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
    if(is.null(pdfFilename)) stop("Please provide the pdf file name to output the figures, eg, myfile.pdf")
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
    pdf(file=pdfFilename,pointsize=pdf.pointsize,height = 8.5, width = 11)
    plot(d,leaflab=leaflab,main=dendrogram.title)
    
    if ( plotProfile & (plotRepresentatives | plotIndivProfile)) {
      par(mfrow = profile.mfrow)
      for (i in gidx ) {
        if (plotIndivProfile){
          pfc.plotProfile(dt2plot[PfCluster$cohGroups$members[[i]],],
                                main = paste("Group:",PfCluster$cohGroups$grpid[i]),
                                ylimit = c(dt.min,dt.max),
                                sub = paste("Coherence Index =",round(PfCluster$cohGroups$cohereIndex[i],digits=3))
                                )
        }
        if (plotRepresentatives) {
          lines(PfCluster$representatives[[i]],col=rep.line.col,lty=rep.line.lty,lwd=rep.line.lwd)
        }
      }
    }
    dev.off()
  }
