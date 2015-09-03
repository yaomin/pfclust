plot.Mvh <-
  function(M,
           base.col='blue',
           base.pch=16,
           base.cex=2,
           xlab='h',
           ylab='M',
           highlight.pch=4,
           highlight.cex=4,
           highlight.col='red')
  ## plot M vs h for diagnostics
  {
    plot(M$h,M$M,
         pch=base.pch,
         xlab=xlab,
         ylab=ylab,
         cex=base.cex,
         col=base.col)
    points(M$d.inf[1],
           M$d.inf[2],
           pch=highlight.pch,
           cex=highlight.cex,
           col=highlight.col)
  }
