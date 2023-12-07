#' Plots for lqmix objects
#'
#' Graphically display component and/or transition probabilities of a fitted model of \code{\link{class}} \code{lqmix}
#'
#' @param x an object of class \link{search_lqmix}
#' @param \dots not used
#'
#' @importFrom grDevices devAskNewPage
#' @importFrom graphics axis par legend par matplot
#' @importFrom diagram plotmat
#'
#' @export



plot.lqmix = function(x,...){

  if(x$model == "TC"){
    par(mar = c(0,2,2,2))
    pg = round(diag(x$pg),2)
    plotmat(pg,relsize=0.7,
            name = as.character(1:x$G), box.col="lightgray",
            box.size = diag(pg)*0.5,
            box.lwd = 1,
            self.cex = 0.8,add = F,
            self.lwd = 0,
            cex.txt = 0.8,
            arr.type = "none",arr.width = 0.4,
            box.prop = 0.5,
            main="Component probabilities")

  }else if(x$model == "TV"){
    par(mar = c(0,2,2,2))
    Gamma = round(t(x$Gamma),2)
    plotmat(Gamma,relsize=0.7,
            name = as.character(1:x$m), box.col="lightgray",
            arr.lwd = Gamma*3,
            box.lwd = 1,
            self.cex = 0.8,
            self.lwd = diag(Gamma)*3,
            cex.txt = 0.8,
            box.size = 0.1,add = F,
            arr.type = "triangle",arr.width = 0.4,
            box.prop = 0.5,main="Transition probabilities")

  } else{
    par(mar = c(0,2,2,2))
    pg = round(diag(x$pg),2)
    plotmat(pg,relsize=0.7,
            name = as.character(1:x$G), box.col="lightgray",
            box.size = diag(pg)*0.5,
            box.lwd = 1,
            self.cex = 0.8,
            self.lwd = 0,
            cex.txt = 0.8,
            arr.type = "none",add = F,arr.width = 0.4,
            box.prop = 0.5,
            main="Component probabilities")

    devAskNewPage(ask = T)
    par(mar = c(0,2,2,2))
    Gamma = round(t(x$Gamma),2)
    plotmat(Gamma,relsize=0.7,
            name = as.character(1:x$m), box.col="lightgray",
            arr.lwd = Gamma*3,
            box.lwd = 1,
            self.cex = 0.8,
            self.lwd = diag(Gamma)*3,
            cex.txt = 0.8,
            box.size = 0.1,
            arr.type = "triangle",
            arr.width = 0.4,
            box.prop = 0.5,main="Transition probabilities")
  }

}
