#' Plots for lqmix objects
#'
#' Graphically display component and/or transition probabilities of a fitted model of \code{\link{class}} \code{lqmix}
#'
#' @param x an object of class \link{search_lqmix}
#' @param \dots not used
#'
#' @importFrom grDevices devAskNewPage
#' @importFrom graphics axis par legend par matplot barplot
#' @importFrom diagram plotmat
#'
#' @export



plot.lqmix = function(x,...){

  if(x$model == "TC"){

    par(mar = c(4,5,3,3))
    pg = round(x$pg,2)
    names(pg) = as.character(1:length(pg))
    barplot(pg, ylim = c(0,1), ylab = "pg", main = "Component probabilities", cex.main = 1)


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
    par(mar = c(4,5,3,3))
    pg = round(x$pg,2)
    names(pg) = as.character(1:length(pg))
    barplot(pg, ylim = c(0,1), ylab = "pg", main = "Component probabilities", cex.main = 1)

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
