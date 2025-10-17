#' Plots for \code{search_lqmix} objects
#'
#' Graphically display model selection criteria and component and/or transition probabilities of the optimal fitted model of \code{\link{class}} \code{\link{search_lqmix}}.
#'
#' @param x an object of class \code{search_lqmix}
#' @param \dots not used
#'
#' @importFrom grDevices devAskNewPage
#' @importFrom graphics axis par legend par matplot barplot
#' @importFrom diagram plotmat
#'
#' @export


plot.search_lqmix = function(x, ...){


    # plot model selection
    par(mar=c(5, 4, 4, 2) + 0.1,mfrow=c(1,1), xpd = TRUE)

  if(is.null(x$mv) & !is.null(x$Gv)){
      model = "TC"
      if(x$method == "bic"){
        bbicv = x$bicv
        matplot(cbind(bbicv), type="b", pch = x$Gv, lty=x$Gv, add = F,col = x$Gv, xaxt="n", xlab="Number of components", ylab = "", main = "Model selection criterion: BIC")
      }else if(x$method == "aic"){
        aaicv = x$aicv
        matplot(cbind(aaicv), type="b", pch = x$Gv, lty=x$Gv, add = F,col = x$Gv, xaxt="n", xlab="Number of components", ylab = "", main = "Model selection criterion: BIC")
      }else{
        llkv = x$lkv
        matplot(cbind(llkv), type="b", pch = x$Gv, lty=x$Gv, add = F,col = x$Gv, xaxt="n", xlab="Number of components", ylab = "", main = "Model selection criterion: logLik")
      }
      axis(side=1,at=1:length(x$Gv), labels = x$Gv)
    } else if(!is.null(x$mv) & is.null(x$Gv)){
    model = "TV"

    if(x$method == "bic"){
      bbicv = x$bicv
      matplot(cbind(bbicv), type="b", pch = x$mv, lty=x$mv, add = F,col = x$mv, xaxt="n", xlab="Number of states", ylab = "", main = "Model selection criterion: BIC")

    }else if(x$method == "aic"){
      aaicv = x$aicv
      matplot(cbind(aaicv), type="b", pch = x$mv, lty=x$mv, add = F,col = x$mv, xaxt="n", xlab="Number of states", ylab = "", main = "Model selection criterion: BIC")
    }else{
      llkv = x$lkv
      matplot(cbind(llkv), type="b", pch = x$mv, lty=x$mv, add = F,col = x$mv, xaxt="n", xlab="Number of states", ylab = "", main = "Model selection criterion: logLik")
    }
    axis(side=1,at=1:length(x$mv), labels = x$mv)

  }else{
    model = "TCTV"
    if(x$method == "bic"){
      bbicv = matrix(x$bicv, length(x$mv), length(x$Gv), byrow = T)
      xx = matrix(rep(1:length(x$mv), each = length(x$Gv)), ncol = length(x$Gv), byrow = T)
      par(mar=c(5, 4, 4, 2) + 0.1,mfrow=c(1,1), xpd = TRUE)

      matplot(xx, cbind(bbicv),
              type="b", pch = x$Gv, lty=x$Gv, add = F,col = x$Gv, xaxt="n", xlab="Number of states", ylab = "",
              main = "Model selection criterion: BIC",
              cex.main = 1.5, cex.axis = 1.2, cex = 1.2)

      legend("topright",horiz = F, legend=c(paste("G=", x$Gv, sep="")),lty=x$Gv,bty="n",
             col = x$Gv, pch = c(x$Gv), cex = 1.2)

    }else if(x$method == "aic"){
      aaicv = matrix(x$aicv, length(x$mv), length(x$Gv), byrow = T)
      xx = matrix(rep(1:length(x$mv), each = length(x$Gv)), ncol = length(x$Gv), byrow = T)
      par(mar=c(5, 4, 4, 2) + 0.1,mfrow=c(1,1), xpd = TRUE)
      matplot(xx, cbind(aaicv), type="b",
              pch = x$Gv, lty=x$Gv, add = F,
              col = x$Gv, xaxt="n",
              xlab="Number of states", ylab = "",
              main = "Model selection criterion: AIC",
              cex.main = 1.5, cex.axis = 1.2, cex = 1.2)

      legend("topright",horiz = F, legend=c(paste("G=", x$Gv, sep="")),lty=x$Gv,bty="n",
             col = x$Gv, pch = c(x$Gv), cex = 1.2)
    }else{
      llkv = matrix(x$lkv, length(x$mv), length(x$Gv), byrow = T)
      xx = matrix(rep(1:length(x$mv), each = length(x$Gv)), ncol = length(x$Gv), byrow = T)
      par(mar=c(5, 4, 4, 2) + 0.1,mfrow=c(1,1), xpd = TRUE)
      matplot(xx, cbind(llkv), type="b",
              pch = x$Gv, lty=x$Gv,add = F,
              col = x$Gv, xaxt="n",
              xlab="Number of states", ylab = "",
              main = "Model selection criterion: logLik",
              cex.main = 1.5, cex.axis = 1.2, cex = 1.2)

      legend("bottomright",horiz = F, legend=c(paste("G=", x$Gv, sep="")),lty=x$Gv,bty="n",
             col = x$Gv, pch = c(x$Gv), cex = 1.2)

    }
    axis(side=1,at=1:length(x$mv), labels = x$mv)

  }


  # plot probs for the latent process
  if(x$optimal$model == "TC"){
    message(cat(paste("Plots for the lqmix optimal model: TC random coefficients with G=", x$optimal$G, " at qtl=", x$optimal$qtl, sep="")))
    devAskNewPage(ask = T)
    par(mar = c(4,5,3,3))
    pg = round(x$optimal$pg,2)
    names(pg) = as.character(1:length(pg))
    barplot(pg, ylim = c(0,1), ylab = "pg",
            main = "Component probabilities",
            cex.main = 1.5)

  }else if(x$optimal$model == "TV"){
    message(cat(paste("Plots for the lqmix optimal model: TV random coefficients with m=", x$optimal$m, " at qtl=", x$optimal$qtl, sep="")))
    devAskNewPage(ask = T)
    par(mar = c(0,2,2,2))
    Gamma = round(t(x$optimal$Gamma),2)
    plotmat(Gamma,relsize=0.9,
              name = as.character(1:x$optimal$m), box.col="lightgray",
              arr.lwd = Gamma*3,
              box.lwd = 1,
              self.cex = 1,
              self.lwd = diag(Gamma)*3,
              cex.txt = 1.2,
              box.size = 0.1,add = F,
              arr.type = "triangle",arr.width = 0.2,
              box.prop = 0.5, cex.main = 1.5, main="Transition probabilities")

  } else{
    message(cat(paste("Plots for the lqmix optimal model: TC and TV random coefficients with m=", x$optimal$m, " and G=", x$optimal$G, " at qtl=", x$optimal$qtl, sep="")))

    devAskNewPage(ask = T)
    par(mar = c(4,5,3,3))
    pg = round(x$optimal$pg,2)
    names(pg) = as.character(1:length(pg))
    barplot(pg, ylim = c(0,1), ylab = "pg",
            main = "Component probabilities",
            cex.main = 1.5)

    devAskNewPage(ask = TRUE)
    par(mar = c(0,2,2,2))
    Gamma = round(t(x$optimal$Gamma),2)
    plotmat(Gamma,relsize=0.9,
            name = as.character(1:x$optimal$m), box.col="lightgray",
            arr.lwd = Gamma*3,
            box.lwd = 1,
            self.cex = 1,
            self.lwd = diag(Gamma)*3,
            cex.txt = 1.2,
            box.size = 0.1,
            arr.type = "triangle",
            arr.width = 0.2,
            box.prop = 0.5,
            cex.main = 1.5,
            cex.axis = 1.2,
            main="Transition probabilities")
  }

}

