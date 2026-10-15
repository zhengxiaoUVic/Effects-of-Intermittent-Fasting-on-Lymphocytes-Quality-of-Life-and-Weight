library(mgcv)
library(gratia)
library(plotrix)

df <- read.csv("~/Downloads/IF_v3.1.csv", stringsAsFactors=T)
dat <- df[!is.na(df$LYMP),] #IF-013 contains one missing value

  ## IF-017 and 019 experienced a steady decrease.
  ## IF-013 and 020 have slowed lymphocyte increase.
  

viz <- function(subj, month_vec=c(0, 1), no_labels=F, panel_text=NULL,
                ylim=NULL, tbl_pos=c(-13, NULL), tbl_pad=2) {
  #' Input -
  #' subj:       a character specifying the target subject
  #' month_vec:  a numeric vector of months at which you want the derivatives to be computed
  #' no_labels:  whether not to display the axis labels (for the common labels)
  #' panel_text: a character for adding the panel letter
  #' ylim:       a numeric vector of length 2, giving y coordinates ranges
  #' tbl_pos:    a numeric vector of length 2, giving coordinates to locate the table (no input no output)
  #' tbl_pad:    a numeric value of padding around table text
  #' Output -    plot a generalized additive model with CIs and derivatives
  
  sub <- subset(dat, ID==subj)
  months_to_der <- month_vec
  phases <- c(tail(subset(sub, Phase==0)[,"Time"],1), #start of fasting
              tail(subset(sub, Phase==1)[,"Time"],1)) #start of posttest
  
  
  fit <- gam(LYMP ~ s(Time, bs="cr", k=6), data=sub, method="REML") #fit a GAM
  pred <- predict(fit, type="response", se.fit=T) #get standard errors
  sub$fit <- pred$fit
  sub$se <- pred$se.fit
  sub$lwr <- sub$fit - 1.96 * sub$se
  sub$upr <- sub$fit + 1.96 * sub$se #compute 95% CIs
  
  der1 <- derivatives(fit, data=data.frame("Time"=months_to_der)) #compute the first derivative
  der2 <- derivatives(fit, order=2, data=data.frame("Time"=months_to_der)) #compute the second derivative
  tbl <- data.frame("Est."=c(der1$derivative, der2$derivative),
                    "Lwr."=c(der1$lower, der2$lower),
                    "Upr."=c(der1$upper, der2$upper))
  rownames(tbl) <- paste0(rep(c("1st@","2nd@"), each=length(months_to_der)), months_to_der)
  tbl <- data.frame(apply(tbl, 2, function(x) formatC(x, format="f", digits=1)))
  xvals <- seq(min(sub$Time), max(sub$Time), length.out=100)
  spline <- predict(fit, newdata=data.frame("Time"=xvals), type="response")
  
  
  if(!no_labels) {
    lab.x <- "Month"
    lab.y <- expression(paste("Lymphocyte Count (\u00d7",10^9/L,")"))
    par(mar=c(5,5,2,2))
  } else {
    lab.x <- lab.y <- ""
  }
  if(is.null(ylim)) {
    ylim <- range(sub$LYMP)+c(-10,10)
  }
  plot(sub$Time, sub$LYMP, ylim=ylim, xlim=c(-15,6),
       main=paste0(subj, ": GAM Fit (edf ", round(summary(fit)$edf, 1), ")"), cex=1.2,
       xlab=lab.x, ylab=lab.y)
  for(i in 1:nrow(sub)) { #add pointwise CIs
    lines(c(sub$Time[i], sub$Time[i]), c(sub$lwr[i], sub$upr[i]), col="blue", lwd=2)
  }
  lines(xvals, spline, col="purple", lwd=2)
  abline(v=phases, col="gray", lty=2)
  text(x=mean(phases), y=100, labels="IF", col="red", cex=1.2)
  addtable2plot(x=tbl_pos[1], y=tbl_pos[2], table=tbl,
                bty="n", ypad=tbl_pad, display.rownames=T)
  mtext(panel_text, side=3, line=1, adj=0, font=2)
}


png(filename="Fig.1B.png", width=750*4, height=700*4, res=300)
par(mfrow=c(2,2), oma=c(2,2,0,0))
viz("IF-013", ylim=c(30,110), tbl_pos=c(-13, 70), no_labels=T)
viz("IF-020", ylim=c(30,110), tbl_pos=c(-13, 70), no_labels=T)
viz("IF-017", ylim=c(30,110), tbl_pos=c(-13, 35), no_labels=T)
viz("IF-019", ylim=c(30,110), tbl_pos=c(-13, 70), no_labels=T)
mtext("Month", side=1, line=0, outer=T, cex=1.2)
mtext(expression(paste("Lymphocyte Count (\u00d7",10^9/L,")")),
      side=2, line=0, outer=T, las=0, cex=1.2)
dev.off()

