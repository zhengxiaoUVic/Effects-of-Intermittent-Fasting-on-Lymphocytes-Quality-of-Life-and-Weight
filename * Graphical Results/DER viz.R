# Error in `.rowNamesDF<-`(x, value = value) : invalid 'row.names' length
# In addition: Warning messages:
# 1: Unknown or uninitialised column: `derivative`. 
# 2: Unknown or uninitialised column: `derivative`. 
# 3: Unknown or uninitialised column: `lower`. 
# 4: Unknown or uninitialised column: `lower`. 
# 5: Unknown or uninitialised column: `upper`. 
# 6: Unknown or uninitialised column: `upper`. 
devtools::install_version("gratia", version="0.8.1")

library(mgcv)
library(gratia) #Errors in v 0.9.0
library(plotrix)

df <- read.csv("~/Downloads/IF_v3.1.csv", stringsAsFactors=T)
df$ID_new <- df$ID
levels(df$ID_new) <- c(paste0("S0",1:9), paste0("S",10:15))
dat <- df[!is.na(df$LYMP),] #IF-013 (S09) contains one missing value

viz <- function(subj, month_vec=c(0, 1), no_labels=F, panel_text=NULL,
                ylim=NULL, tbl_pos=c(-13, NULL), tbl_pad=2,
                line_width=2, point_size=1.2, tbl_size=1, title_size=1.2, label_size=1.2) {
  #' Input -
  #' subj:       a character specifying the target subject
  #' month_vec:  a numeric vector of months at which you want the derivatives to be computed
  #' no_labels:  whether not to display the axis labels (for the common labels)
  #' panel_text: a character for adding the panel letter
  #' ylim:       a numeric vector of length 2, giving y coordinates ranges
  #' tbl_pos:    a numeric vector of length 2, giving coordinates to locate the table (no input no output)
  #' tbl_pad:    a numeric value of padding around table text
  #' ...
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
  } else {
    lab.x <- lab.y <- ""
  }
  if(is.null(ylim)) {
    ylim <- range(sub$LYMP)+c(-10,10)
  }
  plot(sub$Time, sub$LYMP, ylim=ylim, xlim=c(-15,6), cex=point_size, lwd=point_size*.8,
       main=paste0(sub$ID_new[1], ": GAM Fit (edf ", round(summary(fit)$edf, 1), ")"),
       xlab=lab.x, ylab=lab.y, cex.main=title_size, cex.axis=label_size*.6)
  abline(v=phases, col="gray", lty=2, lwd=line_width*.7)
  lines(xvals, spline, col="purple", lwd=line_width)
  for(i in 1:nrow(sub)) { #add pointwise CIs
    lines(c(sub$Time[i], sub$Time[i]), c(sub$lwr[i], sub$upr[i]), col="blue", lwd=line_width)
  }
  text(x=mean(phases), y=ylim[2]-10, labels="IF", col="red", cex=label_size)
  addtable2plot(x=tbl_pos[1], y=tbl_pos[2], table=tbl, cex=tbl_size,
                bty="n", ypad=tbl_pad, display.rownames=T)
  mtext(panel_text, side=3, line=1, adj=0, font=2)
}


viz2 <- function(subj, tbl_pos) {
  viz(subj, ylim=c(0,140), no_labels=T, tbl_pos=tbl_pos,
      line_width=2, point_size=1.5, tbl_size=.9, title_size=1.5, label_size=1.5)
}

{png(filename="Fig.1B.png", width=3270, height=5085, res=300)
par(mfrow=c(5,3), mar=c(2,2,2,2), oma=c(3,6,2,1))
viz2("IF-009", tbl_pos=c(-13, 80))
viz2("IF-017", tbl_pos=c(-13, 10))
viz2("IF-002", tbl_pos=c(-13, 80))
viz2("IF-010", tbl_pos=c(-13, 80))
viz2("IF-019", tbl_pos=c(-13, 80))
viz2("IF-004", tbl_pos=c(-13, 80))
viz2("IF-012", tbl_pos=c(-13, 80))
viz2("IF-020", tbl_pos=c(-13, 80))
viz2("IF-005", tbl_pos=c(-13, 80))
viz2("IF-013", tbl_pos=c(-13, 80))
viz2("IF-021", tbl_pos=c(-13, 80))
viz2("IF-006", tbl_pos=c(-13, 80))
viz2("IF-018", tbl_pos=c(-13, 80))
viz2("IF-022", tbl_pos=c(-13, 80))
viz2("IF-007", tbl_pos=c(-13, 80))
mtext("Month", side=1, line=2, outer=T, cex=1.2)
mtext(expression(paste("Lymphocyte Count (\u00d7",10^9/L,")")),
      side=2, line=1, outer=T, las=0, cex=1.2)
dev.off()}

