library(nlme)

df <- read.csv("~/Downloads/IF_v3.1.csv", stringsAsFactors=T)
lymphocyte <- df[df$Phase==1, -2]
autophagy <- df[!is.na(df$ATG) & df$Phase==1, -1] #Excluded IF-007 and 012 (no records)
fast6 <- paste0("IF-0",c(17,19,20,21,22))
lymphocyte$Fast <- factor(ifelse(lymphocyte$ID %in% fast6, "Six-Month", "Three-Month"))
autophagy$Fast <- factor(ifelse(autophagy$ID %in% fast6, "Six-Month", "Three-Month"))


fit.LYMP <- lme(fixed=LYMP ~ Fast * Time_Within_Phase, #time starts at 0 for better interpretation of the baseline
                random=~Time_Within_Phase | ID, data=lymphocyte, method="ML",
                control=lmeControl(opt="optim", optCtrl=list(maxit=1e5)))

fit.ATG <- lme(fixed=ATG ~ Fast * Time_Within_Phase,
               random=~Time_Within_Phase | ID, data=autophagy, method="ML")


library(ggeffects)
library(ggplot2)

viz2 <- function(fit, lwd=1, pos=c(NULL,NULL,NULL)) {
  #' Input - 
  #' fit:     linear mixed-effects model fit
  #' lwd:     line width
  #' pos:     positioning aesthetics of the model expressions
  #' Output - return a ggplot of the fit
  dat <- ggpredict(fit, terms=c("Time_Within_Phase", "Fast")) #95% confidence intervals
  coef <- round(coef(summary(fit))[,"Value"], 2) #extract fixed effect coefficients
  tex_6 <- paste0("hat(italic(Y)) == ", coef[3], " * italic(x) + ", coef[1]) #fitted base model
  tex_3 <- paste0("hat(italic(Y)) == ", coef[3]+coef[4], " * italic(x) + ", coef[1]+coef[2])
  ggplot(dat, aes(x, predicted)) + 
    geom_line(aes(linetype=group, color=group), linewidth=lwd) +
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
    scale_linetype_manual(values=c("dotted", "solid")) +
    scale_color_manual(values=c("#D74E44", "#5E89BB")) +
    scale_fill_manual(values=c("#D74E44", "#5E89BB")) +
    annotate(geom="text", x=pos[1], y=pos[3], label=tex_6, hjust=0, parse=T, size=5, color="#D74E44") +
    annotate(geom="text", x=pos[1], y=pos[2], label=tex_3, hjust=0, parse=T, size=5, color="#5E89BB") +
    labs(x="Month", color="Fast", fill="Fast", linetype="Fast") +
    theme_classic()
}

p.LYMP <- viz2(fit.LYMP, pos=c(0.5, 100, 95)) +
  labs(y=expression(paste("Lymphocyte Count (\u00d7",10^9/L,")"))) +
  scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand=c(0,0))
p.ATG <- viz2(fit.ATG, pos=c(3, 0.965, 0.88)) +
  labs(y="Relative Net LC3-II Expression") +
  scale_y_continuous(expand=c(0,0)) + coord_cartesian(ylim=c(0,1)) + scale_x_continuous(expand=c(0,0))

library(ggpubr)
fig <- ggarrange(p.LYMP, p.ATG, ncol=2, common.legend=T, legend="bottom", font.label=list(size=2))
fig <- annotate_figure(fig, top=text_grob("Linear Mixed-Effects Regression Models", size=14, face="bold"))

ggsave("~/Desktop/Fig.png", fig, width=2800, height=1600, units="px", device="png", dpi=300)





# png(filename="Fig.png", width=2800, height=1600, res=300)
# par(mfrow=c(1,2), mar=c(4,5,3,1))
# plot(1, xlim=c(0, 5), ylim=c(60, 85), main="", cex.lab=1.2,
#      type="n", xlab="Month", ylab=expression(paste("Lymphocyte Count (\u00d7",10^9/L,")")))
# abline(a=coef.fixed.LYMP[1], b=coef.fixed.LYMP[3], lty=2, lwd=2)
# abline(a=coef.fixed.LYMP[1]+coef.fixed.LYMP[2], b=coef.fixed.LYMP[3]+coef.fixed.LYMP[4], col="red", lwd=2)
# legend("bottomright", c("6-Month", "3-Month"),
#        col=c("black","red"), lty=c(2,1), lwd=2, bty="n")
# text(1,83,bquote(hat(italic(Y)) ==
#                      .(format(coef.fixed.LYMP[1],digits=4,nsmall=2))
#                    + .(round(coef.fixed.LYMP[3],2)) * italic(x)))
# text(1,85,bquote(hat(italic(Y)) == .(round(coef.fixed.LYMP[1]+coef.fixed.LYMP[2],2)) 
#                    + .(round(coef.fixed.LYMP[3]+coef.fixed.LYMP[4],2)) * italic(x)), col="red")
# 
# 
# plot(1, xlim=c(0, 5), ylim=c(0, 0.7), main="", cex.lab=1.2,
#      type="n", xlab="Month", ylab="Autophagy")
# abline(a=coef.fixed.ATG[1], b=coef.fixed.ATG[3], lty=2, lwd=2)
# abline(a=coef.fixed.ATG[1]+coef.fixed.ATG[2], b=coef.fixed.ATG[3]+coef.fixed.ATG[4], col="red", lwd=2)
# legend("bottomleft", c("6-Month", "3-Month"),
#        col=c("black","red"), lty=c(2,1), lwd=2, bty="n")
# text(4,0.65,bquote(hat(italic(Y)) == 
#                      .(round(coef.fixed.ATG[1],2)) 
#                  - .(abs(round(coef.fixed.ATG[3],2))) * italic(x)))
# text(4,0.7,bquote(hat(italic(Y)) == .(round(coef.fixed.ATG[1]+coef.fixed.ATG[2],2)) 
#                  - .(abs(round(coef.fixed.ATG[3]+coef.fixed.ATG[4],2))) * italic(x)), col="red")
# mtext("Linear Mixed-Effects Regression Models",
#       side=3, line=-2, outer=T, cex=1.5, font=2)
# dev.off()

