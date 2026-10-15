## Data
df <- read.csv("~/Downloads/IF_v3.1.csv", stringsAsFactors=T)

df$Phase <- factor(df$Phase)
df$ID_old <- df$ID
levels(df$ID) <- c(paste0("IF0",1:9), paste0("IF",10:15))
dat <- df[!is.na(df$LYMP),] #IF-013 (S09) contains one missing value


## Figure 2A: Raw Data
png(filename="Fig.2A*.png", width=2500, height=1800, res=300)

library(lattice) #v 0.23-1

xyplot(LYMP ~ Time | ID, data=df, group=Phase, type="p", layout=c(5,3), cex=1.2, pch=19,
       index.cond=list(c(6,7,8,9,11,
                         10,12,13,14,15,
                         1,2,3,4,5)), #change panel order
       scales=list(x=list(labels=c(-15,-5,0,5) * 30, at=c(-15,-5,0,5), cex=1.2),
                   y=list(at=c(20,50,80,110), cex=1.2)),
       xlab=list(label="Day Increments", cex=1.2), 
       ylab=list(label=expression(paste("Lymphocyte Count (\u00d7",10^9/L,")")), cex=1.2),
       ylab.right=list(label=" Top: \"3+1\"            Middle: \"6\"          Bottom: \"3+3\"", cex=1.2, rot=270))
dev.off()


## Linear Mixed-Effects (LME) Model
library(nlme) #v 3.1-171

fit <- lme(fixed=LYMP ~ Phase * Time, 
           random=~Time | ID, #random intercept and slope for each subject
           data=df[!is.na(df$LYMP),], method="REML")
summary(fit)


library(dplyr) #v 1.2.1

df <- df %>%
  full_join(data.frame(expand.grid("ID"=levels(df$ID), "Time"=c(-3,-2,-1)), "Phase"="0"),
            by=c("ID","Time","Phase")) #add pre-intervention time points for predictions

df$pred <- predict(fit, newdata=df) #add predictions from the LME model
pooled_mean <- mean(subset(df, Time==0)$pred) #mean predicted baseline

df <- df %>%
  group_by(ID) %>%
  mutate(Baseline_pred = pred[Time==0]) %>%
  mutate(pred_centered = (pred - Baseline_pred + pooled_mean) / pooled_mean) %>% #create normalized predictions
  mutate(Seg = factor(if_else(Time==0, "1", Phase))) %>% #create segments for plotting
  select(-Baseline_pred)

last_LYMP <- df %>%
  group_by(ID) %>%
  slice_max(order_by=Time) #labeling use in figure, geom_text_repel()

new_rows <- df %>%
  group_by(ID) %>%
  filter(Phase==1) %>%
  slice(n()) %>% #duplicate each subject's last entry at the end of IF
  mutate(Seg="2") %>% #modify the segment factor
  ungroup()
new_rows2 <- df %>%
  group_by(ID) %>%
  filter(Time==0) %>%
  slice(n()) %>% #duplicate each subject's first entry at the start of IF
  mutate(Seg="0") %>% #modify the segment factor
  ungroup()
df <- rbind(df, new_rows, new_rows2)


## Figure 2B: LME Model Fit
library(ggplot2) #v 4.0.3 (many deprecated theme controls)
library(ggrepel) #v 0.9.8

fig_3Phases <- ggplot(subset(df, Time >= -3),
                      aes(x=Time, y=pred_centered, group=interaction(ID, Seg), color=Seg, linetype=Seg)) +
  geom_line() +
  labs(x="\nDay Increments", y="Normalized Lymphocyte Prediction\n",
       color="", linetype="") +
  scale_x_continuous(limit=c(-3.5,6.5), breaks=-3:6, expand=c(0,0), labels=c(-3:6)*30) +
  scale_color_manual(values=c("#3071ad","#d79f36","#469b76"),
                     labels=c("0"="Baseline", "1"="TRE", "2"="Post-TRE")) +
  scale_linetype_manual(values=c("longdash","solid","dashed"),
                        labels=c("0"="Baseline", "1"="TRE", "2"="Post-TRE")) +
  geom_text_repel(data=last_LYMP, aes(x=Time, y=pred_centered, label=ID),
                  vjust=0.8, hjust=-0.4, segment.color=NA, size=4, col="red") +
  geom_text_repel(data=subset(df, Time==-3), aes(x=Time, y=pred_centered, label=ID),
                  vjust=0.5, hjust=0.4, segment.color=NA, size=4, col="red") +
  geom_hline(yintercept=1, linetype="dotted", color="gray60") +
  theme_classic() +
  theme(legend.position=c(.12, .85), 
        plot.margin=margin(t=25, r=10, l=10),
        plot.title=element_text(size=20),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        legend.title=element_text(size=15),
        legend.text=element_text(size=12))

ggsave("~/Documents/IF/Fig.2B*.png", fig_3Phases, width=2500, height=1200, units="px", device="png", dpi=300)


## Figure 2, Panels A & B

library(magick) #v 2.9.1, macOS needs to install XQuartz, https://www.xquartz.org/

fig2a <- image_read("~/Documents/IF/Fig.2A*.png")
fig2b <- image_read("~/Documents/IF/Fig.2B*.png")
fig2ab <- image_append(c(fig2a, fig2b), stack=T)
fig2ab <- image_annotate(fig2ab, "A", size=80, weight=650)
fig2ab <- image_annotate(fig2ab, "B", size=80, weight=650, location="+0+1700")
image_write(fig2ab, path="~/Documents/IF/Figure 2AB.png", format="png")


library(mgcv) #v 1.9-4
library(gratia) #v 0.11.2
library(plotrix) #v 3.8-14

## Figure 2C: Generalized Additive Models (GAMs)

viz <- function(subj, month_vec=c(0, 1), no_labels=F, panel_text=NULL,
                ylim=NULL, tbl_pos=c(-13, NULL), tbl_pad=1.5,
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
  
  sub <- subset(dat, ID_old==subj)
  phases <- c(tail(subset(sub, Phase==0)[,"Time"],1), #start of fasting
              tail(subset(sub, Phase==1)[,"Time"],1)) #start of posttest
  
  
  fit <- gam(LYMP ~ s(Time, bs="cr", k=6), data=sub, method="REML") #fit a GAM
  pred <- predict(fit, type="response", se.fit=T) #get standard errors
  sub$fit <- pred$fit
  sub$se <- pred$se.fit
  sub$lwr <- sub$fit - 1.96 * sub$se
  sub$upr <- sub$fit + 1.96 * sub$se #compute 95% CIs
  
  der1 <- derivatives(fit, data=data.frame("Time"=month_vec)) #compute the first derivative
  der2 <- derivatives(fit, order=2, data=data.frame("Time"=month_vec)) #compute the second derivative
  tbl <- data.frame("Est."=c(der1$.derivative, der2$.derivative),
                    "Lwr."=c(der1$.lower_ci, der2$.lower_ci),
                    "Upr."=c(der1$.upper_ci, der2$.upper_ci))
  rownames(tbl) <- paste0(rep(c("1st@","2nd@"), each=length(month_vec)), month_vec)
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
       main=paste0(sub$ID[1], " (edf ", round(summary(fit)$edf, 1), ")"),
       xlab=lab.x, ylab=lab.y, cex.main=title_size, cex.lab=title_size, cex.axis=label_size)
  abline(v=phases, col="gray", lty=2, lwd=line_width*.7)
  lines(xvals, spline, col="purple", lwd=line_width)
  for(i in 1:nrow(sub)) { #add pointwise CIs
    lines(c(sub$Time[i], sub$Time[i]), c(sub$lwr[i], sub$upr[i]), col="blue", lwd=line_width)
  }
  text(x=mean(phases), y=ylim[2]-10, labels="TRE", col="red", cex=label_size)
  addtable2plot(x=tbl_pos[1], y=tbl_pos[2], table=tbl[1:2,], cex=tbl_size,
                bty="n", ypad=tbl_pad, display.rownames=T)
  mtext(panel_text, side=3, line=1, adj=0, font=2)
}


viz2 <- function(subj, tbl_pos) {
  viz(subj, ylim=c(0,140), no_labels=T, tbl_pos=tbl_pos,
      line_width=2.5, point_size=2, tbl_size=2, title_size=2, label_size=2)
}

 {png(filename="Fig.2C1*.png", width=4000, height=3750, res=300)
  par(mfrow=c(3,3), mar=c(2.5,2.5,2.5,2.5), oma=c(3,6,2,1))
  viz2("IF-009", tbl_pos=c(-14, 80))
  viz2("IF-017", tbl_pos=c(-14, 5))
  viz2("IF-002", tbl_pos=c(-14, 5))
  viz2("IF-010", tbl_pos=c(-14, 80))
  viz2("IF-019", tbl_pos=c(-14, 80))
  viz2("IF-004", tbl_pos=c(-14, 80))
  viz2("IF-012", tbl_pos=c(-14, 80))
  viz2("IF-020", tbl_pos=c(-14, 80))
  viz2("IF-005", tbl_pos=c(-14, 80))
  mtext("Month", side=1, line=2, outer=T, cex=2)
  mtext(expression(paste("Lymphocyte Count (\u00d7",10^9/L,")")),
        side=2, line=1, outer=T, las=0, cex=2)
  dev.off()}

 {png(filename="Fig.2C2*.png", width=4000, height=2500, res=300)
  par(mfrow=c(2,3), mar=c(2.5,2.5,2.5,2.5), oma=c(3,6,2,1))
  viz2("IF-013", tbl_pos=c(-10, 5))
  viz2("IF-021", tbl_pos=c(-14, 5))
  viz2("IF-006", tbl_pos=c(-14, 80))
  viz2("IF-018", tbl_pos=c(-14, 5))
  viz2("IF-022", tbl_pos=c(-14, 80))
  viz2("IF-007", tbl_pos=c(-14, 80))
  mtext("Month", side=1, line=2, outer=T, cex=2)
  mtext(expression(paste("Lymphocyte Count (\u00d7",10^9/L,")")),
        side=2, line=1, outer=T, las=0, cex=2)
  dev.off()}





dat_GAM <- data.frame()
tbl_DER <- data.frame()
for (subj in unique(dat$ID)) {
  sub <- subset(dat, ID==subj)
  
  fit <- gam(LYMP ~ s(Time, bs="cr", k=6), data=sub, method="REML") #fit a GAM
  pred <- predict(fit, type="response", se.fit=T) #get standard errors
  sub$fit <- pred$fit
  sub$se <- pred$se.fit
  sub$lwr <- sub$fit - 1.96 * sub$se
  sub$upr <- sub$fit + 1.96 * sub$se #compute 95% CIs
  dat_GAM <- rbind(dat_GAM, sub)
  
  der1 <- derivatives(fit, data=data.frame("Time"=c(0, 1))) #compute the first derivative
  der2 <- derivatives(fit, order=2, data=data.frame("Time"=c(0, 1))) #compute the second derivative
  tbl <- data.frame("ID"=subj,
                    "Time"=c(0, 1, 0, 1),
                    "Der."=c(1, 1, 2, 2),
                    "Est."=c(der1$derivative, der2$derivative),
                    "Lwr."=c(der1$lower, der2$lower),
                    "Upr."=c(der1$upper, der2$upper))
  tbl_DER <- rbind(tbl_DER, tbl)
}
addWorksheet(wb, "Figure 2C - GAMs")
writeData(wb, 4, dat_GAM) # fourth tab
addWorksheet(wb, "Figure 2C - Derivatives")
writeData(wb, 5, tbl_DER) # fifth tab

saveWorkbook(wb, "~/Documents/IF/singleXLS.xlsx", overwrite=T) # save the workbook to an Excel file (52 KB)

