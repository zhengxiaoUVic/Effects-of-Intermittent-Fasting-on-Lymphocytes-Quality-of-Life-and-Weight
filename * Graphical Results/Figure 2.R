## Data
df <- read.csv("~/Downloads/IF_v3.1.csv", stringsAsFactors=T)

library(openxlsx) # v 4.2.7
wb <- createWorkbook() # create a new workbook
addWorksheet(wb, "Raw Data")
writeData(wb, 1, df) # add the data frame to a separate sheet in the workbook

df$Phase <- factor(df$Phase)
df$ID_old <- df$ID
levels(df$ID) <- c(paste0("S0",1:9), paste0("S",10:15))
dat <- df[!is.na(df$LYMP),] #IF-013 (S09) contains one missing value

addWorksheet(wb, "Figure 2A")
writeData(wb, 2, dat) # second tab


## Figure 2A: Raw Data
png(filename="Fig.2A.png", width=2800, height=3485, res=300)

library(lattice) #v 0.22-6

xyplot(LYMP ~ Time | ID, data=df, group=Phase, type="p", layout=c(3,5), cex=1.5, pch=19,
       index.cond=list(c(11,15,5,
                         9,14,4,
                         8,13,3,
                         7,12,2,
                         6,10,1)), #change panel order
       xlab=list(label="Month", cex=1.2), 
       ylab=list(label=expression(paste("Lymphocyte Count (\u00d7",10^9/L,")")), cex=1.2),
       sub=list(label="   Left: \"3+3\"                                  Center: \"6\"                                  Right: \"3+1\"", cex=1.2))
dev.off()


## Linear Mixed-Effects (LME) Model
library(nlme) #v 3.1-166

fit <- lme(fixed=LYMP ~ Phase * Time, 
           random=~Time | ID, #random intercept and slope for each subject
           data=df[!is.na(df$LYMP),], method="REML")
summary(fit)


library(dplyr) #v 1.1.4

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

addWorksheet(wb, "Figure 2B")
writeData(wb, 3, df) # third tab


## Figure 2B: LME Model Fit
library(ggplot2) #v 3.5.1 (many deprecated theme controls)
library(ggrepel) #v 0.9.5

fig_3Phases <- ggplot(subset(df, Time >= -3),
                      aes(x=Time, y=pred_centered, group=interaction(ID, Seg), color=Seg, linetype=Seg)) +
  geom_line() +
  labs(x="\nMonth", y="Normalized Lymphocyte Prediction\n",
       title="Linear Mixed-Effects Model", color="Treatment", linetype="Treatment") +
  scale_x_continuous(limit=c(-3.5,6.5), breaks=-3:6, expand=c(0,0)) +
  scale_color_manual(values=c("#3071ad","#d79f36","#469b76"),
                     labels=c("0"="Pre-Intervention", "1"="Intermittent Fasting", "2"="Follow-Up")) +
  scale_linetype_manual(values=c("longdash","solid","dashed"),
                        labels=c("0"="Pre-Intervention", "1"="Intermittent Fasting", "2"="Follow-Up")) +
  geom_text_repel(data=last_LYMP, aes(x=Time, y=pred_centered, label=ID),
                  vjust=0.8, hjust=-0.4, segment.color=NA, size=3, col="red") +
  geom_text_repel(data=subset(df, Time==-3), aes(x=Time, y=pred_centered, label=ID),
                  vjust=0.5, hjust=0.4, segment.color=NA, size=3, col="red") +
  geom_hline(yintercept=1, linetype="dashed", color="gray60") +
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

ggsave("~/Documents/IF/Fig.2B.png", fig_3Phases, width=2800, height=1600, units="px", device="png", dpi=300)


library(mgcv) #v 1.9-1
# devtools::install_version("gratia", version="0.8.1")
library(gratia) #Errors in v 0.9.0+
library(plotrix) #v 3.8-4

## Figure 2C: Generalized Additive Models (GAMs)

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
  tbl <- data.frame("Est."=c(der1$derivative, der2$derivative),
                    "Lwr."=c(der1$lower, der2$lower),
                    "Upr."=c(der1$upper, der2$upper))
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
       main=paste0(sub$ID[1], ": GAM Fit (edf ", round(summary(fit)$edf, 1), ")"),
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

{png(filename="Fig.2C.png", width=3270, height=5085, res=300)
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


## Figure 2

library(magick) #v 2.8.4, macOS needs to install XQuartz, https://www.xquartz.org/

fig2a <- image_read("~/Documents/IF/Fig.2A.png")
fig2b <- image_read("~/Documents/IF/Fig.2B.png")
fig2c <- image_read("~/Documents/IF/Fig.2C.png")

fig2ab <- image_append(c(fig2a, fig2b), stack=T)
fig2ab <- image_annotate(fig2ab, "A", size=100, weight=700)
fig2ab <- image_annotate(fig2ab, "B", size=100, weight=700, location="+0+3550")
fig2 <- image_append(c(fig2ab, fig2c))
fig2 <- image_annotate(fig2, "C", size=100, weight=700, location="+2900")
image_write(fig2, path="~/Documents/IF/Figure 2.png", format="png")


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

