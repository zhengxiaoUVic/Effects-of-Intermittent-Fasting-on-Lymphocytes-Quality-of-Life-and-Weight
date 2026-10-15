## Data
df <- read.csv("~/Downloads/IF_v3.1.csv", stringsAsFactors=T)
df$Phase <- factor(df$Phase)
df$ID_old <- df$ID
levels(df$ID) <- c(paste0("S0",1:9), paste0("S",10:15))


## Generalized Additive Models (GAMs)
library(mgcv) #v 1.9-1

N <- 100 #number of time points
nTime <- nPred <- nPhase <- numeric(N); i <- 0
for (subj in unique(df$ID)) {
  sub <- subset(df, ID==subj)
  fit <- gam(LYMP ~ s(Time, bs="cr", k=6), data=sub, method="REML") #fit a GAM
  xvals <- sort(c(seq(min(sub$Time), max(sub$Time), length.out=N-2), 1E-3, 3)) #time points
  nPred[i*N+c(1:N)] <- predict(fit, newdata=data.frame("Time"=xvals), type="response") #add predictions
  nTime[i*N+c(1:N)] <- xvals
  if (tail(sub$Phase,1)==1) { #for those without posttest
    nPhase[i*N+c(1:N)] <- (xvals >= 0)
  } else {
    nPhase[i*N+c(1:N)] <- (xvals >= 3) + (xvals >= 0)
  }
  i <- i + 1
}
dat <- data.frame("ID"=rep(unique(df$ID), each=N),
                  "Time"=nTime,
                  "Phase"=factor(nPhase),
                  "pred"=nPred)
pooled_mean <- mean(subset(dat, Time==1E-3)$pred) #mean baseline

library(dplyr) #v 1.1.4

baseline <- dat %>%
  filter(Time==1E-3) %>%
  select(ID, pred) %>%
  rename(baseline_pred=pred)

dat <- dat %>%
  left_join(baseline, by="ID") %>%
  mutate(pred_centered = (pred - baseline_pred + pooled_mean) / pooled_mean) %>% #create normalized predictions
  select(-baseline_pred)

last_LYMP <- dat %>%
  group_by(ID) %>%
  slice_max(order_by=Time) #labeling use in figure, geom_text_repel()

new_rows <- dat %>%
  filter(Phase==2 & Time==3) %>% #duplicate each subject's last entry at the end of IF == Month 3
  mutate(Phase="1") #modify the phase factor for plotting
dat <- rbind(dat, new_rows)


## Visualization
library(ggplot2) #v 3.5.0 (many deprecated theme controls)
library(ggrepel) #v 0.9.5

fig <- ggplot(subset(dat, Time >= 0),
              aes(x=Time, y=pred_centered, group=interaction(ID, Phase), color=Phase, linetype=Phase)) +
  geom_hline(yintercept=1, linetype="dashed", color="gray80") +
  geom_line() +
  labs(x="Month", y="Normalized Lymphocyte Count",
       title="Generalized Additive Models", color="Treatment", linetype="Treatment") +
  scale_x_continuous(limit=c(0,6.5), breaks=c(0,1,2,3,4,5,6), expand=c(0,0)) +
  scale_color_manual(values=c("black","dodgerblue"),
                     labels=c("1"="Intermittent Fasting", "2"="Follow-Up")) +
  scale_linetype_manual(values=c("solid","dashed"),
                        labels=c("1"="Intermittent Fasting", "2"="Follow-Up")) +
  geom_text_repel(data=last_LYMP, aes(x=Time, y=pred_centered, label=ID),
                  vjust=-0.1, hjust=-0.4, segment.color=NA, size=3, col="red") +
  theme_classic() +
  theme(legend.position=c(.15, .9))

ggsave("~/Desktop/Figure_GAM.png", fig, width=2200, height=1600, units="px", device="png", dpi=300)

