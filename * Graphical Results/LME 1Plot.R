## Data
df <- read.csv("~/Downloads/IF_v3.1.csv", stringsAsFactors=T)
df$Phase <- factor(df$Phase)
df$ID_old <- df$ID
levels(df$ID) <- c(paste0("S0",1:9), paste0("S",10:15))


## Linear Mixed-Effects (LME) Model
library(nlme) #v 3.1-165

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
pooled_mean_raw <- mean(subset(df, Time==0)$LYMP) #mean raw baseline

df <- df %>%
  group_by(ID) %>%
  mutate(Baseline = LYMP[Time==0]) %>%
  mutate(Baseline_pred = pred[Time==0]) %>%
  mutate(raw_centered = (LYMP - Baseline + pooled_mean_raw) / pooled_mean_raw) %>% #create normalized raw
  mutate(pred_centered = (pred - Baseline_pred + pooled_mean) / pooled_mean) %>% #create normalized predictions
  mutate(LYMP_normalized = LYMP / Baseline) %>% #create self-normalized raw
  mutate(Seg = factor(if_else(Time==0, "1", Phase))) %>% #create segments for plotting
  select(-c(Baseline, Baseline_pred))

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


## Visualization
library(ggplot2) #v 3.5.1 (many deprecated theme controls)
library(ggrepel) #v 0.9.5

# LME fit for the intermittent fasting and follow-up (Figure 2A)
fig_2Phases <- ggplot(subset(df, Seg != 0),
              aes(x=Time, y=pred_centered, group=interaction(ID, Seg), color=Seg, linetype=Seg)) +
  geom_line() +
  labs(x="Month", y="Normalized Lymphocyte Prediction",
       title="Linear Mixed-Effects Model", color="Treatment", linetype="Treatment") +
  scale_x_continuous(limit=c(0,6.5), breaks=c(0,1,2,3,4,5,6), expand=c(0,0)) +
  scale_color_manual(values=c("black","dodgerblue"),
                     labels=c("1"="Intermittent Fasting", "2"="Follow-Up")) +
  scale_linetype_manual(values=c("solid","dashed"),
                        labels=c("1"="Intermittent Fasting", "2"="Follow-Up")) +
  geom_text_repel(data=last_LYMP, aes(x=Time, y=pred_centered, label=ID),
                  vjust=0.8, hjust=-0.4, segment.color=NA, size=3, col="red") +
  geom_hline(yintercept=1, linetype="dashed", color="gray60") +
  theme_classic() +
  theme(legend.position=c(.15, .9))

ggsave("~/Desktop/Figure_LMER_2Phases.png", fig_2Phases, width=2200, height=1600, units="px", device="png", dpi=300)


# LME fit for the pre-intervention, intermittent fasting and follow-up
fig_3Phases <- ggplot(subset(df, Time >= -3),
              aes(x=Time, y=pred_centered, group=interaction(ID, Seg), color=Seg, linetype=Seg)) +
  geom_line() +
  labs(x="Month", y="Normalized Lymphocyte Prediction",
       title="Linear Mixed-Effects Model", color="Treatment", linetype="Treatment") +
  scale_x_continuous(limit=c(-3.5,6.5), breaks=-3:6, expand=c(0,0)) +
  scale_color_manual(values=c("dodgerblue","black","dodgerblue2"),
                     labels=c("0"="Pre-Intervention", "1"="Intermittent Fasting", "2"="Follow-Up")) +
  scale_linetype_manual(values=c("longdash","solid","dashed"),
                        labels=c("0"="Pre-Intervention", "1"="Intermittent Fasting", "2"="Follow-Up")) +
  geom_text_repel(data=last_LYMP, aes(x=Time, y=pred_centered, label=ID),
                  vjust=0.8, hjust=-0.4, segment.color=NA, size=3, col="red") +
  geom_text_repel(data=subset(df, Time==-3), aes(x=Time, y=pred_centered, label=ID),
                  vjust=0.5, hjust=0.4, segment.color=NA, size=3, col="red") +
  geom_hline(yintercept=1, linetype="dashed", color="gray60") +
  theme_classic() +
  theme(legend.position=c(.12, .85))

ggsave("~/Desktop/Figure_LMER_3Phases.png", fig_3Phases, width=2550, height=1600, units="px", device="png", dpi=300)


# raw lymphocytes for the intermittent fasting and follow-up
fig_2Phases_raw <- ggplot(subset(df, Seg != 0),
                      aes(x=Time, y=raw_centered, group=interaction(ID, Seg), color=Seg, linetype=Seg)) +
  geom_line() +
  labs(x="Month", y="Normalized Lymphocyte Count",
       title="Raw Data", color="Treatment", linetype="Treatment") +
  scale_x_continuous(limit=c(0,6.5), breaks=c(0,1,2,3,4,5,6), expand=c(0,0)) +
  scale_color_manual(values=c("black","dodgerblue"),
                     labels=c("1"="Intermittent Fasting", "2"="Follow-Up")) +
  scale_linetype_manual(values=c("solid","dashed"),
                        labels=c("1"="Intermittent Fasting", "2"="Follow-Up")) +
  geom_text_repel(data=last_LYMP, aes(x=Time, y=raw_centered, label=ID),
                  vjust=0.8, hjust=-0.4, segment.color=NA, size=3, col="red") +
  geom_hline(yintercept=1, linetype="dashed", color="gray60") +
  theme_classic() +
  theme(legend.position=c(.15, .9))

ggsave("~/Desktop/Figure_raw_2Phases.png", fig_2Phases_raw, width=2200, height=1600, units="px", device="png", dpi=300)

library(plotly) #v 4.10.4
ggplotly(ggplot(subset(df, Seg != 0), aes(x=Time, y=raw_centered, color=ID))+geom_line()) #IF-013 contains one missing value at Month 5

# Note: the three figures above used the subject-centering transformation


# self-normalized raw lymphocytes for the intermittent fasting and follow-up
fig_raw <- ggplot(subset(df, Seg != 0),
                  aes(x=Time, y=LYMP_normalized, group=interaction(ID, Seg), color=Seg, linetype=Seg)) +
  geom_line() +
  labs(x="Month", y="Relative Lymphocyte Count",
       title="Self-Normalized Raw Data", color="Treatment", linetype="Treatment") +
  scale_x_continuous(limit=c(0,6.5), breaks=c(0,1,2,3,4,5,6), expand=c(0,0)) +
  scale_color_manual(values=c("black","dodgerblue"),
                     labels=c("1"="Intermittent Fasting", "2"="Follow-Up")) +
  scale_linetype_manual(values=c("solid","dashed"),
                        labels=c("1"="Intermittent Fasting", "2"="Follow-Up")) +
  geom_text_repel(data=last_LYMP, aes(x=Time, y=LYMP_normalized, label=ID),
                  vjust=0.8, hjust=-0.4, segment.color=NA, size=3, col="red") +
  geom_hline(yintercept=1, linetype="dashed", color="gray60") +
  theme_classic() +
  theme(legend.position=c(.15, .9))

ggsave("~/Desktop/Figure_self.raw.png", fig_raw, width=2200, height=1600, units="px", device="png", dpi=300)

