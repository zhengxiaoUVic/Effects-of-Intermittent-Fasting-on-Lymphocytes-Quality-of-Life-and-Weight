######################
## Data Preparation ##
######################
df <- read.csv("~/Downloads/IF_v3.1.csv", stringsAsFactors=T)
df$Fast <- factor(ifelse(df$ID %in% paste0("IF-0",c(17,19,20,21,22)), 
                         "Six", "Three")) #add a dummy variable for fasting groups
df$Phase <- factor(df$Phase)
dat <- df[!is.na(df$LYMP) & df$Phase!="2",] #IF-013 contains one missing value

######################
##     Modeling     ##
######################
library(nlme)

fit <- lme(fixed=LYMP ~ Phase * Fast * Time,
           random=~Time | ID, data=dat, method="ML")
#can also test `random=~Time * Phase | ID`
#         and `random=~Time * Phase * Fast | ID`
#         with `control=lmeControl(opt="optim", optCtrl=list(maxit=1e5))`

summary(fit) #p-values
coef(fit) #coefficients
intervals(fit) #95% confidence intervals

