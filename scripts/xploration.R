## Time dependent survival model
library(tidyverse)
library(ggplot2)

library(RColorBrewer)

library(psych)

d <- read_csv("../data/parameters_v3.csv")

## Polyfunctionality
polyfunc.cols <- colnames(d) %>% str_subset(., "ADNKD|ADCP|ADNP|ADCD")
polyfunc.clustersize  <- c(1,2,3,4,5,6)

set.seed(11258)
polyfunc.withinss  <- polyfunc.clustersize %>% map_dbl(function(x){
    k <- d %>% select(polyfunc.cols) %>% kmeans(., x, nstart = 20)
    k$tot.withinss
})

plot(polyfunc.withinss)

polyfunc.kmeans <- d %>% select(polyfunc.cols) %>% kmeans(., 3, nstart = 20)

polyfunc.pca <- d %>% select(polyfunc.cols) %>% prcomp(., center = TRUE, scale. = TRUE)
polyfunc.pcadf <- data.frame(polyfunc.pca$x[,"PC1"], polyfunc.pca$x[,"PC2"], d$Survival, d$`Time of death`, polyfunc.kmeans$cluster)
colnames(polyfunc.pcadf) <- c("PC1", "PC2", "Survival", "Time of death", "cluster")
polyfunc.pcadf <- polyfunc.pcadf %>% mutate(cluster=as.factor(cluster))


pdf("ief_clustering.pdf")
ggplot(polyfunc.pcadf, aes(PC1, PC2, color=cluster, shape=Survival)) + geom_point()
dev.off()

## Pairwise correlations
d.cor <- d %>% select_if(~!any(is.na(.))) %>% select(-NHP, -`Vaccine Group`) %>% select_if(is.numeric) %>% data.matrix %>% corr.test(., method="kendall", adjust="holm")

## pvalues
d.cor.pval <- d.cor$p %>% as_tibble %>% mutate(name=colnames(.)) %>% gather(variable, value, -name)

## Kendall's tau value
d.cor.r <- d.cor$r %>% as_tibble %>% mutate(name=colnames(.)) %>% gather(variable, value, -name) %>% mutate(p.value = d.cor.pval$value)
d.cor.r %>% filter(name %in% polyfunc.cols & variable %in% polyfunc.cols & name != variable) %>% arrange(-desc(value))

## pvalues vs value of Kendall's Tau
d.cor.limit  <- 0.7
d.cor.r %>% filter((value >= d.cor.limit | value <= -d.cor.limit) & p.value <=0.05 & variable!=name)

mypltt<-brewer.pal(7,"RdBu")

pdf("../plots/cor.pdf", w = 25, h = 20)
d.cor.r %>% filter((p.value<=0.05 & (value >= d.cor.limit | value <= -d.cor.limit)) & (variable !=name)) %>% mutate(variable = fct_drop(variable), name=fct_drop(name)) %>% ggplot(aes(variable, name, fill=value)) + geom_tile() + geom_label(aes(label=round(value, 1))) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_gradientn(colors=mypltt)
dev.off()

#Remove min and max values
## d  <- d %>% select(names(.)[!str_detect(names(.), "min|max")])

## Extract colnames done before infection
exo.var  <- names(d)[!str_detect(names(d), "d0|d3|d6|d9|d12|d15|d22|d28")]

## ## Extract time and long by time
## d.lng <- d %>% gather(variable, value, -exo.var) %>% mutate(time = ifelse(variable != "Time of death", str_extract(variable, "(?<=d)[0-9]+"),  value)) %>% mutate(., time=as.numeric(time)) %>% mutate(., variable = ifelse(variable != "Time of death", str_extract(variable, ".+(?= d[0-9]+)") %>% str_trim(.), variable)) %>% spread(variable, value)

d.lng <- d %>% gather(variable, value, -exo.var) %>% mutate(time = str_extract(variable, "(?<=d)[0-9]+")) %>% mutate(., time=as.numeric(time)) %>% mutate(., variable = str_extract(variable, ".+(?= d[0-9]+)") %>% str_trim(.)) %>% spread(variable, value)

#######################
## Survival analysis ##
#######################

library(survival)
library(survminer)

## Notes:
## NHP 7 died on day 9
## NHP 20 euthanized on day 28 but survived infection

## Setup time intervals
time.measured.points <- d.lng$time %>% unique %>% sort
time.start.points <- c(d$`Time of death`, d.lng$time) %>% unique %>% sort
time.end.points <- lead(time.start.points, default=Inf)

d.time <- map2_dfr(time.start.points, time.end.points, function(x,y){
    d.interval <- d.lng %>% filter(time == x)
    if(d.interval %>% count == 0){
        time.intersect <- intersect(1:time.start.points[which(time.start.points == 7)], time.measured.points)
        d.lng %>% filter( time== time.intersect %>% last) %>% mutate(tstart =x , tend = y) %>% data.frame
    } else {
        d.interval %>% mutate(tstart =x , tend = y) %>% data.frame
    }
}) %>% as_tibble

d.time %>% filter(NHP == 3) %>% select(tstart, tend, time, Time.of.death)

d.time  <- d.time %>% mutate(death = 
d.time %>% select(tstart, tend, Time.of.death) %>% pmap_dbl(function(tstart, tend, Time.of.death){
    if(tend < Time.of.death){                        #Before death death == 0
        0
    } else if (tend == Time.of.death){               #On time of death if day is  not 28 then 1 else 0. All mice survived on day 28
        if(tend == 28){
            0
        } else {
            1
        }
    } else {                            #If time > time of death, death is 1
        NA
    }
})) %>% filter(!is.na(death))

d.time %>% filter(NHP == 2) %>% select(tstart, tend, time, Time.of.death, death)

write_csv(d.time, "../data/parameters_with_time_intervals.csv")

d.time %>% filter(death == 1) %>% count(NHP)

## Create survival variable

## Temporary
d.time$Group <- as.factor(d.time$Vaccine.Group)
levels(d.time$Group) <- c(1:6)

fit <- survfit(Surv(tstart, tend, death == 1) ~ Group, data = d.time)
png("../plots/surv.png")
ggsurvplot(fit, data = d.time, censor.shape="x", censor.size = 4, risk.table=TRUE, ncensor.plot=TRUE)
dev.off()

d.cox <- coxph(Surv(tstart, tend, death) ~ cluster(NHP) + strata(Vaccine.Group) + (Body.temperture + Kidney.BUN  + Kidney.Creatinine + Liver.ALB  + Liver.ALP  + Liver.ALT + Liver.AST + Weight), d.time)

zp <- cox.zph(d.cox)

## BC/LSFHN/EboGP
## HPIV1/EboGP

######################
## Outcome analysis ##
######################

## Select only values that are not measured over time
d.outcome <- d %>% mutate(outcome=ifelse(`Time of death` < 28, 1, 0))

library(rstanarm)
library(projpred)

d.outcome  <- d.outcome %>% select(names(.)[str_detect(names(.), "min|max") | names(.) %in% c("Vaccine Group", "outcome")])

d.outcome <- d.outcome %>% mutate_at(str_subset(names(d.outcome), "max"), scale)

d.outcome %>% select_if(~!any(is.na(.))) %>% select(-`Vaccine Group`) %>% select_if(is.numeric) %>% data.matrix %>% corr.test(., method="kendall", adjust="holm")

## Define prior
SEED <- 112358
fit_partialpool <- 
  stan_glm( outcome ~ `Clinical Score max` + `max Viremia PFU (log10)` + `max Viremia RT-PCR` + `Body temperture max` + `Weight max (%)` + `Liver AST max` + `Liver ALT max` + `Liver ALB max` + `Liver ALP max` + `Kidney Creatinine max` + `Kidney BUN max` + `Body temperture min` + `Weight min (%)` + `Liver AST min` + `Liver ALT min` + `Liver ALB min` + `Liver ALP min` + `Kidney Creatinine min` + `Kidney BUN min`, data = d.outcome, 
             family = binomial("logit"), seed = SEED, iter=6000, warmup=500, chains=4)

pdf("../plots/GLM.pdf")
plot(fit_partialpool)
dev.off()

launch_shinystan(fit_partialpool)

posterior_vs_prior(fit_partialpool)

## Random Forest

d.outcome <- d.outcome %>% mutate(
                  outcome = as.factor(d.outcome$outcome)
                  )
levels(d.outcome$outcome) <- c("alive","dead")

library(randomForest)
library(caret)

rfControl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats= 10,
  verboseIter = FALSE,
  returnData = FALSE,
  allowParallel = TRUE,
  summaryFunction = twoClassSummary,
  classProbs = TRUE,
  savePredictions = T,
  search="grid"
)
set.seed(11258)
mtry <- sqrt(21)
rf_random <- train(outcome ~ ., data=data.frame(d.outcome), method="rf", tuneLength=15, trControl=rfControl)
print(rf_random)
plot(rf_random)

pdf("../plots/var_imp.pdf")
v <- varImp(rf_random$finalModel)
tibble(gini=v, feature=rownames(v)) %>% ggplot(aes(x=reorder(feature, gini$Overall), y=gini$Overall)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
