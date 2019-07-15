## Time dependent survival model
library(tidyverse)
library(ggplot2)

library(RColorBrewer)

library(psych)
library(Rtsne)

mypltt<-brewer.pal(7,"RdBu")

d <- read_csv("../data/parameters_v3.csv")

## Filtering features
## Removing clinical features measured over time
to.exclude <- c("Body temperture", "Weight", "Liver AST", "Liver ALT", "Liver ALB", "Liver ALP", "Kidney Creatinine", "Kidney BUN")
to.exclude <- str_which(colnames(d),paste(to.exclude, collapse="|"))
d <- d %>% select(colnames(d)[-to.exclude])
## Check variance
multiple.measures <- d %>% select(str_which(colnames(d), "d26|d54")) %>% gather %>% select(key) %>% mutate(name=str_replace(key,"d[0-9]+", "")) %>% count(name) %>% filter(n ==44) %>% select(name) %>% as_vector %>% paste(collapse="|")

pdf("../plots/hist_d54_d26.pdf", w = 10, h = 75)
d %>% filter(`Vaccine Group`!="Control") %>% select(str_which(colnames(.), multiple.measures)) %>% mutate_all(scale)  %>% gather %>%  mutate(
                                                                                                                                  name = str_replace(key, " d[0-9]+", ""),
                                                                                                                                  day = str_extract(key, "d[0-9]+")
                                                                                                                              ) %>% ggplot(aes(value, fill=day)) + geom_histogram() +  facet_grid(name ~.)
dev.off()

pdf("../plots/hist_d54_d26.pdf", w = 10, h = 75)
d %>% filter(`Vaccine Group`!="Control") %>% select(str_which(colnames(.), multiple.measures)) %>% mutate_all(scale) %>% summarise_all(var)  %>% gather %>%  mutate(
                                                                                                                                  name = str_replace(key, " d[0-9]+", ""),
                                                                                                                                  day = str_extract(key, "d[0-9]+")
                                                                                                                              ) %>% ggplot(aes(day, value)) + geom_col() +  facet_grid(name ~.)
dev.off()

## Polyfunctionality
polyfunc.cols <- colnames(d) %>% str_subset(., "ADNKD|ADCP|ADNP|ADCD")
polyfunc.clustersize  <- c(1,2,3,4,5,6)
d.polyfunc <- d %>% select(polyfunc.cols) %>% mutate_all(scale)

set.seed(11258)
polyfunc.withinss  <- polyfunc.clustersize %>% map_dbl(function(x){
    k <- d.polyfunc %>% kmeans(., x, nstart = 2)
    k$tot.withinss
})

plot(polyfunc.withinss)

polyfunc.kmeans <- d.polyfunc %>% kmeans(., 5, nstart = 20)

polyfunc.pca <- d.polyfunc %>% prcomp(., center = TRUE, scale. = TRUE)
polyfunc.pcadf <- data.frame(polyfunc.pca$x[,"PC1"], polyfunc.pca$x[,"PC2"], d$Survival, d$`Time of death`, polyfunc.kmeans$cluster)
colnames(polyfunc.pcadf) <- c("PC1", "PC2", "Survival", "Time of death", "cluster")
polyfunc.pcadf <- polyfunc.pcadf %>% mutate(cluster=as.factor(cluster))

pdf("../plots/ief_clustering.pdf")
ggplot(polyfunc.pcadf, aes(PC1, PC2, color=cluster, shape=Survival)) + geom_point()
dev.off()

polyfunc.tsne <- d.polyfunc %>% Rtsne(perplexity=5)
polyfunc.tsnedf <- data.frame(polyfunc.tsne$Y[,1], polyfunc.tsne$Y[,2], d$Survival, d$`Time of death`, polyfunc.kmeans$cluster)
colnames(polyfunc.tsnedf) <- c("comp1", "comp2", "Survival", "Time of death", "cluster")
polyfunc.tsnedf <- polyfunc.tsnedf %>% mutate(cluster=as.factor(cluster))

pdf("../plots/ief_clustering_tsne.pdf")
ggplot(polyfunc.tsnedf, aes(comp1, comp2, color=cluster, shape=Survival)) + geom_point()
dev.off()

## Cluster differences
d.cluster <- d %>% mutate(cluster = as.factor(polyfunc.kmeans$cluster)) %>% select(c(polyfunc.cols, "cluster"))

d.cluster.pval <- d.cluster %>% select(polyfunc.cols) %>% names(.) %>% map(~kruskal.test(get(.) ~ d.cluster$cluster, data=d.cluster)$p.value)
tibble(name=polyfunc.cols, pval=p.adjust(d.cluster.pval)) %>% filter(pval <= 0.05)

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

pdf("../plots/cor.pdf", w = 25, h = 20)
d.cor.r %>% filter((p.value<=0.05 & (value >= d.cor.limit | value <= -d.cor.limit)) & (variable !=name)) %>% mutate(variable = fct_drop(variable), name=fct_drop(name)) %>% ggplot(aes(variable, name, fill=value)) + geom_tile() + geom_label(aes(label=round(value, 1))) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_gradientn(colors=mypltt)
dev.off()

#Remove min and max values
## d  <- d %>% select(names(.)[!str_detect(names(.), "min|max")])

## Extract colnames done before infection
exo.var  <- names(d)[!str_detect(names(d), "d0|d3|d6|d9|d12|d15|d22|d28")]

#######################
## Survival analysis ##
#######################

library(survival)
library(survminer)

## Notes:
## NHP 7 died on day 9
## NHP 20 euthanized on day 28 but survived infection

d <- d %>% mutate(
               status = ifelse(d$Survival == "No", 1, 0)
           )

d$`Vaccine Group` <- as.factor(d$`Vaccine Group`)

names(d) <- gsub(" ", ".", names(d))

model <- coxph(Surv(Time.of.death, status) ~ Vaccine.Group, data = d)

fit <- survfit(Surv(Time.of.death, status) ~ Vaccine.Group, data = d)
pdf("../plots/surv.pdf")
ggsurvplot(fit, data = d, pval = TRUE, risk.table = TRUE, risk.table.col = "strata", ggtheme=theme_bw())
dev.off()

######################
## Outcome analysis ##
######################

library(rstanarm)
options(mc.cores = parallel::detectCores())
library(projpred)
library(caret)

zero.var <- d %>% select_if(is.numeric) %>% nearZeroVar(., names=TRUE)
fixed.vars <- names(d)[!str_detect(names(d), "d0|d3|d6|d9|d12|d15|d22|d28|min|max")]

## Select only values that are not measured over time
d.outcome  <- d.outcome %>% select(names(.)[(str_detect(names(.), "min|max") | names(.) %in% fixed.vars) & (!names(.) %in% zero.var )])

## Scale features
d.outcome <- d.outcome %>% mutate_if(is.numeric, scale)

d.outcome <- d.outcome %>% mutate(outcome=ifelse(d$`Time of death` < 28, 1, 0))

## Identify correlations
d.cor <- d.outcome %>% select_if(~!any(is.na(.))) %>% select(-`Vaccine Group`) %>% select_if(is.numeric) %>% data.matrix %>% corr.test(., method="kendall", adjust="holm")

## pvalues
d.cor.pval <- d.cor$p %>% as_tibble %>% mutate(name=colnames(.)) %>% gather(variable, value, -name)

## Kendall's tau value
d.cor.r <- d.cor$r %>% as_tibble %>% mutate(name=colnames(.)) %>% gather(variable, value, -name) %>% mutate(p.value = d.cor.pval$value)

d.cor.limit <- 0.7

d.cor.r %>% filter((p.value<=0.05 & (value >= d.cor.limit | value <= -d.cor.limit)) & (variable !=name)) %>% mutate(variable = fct_drop(variable), name=fct_drop(name)) %>% ggplot(aes(variable, name, fill=value)) + geom_tile() + geom_label(aes(label=round(value, 1))) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_gradientn(colors=mypltt)

## Define prior
SEED <- 112358
fit_partialpool <- stan_glmer( outcome ~ (1 + `Clinical Score max` + `max Viremia PFU (log10)` + `max Viremia RT-PCR` + `Body temperture max` + `Weight max (%)` + `Liver AST max` + `Liver ALT max` + `Liver ALB max` + `Liver ALP max` + `Kidney Creatinine max` + `Kidney BUN max` + `Body temperture min` + `Weight min (%)` + `Liver AST min` + `Liver ALT min` + `Liver ALB min` + `Liver ALP min` + `Kidney Creatinine min` + `Kidney BUN min` |`Vaccine Group`) + `Clinical Score max` + `max Viremia PFU (log10)` + `max Viremia RT-PCR` + `Body temperture max` + `Weight max (%)` + `Liver AST max` + `Liver ALT max` + `Liver ALB max` + `Liver ALP max` + `Kidney Creatinine max` + `Kidney BUN max` + `Body temperture min` + `Weight min (%)` + `Liver AST min` + `Liver ALT min` + `Liver ALB min` + `Liver ALP min` + `Kidney Creatinine min` + `Kidney BUN min`, data = d.outcome, family = binomial("logit"), seed = SEED, iter=6000, warmup=500, chains=4)

pdf("../plots/GLM.pdf")
plot(fit_partialpool)
dev.off()

loo(fit_partialpool)

draws <- as.matrix(fit_partialpool)
alphas  <- sweep(draws[,-1], 1, draws[,1],"+")

partialpool <- summary_stats(alphas)


partialpool <- partialpool[-nrow(partialpool),]


launch_shinystan(fit_partialpool)

posterior_vs_prior(fit_partialpool)

## Random Forest

d <- d %>% mutate(
                  outcome = as.factor(d$status)
                  )
levels(d$outcome) <- c("alive","dead")

library(doMC)
registerDoMC(cores = 16)
library(caret)
library(randomForest)

d.subset <- d %>% filter(Vaccine.Group!="Control") %>% select(-matches("temperture|kidney|liver|shedding|viremia|weight|clinical.score")) %>% select(-c("NHP", "Survival", "status", "Time.of.death", "Vaccine.Group")) 

d.dummy <- dummyVars(~Vaccine.Group , d)
d.dummy <- predict(d.dummy, d)

repeats <- 10

rfControl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats= repeats,
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
rf_random <- train(outcome ~ ., data=d.subset, method="rf", tuneLength=15, trControl=rfControl)
print(rf_random)
plot(rf_random)

pdf("../plots/var_imp.pdf")
v <- varImp(rf_random$finalModel)
tibble(gini=v, feature=rownames(v)) %>% ggplot(aes(x=reorder(feature, gini$Overall), y=gini$Overall)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + coord_flip()
dev.off()

predRf <- rf_random$pred[rf_random$pred$mtry == 2,]

getRocCV <- function(pred){
    temp.obs <- c()
    temp.pred <- c()
    for(i in seq(1, repeats)){
        temp <- paste("Rep", sprintf("%02d",i), sep="")
        temp.p <- pred[grepl(temp, pred[,"Resample"]), "dead"]
        temp.o <- pred[grepl(temp, pred[,"Resample"]), "obs"]
        temp.pred[[length(temp.pred)+1]] <- temp.p
        temp.obs[[length(temp.obs)+1]] <- temp.o
    }
    return(list(temp.pred, temp.obs));
}

getAvgRoc <- function(roc){
    m <- max(sapply(roc@x.values, length))
    resx <- sapply(roc@x.values, function(x){
        x <- c(x, rep(NA, m-length(x)));
    });
    resy <- sapply(roc@y.values, function(x){
        x <- c(x, rep(NA, m-length(x)));
    });
    roc.df <- data.frame(rowMeans(as.data.frame(resx), na.rm=T), rowMeans(as.data.frame(resy), na.rm=T))
    colnames(roc.df) <- c("meanx", "meany")
    return(roc.df)
}

library(ROCR)
predRoc <- getRocCV(predRf)
p <- prediction(predRoc[1][[1]], predRoc[2][[1]])
rocRf <- performance(p, "tpr", "fpr")
rocRf.avg <- getAvgRoc(rocRf)

auc <- performance(p, "auc")
auc <- mean(unlist(auc@y.values))

acc <- performance(p, "acc")
acc <- mean(unlist(acc@y.values))

roc.x <- do.call(c, rocRf@x.values)
roc.y <- do.call(c, rocRf@y.values)

npoints <- lapply(rocRf@x.values, length)

roc.fold <- sapply(seq(1:length(npoints)), function(x){
    rep(x, npoints[[x]])
}) %>% unlist


roc.df <- tibble(x=roc.x, y=roc.y, fold=roc.fold)
roc.avg.df <- roc.df %>% group_by(x) %>% summarise(maxy = max(y), miny=min(y), meany=mean(y), conf25=quantile(y,probs=c(.025,.975))[[1]], conf975=quantile(y,probs=c(.025,.975))[[2]])

pdf("../plots/roc.df")
ggplot(roc.avg.df) + geom_line(aes(x = x, y = meany), alpha=1, size=0.5) + geom_ribbon(aes(ymin = conf25, ymax=conf975, x =x), fill = "grey70", alpha =0.5) + geom_abline(color="#707070") +xlab("FPR") +ylab("TPR") + ggtitle(paste("ROC Curves for Random Forest model with 10-fold CV repeated", repeats, "times"))+ coord_cartesian(xlim = c(0, 1), ylim=c(0,1)) + theme_bw() + annotate("text", label=paste0("AUC: ", round(auc,2)), x=0.1,y=1) + annotate("text", label=paste0("ACC: ", round(acc,2)), x=0.1,y=0.95)
dev.off()
