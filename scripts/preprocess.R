library(tidyverse)
library(psych)
library(Rtsne)
library(RColorBrewer)
library(survival)
library(Cyclops)
library(BrokenAdaptiveRidge)

df <- read_csv("data/parameters_v3.csv")
names(df) <- make.names(names(df)) %>% str_replace_all("(\\.){2,5}", ".")
iga <- read_csv("data/IgA_peptide_array.csv")
igg <- read_csv("data/IgG_peptide_array.csv")

pep.arr <- inner_join(iga, igg, "NHP", suffix=c(".iga", ".igg"))
df.merged <- inner_join(pep.arr, df, "NHP")

## Exclude clinical attributes and scale
df.subset <- df.merged %>% select(-matches("temperture|kidney|liver|shedding|viremia|weight|clinical.score")) %>% mutate(
                                                                                                                      NHP=as.factor(NHP),
                                                                                                                      BC.LSFcHN.EboGP = ifelse(Vaccine.Group == "BC/LSFcHN/EboGP", 1, 0),
                                                                                                                      BC.LSFHN.EboGP = ifelse(Vaccine.Group == "BC/LSFHN/EboGP", 1, 0),
                                                                                                                      HPIV3..FHN.EboGP = ifelse(Vaccine.Group == "HPIV3/?FHN/EboGP", 1, 0),
                                                                                                                      HPIV1.EboGP = ifelse(Vaccine.Group == "HPIV1/EboGP", 1, 0),
                                                                                                                      HPIV3.EboGP = ifelse(Vaccine.Group == "HPIV3/EboGP", 1, 0)

                                                                                                                       ) %>% mutate_if(is.numeric, scale)
df.subset <- df.subset %>% mutate(
                               Survival = ifelse(Survival == "Yes", 0, 1),
                               Time.of.death = df$Time.of.death[1:20] # remove Control
                           )

## Peptide array: IgA and IgG
df.subset %>% select(matches("Peptide_[0-9]+.iga|Vaccine.Group|Survival")) %>% gather(var, val, -Vaccine.Group, -Survival) %>% ggplot(aes(var, val)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_grid(Vaccine.Group + Survival ~ .)

getSignificantPairwiseCorrelations <- function(df){
    df <- corr.test(ief, method="kendall", adjust="holm")
    df.r <- df$r %>% data.frame %>% mutate(var1 = rownames(.)) %>% gather(var2, value, -var1)
    df.p  <- df$p %>% data.frame %>% mutate(var1 = rownames(.)) %>% gather(var2, value, -var1)
    df.merged <- inner_join(df.r, df.p, c("var1", "var2"), suffix=c(".r", ".p"))
    return(df.merged);
}

## Create polyfunctionality from Immune effectors functions
ief <- df.subset %>% select(matches("ADCD|ADCP|ADNKD|ADNP"))
## Check correlations between ief
ief.corr <- getSignificantPairwiseCorrelations(ief)
## Plot heatmap
ief.corr %>% filter(value.r >= 0.5 & var1!=var2) %>% ggplot(aes(var1, var2)) + geom_tile(aes(fill=value.r)) + geom_text(aes(label=round(value.r, 2))) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

## k means
set.seed(112358)
ief.kmeans <- kmeans(ief, 3)

## Plot clusters
set.seed(112358)
ief.tsne <- ief %>% Rtsne(perplexity=5)

ief.kmeans.tsne <- data.frame(
    x = ief.tsne$Y[,1],
    y = ief.tsne$Y[,2],
    cluster = as.factor(ief.kmeans$cluster),
    survival = as.factor(df.subset$Survival),
    vaccine = df.subset$Vaccine.Group
) %>% cbind(., ief)

ggplot(ief.kmeans.tsne, aes(x,y, color = cluster)) + geom_point(aes(shape=as.factor(survival)), size = 5) + geom_text(aes(label = vaccine), hjust = 0, vjust = 1) + scale_colour_brewer(palette = "Set1") + scale_x_continuous(limits = c(-20, 20))

pdf("plots/polyfunc_clusters_differences.pdf", w=10, h = 70)
ief.kmeans.tsne %>% select(-c(survival, x, y, vaccine)) %>% gather(var, val, -cluster) %>% ggplot(aes(x=cluster, y=val)) + geom_boxplot() + facet_grid(var ~ .)
dev.off()

## Survival models
all_cox_formula <- as.formula(paste("Surv(Time.of.death, Survival) ~", paste0(
                                                                           names(df.subset) %>% str_subset("Survival|Time|NHP|Vaccine.Group", negate=TRUE),
  collapse=" + ")))

cyclopsDataCox <- createCyclopsData(all_cox_formula, data = df.subset, modelType = "cox")

reg_fit_cox <- fitCyclopsModel(cyclopsDataCox, prior = createPrior("laplace", variance = 10), forceNewObject = TRUE)
coef(reg_fit_cox)[which(coef(reg_fit_cox) != 0)]

## Checking confidence intervals
confint(reg_fit_cox, names(coef(reg_fit_cox)[which(coef(reg_fit_cox) != 0)]), overrideNoRegularization = TRUE)

bar_fit_cox <- fitCyclopsModel(cyclopsDataCox,
                               prior = BrokenAdaptiveRidge::createBarPrior(
                                 penalty = "bic", initialRidgeVariance = 1000, fitBestSubset = TRUE),
                               forceNewObject = TRUE)
coef(bar_fit_cox)[which(coef(bar_fit_cox) != 0)]

model.var <- lapply(seq(10, 3200, 100), function(x){
    bar_fit_cox <- fitCyclopsModel(cyclopsDataCox,
                                   prior = BrokenAdaptiveRidge::createBarPrior(
                                                                    penalty = "bic", initialRidgeVariance = x, fitBestSubset = TRUE),
                                   forceNewObject = TRUE)
    print(x)
    coef(bar_fit_cox)[which(coef(bar_fit_cox) != 0)]
})

## Only vaccine groups
vaccine.groups <- c("BC.LSFcHN.EboGP", "BC.LSFHN.EboGP", "HPIV3..FHN.EboGP", "HPIV1.EboGP", "HPIV3.EboGP")

interaction.terms <- apply(expand.grid(vaccine.groups, vaccine.groups), 1, paste0, collapse=" * ") %>% paste0(., collapse=" + ")
groups_formula_cox <- as.formula(paste0(c("Surv(Time.of.death, Survival) ~ BC.LSFcHN.EboGP + BC.LSFHN.EboGP + HPIV3..FHN.EboGP + HPIV1.EboGP + HPIV3.EboGP", interaction.terms), collapse=" + "))

cyclopsDataCox <- createCyclopsData(groups_formula_cox, data = df.subset, modelType = "cox")

reg_fit_cox <- fitCyclopsModel(cyclopsDataCox, prior = createPrior("laplace", variance = 10), forceNewObject = TRUE)
coef(reg_fit_cox)[which(coef(reg_fit_cox) != 0)]

## Checking confidence intervals
confint(reg_fit_cox, names(coef(reg_fit_cox)[which(coef(reg_fit_cox) != 0)]), overrideNoRegularization = TRUE)

bar_fit_Lr <- fitCyclopsModel(cyclopsDataCox,
                               prior = BrokenAdaptiveRidge::createBarPrior(
                                 penalty = "bic", initialRidgeVariance = 1000, fitBestSubset = TRUE),
                               forceNewObject = TRUE)
coef(bar_fit_Lr)[which(coef(bar_fit_Lr) != 0)]



## LR models
all_lr_formula <- as.formula(paste("Survival ~", paste0(
                                                                           names(df.subset) %>% str_subset("Survival|Time|NHP|Vaccine.Group", negate=TRUE),
  collapse=" + ")))

cyclopsDataLr <- createCyclopsData(all_lr_formula, data = df.subset, modelType = "lr")

reg_fit_Lr <- fitCyclopsModel(cyclopsDataLr, prior = createPrior("laplace", variance = 10, exclude = "(Intercept)"), forceNewObject = TRUE)
coef(reg_fit_Lr)[which(coef(reg_fit_Lr) != 0)]
## Checking confidence intervals
confint(reg_fit_Lr, names(coef(reg_fit_Lr)[which(coef(reg_fit_Lr) != 0)]), overrideNoRegularization = TRUE)

bar_fit_Lr <- fitCyclopsModel(cyclopsDataLr,
                               prior = BrokenAdaptiveRidge::createBarPrior(
                                 penalty = "bic", initialRidgeVariance = 1000, fitBestSubset = TRUE, exclude = "(Intercept)"),
                               forceNewObject = TRUE)
coef(bar_fit_Lr)[which(coef(bar_fit_Lr) != 0)]

model.var <- lapply(seq(10, 1000, 10), function(x){
    bar_fit_Lr <- fitCyclopsModel(cyclopsDataLr,
                                   prior = BrokenAdaptiveRidge::createBarPrior(
                                                                    penalty = "bic", initialRidgeVariance = x, fitBestSubset = TRUE),
                                   forceNewObject = TRUE)
    print(x)
    coef(bar_fit_Lr)[which(coef(bar_fit_Lr) != 0)]
})

