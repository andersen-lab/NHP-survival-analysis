library(tidyverse)
library(psych)
library(Rtsne)
library(RColorBrewer)
library(survival)
library(Cyclops)
library(BrokenAdaptiveRidge)
library(glmnet)
library(network)
library(sna)
library(ggnetwork)

df <- read_csv("data/parameters_v5-3_zero.csv")
names(df) <- make.names(names(df)) %>% str_replace_all("(\\.){2,5}", ".")
iga <- read_csv("data/IgA_peptide_array.csv")
igg <- read_csv("data/IgG_peptide_array.csv")

## Peptide array to domains
GP1.base <- c(6:18, 21:27, 37:42, 42:48)
GP1.head <- c(15:24, 24:40, 40:44, 51:57)
Rbd <- c(24:40)
glycan.cap <- c(54:79)
mucin <- c(76:116)
fusion.loop <- c(125:139)
hr1 <- c(137:150)
hr2 <- c(147:158)
mper <- c(156:163)

GP.regions <- list(
    GP1.base, GP1.head, Rbd, glycan.cap, mucin, fusion.loop, hr1, hr2, mper
)
domain.names <- as.character(substitute(c(GP1.base, GP1.head, Rbd, glycan.cap, mucin, fusion.loop, hr1, hr2, mper)))[2:10] %>% str_c("Peptide.", .)

## Take average of val for each domain
domain.mean <- sapply(GP.regions, function(pep){
    igg %>% select(c(NHP, str_c("Peptide_", str_pad(pep, 3, pad="0")))) %>% gather(var, val, -NHP) %>% spread(NHP, val) %>% summarise_all(mean) %>% select(-var) %>% as_vector
}) %>% data.frame
colnames(domain.mean) <- domain.names

pep.arr <- domain.mean %>% add_rownames("NHP") %>% mutate(NHP=as.integer(NHP))

pep.arr <- inner_join(iga, igg, "NHP", suffix=c(".iga", ".igg"))
df.merged <- inner_join(pep.arr, df, "NHP")

## Exclude clinical attributes and scale
df.subset <- df.merged %>% select(-matches("temperture|kidney|liver|shedding|viremia|weight")) %>% mutate(
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
df.subset %>% select(matches("Peptide.|Vaccine.Group|Survival")) %>% gather(var, val, -Vaccine.Group, -Survival) %>% ggplot(aes(var, val)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_grid(Vaccine.Group + Survival ~ .)

## Hierarchical clustering
l <- df.subset %>% select(-matches("Peptide")) %>% bind_cols(df.subset.peptide) %>% select_if(is.numeric) %>% corr.test(., method="kendall", adjust="holm")
df.subset.hclust <- hclust(dist(l$r), method="centroid")

getSignificantPairwiseCorrelations <- function(df){
    df <- df %>% select_if(is.numeric)
    df <- corr.test(df, method="kendall", adjust="BH")
    df$r[lower.tri(df$r)] <- NA
    df$p[lower.tri(df$p)] <- NA
    df.r <- df$r %>% data.frame %>% mutate(var1 = rownames(.)) %>% gather(var2, value, -var1) %>% drop_na(value)
    df.p  <- df$p %>% data.frame %>% mutate(var1 = rownames(.)) %>% gather(var2, value, -var1) %>% drop_na(value)
    df.merged <- inner_join(df.r, df.p, c("var1", "var2"), suffix=c(".r", ".p"))
    return(df.merged);
}

## Pairwise correlations
## Subset features to test
iga.peptides <- c(79, 89, 111, 122, 123) %>% str_pad(3, pad="0") %>% paste0("Peptide_", ., ".iga")
igg.peptides <- c(79, 81, 109, 111, 118, 119, 122) %>% str_pad(3, pad="0") %>% paste0("Peptide_", ., ".iga")
pep.breadth <- colnames(df.subset) %>% str_subset("Peptide.Breadth")

df.subset.infectivity <- colnames(df.subset) %>%str_subset("Restoration") %>% str_subset("\\.2ug")
df.subset.ief <- colnames(df.subset) %>% str_subset("ADCD|ADNP|ADCP") %>% str_subset("\\.004|\\.0004|\\.0001") %>% str_subset("d54")

covariates.to.include <- df.subset %>% select(matches("Mab|GP.binding|GPcl.binding|sGP.binding|Serum.Neut")) %>% colnames
covariates.to.include <- c(covariates.to.include, "Proportion.binding.to.RBD.GC.GP2.d54.GPmuc.immobilized.", "Proportion.binding.to.RBD.d54.GPcl.immobilized.vs.sGP.", "Serum.IgA.d54", "Serum.IgG.d54")
source.nodes.to.include <- c(df.subset.ief, "Proportion.binding.to.RBD.d54.GPcl.immobilized.vs.sGP.", "GP.binding.BLI.d54", "GP.binding.BLI.d26")

## igg.peptides, iga.peptides
df.subset.selected <- df.subset %>% select(c(covariates.to.include, df.subset.infectivity, df.subset.ief, pep.breadth))

## select(-matches("Peptide"))
df.merged.corr <- df.subset %>% select(-matches("Peptide|Restoration|ADCD|ADNP|ADCP|Mab.Comp|ADNKD|Proportion|binding|Serum|Restoration|Mucosal")) %>% bind_cols(df.subset.selected)  %>% mutate(
                                                                                                                                 Survival = 1 - Survival
                                                                                                                             ) %>%  getSignificantPairwiseCorrelations

## Map to categories
df.categories <- data.frame(name=unique(c(df.merged.corr$var1, df.merged.corr$var2)))

getGroup  <- function(x){
    if(str_detect(x, "IgG|IgA|IgM")){
        "IgG/IgA/IgM";
    } else if(str_detect(x, "ADCD|ADCP|ADNKD|ADNP")){
        "Immune effector function"
    } else if(str_detect(x, "binding|Mab|Proportion")){
        "Binding assay";
    } else if(str_detect(x, "BC.LSFcHN.EboGP|BC.LSFHN.EboGP|HPIV3..FHN.EboGP|HPIV1.EboGP|HPIV3.EboGP")){
        "Vaccine group";
    } else if(str_detect(x, "Restoration")){
        "Restoration of infectivity";
    } else if(str_detect(x, "Neut")){
        "Neutralization assay";
    } else{
        "Clinical";
    }
}

df.categories <- df.categories %>% mutate(group = df.categories$name %>% map_chr(getGroup))
df.categories %>% write_csv("data/node_groups.csv")

## Write to file
df.merged.corr %>% filter(((var1 != var2 & value.p <= 0.05) | (var1 %in% covariates.to.include) | (var2 %in% covariates.to.include)) & abs(value.r) >= 0.6) %>% mutate(id = c(1:nrow(.))) %>% write_csv("data/pairwise_correlations.csv")

df.merged.corr %>% filter(var1 != var2 & abs(value.r) >= 0.5 & value.p <= 0.05) %>% write_csv("data/pairwise_correlations_network.csv")

df.merged.corr %>% filter(value.r >= 0.6 & var1!=var2 & value.p <= 0.05) %>% ggplot(aes(var1, var2)) + geom_tile(aes(fill=value.r)) + geom_text(aes(label=round(value.r, 2))) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/pairwise_univariate.pdf", w = 30, h = 30)

df.merged.corr %>% filter(value.r >= 0.6 & var1!=var2) %>% ggplot(aes(var1, var2)) + geom_tile(aes(fill=value.r)) + geom_text(aes(label=round(value.r, 2))) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/pairwise_univariate_all.pdf", w = 30, h = 30)

## Plot network
df.corr <- df.merged.corr %>% filter(var1 != var2 & abs(value.r) >= 0.5 & (var1 %in% c("survival.index", df.subset.ief) | (var2 %in% c("survival.index", df.subset.ief))))
val <- df.corr %>% select(value.r) %>% as_vector
pval <- df.corr %>% select(value.p) %>% as_vector
abs.val <- df.corr %>% select(value.r) %>% abs %>% as_vector
abs.val <- abs.val * 10

nodes <- unique(c(df.corr$var1, df.corr$var2))
num.nodes <- length(nodes)

n <- network.initialize(num.nodes)
network.vertex.names(n) <- nodes
n[as.matrix(df.corr %>% select(var1, var2))] <- 1
set.edge.attribute(n, "abs.val", abs.val)
set.edge.attribute(n, "val", val)
set.edge.attribute(n, "pval", pval)
set.vertex.attribute(n, "group", df.categories$group)

n <- ggnetwork(n, layout = "fruchtermanreingold", cell.jitter = 0, weights = "abs.val")

ggplot(n, aes(x, y, xend = xend, yend = yend)) +
  geom_edges(aes(linetype = pval <= 0.05, color = val >= 0)) +
    geom_nodes(aes(color = group)) +
    geom_nodelabel_repel(aes(color = group, label = vertex.names), fontface = "bold", box.padding = unit(1, "lines")) +
    theme_blank()
ggsave("plots/force_network.svg", w= 10, h = 10, device="svg")

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
                                                                           names(df.subset) %>% str_subset("Survival|Time|NHP|Vaccine.Group|survival.index|Clinical.Score.max", negate=TRUE),
  collapse=" + ")))

cyclopsDataCox <- createCyclopsData(all_cox_formula, data = df.subset, modelType = "cox")

reg_fit_cox <- fitCyclopsModel(cyclopsDataCox, prior = createPrior("laplace", variance = 10), forceNewObject = TRUE)

## Checking confidence intervals
coef.cox <- coef(reg_fit_cox)[which(coef(reg_fit_cox) != 0)]
ci.cox <- confint(reg_fit_cox, names(coef.cox), overrideNoRegularization = TRUE) %>% data.frame %>% rownames_to_column %>% mutate(coef=c(coef.cox))

ci.cox %>% ggplot(aes(x=as.factor(rowname), y=coef)) + geom_point() + geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..),width=0.1) + geom_hline(yintercept=0) + coord_flip() + ylab("Coefficient") + xlab("Variable") + theme_bw()
ggsave("plots/cox_coef.pdf", w= 6, h = 6)


bar_fit_cox <- fitCyclopsModel(cyclopsDataCox,
                               prior = BrokenAdaptiveRidge::createBarPrior(
                                 penalty = "bic", initialRidgeVariance = 1000, fitBestSubset = TRUE),
                               forceNewObject = TRUE)
coef(bar_fit_cox)[which(coef(bar_fit_cox) != 0)]


## Only vaccine groups
vaccine.groups <- c("BC.LSFcHN.EboGP", "BC.LSFHN.EboGP", "HPIV3..FHN.EboGP", "HPIV1.EboGP", "HPIV3.EboGP")

interaction.terms <- apply(expand.grid(vaccine.groups, vaccine.groups), 1, paste0, collapse=" * ") %>% paste0(., collapse=" + ")

groups_formula_cox <- as.formula(paste0(c("Surv(Time.of.death, Survival) ~ BC.LSFcHN.EboGP + BC.LSFHN.EboGP + HPIV3..FHN.EboGP + HPIV1.EboGP + HPIV3.EboGP", interaction.terms), collapse=" + "))

cyclopsDataCox <- createCyclopsData(groups_formula_cox, data = df.subset, modelType = "cox")

reg_fit_cox <- fitCyclopsModel(cyclopsDataCox, prior = createPrior("laplace", variance = 10), forceNewObject = TRUE)
coef.cox <- coef(reg_fit_cox)[which(coef(reg_fit_cox) != 0)]

ci.cox <- confint(reg_fit_cox, names(coef.cox), overrideNoRegularization = TRUE) %>% data.frame %>% rownames_to_column %>% mutate(coef=c(coef.cox))

pdf("plots/cox_coef_vg.pdf")
ci.cox %>% ggplot(aes(x=as.factor(rowname), y=coef)) + geom_point() + geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..),width=0.1) + geom_hline(yintercept=0) + coord_flip() + ylab("Coefficient") + xlab("Variable") + theme_bw()
dev.off()

## Checking confidence intervals
confint(reg_fit_cox, names(coef(reg_fit_cox)[which(coef(reg_fit_cox) != 0)]), overrideNoRegularization = TRUE)


## LR models
all_lr_formula <- as.formula(paste("Survival ~", paste0(
                                                                           names(df.subset) %>% str_subset("Survival|survival.index|Time|NHP|Vaccine.Group", negate=TRUE),
  collapse=" + ")))

cyclopsDataLr <- createCyclopsData(all_lr_formula, data = df.subset, modelType = "lr")

reg_fit_Lr <- fitCyclopsModel(cyclopsDataLr, prior = createPrior("laplace", variance = 10, exclude = "(Intercept)"), forceNewObject = TRUE)
## Checking confidence intervals

coef.lr <- coef(reg_fit_Lr)[which(coef(reg_fit_Lr) != 0)]
ci.lr <- confint(reg_fit_Lr, names(coef.lr), overrideNoRegularization = TRUE) %>% data.frame %>% rownames_to_column %>% mutate(coef=c(coef.lr))

pdf("plots/lr_coef.pdf")
ci.lr %>% ggplot(aes(x=as.factor(rowname), y=coef)) + geom_point() + geom_errorbar(aes(ymin=X2.5..,ymax=X97.5..),width=0.1) + geom_hline(yintercept=0) + coord_flip() + ylab("Coefficient") + xlab("Variable") + theme_bw()
dev.off()

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

## Lasso survival index
y <- df.subset$survival.index
x <- df.subset %>% select(-c(Survival,survival.index,Time.of.death,NHP,Vaccine.Group)) %>% select_if(is.numeric)
fit <- glmnet(x=as.matrix(x), y=as_vector(y), alpha = 1)

coef(fit, s = 0.008)

cvfit <- cv.glmnet(as.matrix(x), as_vector(y), type.measure = "mse", nfolds = 20)

rmse <- cvfit$cvm[which(cvfit$lambda == cvfit$lambda.min)]

coef.glmnet <- coef(cvfit, s = "lambda.min") %>% as.matrix %>% data.frame %>% rownames_to_column %>% filter(X1!=0)

pdf("plots/coef_glmnet.pdf")
coef.glmnet %>% ggplot(aes(x=as.factor(rowname), y=X1)) + geom_point() + geom_hline(yintercept=0) + coord_flip() + ylab("Coefficient") + xlab("Variable") + theme_bw()
dev.off()

plot(cvfit)


## Peptide 109.igg and 108.iga
pdf("plots/peptide_109.pdf")
df.subset %>% select(c(Survival, `Peptide_109.igg`)) %>% mutate(Survival = as.factor(Survival)) %>% ggplot(aes(Survival, `Peptide_109.igg`)) + geom_boxplot()
dev.off()

pdf("plots/peptide_108.pdf")
df.subset %>% select(c(Survival, `Peptide_108.iga`)) %>% mutate(Survival = as.factor(Survival)) %>% ggplot(aes(Survival, `Peptide_108.iga`)) + geom_boxplot()
dev.off()


df.subset %>% select(c(Survival, `Peptide_109.igg`)) %>% wilcox.test(`Peptide_109.igg` ~ Survival, .)

df.subset %>% select(c(Survival, `Peptide_108.iga`)) %>% wilcox.test(`Peptide_108.iga` ~ Survival, .)

