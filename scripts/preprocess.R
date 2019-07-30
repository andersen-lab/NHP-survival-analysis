library(tidyverse)
library(psych)
library(Rtsne)
library(RColorBrewer)

df <- read_csv("data/parameters_v3.csv")
names(df) <- make.names(names(df)) %>% str_replace_all("(\\.){2,5}", ".")
iga <- read_csv("data/IgA_peptide_array.csv")
igg <- read_csv("data/IgG_peptide_array.csv")

pep.arr <- inner_join(iga, igg, "NHP", suffix=c(".iga", ".igg"))
df.merged <- inner_join(pep.arr, df, "NHP")

## Exclude clinical attributes and scale
df.subset <- df.merged %>% select(-matches("temperture|kidney|liver|shedding|viremia|weight|clinical.score|time"))
df.subset <- df.subset %>% mutate(
                  Survival = ifelse(Survival == "Yes", 1, 0)
                  ) %>% mutate_if(is.numeric, scale)


getSignificantPairwiseCorrelations <- function(df){
    df <- corr.test(ief, method="kendall", adjust="holm")
    df.r <- df$r %>% data.frame %>% mutate(var1 = rownames(.)) %>% gather(var2, value, -var1)
    df.p  <- df$p %>% data.frame %>% mutate(var1 = rownames(.)) %>% gather(var2, value, -var1)
    df.merged <- inner_join(df.r, df.p, c("var1", "var2"), suffix=c(".r", ".p"))
    return(df.merged);
}

## Create polyfunctionality form Immune effectors functions
ief <- df.subset %>% select(matches("ADCD|ADCP|ADNKD|ADNP"))
## Check correlations between ief
ief.corr <- getSignificantPairwiseCorrelations(ief)
## Plot heatmap
ief.corr %>% filter(value.r >= 0.5 & var1!=var2) %>% ggplot(aes(var1, var2)) + geom_tile(aes(fill=value.r)) + geom_text(aes(label=round(value.r, 2))) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
## k means
set.seed(112358)
ief.kmeans <- kmeans(ief, 5)

## Plot clusters
set.seed(112358)
ief.tsne <- ief %>% Rtsne(perplexity=5)

ief.kmeans.tsne <- data.frame(
    x = ief.tsne$Y[,1],
    y = ief.tsne$Y[,2],
    cluster = as.factor(ief.kmeans$cluster)
) %>% cbind(., ief)

ggplot(ief.kmeans.tsne, aes(x,y, color = cluster)) + geom_point() + scale_colour_brewer(palette = "Set1")

pdf("plots/polyfunc_clusters_differences.pdf", w=10, h = 70)
ief.kmeans.tsne %>% gather(var, val, -cluster) %>% ggplot(aes(x=cluster, y=val)) + geom_boxplot() + facet_grid(var ~ .)
dev.off()
