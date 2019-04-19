library(tidyverse)
library(ggplot2)
library(RColorBrewer)

library(Rtsne)
library(dummies)

d  <- read_csv("../data/parameters.csv")
peptide  <- read_csv("../data/IgG_peptide_array.csv")

## Select paramters that haven't been measured at 8 different time points
t <- grepl("d[0-9]+", names(d))
time.dep <- unique(str_extract(names(d)[t], ".+(?= d[0-9]+)"))

## Select only those with less than 8 timepoints
t <- sapply(time.dep, function(x){
    return(sum(str_detect(names(d), paste(x," .*d[0-9]+",sep="")))==8)
})

time.dep.str <- paste(time.dep[t], collapse="|")

d.time.indep <- d %>% select(-matches(time.dep.str))

## Combine peptide and time independent paramters
d.time.indep <- inner_join(d.time.indep, peptide, "NHP")

## Pairwise correlations
d.time.indep.cor <- d.time.indep %>% select_if(is.numeric) %>% select(-NHP) %>% data.matrix %>% cor(., method="kendall")

pdf("../plots/cor.pdf", w = 20, h = 20)
mypalette<-brewer.pal(7,"RdBu")
d.time.indep.cor %>% as_tibble %>% mutate(name=colnames(.), abs.sum=map_dbl(., function(x){
    sum(abs(x))
})) %>% gather(variable, value, -name, -abs.sum) %>% filter((value >= 0.8 | value <= -0.8) & (variable !=name)) %>% mutate(variable = fct_drop(variable), name=fct_drop(name)) %>% ggplot(aes(variable, name, fill=value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_gradientn(colors=mypalette)
dev.off()

## Peptide analysis
peptide.kmeans <- map(c(1:5), function(x){
    peptide %>% select(-NHP) %>% data.matrix %>% kmeans(., centers = x)
})

peptide.kmeans.withinss <- peptide.kmeans %>% map(function(x){
    x$tot.withinss
}) %>% as_vector

plot(peptide.kmeans.withinss)

## 3 clusters
peptide <- peptide %>% mutate(cluster = as.factor(peptide.kmeans[[3]]$cluster))

peptide.tsne <- peptide  %>% select(-NHP) %>% data.matrix %>% Rtsne(perplexity = 2)

peptide %>% mutate(x=peptide.tsne$Y[,1], y=peptide.tsne$Y[,2]) %>% ggplot(., aes(x, y, color=cluster)) + geom_point()

