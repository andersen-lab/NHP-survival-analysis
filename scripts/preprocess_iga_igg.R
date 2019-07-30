library(tidyverse)
library(Rtsne)

iga <- read_csv("data/IgA_peptide_array.csv") %>% mutate_all(scale)
igg <- read_csv("data/IgG_peptide_array.csv") %>% mutate_all(scale)

set.seed(11258)
iga.clustersize <- c(1:20)
iga.withinss  <- iga.clustersize %>% map_dbl(function(x){
    k <- iga %>% kmeans(., x, nstart = 2)
    k$tot.withinss
})

plot(iga.withinss)
lines(iga.withinss)

set.seed(11258)
igg.clustersize <- c(1:15)
igg.withinss  <- igg.clustersize %>% map_dbl(function(x){
    k <- igg %>% kmeans(., x, nstart = 2)
    k$tot.withinss
})

plot(igg.withinss)
lines(igg.withinss)


