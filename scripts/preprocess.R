library(tidyverse)

df <- read_csv("data/parameters_v3.csv")
names(df) <- make.names(names(df)) %>% str_replace_all("(\\.){2,5}", ".")

