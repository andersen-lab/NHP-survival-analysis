## Time independent variables

library(tidyr)
library(stringr)
library(readr)
library(dplyr)
library(ggplot2)

df <- read_csv("../data/parameters.csv")
pep <- read_csv("../data/IgG_peptide_array.csv")

## Extract colnames with time specified
t <- grepl("d[0-9]+", names(df))
before.infection <- unique(str_extract(names(df)[t], ".+(?= d[0-9]+)"))
## Select only those before first infection point d56(?) i.e., less than 8 timepoints
t <- sapply(before.infection, function(x){
    return(sum(str_detect(names(df), paste(x,".*d[0-9]+",sep=""))) < 8)
})

before.infection.str <- paste(before.infection[t], collapse="|")

## Exclude for gather()
exclude.var <- names(df)[1:7]
exclude.var.str <- paste(names(df)[1:7], collapse = "|")
exclude.var.str <- str_replace_all(exclude.var.str, "\\(", "\\\\(") %>% str_replace_all(., "\\)", "\\\\)") #Escape ()

df.before.infection <- df %>% select_if((grepl(before.infection.str, names(.)) & grepl("d[0-9]+", names(.))) | grepl(exclude.var.str, names(.)) ) %>% gather(., variable, value, -exclude.var) %>% mutate(., time = str_extract(variable, "(?<=d)[0-9]+")) %>% mutate(., time=as.numeric(time)) %>% mutate(., variable = str_extract(variable, ".+(?= d[0-9]+)") %>% str_trim(.))

write_csv(df.before.infection, "../data/before_infection_parameters.csv")
