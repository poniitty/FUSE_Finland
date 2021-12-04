library(tidyverse)
library(data.table)

# Trait names
trts <- fread("data/TRY_trait_ids.txt") %>% 
  select(-V6)

trts %>% filter(AccSpecNum >= 1000) %>% 
  mutate(cums = cumsum(ObsNum)) %>% 
  mutate(class = cut(cums, breaks = seq(0, 10000000, by = 2000000), labels = 1:5)) -> trts

write.table(paste(trts %>% filter(class == 1) %>% pull(TraitID), collapse = ", "),
            "data/selected_traits/traits_1.txt",col.names = FALSE, row.names = FALSE)
write.table(paste(trts %>% filter(class == 2) %>% pull(TraitID), collapse = ", "),
            "data/selected_traits/traits_2.txt",col.names = FALSE, row.names = FALSE)
write.table(paste(trts %>% filter(class == 3) %>% pull(TraitID), collapse = ", "),
            "data/selected_traits/traits_3.txt",col.names = FALSE, row.names = FALSE)
write.table(paste(trts %>% filter(class == 4) %>% pull(TraitID), collapse = ", "),
            "data/selected_traits/traits_4.txt",col.names = FALSE, row.names = FALSE)
write.table(paste(trts %>% filter(class == 5) %>% pull(TraitID), collapse = ", "),
            "data/selected_traits/traits_5.txt",col.names = FALSE, row.names = FALSE)


#
try <- fread("data/TRY_species_names.txt")
fnames <- fread("data/Finnish_vascular_plant_list.tsv", encoding = "UTF-8") %>% 
  filter(`Taksonominen taso` == "laji")
anames <- fread("data/NamesFromDatabasesFINAL.csv", encoding = "Latin-1")
anames[grepl("Isoetes",anames$final),]

snames <- unique(na.omit(c(fnames$`Tieteellinen nimi`, fnames$Synonyymit,
                 anames$OldName, anames$final)))

gnames <- sort(unique(unlist(lapply(snames, function(x) str_split(x, " ")[[1]][1]))))

try %>% mutate(genus = unlist(lapply(AccSpeciesName, function(x) str_split(x, " ")[[1]][1]))) %>% 
  filter(genus %in% gnames) -> try

