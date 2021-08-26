options(warn=-1)
setwd("~/Github/school-markets-mexico/scripts")
source('./utils.R')

"
This script computes commuting zones at high school levels using 5 km, 8km, 
10 km, and 15 km buffers. 

Last update: 2021-03-12
"

# Load the data
db <- readRDS("../data/master_db_filtered.rds")

# Preparing the data for computing buffers

# coordinates
db$latitud <- as.numeric(db$latitud)
db$longitud <- as.numeric(db$longitud)

# schools with migrations data
dic_sec <- readRDS("../data/agregados/agregado_dist_sec.rds")
# unique schools using orig
unique_orig_sec <- dic_prim %>% select(cct_o) %>% unique() %>% rename(cct = cct_o)
# unique schools using dest
unique_dest_sec <- dic_prim %>% select(cct_d) %>% unique() %>% rename(cct = cct_d)
# db primaria
db_sec <- dplyr::left_join(unique_orig_sec, db) %>% na.omit()


# calculate buffers fr vaious radious and save to disk
for (buff in c(5000, 8000, 10000, 15000)){
  ## see docstrings in ../data/utils.R
  calc_buffers(db_sec, buff, save=TRUE)  
}