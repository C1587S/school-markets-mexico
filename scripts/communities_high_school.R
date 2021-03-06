options(warn=-1)
setwd("~/Github/school-markets-mexico/scripts")
source('../scripts/utils.R')
library(leaflet)
library(mapview)
library(htmltools)
library(crayon)

"
This script computes school markets at high school levels using 5 km and 10 km
buffers. 

Last update: 2021-06-05
"

#  buffer size
type_list <- c("5kms", "10kms")
# algorithms list
algo_list <- c("fg", "wt", "lp", "cl")  # "le"
#  Secundaria
df <- readRDS("../data/schools/agregado_dist_sec.rds")
# additionals
save <- TRUE
valid_buffers_comm <- c() # stores buffers with at least one connection

for (type in type_list){
  cat(yellow("***********************", "\n"))
  cat(yellow("*** Kilómetros", type,"***\n"))
  cat(yellow("***********************", "\n"))
  
  valid_buffers_comm <- c()
  
# Read buffer data
  if (type=="5kms"){
    dfbuffers <- readRDS("../data/buffers/df_buffers_5kms_sec.rds")$df_buffers  
  } else if(type=="10kms"){
    dfbuffers <- readRDS("../data/buffers/df_buffers_10kms_sec.rds")$df_buffers 
  }
  # buffer list
  buff_list <- dfbuffers %>% 
                  group_by(buffer) %>% 
                  summarise(n=n()) %>%
                  filter(n>100) %>%  
                  arrange(n) %>%  
                  select(buffer) %>% 
                  pull()
  # iterate over buffer sizes
  i <- 1
  for (buff in buff_list){
    cat(blue("***********************", "\n"))
    cat(blue("*** Buffer", buff,"***\n"))
    cat(blue("***********************", "\n"))
    
    current_group <- buff
    # nodes
    nodos <- dfbuffers %>% 
            filter(buffer==current_group) %>% 
            rename(grupo=buffer, name=cct, lat=latitud, lon=longitud)
    # relationships
    relations <- get_relations(nodos, df)
    
    select_relations <- get_select_relations(nodos, current_group, relations)
    select_nodos <- get_select_nodos(select_relations, current_group)
    # adds buffers with mobility
    if (nrow(select_relations)>0){ # filter if not enough migrations
      valid_buffers_comm[i] <-current_group
    } else {
      next
    }
    # compute network
    school_network <- graph_from_data_frame(select_relations, directed=FALSE,
                                            vertices=select_nodos)
    
    tbl_centrality <- get_centrality_stats(school_network, current_group, 
                                           save=save)
    
    # iterate over algorithms for community detection
    for (algo in algo_list){
      cat(red("***********************", "\n"))
      cat(red("*** Algorithm", algo,"***\n"))
      cat(red("***********************", "\n")) 
      
      algorithm <- algo
      fc <- c()
      if (algorithm=="fg"){
        try(fc <- cluster_fast_greedy(school_network) )
      } else if (algorithm=="wt"){
        try(fc <- cluster_walktrap(school_network))
      } else if (algorithm=="lp"){
        try(fc <- cluster_label_prop(school_network))
      } else if (algorithm=="le"){
        try(fc <- cluster_leading_eigen(school_network))
      } else if (algorithm=="cl"){
        try(fc <- cluster_louvain(school_network))
      }
      
      algorithm <- str_replace(algorithm(fc), " ", "_")
      select_nodos <- save_subgroups(fc, select_nodos, algorithm, current_group,
                                     save=save, buffer=current_group, type=type)
      
      # tables with community members
      mrkt_members_tbl <- tbls_mrkt_members(select_nodos, current_group,
                                            algorithm, save=save, type=type)
      # save cluster results
      folder <- str_c("../data/results/school_markets/primaria", type)
      dir.create(file.path(folder)) # create directory if not exists
      path_nwk <- str_c(folder, "/buffers_with_comms.rds")
      save_network(fc, type, current_group, algorithm)
    }
    i <- i + 1
  }
  # save matrix with community members
  path <- str_c(folder, "/buffers_with_comms.csv")
  write_csv(data.frame("buffers_with_communities"=valid_buffers_comm), path)
}
