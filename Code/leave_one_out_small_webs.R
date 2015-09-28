library(vegan)
library(phytools)
library(geiger)
library(igraph)
library(magrittr)

source("RDPG_functions.R")

list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE,
  full.names=FALSE, ignore.case=FALSE) {
  # use full.names=TRUE to pass to file.info
  all <- list.files(path, pattern, all.dirs,
           full.names=TRUE, recursive=FALSE, ignore.case)
  dirs <- all[file.info(all)$isdir]
  # determine whether to return full names or just dir names
  if(isTRUE(full.names))
    return(dirs)
  else
    return(basename(dirs))
}

max_rank <- 15

for(dir_name in list.dirs()){

do_all(dir_name, max_rank)

}
