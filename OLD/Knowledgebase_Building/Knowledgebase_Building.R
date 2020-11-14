#########################################################################################
# Loading Library
#########################################################################################
library(igraph)
library(biomaRt)
library(igraph)
library(dplyr)
library(biomaRt)
library(data.table)
library(dnet)
library(tidyverse)
library(pheatmap)

#########################################################################################
# Loading data
#########################################################################################
omni_original <- read.delim("data_for_markdown/interactions_omnipathdb.txt")
STRING <- read.delim("data_for_markdown/9606.protein.links.v10.5.txt",sep = " ")

#########################################################################################
# Query annotation
#########################################################################################
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
uniprot2name <- getBM(attributes = c("uniprot_gn","hgnc_symbol"),mart = ensembl)
ENSP2name <- getBM(attributes = c("ensembl_peptide_id","hgnc_symbol"),mart = ensembl)


#########################################################################################
# Tidy Omnipath
#########################################################################################
# annotating omnipath's source name
omni_original<- merge(x = omni_original,y = uniprot2name, by.x = "source",by.y = "uniprot_gn")
colnames(omni_original) <- c("source_name", "target_name", "is_directed", "is_stimulation", "is_inhibition", "dip_url")

# annotating omnipath's target name
omni_original<- merge(x = omni_original,y = uniprot2name, by.x = "target_name",by.y = "uniprot_gn"
omni_original <- omni_original[,c(7,8,3:5)]
colnames(omni_original) <-c("source_name", "target_name", "is_directed", "is_stimulation", "is_inhibition")

# Only use directed edges
omni <- omni_original[which(omni_original$is_directed == 1),]

# Only use edges have single direction
omni <- omni[which(omni$is_stimulation + omni$is_inhibition == 1),]

# final clean
omni$sign <- omni$is_stimulation - omni$is_inhibition
omni$interaction.type <- "PPI"
omni <- omni[,-(4:5)]
rm(omni_original,uniprot2name)

#########################################################################################
# Tidy STRING
#########################################################################################
# get rid of mess
STRING$protein1 <- sub("9606.", "", STRING$protein1)
STRING$protein2 <-sub("9606.", "", STRING$protein2)

# Annotating first column
STRING<- merge(x = STRING,y = ENSP2name, by.x = "protein1",by.y = "ensembl_peptide_id")
STRING <- STRING[,c(4,2,3)]
colnames(STRING) <- c("protein1","protein2","combined_score")

# Annotating second column
STRING<- merge(x = STRING,y = ENSP2name, by.x = "protein2",by.y = "ensembl_peptide_id")
STRING <- STRING[,c(2,4,3)]
colnames(STRING) <- c("protein1","protein2","weight")

# calculate weight
STRING$weight <- STRING$weight/1000

omnipath <- merge(x = omni,y = STRING,by.x = c("source_name","target_name"),by.y = c("protein1","protein2"))
omnipath <- omnipath[-(which(omnipath$source_name == "" | omnipath$target_name == "")),]

# remove isolated loop
omnipath <- omnipath[-c(2965,10991,7044,7725,8114,10010,5055,5057),]
rm(STRING,omni,ENSP2name,ensembl)

#########################################################################################
# Mearging 2 dataset
#########################################################################################
omnipath <- merge(x = omni,y = STRING,by.x = c("source_name","target_name"),by.y = c("protein1","protein2"))
omnipath <- omnipath[-(which(omnipath$source_name == "" | omnipath$target_name == "")),]

rm(STRING,omni,ENSP2name,ensembl)

omniNet_messy <- graph_from_data_frame(d = omnipath,directed = T)
omniNet_messy <- igraph::simplify(graph = omniNet_messy,remove.multiple = T,
                                  remove.loops = T,
                                  edge.attr.comb = list(weight = function(x) mean(x),
                                                        sign = function(x) sign(sum(x)),"ignore"))
# delete isolated edges
All_edges <- E(omniNet_messy) 
k = 0
isolated_edges <- c()
for (i in 1:length(All_edges)) {
  edgeHead <- head_of(graph = omniNet_messy,es = All_edges[i])
  edgeTail <- tail_of(graph = omniNet_messy,es = All_edges[i])
  if (unname(degree(v = edgeHead, graph = omniNet_messy, mode = "all")) == 1 & unname(degree(v = edgeTail,graph = omniNet_messy,mode = "all")) == 1) {
    isolated_edges <- c(isolated_edges,i)
  }
}
# delete isolated edges
omniNet <- delete.edges(omniNet_messy,All_edges[isolated_edges])
# get all node name
nodeNames <- V(omniNet)$name

### Clean up
rm(omniNet_messy,edgeHead,edgeTail,All_edges)
