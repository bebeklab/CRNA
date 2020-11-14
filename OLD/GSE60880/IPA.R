library(igraph)
library(dplyr)
library(data.table)
library(dnet)
library(tidyverse)
library(pheatmap)
library(predictionet)
library(ggplot2)

#########################################################################################
#Calculte p-value
#########################################################################################
p <- function(r,omniNet,omnipath,miarray) {
  if (r %in% V(omniNet)$name) {
    Vrg <- unique(omnipath$target_name)
    D <- unique(rownames(miarray))
    R <- unique(neighbors(omniNet,r,"out")$name)
    n <- length(Vrg)
    O <- unique(R[which(R %in% D)])
    a <- length(O)
    b <- length(unique(D[which(D %in% Vrg)])) - a
    c <- length(R) - a
    d <- n - a - b - c
    p.value <- 0
    
    for (k in 0:min(c,b)) {
      p.value <- p.value +
        fisher.test(as.table(matrix(base::c((a+k),(c-k),(b-k),(d+k)),
                                    nrow = 2)),alternative = "greater")$p.value} 
    return(p.value)}
  else {return(NA)}
}


#########################################################################################
#Calculte Z-score
#########################################################################################

z <- function(r,omniNet_directed,miarray)  {
  if(r %in% V(omniNet_directed)$name){
    D <- unique(rownames(miarray))
    OO <- unique(neighbors(omniNet_directed,r,"out")$name[(neighbors(omniNet_directed,r,"out")$name) %in% D])
    if (length(OO)>0) {
      numerator <- 0
      denominator <- 0
      for (i in 1:length(OO)) {
        v <- OO[i]
        wR <- max(E(omniNet_directed)[which(V(omniNet_directed)$name ==r) %->% which(V(omniNet_directed)$name == v )]$weight)
        sR <- max(E(omniNet_directed)[which(V(omniNet_directed)$name ==r) %->% which(V(omniNet_directed)$name == v )]$sign)
        sD <- sign(miarray[which(rownames(miarray) == v),"logFC"])
        numerator <- numerator + wR*sR*sD
        denominator <- denominator + wR^2
      }
      z.score <- as.numeric(numerator/sqrt(denominator))
      return(z.score)                 
    }
    else return(NA)}
  else return(NA)     
}


#########################################################################################
# p-value for 2 nodes in BFS
#########################################################################################

p1 <- function(r1,r2,omniNet,miarray) {
  if ( length(E(omniNet)[r1 %->% r2])>0) {
    D <- unique(rownames(miarray))   # DEGs
    R1 <- unique(neighbors(omniNet,r1,"out")$name[neighbors(omniNet,r1,"out")$name %in% D])  # regulated genes by r1 in network
    O1 <- unique(R1[which(R1 %in% D)])  # regulated genes by r1 in both DEGs and network
    R2 <- unique(neighbors(omniNet,r2,"out")$name[neighbors(omniNet,r2,"out")$name %in% D])  # regulated genes by r1 in network
    O2 <- unique(R2[which(R2 %in% D)])  # regulated genes by r1 in both DEGs and network
    Rrg <- c(R1,R2) # which genes are regulated by r1 or r2 in DEGs and network
    n <- length(Rrg) #universe
    a <- length(O1[O1 %in% O2]) # Overlap
    b <- length(O1) - a # only r1
    c <- length(O2) - a # only r2
    d <- n - a - b - c
    p.value <- 0
    for (k in 0:min(b,c)) {
      p.value <- p.value +  fisher.test(as.table(matrix(base::c((a+k),(c-k),(b-k),(d+k)),nrow = 2)),alternative = "greater")$p.value
    } 
    return(p.value)
  }else{return(NA)}
}


#########################################################################################
# BFS in IPA
#########################################################################################
ipa_bfs <- function(r1,omniNet,miarray,IPAresultPZ){
  passedNodes <- as.vector(r1)
  passedEdges <- c()
  IPAresultPZtemp_r1 <- IPAresultPZ
  ### Find r2
  possible_r2 <- unique(neighbors(omniNet,r1,
                                  "out")$name[neighbors(omniNet,r1,
                                                        "out")$name %in% IPAresultPZ$SYMBOL]) 
  
  possible_r2 <- possible_r2[!possible_r2 %in% passedNodes]      ### That's r2
  
  ###if r1 has childs(r2)
  if (length(possible_r2)>0) {  ### r1 has r2
    
    ###Order r2, select 3 most significant r2
    IPAresultPZtemp_r2 <- IPAresultPZ[which(IPAresultPZ$SYMBOL %in% possible_r2),]
    for (i in 1:length(possible_r2)) {
      r2 <- IPAresultPZtemp_r2[,"SYMBOL"]
      IPAresultPZtemp_r2[i,"p1"] <- p1(r1,r2,omniNet,miarray)
    }
    
    ### order r2
    IPAresultPZtemp_r2 <- IPAresultPZtemp_r2[order(IPAresultPZtemp_r2$p1),]
    
    ### if has more than 3 r2, choose top 3
    # if(nrow(IPAresultPZtemp_r2)>=3){ ### if has more than 3 r2, choose top 3
    # IPAresultPZtemp_r2 <- IPAresultPZtemp_r2[1:3,]}
    passedNodes <- c(passedNodes,IPAresultPZtemp_r2$SYMBOL) ### record r2
    for (x in 1:nrow(IPAresultPZtemp_r2)) {
      passedEdges <- c(passedEdges, get.edge.ids(omniNet,c(r1,IPAresultPZtemp_r2$SYMBOL[x])))
    }
    
    ### for each r2, find r3
    for (j in 1:nrow(IPAresultPZtemp_r2)) {
      r2 <- IPAresultPZtemp_r2$SYMBOL[j]
      possible_r3 <- neighbors(omniNet,IPAresultPZtemp_r2$SYMBOL[j],
                               "out")$name[neighbors(omniNet,IPAresultPZtemp_r2$SYMBOL[j],
                                                     "out")$name %in% IPAresultPZ$SYMBOL] 
      
      possible_r3 <- possible_r3[!possible_r3 %in% passedNodes]      ### That's r3
      
      ###if r3 exist
      if (length(possible_r3)>0) {  ### r2 has r3
        
        ###calculate p1 for r3
        IPAresultPZtemp_r3 <- IPAresultPZ[which(IPAresultPZ$SYMBOL %in% possible_r3),]
        for (k in 1:length(possible_r3)) {
          r3 <- IPAresultPZtemp_r3[k,"SYMBOL"]
          IPAresultPZtemp_r3[k,"p1"] <- p1(r2,r3,omniNet,miarray)
        }
        
        ### Order r3
        IPAresultPZtemp_r3 <- IPAresultPZtemp_r3[order(IPAresultPZtemp_r3$p1),]
        
        ###If has more than 3, select top 3
        # if(nrow(IPAresultPZtemp_r3)>=3){
        # IPAresultPZtemp_r3 <- IPAresultPZtemp_r3[1:3,]}
        
        ### record r3
        passedNodes <- c(passedNodes,IPAresultPZtemp_r3$SYMBOL) 
        for (x in 1:nrow(IPAresultPZtemp_r3)) {
          passedEdges <- c(passedEdges, get.edge.ids(omniNet,c(r2,IPAresultPZtemp_r3$SYMBOL[x])))
        }
        
        ### for each r3, find r4
        for (l in 1:nrow(IPAresultPZtemp_r3)) {
          r3 <- IPAresultPZtemp_r3$SYMBOL[l]
          possible_r4 <- neighbors(omniNet,IPAresultPZtemp_r3$SYMBOL[l],
                                   "out")$name[neighbors(omniNet,IPAresultPZtemp_r3$SYMBOL[l],
                                                         "out")$name %in% IPAresultPZ$SYMBOL] 
          
          possible_r4 <- possible_r4[!possible_r4 %in% passedNodes]      ### That's r4
          
          ###if r4 exist
          if (length(possible_r4)>0) {  ### r3 has r4
            
            ###calculate p1 for each r4
            IPAresultPZtemp_r4 <- IPAresultPZ[which(IPAresultPZ$SYMBOL %in% possible_r4),]
            for (m in 1:length(possible_r4)) {
              r4 <- IPAresultPZtemp_r4[m,"SYMBOL"]
              IPAresultPZtemp_r4[m,"p1"] <- p1(r3,r4,omniNet,miarray)
            }
            
            ### Order r4
            IPAresultPZtemp_r4 <- IPAresultPZtemp_r4[order(IPAresultPZtemp_r4$p1),]
            
            
            ###If has more than 3, select top 3
            # if(nrow(IPAresultPZtemp_r4)>=3){
            # IPAresultPZtemp_r4 <- IPAresultPZtemp_r4[1:3,]
            # }
            
            ### record r4
            passedNodes <- c(passedNodes,IPAresultPZtemp_r4$SYMBOL) 
            for (x in 1:nrow(IPAresultPZtemp_r4)) {
              passedEdges <- c(passedEdges, get.edge.ids(omniNet,c(r3,IPAresultPZtemp_r4$SYMBOL[x])))
            }
          }#if r4 exist 
        }# for each r3, find r4
      }#if r3 exist
    }# for each r2, find r3
  }# if r2 exist 
  return(subgraph.edges(omniNet,eids = passedEdges))
}



#########################################################################################
# Loading data
#########################################################################################
omnipath <- read.csv("data_for_markdown/omnipath_final.csv",row.names = 1)
DEGList <- read.csv("data_for_markdown/GSE60880_DEGList.csv",row.names = 1)
omniNet<- graph_from_data_frame(d = omnipath,directed = T)

#########################################################################################
# IPA p-value
#########################################################################################
miarray <- as.data.frame(DEGList[which(abs(DEGList$logFC) >= log2(2)),])
rownames(miarray) <- miarray$SYMBOL

IPAresult_FC2 <- as.data.frame(matrix(NA,ncol = 5,nrow = length(nodeNames)))
IPAresult_FC2[,1] <- nodeNames
colnames(IPAresult_FC2) <- c("SYMBOL","pValue","Zscore","mu","p1")

for (i in 1:length(nodeNames)) {
  IPAresult_FC2[i,2] <-  p(nodeNames[i],omniNet,omnipath,miarray)
}

## p-Value 0.05 ########################################################################################################################

IPAresult_FC2P0.05<- IPAresult_FC2[which(IPAresult_FC2$pValue <=0.05),]
print(paste("When using |FC| >= 2 as the threshold,  There are",nrow(IPAresult_FC2P0.05),
            "genes has significant overlapping p-value(0.05)"))


#########################################################################################
# IPA  Z-score
#########################################################################################

for (i in 1:nrow(IPAresult_FC2P0.05)) {
  IPAresult_FC2P0.05[i,3] <- z(IPAresult_FC2P0.05[i,1],omniNet,miarray)
}

ggplot(data = IPAresult_FC2P0.05, aes(x =Zscore)) +
  geom_histogram(color = "white", bins = 80,fill = "royalblue") +
  labs(title = "Distribution of IPA Activation Z-score",
       x = "Activation Z-score",y = "Frequency")

IPAresult_FC2P0.05 <- IPAresult_FC2P0.05[order(abs(IPAresult_FC2P0.05$Zscore),decreasing = T),]
IPAresult_FC2P0.05_Z2 <- IPAresult_FC2P0.05[which(abs(IPAresult_FC2P0.05$Zscore)>=2),]
IPAresult_FC2P0.05_Z1 <- IPAresult_FC2P0.05[which(abs(IPAresult_FC2P0.05$Zscore)>=1),]
IPAresult_FC2P0.05_Z0.5 <- IPAresult_FC2P0.05[which(abs(IPAresult_FC2P0.05$Zscore)>=0.5),]
IPAresult_FC2P0.05_Z0.25 <- IPAresult_FC2P0.05[which(abs(IPAresult_FC2P0.05$Zscore)>=0.25),]

print(paste("Among significant genes,",nrow(IPAresult_FC2P0.05_Z2),
            "genes' |Z| >= 2,",nrow(IPAresult_FC2P0.05_Z1),
            "genes' |Z| >= 1,",nrow(IPAresult_FC2P0.05_Z0.5),
            "genes' |Z| >= 0.5,",nrow(IPAresult_FC2P0.05_Z0.25),
            "genes' |Z| >= 0.25"))

# get FC2 P0.05 Z0.5 plot 
IPAresult_FC2P0.05_Z2_graph <-  induced.subgraph(omniNet, IPAresult_FC2P0.05_Z2$SYMBOL)
IPAresult_FC2P0.05_Z2_graph <- induced.subgraph(IPAresult_FC2P0.05_Z2_graph,
                                                  names(which(degree(IPAresult_FC2P0.05_Z2_graph,
                                                                     V(IPAresult_FC2P0.05_Z2_graph)) >0)))
IPAresult_FC2P0.05_Z2_edgeList <- igraph::as_data_frame(IPAresult_FC2P0.05_Z2_graph)

print(paste("When |FC| >= 2, p-value >= 0.05, Z >= 2,  There are",nrow(IPAresult_FC2P0.05_Z2_edgeList),
            "edges,",length(unique(c(IPAresult_FC2P0.05_Z2_edgeList$from,IPAresult_FC2P0.05_Z2_edgeList$to))),
            "nodes in the network"))

# plot the network

plot(IPAresult_FC2P0.05_Z2_graph, edge.arrow.size=1, edge.curved=0,

     vertex.color="orange", vertex.frame.color="#555555",

     vertex.label=V(IPAresult_FC2P0.05_Z2_graph)$name, vertex.label.color="black",

     vertex.label.cex=.5,vertex.size=7)
