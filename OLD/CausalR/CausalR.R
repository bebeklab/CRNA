#########################################################################################
# Loading library
#########################################################################################
library(CausalR)

#########################################################################################
# Loading data
#########################################################################################
# omnipath <- read.csv("data_for_markdown/omnipath_final.csv",row.names = 1)
# DEGList <- read.csv("data_for_markdown/GSE11352_DEGList.csv",row.names = 1)

#########################################################################################
# Building graph
#########################################################################################

# causalsss <- omnipath[,c(1,3,2)]
# causalsss[which(causalsss$sign == 1),2] <- "Activates"
# causalsss[which(causalsss$sign == -1),2] <- "Inhibits"
# write.table(causalsss,"data_for_markdown/causalR_network.sif",quote = F,sep = " ")

# miarray <- as.data.frame(DEGList[which(abs(DEGList$logFC) >= log2(2)),])
# causal_exprs <- miarray[,c(1,2)]
# causal_exprs[which(causal_exprs$logFC > 0),2] <- 1
# causal_exprs[which(causal_exprs$logFC < 0),2] <- -1
# write.table(causal_exprs,"data_for_markdown/GSE1352_causal_exprs.txt",quote = F,sep = " ")


omniNet <- CreateCCG("data_for_markdown/CausalR_omnipath.sif")


GSE60880_miarray <- ReadExperimentalData("data_for_markdown/GSE60880_causal_exprs.txt",network = omniNet)

GSE11352_miarray <- ReadExperimentalData("data_for_markdown/GSE11352_causal_exprs.txt",network = omniNet)


runSCANR(network =omniNet,experimentalData =GSE11352_miarray, writeNetworkFiles = "all",outputDir = "CausalR/GSE11352/" )
runSCANR(network =omniNet,experimentalData =GSE60880_miarray, writeNetworkFiles = "all",outputDir = "CausalR/GSE60880/" )

