# # You can change the simulation design in line 283.
# ##########################################################################################################################################
# # Graph Simulations
# ##########################################################################################################################################

#=========================================================================================#
# Set some global options
#=========================================================================================#
options(contrasts=c("contr.treatment", "contr.poly"))
options(max.print=1e6)
options(digits=12)
options(echo=TRUE)
#options(warn=-1)
#options(error=utils::recover)
# install.packages("Rcpp",repos = "http://cran.us.r-project.org")
library("Rcpp")
#=========================================================================================#
# Loading Package Libraries
#=========================================================================================#
library("parallel")
library("limma")
library("Matrix")
library("abind")
library("igraph")
library("glmnet")
library("testit")
library("pheatmap")
library("grid")
library("gridExtra")

library("magrittr")

#=========================================================================================#
# Loading Additional Libraries for Simulations
#=========================================================================================#
library("Biobase")
library("ROCR")
library("CausalR")
library("foreach")
library("doParallel")
library("doRNG")
library("CRNA")
# #==========================================================================================#
# # Set and check the working directory
# #==========================================================================================#
if ((.Platform$OS.type == "windows") || (.Platform$OS.type == "unix")) {
  HOME.path <- Sys.getenv("HOME")
  setwd(dir=file.path(HOME.path, "RESULTS/METHODS/CRNA/SIMS", fsep=.Platform$file.sep))
} else {
  warning("OS not recognized \n")
}
getwd()

# #==========================================================================================#
# # Set some path shortcuts
# #==========================================================================================#
CODE.path <- "CODES"
RESULT.path <- "RESULTS/METHODS/CRNA/REAL"
DATA.path <- "DATA"

#==========================================================================================#
#Saving the workspace
#==========================================================================================#
save.image(file=file.path(HOME.path, RESULT.path, "CRNA_SIMS_ascending_largekb_all_2408.RData", fsep=.Platform$file.sep))

#==========================================================================================#
#Load the workspace
#==========================================================================================#
load(file=file.path(HOME.path, RESULT.path, "CRNA_SIMS_ascending_largekb_all_2408.RData", fsep=.Platform$file.sep))

#==========================================================================================#
# set and check the working directory
#==========================================================================================#
if ((.Platform$OS.type == "windows") || (.Platform$OS.type == "unix")) {
  HOME.path <- Sys.getenv("HOME")
  setwd(dir=file.path(HOME.path, "RESULTS/METHODS/CRNA/SIMS", fsep=.Platform$file.sep))
} else {
  warning("OS not recognized \n")
}
getwd()

#==========================================================================================#
# Set some path shortcuts
#==========================================================================================#
CODE.path <- "CODES"
RESULT.path <- "RESULTS/METHODS/CRNA/REAL"
DATA.path <- "DATA"

#==========================================================================================#
# Source some R procedure files
#==========================================================================================#
source(file=file.path("sim_functions.r", fsep=.Platform$file.sep))

if (exists(".Random.seed")) rm(.Random.seed)
RNGkind(kind="L'Ecuyer-CMRG")
#==========================================================================================#
# Functions, constants, and parameters
#==========================================================================================#
fdrmin <- 0.01
fdropt <- 0.05
fdrmax <- 0.10

#=========================================================================================#
# Cluster configuration
#=========================================================================================#
if (require("parallel")) {
  print("'parallel' is attached correctly \n")
} else {
  stop("'parallel' must be attached first \n")
}

if (.Platform$OS.type == "windows") {
  # Windows OS
  cpus <- 39
  conf <- list("spec"=rep("localhost", cpus),
              "type"="SOCKET",
              "homo"=TRUE,
              "verbose"=TRUE,
              "outfile"=file.path(getwd(), "output_SOCK.txt", fsep=.Platform$file.sep))
} else if (.Platform$OS.type == "unix") {
  # Linux or Mac OS
  argv <- commandArgs(trailingOnly=TRUE)  
  if (is.empty(argv)) {
    # Mac OS : No argument "argv" in this case  
    cpus <- detectCores(logical = TRUE)
    conf <- list("spec"=rep("localhost", cpus),
                 "type"="SOCKET",
                 "homo"=TRUE,
                 "verbose"=TRUE,
                 "outfile"=file.path(getwd(), "output_SOCK.txt", fsep=.Platform$file.sep))
  } else {
    # Linux OS : Retrieving arguments "type" and "cpus" from the SLURM script
    type <- as.character(argv[1])
    cpus <- 39
    if (type == "SOCKET") {
      conf <- list("spec"=rep("localhost", cpus),
                  "type"="SOCKET",
                  "homo"=TRUE,
                  "verbose"=TRUE,
                  "outfile"=file.path(getwd(), "output_SOCK.txt", fsep=.Platform$file.sep))
    } else if (type == "MPI") {
      if (require("Rmpi")) {
          print("Rmpi is loaded correctly \n")
      } else {
          stop("Rmpi must be installed first to use MPI\n")
      }
      conf <- list("spec"=cpus,
                  "type"="MPI",
                  "homo"=TRUE,
                  "verbose"=TRUE,
                  "outfile"=file.path(getwd(), "output_MPI.txt", fsep=.Platform$file.sep))
    } else {
      stop("Unrecognized cluster type: you must specify a \"SOCKET\" or \"MPI\" cluster type\n")
    }
  }
} else {
  stop("Unrecognized platform: you must specify a \"windows\" or \"unix\" platform type\n")
}

cat("Cluster configuration:\n")
print(conf)

#=========================================================================================#
# Cluster configuration
#=========================================================================================#
#srun -N 1 -n 4 -p class --x11 --pty /bin/bash
#srun -N 1 -n 4 --x11 --pty /bin/bash
#
numCores <- 39   # nsizes below
registerDoParallel(numCores)

#==========================================================================================#
#Saving the workspace
#==========================================================================================#
save.image(file=file.path(HOME.path, RESULT.path, "CRNA_SIMS_ascending_largekb_all_2408.RData", fsep=.Platform$file.sep))

#==========================================================================================#
# Definition of the factors:
# Factor #1: Treatment Fator
#==========================================================================================#
TF <- factor(c(rep(1,50), rep(2,50)), ordered=FALSE, labels=c("C", "T"))
TF.lev <- levels(TF)
TF.ng <- nlevels(TF)
TF.def <- vector(mode="list", length=TF.ng)
TF.tab <- numeric(TF.ng)
for (g in 1:TF.ng) {
  TF.def[[g]] <- which(TF == TF.lev[g])
  TF.tab[g] <- length(TF.def[[g]])
}

rm(g)

# #==========================================================================================#
# # Sample Names
# #==========================================================================================#
sample.names <- c(paste(TF.lev[1], 1:50, sep=""), 
                  paste(TF.lev[2], 1:50, sep=""))

#==========================================================================================#
# Initializations - Constants - Parameters
#==========================================================================================#
seed <- 2408

nnodes <- c(5000, 10000, 15000, 20000)     # |V| initial nodes
nsizes <- length(nnodes)
ntdegs <- seq(from=30, by=30, length=nsizes) #30 5k 120 180 240 

#0.01 0.05 0.95 0.99
ndepths <- 4
depth.levels <- 1:ndepths

# # m = c(8,12)
# # p = c(0.5,0.8)

mu <- c(0.5, 1, 1.5, 2)
sigma2 <- c(0.5, 1, 1.5, 2)

R <- 100
B <- 50
FR <- 0.1
FC <- 2

alpha.set <- seq(0.10, 0.90, by=0.10)
alpha.set <- c(alpha.set, 0.01, 0.05, 0.95, 0.99)
fdr.set <- sort(c(10^seq(-2, -1, by=1), 5* 10^(-2)))

fc.set <- seq(from=1, to=3, by=0.5)
fr.set <- seq(from=0.0, to=0.2, by=0.02)
rep.set <- seq(from=10, to=100, by=10)

# #==========================================================================================#
# # Generation of Knowldge Base Graph
# #==========================================================================================#
set.seed(seed)

result <- foreach (vindex=1:nsizes, .options.RNG=2408) %dorng% {

  #---------------------------- Generating Initial Graph ------------------------#
  # Generating scale-free PPI graph according to the Barabasi-Albert model with |V| nodes 
  PPI.seed <- igraph::sample_pa(n=nnodes[vindex], directed=TRUE, m=12, zero.appeal=.3) %>%
    igraph::rewire(igraph::each_edge(p=0.8, loops=FALSE, multiple=TRUE))

  # randomly assign weights to edges
  igraph::E(PPI.seed)$weight <- sample(x=100, size=igraph::ecount(PPI.seed), replace=TRUE)/100

  # Filtering edges that have low weights to reduce connectivity
  PPI.seed.ig <- PPI.seed - igraph::edge(igraph::E(PPI.seed)[which(igraph::E(PPI.seed)$weight < 0.7)])

  # Edge sign for directed graph
  igraph::E(PPI.seed.ig)$sign <- 1

  # Vertex name
  igraph::V(PPI.seed.ig)$name <- as.character(1:igraph::vcount(PPI.seed.ig))

  # Removing singleton nodes
  PPI.seed.ig <- igraph::graph_from_data_frame(igraph::as_data_frame(PPI.seed.ig))

  #-------------------------- Generating Regulatory Graph -----------------------#
  # Randomly selecting true DEGs from previous graph
  tdegs <- as.character(sample(igraph::V(PPI.seed.ig)$name, ntdegs[vindex], FALSE))

  # Rebuilding the graph to include all the true DEGS to true DEGs paths including non true DEGs in between
  edges.set <- igraph::shortest_paths(PPI.seed.ig, from=tdegs[1], to=tdegs[-1], output="epath", mode="out")$epath[[1]]
  for (j in 1:length(tdegs)) {
    tdegs.sp <- igraph::shortest_paths(PPI.seed.ig, from=tdegs[j], to=tdegs[-j], output="epath", mode="out")$epath
    for (k in 1:length(tdegs.sp)) {
      edges.set <- c(edges.set, tdegs.sp[[k]])
    }
  }
  PPI.seed.ig <- igraph::subgraph.edges(PPI.seed.ig, eids=unique(edges.set), delete.vertices=TRUE)

  # Updating true DEGs
  tdegs <- as.character(intersect(igraph::V(PPI.seed.ig)$name, tdegs))

  # Simulated Regulatory Network and Expression values at each depth
  # For a given simulation design (ascending, descending, none)
  sim.obj <- sim(depth.levels=depth.levels,
                 tdegs=tdegs,
                 ppi.seed.ig=PPI.seed.ig,
                 meanlog=mu[1],
                 sdlog=sigma2[3],
                 factor=TF,
                 descent= FALSE,
                 seed=seed)

  X <- vector(mode="list", length=ndepths)
  TREGs <- vector(mode="list", length=ndepths)
  TDEGs <- vector(mode="list", length=ndepths)
  PPI.ig <- vector(mode="list", length=ndepths)
  for (depth in depth.levels) {
    X[[depth]] <- sim.obj$express[[depth]]
    TREGs[[depth]] <- sim.obj$tregs[[depth]]
    TDEGs[[depth]] <- sim.obj$tdegs[[depth]]
    PPI.ig[[depth]] <- sim.obj$ppi.ig[[depth]]
  }

  # Intersection of expression matrix genes and graph nodes at each depth
  # Updating final Expression Data matrix
  for (depth in depth.levels) {
    nodes_names <- igraph::V(PPI.ig[[depth]])$name
    inter_names <- intersect(x=rownames(X[[depth]]), y=nodes_names)
    X[[depth]] <- X[[depth]][inter_names,]
    colnames(X[[depth]]) <- sample.names
  }
  n <- ncol(X[[1]])
  p <- nrow(X[[1]])

  # Updating final Knowledge Base graph
  for (depth in depth.levels) {
    PPI.ig[[depth]] <- delete_vertices(graph=PPI.ig[[depth]], v=setdiff(x=nodes_names, y=inter_names))
  }
  p <- igraph::vcount(PPI.ig[[1]])
  e <- igraph::ecount(PPI.ig[[1]])

  size.set <- ceiling(seq(from=p/100, to=p, length=100))

  # Estimated simulated DEGs at each depth
  # For a given Simulation design (ascending, descending, fixed)
  # For a range of DE FDR values 
  DEGs <- vector(mode="list", length=ndepths)
  d <- vector(mode="list", length=ndepths)
  for (depth in depth.levels) {
    degfdr.obj <- deg.fdr(fit=deg(express=X[[depth]], factor=TF), fdr=fdr.set)
    DEGs[[depth]] <- vector(mode="list", length=length(fdr.set))
    d[[depth]] <- vector(mode="list", length=length(fdr.set))
    for (f in 1:length(fdr.set)) {
      names(degfdr.obj$d[[f]]) <- paste("DE FDR=", format(x=fdr.set[f], digits=3, scientific=TRUE), sep="")
      DEGs[[depth]][[f]] <- degfdr.obj$degs[[f]]
      d[[depth]][[f]] <- degfdr.obj$d[[f]]
    }
  }

  list("PPI.seed"=PPI.seed, "PPI.ig"=PPI.ig, 
      "X"=X, "TREGs"=TREGs, "TDEGs"=TDEGs, "DEGs"=DEGs, 
      "d"=d, "p"=p, "e"=e, "n"=n, "size.set"=size.set)

}
save.image(file=file.path(HOME.path, RESULT.path, "CRNA_SIMS_ascending_largekb_all_2408.RData", fsep=.Platform$file.sep))

for (i in 1:nsizes) {
  cat("\n KB size: ", nnodes[i], "\n\n")
  cat("|V|: ", igraph::vcount(result[[i]]$PPI.seed), "\n", sep="")
  cat("|E|: ", igraph::ecount(result[[i]]$PPI.seed), "\n", sep="")
  cat("p: ", result[[i]]$p, "\n", sep="")
  cat("e: ", result[[i]]$e, "\n", sep="")
  for (depth in depth.levels) {
    cat("\n Depth: ", depth, "\n", sep="")
    dist <- igraph::distances(
      graph = result[[i]]$PPI.ig[[depth]],
      v = igraph::V(result[[i]]$PPI.ig[[depth]])$name,
      to = igraph::V(result[[i]]$PPI.ig[[depth]])$name,
      mode = "out",
      weights = NULL)[1,]
    dist <- dist[!is.infinite(dist)]
    cat("Median degree (\\nu): ", median(igraph::degree(result[[i]]$PPI.ig[[depth]])), "\n", sep="")
    cat("Median path length (\\lambda): ", median(dist), "\n", sep="")
    cat("e/p: ", result[[i]]$e/result[[i]]$p, "\n", sep="")
  }
}

rm(i, j, k, f, n, p, e, d, depth, PPI.seed, PPI.ig, PPI.seed.ig, tdegs.sp, edges.set, tdegs, 
  X, TREGs, TDEGs, size.set, dist,
  sim.obj, inter_names, nodes_names, degfdr.obj)

# ==========================================================================================#
# Running cof complete analysis for each size of Knowledge Base
# For each DE FDR
# For each depth
# For a given Simulation design (ascending, descending, fixed)
# ==========================================================================================#
save.image(file=file.path(HOME.path, RESULT.path, "CRNA_SIMS_ascending_largekb_all_2408.RData", fsep=.Platform$file.sep))

alpha.cv <- vector(mode="list", length=nsizes)
size.cv <- vector(mode="list", length=nsizes)
amin <- vector(mode="list", length=nsizes)
amin.id <- vector(mode="list", length=nsizes)
qmin <- vector(mode="list", length=nsizes)
qmin.id <- vector(mode="list", length=nsizes)
CRNA.output <- vector(mode="list", length=nsizes)
IPA.output <- vector(mode="list", length=nsizes)
CAUSALR.output <- vector(mode="list", length=nsizes)

for (i in 1:nsizes) {
  cat("\n KB size: ", nnodes[i], "\n\n")

#==========================================================================================#
# Tuning of RWR parameters q (model size) and alpha (restart probability)
# For each DE FDR
# For each depth
# For a given Simulation design (ascending, descending, fixed)
#==========================================================================================#

#==============================================================================#
# Computation of OOB PE as a function of RWR parameters for each DE FDR

  alpha.cv[[i]] <- vector(mode="list", length=ndepths)
  size.cv[[i]] <- vector(mode="list", length=ndepths)

  for (depth in depth.levels) {
    cat("Depth: ", depth, "\n")
    crna.fit <- crna.tuning.sl(graph=result[[i]]$PPI.ig[[depth]],
                              express=result[[i]]$X[[depth]],
                              factors=list(TF),
                              B=B,
                              alpha=alpha.set,
                              size=result[[i]]$size.set,
                              fdr=fdr.set[1],
                              lag=2, span=0.40, degree=2, family="gaussian",
                              conf=conf,
                              parallel=TRUE,
                              seed=seed)
    alpha.cv[[i]][[depth]] <- crna.fit$alpha.cv
    size.cv[[i]][[depth]] <- crna.fit$size.cv
    save.image(file=file.path(HOME.path, RESULT.path, "CRNA_SIMS_ascending_largekb_all_2408.RData", fsep=.Platform$file.sep))
  }
  q <- 1   # fixed size qm of maximum variance

  amin[[i]] <- vector(mode="list", length=ndepths)
  amin.id[[i]] <- vector(mode="list", length=ndepths)
  qmin[[i]] <- vector(mode="list", length=ndepths)
  qmin.id[[i]] <- vector(mode="list", length=ndepths)
  for (depth in c(1,4)) {
    cat("Depth: ", depth, "\n")
    amin[[i]][[depth]] <- vector(mode="list", length=length(fdr.set))
    amin.id[[i]][[depth]] <- vector(mode="list", length=length(fdr.set))
    qmin[[i]][[depth]] <- vector(mode="list", length=length(fdr.set))
    qmin.id[[i]][[depth]] <- vector(mode="list", length=length(fdr.set))
    for (f in 1:1) {
      cat("DE FDR: ", fdr.set[f], "\n")
      amin[[i]][[depth]][[f]] <- alpha.cv[[i]][[depth]]$amin[q,f]
      amin.id[[i]][[depth]][[f]] <- alpha.cv[[i]][[depth]]$amin.id[q,f]
      qmin[[i]][[depth]][[f]] <- size.cv[[i]][[depth]]$qmin[amin.id[[i]][[depth]][[f]],f]
      qmin.id[[i]][[depth]][[f]] <- size.cv[[i]][[depth]]$qmin.id[amin.id[[i]][[depth]][[f]],f]
    }
  }
  
  #==========================================================================================#
  # CRNA
  # For fixed graph size q, fixed restart probability alpha, and fixed fold change
  # For fixed number of replications, and fixed frequency
  #==========================================================================================#
  CRNA.output[[i]] <- vector(mode="list", length=ndepths)
  save.image(file=file.path(HOME.path, RESULT.path, "CRNA_SIMS_ascending_largekb_all_2408.RData", fsep=.Platform$file.sep))
  for (depth in depth.levels) {
    cat("Depth: ", depth, "\n")
    CRNA.output[[i]][[depth]] <- vector(mode="list", length=length(fdr.set))
    for (f in 1:(length(fdr.set)-2)) {
      cat("DE FDR: ", fdr.set[f], "\n")
      CRNA.output[[i]][[depth]][[f]]  <- CRNA::crna.sl(graph=result[[i]]$PPI.ig[[depth]],
                                                express=result[[i]]$X[[depth]],
                                                factors=list(TF),
                                                alpha=0.7,
                                                size=400,
                                                FC=FC,
                                                FR=FR,
                                                R=R,
                                                fdr=fdr.set[f],
                                                conf=conf,
                                                parallel=TRUE,
                                                seed=seed)
    }
    save.image(file=file.path(HOME.path, RESULT.path, "CRNA_SIMS_ascending_largekb_all_2408.RData", fsep=.Platform$file.sep))
  }
  
  #==========================================================================================#
  # IPA
  # For fixed number of replications, and fixed frequency
  #==========================================================================================#
  IPA.output[[i]] <- vector(mode="list", length=ndepths)
  
  for (depth in depth.levels) {
    cat("Depth: ", depth, "\n")
    IPA.output[[i]][[depth]] <- vector(mode="list", length=length(fdr.set))
    for (f in 1:(length(fdr.set)-2)) {
      cat("DE FDR: ", fdr.set[f], "\n")
      IPA.output[[i]][[depth]][[f]] <- CRNA::ipa(graph=result[[i]]$PPI.ig[[depth]],
                                          express=result[[i]]$X[[depth]],
                                          factors=list(TF),
                                          FR=FR,
                                          R=R,
                                          fdr=fdr.set[f],
                                          conf=conf,
                                          parallel=TRUE,
                                          seed=seed)
    }
    save.image(file=file.path(HOME.path, RESULT.path, "CRNA_SIMS_ascending_largekb_all_2408.RData", fsep=.Platform$file.sep))
  }
  
  #==========================================================================================#
  # CAUSALR
  # For fixed number of replications, and fixed frequency
  #==========================================================================================#
  CAUSALR.output[[i]] <- vector(mode="list", length=ndepths)
  
  for (depth in depth.levels) {
    cat("Depth: ", depth, "\n")
    CAUSALR.output[[i]][[depth]] <- vector(mode="list", length=length(fdr.set))
    for (f in 1:(length(fdr.set)-2)) {
      cat("DE FDR: ", fdr.set[f], "\n")
      CAUSALR.output[[i]][[depth]][[f]] <- CRNA::causal(graph=result[[i]]$PPI.ig[[depth]],
                                                  express=result[[i]]$X[[depth]],
                                                  factors=list(TF),
                                                  FR=FR,
                                                  R=R,
                                                  fdr=fdr.set[f],
                                                  conf=conf,
                                                  parallel=TRUE,
                                                  seed=seed)
    }
    save.image(file=file.path(HOME.path, RESULT.path, "CRNA_SIMS_ascending_largekb_all_2408.RData", fsep=.Platform$file.sep))
 }
    

}

rm(q, f, i, depth, crna.fit)
save.image(file=file.path(HOME.path, RESULT.path, "CRNA_SIMS_ascending_largekb_all_2408.RData", fsep=.Platform$file.sep))

#==========================================================================================#
# Running of complete analysis for each size of Kwoledge Base
# Comparative Method Prediction Performances
# For fixed DE FDR, fixed number of replications, and fixed frequency
#==========================================================================================#
CRNA.output.perfo <- vector(mode="list", length=nsizes)
IPA.output.perfo <- vector(mode="list", length=nsizes)
CAUSALR.output.perfo <- vector(mode="list", length=nsizes)

for (i in 2:nsizes) {
  cat("\n KB size: ", nnodes[i], "\n\n")
  
  CRNA.output.perfo[[i]] <- vector(mode="list", length=ndepths)
  IPA.output.perfo[[i]] <- vector(mode="list", length=ndepths)
  CAUSALR.output.perfo[[i]] <- vector(mode="list", length=ndepths)
  
  for (depth in depth.levels) {
    cat("\n Depth: ", depth, "\n")
    CRNA.output.perfo[[i]][[depth]] <- vector(mode="list", length=length(fdr.set))
    IPA.output.perfo[[i]][[depth]] <- vector(mode="list", length=length(fdr.set))
    CAUSALR.output.perfo[[i]][[depth]] <- vector(mode="list", length=length(fdr.set))
    for (f in 1:(length(fdr.set)-2)) {
      cat("DE FDR: ", fdr.set[f], "\n")
      obj <- obspred(X=result[[i]]$X[[depth]], tdegs=result[[i]]$TDEGs[[depth]], tregs=result[[i]]$TREGs[[depth]], preds=CRNA.output[[i]][[depth]][[f]]$nodes)
      CRNA.output.perfo[[i]][[depth]][[f]] <- perfo.metrics(x.obs=obj$obs, x.pred=obj$pred, a=1)
      cat("CRNA AUC: ", CRNA.output.perfo[[i]][[depth]][[f]]$`auc`, "\n")
      cat("CRNA FDR: ", CRNA.output.perfo[[i]][[depth]][[f]]$`fdr`, "\n")
      obj <- obspred(X=result[[i]]$X[[depth]], tdegs=result[[i]]$TDEGs[[depth]], tregs=result[[i]]$TREGs[[depth]], preds=IPA.output[[i]][[depth]][[f]]$nodes)
      IPA.output.perfo[[i]][[depth]][[f]] <- perfo.metrics(x.obs=obj$obs, x.pred=obj$pred, a=1)
      cat("IPA AUC: ", IPA.output.perfo[[i]][[depth]][[f]]$`auc`, "\n")
      cat("IPA FDR: ", IPA.output.perfo[[i]][[depth]][[f]]$`fdr`, "\n")
      obj <- obspred(X=result[[i]]$X[[depth]], tdegs=result[[i]]$TDEGs[[depth]], tregs=result[[i]]$TREGs[[depth]], preds=CAUSALR.output[[i]][[depth]][[f]]$nodes)
      CAUSALR.output.perfo[[i]][[depth]][[f]] <- perfo.metrics(x.obs=obj$obs, x.pred=obj$pred, a=1)
      cat("CAUSALR AUC: ", CAUSALR.output.perfo[[i]][[depth]][[f]]$`auc`, "\n")
      cat("CAUSALR FDR: ", CAUSALR.output.perfo[[i]][[depth]][[f]]$`fdr`, "\n")
    }
  }
  
  for (depth in depth.levels) {
    cat("\n Depth: ", depth, "\n")
    for (f in 1:(length(fdr.set)-2)) {
      cat("DE FDR: ", fdr.set[f], "\n")
      
      auc.m <- rep(NA, R)
      fdr.m <- rep(NA, R)
      for (r in 1:R) {
        obj <- obspred(X=result[[i]]$X[[depth]], tdegs=result[[i]]$TDEGs[[depth]], tregs=result[[i]]$TREGs[[depth]], preds=CRNA.output[[i]][[depth]][[f]]$nodes.rep[[r]])
        obj.perfo <- perfo.metrics(x.obs=obj$obs, x.pred=obj$pred, a=1)
        auc.m[r] <- obj.perfo$`auc`
        fdr.m[r] <- obj.perfo$`fdr`
      }
      auc.se <- sd(auc.m, na.rm=TRUE)
      fdr.se <- sd(fdr.m, na.rm=TRUE)
      cat("CRNA AUC se: ", round(auc.se,2), "\n")
      cat("CRNA FDR se: ", round(fdr.se,2), "\n")

      auc.m <- rep(NA, R)
      fdr.m <- rep(NA, R)
      for (r in 1:R) {
        obj <- obspred(X=result[[i]]$X[[depth]], tdegs=result[[i]]$TDEGs[[depth]], tregs=result[[i]]$TREGs[[depth]], preds=IPA.output[[i]][[depth]][[f]]$nodes.rep[[r]])
        obj.perfo <- perfo.metrics(x.obs=obj$obs, x.pred=obj$pred, a=1)
        auc.m[r] <- obj.perfo$`auc`
        fdr.m[r] <- obj.perfo$`fdr`
      }
      auc.se <- sd(auc.m, na.rm=TRUE)
      fdr.se <- sd(fdr.m, na.rm=TRUE)
      cat("IPA AUC se: ", round(auc.se,2), "\n")
      cat("IPA FDR se: ", round(fdr.se,2), "\n")

      auc.m <- rep(NA, R)
      fdr.m <- rep(NA, R)
      for (r in 1:R) {
        obj <- obspred(X=result[[i]]$X[[depth]], tdegs=result[[i]]$TDEGs[[depth]], tregs=result[[i]]$TREGs[[depth]], preds=CAUSALR.output[[i]][[depth]][[f]]$nodes.rep[[r]])
        obj.perfo <- perfo.metrics(x.obs=obj$obs, x.pred=obj$pred, a=1)
        auc.m[r] <- obj.perfo$`auc`
        fdr.m[r] <- obj.perfo$`fdr`
      }
      auc.se <- sd(auc.m, na.rm=TRUE)
      fdr.se <- sd(fdr.m, na.rm=TRUE)
      cat("CAUSALR AUC se: ", round(auc.se,2), "\n")
      cat("CAUSALR FDR se: ", round(fdr.se,2), "\n")
    }
  }
}

rm(depth, f, r, i, obj, obj.perfo, 
   tpr.m, tnr.m, auc.m, fdr.m, for.m,
   tpr.se, tnr.se, auc.se, fdr.se, for.se)

#==========================================================================================#
# Save the workspace
#==========================================================================================#
save.image(file=file.path(HOME.path, RESULT.path, "CRNA_SIMS_ascending_largekb_all_2408.RData", fsep=.Platform$file.sep))