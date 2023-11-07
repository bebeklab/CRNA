#===============================================================================================================================#
#==============#
# Usage         :
#==============#
#                   sim(depth.levels,
#                       tdegs,
#                       ppi.seed.ig,
#                       meanlog=0.5,
#                       sdlog=0.5,
#                       factor,
#                       descent=TRUE,
#                       seed=NULL)
#
#==============#
# Description   :
#==============#
#                   Return the simulated regulatory network, and expression values at each depth
#
#==============#
# Arguments     :
#==============#
#
#==============#
# Values        :
#==============#
#
#===============================================================================================================================#

sim <- function(depth.levels,
                tdegs,
                ppi.seed.ig,
                meanlog=0.5,
                sdlog=0.5,
                factor,
                descent=TRUE,
                seed=NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Definition of the treatment factor
  factor.lev <- levels(factor)
  factor.ng <- nlevels(factor)
  factor.def <- vector(mode="list", length=factor.ng)
  factor.tab <- numeric(factor.ng)
  for (g in 1:factor.ng) {
    factor.def[[g]] <- which(factor == factor.lev[g])
    factor.tab[g] <- length(factor.def[[g]])
  }
  n <- length(factor)
  
  # Initializations
  tdegs.id <- tdegs
  max.depth <- max(depth.levels, na.rm=TRUE)
  depths <- length(depth.levels)
  incr <- round((depths+1)^(1/(depths+1)), 1)
  sample.names <- c(paste0(levels(factor)[1], 1:(n/2)), paste0(levels(factor)[2], 1:(n/2)))
  ppi.ig <- vector(mode="list", length=depths)
  ppi <- vector(mode="list", length=depths)
  express <- vector(mode="list", length=depths)
  degs <- vector(mode="list", length=depths)
  tregs <- vector(mode="list", length=depths)
  tdegs <- vector(mode="list", length=depths)
  
  for (depth in depth.levels) {
    
    each_level_regs <- vector(mode="list", length=depth)
    
    # Defining PPI.ig at each depth
    ppi.ig[[depth]] <- ppi.seed.ig
    
    # True DEGs nodes at each depth
    tdegs[[depth]] <- tdegs.id
    
    # True regulators (REGS) nodes at each depth
    tregs[[depth]] <- as.character(unique(c(tregs[[depth]], names(unlist(lapply(tdegs[[depth]], function(node) igraph::neighborhood(ppi.ig[[depth]], node=node, mode="in", order=depth, mindist=1)))))))
    tregs[[depth]] <- setdiff(tregs[[depth]], tdegs[[depth]])
    
    # Assigning expression values to true REGS nodes
    remaining_regs <- tregs[[depth]]
    
    if (!is.na(descent)) {
      regs_exprs <- matrix(data=NA, nrow=0, ncol=n)
      if (descent) {
        for (x in 1:depth) {
          each_level_regs[[x]] <- names(unlist(lapply(tdegs[[depth]], function(node) igraph::neighborhood(ppi.ig[[depth]], node=node, mode="in", order=x, mindist=x))))
          each_level_regs[[x]] <- intersect(remaining_regs, each_level_regs[[x]])
          remaining_regs <- setdiff(remaining_regs, each_level_regs[[x]])
          regs_exprs_part <- matrix(NA, nrow=length(each_level_regs[[x]]), ncol=n)
          rownames(regs_exprs_part) <- paste0("r_", each_level_regs[[x]])
          for (i in 1:(n/2)) {
            regs_exprs_part[,i] <-  rlnorm(length(each_level_regs[[x]]), meanlog=meanlog, sdlog=sdlog)
            regs_exprs_part[,i+(n/2)] <- rlnorm(length(each_level_regs[[x]]), meanlog=meanlog + incr*(depth - x + 1), sdlog=sdlog)
          }
          regs_exprs <- rbind(regs_exprs, regs_exprs_part)
        }
      } else {
        for (x in 1:depth) {
          each_level_regs[[x]] <- names(unlist(lapply(tdegs[[depth]], function(node) igraph::neighborhood(ppi.ig[[depth]], node=node, mode="in", order=x, mindist=x))))
          each_level_regs[[x]] <- intersect(remaining_regs, each_level_regs[[x]])
          remaining_regs <- setdiff(remaining_regs, each_level_regs[[x]])
          regs_exprs_part <- matrix(NA, nrow=length(each_level_regs[[x]]), ncol=n)
          rownames(regs_exprs_part) <- paste0("r_", each_level_regs[[x]])
          for (i in 1:(n/2)) {
            regs_exprs_part[,i] <- rlnorm(length(each_level_regs[[x]]), meanlog=meanlog, sdlog=sdlog)
            regs_exprs_part[,i+(n/2)] <- rlnorm(length(each_level_regs[[x]]), meanlog=meanlog + incr*x, sdlog=sdlog)
          }
          regs_exprs <- rbind(regs_exprs, regs_exprs_part)
        }
      }
    } else {
      regs_exprs <- matrix(data=NA, nrow=length(tregs[[depth]]), ncol=n)
      rownames(regs_exprs) <- paste0("r_", tregs[[depth]])
      for (i in 1:(n/2)) {
        regs_exprs[,i] <- rlnorm(length(tregs[[depth]]), meanlog=meanlog, sdlog=sdlog)
        regs_exprs[,i+(n/2)] <- rlnorm(length(tregs[[depth]]), meanlog=meanlog + incr, sdlog=sdlog)
      }
    }
    colnames(regs_exprs) <- sample.names
    
    # Assigning expression values to true TDEGs nodes
    tdegs_exprs <- matrix(data=NA, nrow=length(tdegs[[depth]]), ncol=n)
    for (i in 1:(n/2)) {
      tdegs_exprs[,i] <-       rlnorm(length(tdegs[[depth]]), meanlog=meanlog, sdlog=sdlog)
      tdegs_exprs[,i+(n/2)] <- rlnorm(length(tdegs[[depth]]), meanlog=meanlog + incr*(depth + 1), sdlog=sdlog)
    }
    colnames(tdegs_exprs) <- sample.names
    rownames(tdegs_exprs) <- paste0("d_", tdegs[[depth]])
    
    # Assigning unchanged expression values to true background nodes
    bkgs <- setdiff(as.character(V(ppi.ig[[depth]])$name), c(tregs[[depth]], tdegs[[depth]]))
    bkgs_exprs <- matrix(NA, nrow=length(bkgs), ncol=n)
    for (i in 1:(n/2)) {
      bkgs_exprs[,i]       <- rlnorm(length(bkgs), meanlog=meanlog, sdlog=sdlog)
      bkgs_exprs[,i+(n/2)] <- rlnorm(length(bkgs), meanlog=meanlog, sdlog=sdlog)
    }
    colnames(bkgs_exprs) <- sample.names
    rownames(bkgs_exprs) <- bkgs
    
    # Combining expression values at each depth
    express[[depth]] <- rbind(bkgs_exprs, tdegs_exprs, regs_exprs)
    
    # Updating true TDEGs and REGS' names in PPI.ig at each depth
    V(ppi.ig[[depth]])$name[which(as.character(V(ppi.ig[[depth]])$name) %in% tregs[[depth]])] <- c(paste0("r_", V(ppi.ig[[depth]])$name[which(as.character(V(ppi.ig[[depth]])$name) %in% tregs[[depth]])]))
    V(ppi.ig[[depth]])$name[which(as.character(V(ppi.ig[[depth]])$name) %in% tdegs[[depth]])] <- c(paste0("d_", V(ppi.ig[[depth]])$name[which(as.character(V(ppi.ig[[depth]])$name) %in% tdegs[[depth]])]))
    
    # Adjacency matrix at each depth
    ppi[[depth]] <- igraph::as_adjacency_matrix(ppi.ig[[depth]], attr="weight")
    
  }
  
  # Renaming the nodes at each depth
  for (depth in depth.levels) {
    tregs[[depth]] <- paste0("r_", tregs[[depth]])
    tdegs[[depth]] <- paste0("d_", tdegs[[depth]])
  }
  
  return(list("ppi.ig"=ppi.ig, "ppi"=ppi, "express"=express, "tregs"=tregs, "tdegs"=tdegs))
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    ln.param (mu, sigma2, slow="mu")
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

ln.param <- function (mu, sigma2, slow="mu") {
  
  if (slow=="mu") {
    
    comb <- expand.grid(sigma2, mu)
    ncomb <- nrow(comb)
    lnmean <- numeric(ncomb)
    lnvar <- numeric(ncomb)
    for (i in 1:ncomb) {
      lnmean[i] <- exp(comb[i,2] + comb[i,1]/2)
      lnvar[i] <- exp(comb[i,1]-1)*exp(2*comb[i,2]+comb[i,1])
    }
    
    return(cbind("mu"=comb[,2], "sigma2"=comb[,1], "mean"=lnmean, "var"=lnvar))
    
  } else {
    
    comb <- expand.grid(mu, sigma2)
    ncomb <- nrow(comb)
    lnmean <- numeric(ncomb)
    lnvar <- numeric(ncomb)
    for (i in 1:ncomb) {
      lnmean[i] <- exp(comb[i,1] + comb[i,2]/2)
      lnvar[i] <- (exp(comb[i,2])-1)*exp(2*comb[i,1]+comb[i,2])
    }
    
    return(cbind("mu"=comb[,1], "sigma2"=comb[,2], "mean"=lnmean, "var"=lnvar))
    
  }
  
}
#===============================================================================================================================#




#===============================================================================================================================#
#==============#
# Usage         :
#==============#
#                   obspred(X, tdegs, tregs, preds)
#
#==============#
# Description   :
#==============#
#                   Return the vectors of observed and predicted DEGS and REGs for use with the ROCR::prediction() function.
#                   Uses inputs from simulated data.
#
#==============#
# Arguments     :
#==============#
#
#==============#
# Values        :
#==============#
#
#===============================================================================================================================#

obspred <- function(X,
                    tdegs,
                    tregs,
                    preds) {
  
  sorted_nodes_names <- sort(rownames(X))
  regs.obs <- as.numeric(sorted_nodes_names %in% c(tdegs,tregs))
  w <- pmatch(x=sorted_nodes_names, table=preds)
  na <- is.na(w)
  w[na] <- 0
  w[!na] <- 1
  regs.pred <- w
  
  return(list("obs"=regs.obs, "pred"=regs.pred))
  
}
#===============================================================================================================================#




#===============================================================================================================================#
#==============#
# Usage         :
#==============#
#                   perfo.metrics(x.obs, x.pred, a=1)
#
#==============#
# Description   :
#==============#
#                   Prediction accuracy measures:
#                   pred <- prediction(predictions=predictions=x.pred, labels=x.obs)
#                   perf <- performance(prediction.obj=pred, x.measure="spec", measure="sens")
#                   spec <- perf@"x.values"[[1]][2]   # Same as tnr or 1-fpr
#                   sens <- perf@"y.values"[[1]][2]   # Same as tpr or 1-fnr
#                   perf <- performance(prediction.obj=pred, x.measure="tnr", measure="tpr")
#                   tnr <- perf@"x.values"[[1]][2]
#                   tpr <- perf@"y.values"[[1]][2]
#                   perf <- performance(prediction.obj=pred, x.measure="fpr", measure="fnr")
#                   fpr <- perf@"x.values"[[1]][2]
#                   fnr <- perf@"y.values"[[1]][2]
#
#==============#
# Arguments     :
#==============#
#
#==============#
# Values        :
#==============#
#
#===============================================================================================================================#

perfo.metrics <- function(x.obs,
                          x.pred,
                          a=1) {
  
  if (!any(is.na(x.obs)) && !any(is.na(x.obs)) && !is.empty(x.obs) && !any(is.na(x.pred)) && !any(is.na(x.pred)) && !is.empty(x.pred)) {
    pred <- prediction(predictions=x.pred, labels=x.obs)
    perf <- performance(prediction.obj=pred, measure="auc")
    auc <- perf@"y.values"[[1]][1]
    perf <- performance(pred,measure="f", alpha=a)
    f <- perf@"y.values"[[1]][1]
    perf <- performance(prediction.obj=pred, x.measure="tnr", measure="tpr")
    tnr <- perf@"x.values"[[1]][2]
    tpr <- perf@"y.values"[[1]][2]
    perf <- performance(prediction.obj=pred, x.measure="npv", measure="ppv")
    npv <- perf@"x.values"[[1]][2]
    ppv <- perf@"y.values"[[1]][2]
  } else {
    auc <- NA
    f <- NA
    tnr <- NA
    tpr <- NA
    npv <- NA
    ppv <- NA
  }
  
  return(list("auc"=auc, "f"=f, "tnr"=tnr, "tpr"=tpr, "fdr"=1-ppv, "for"=1-npv))
  
}
#===============================================================================================================================#




#===============================================================================================================================#
#==============#
# Usage         :
#==============#
#                   causal(graph,
#                          express,
#                          factors,
#                          FR=0.1,
#                          R=100,
#                          fdr=10^(-2),
#                          conf=NULL,
#                          parallel=FALSE,
#                          seed=NULL)
#
#==============#
# Description   :
#==============#
#                   causalR function at fixed FDR.
#
#==============#
# Arguments     :
#==============#
#
#==============#
# Values        :
#==============#
#
#===============================================================================================================================#

causal <- function (graph,
                    express,
                    factors,
                    FR=0.1,
                    R=100,
                    fdr=10^(-2),
                    conf=NULL,
                    parallel=FALSE,
                    seed=NULL) {
  
  seed <- seed[1]
  nboot <- 1:R
  
  factor1 <- factors[[1]]
  factor1.lev <- levels(factor1)
  factor1.ng <- nlevels(factor1)
  factor1.def <- vector(mode="list", length=factor1.ng)
  factor1.tab <- numeric(factor1.ng)
  for (g in 1:factor1.ng) {
    factor1.def[[g]] <- which(factor1 == factor1.lev[g])
    factor1.tab[g] <- length(factor1.def[[g]])
  }
  
  if (!parallel) {
    
    if (!is.null(seed)) {
      seed <- (0:(R-1)) + seed
    }
    #debug(causal.boot)
    causal.list <- causal.boot(graph=graph,
                               express=express,
                               factor=factor1,
                               nboot=nboot,
                               fdr=fdr,
                               parallel=parallel,
                               seed=seed)
    
  } else {
    
    # Parallel backend registration
    # To be used with the 'parallel' package:
    if (conf$type == "SOCKET") {
      clus <- parallel::makeCluster(spec=conf$spec,
                                    type="PSOCK",
                                    homogeneous=conf$homo,
                                    outfile=conf$outfile,
                                    verbose=conf$verbose)
    } else if (conf$type == "MPI") {
      clus <- parallel::makeCluster(spec=conf$spec,
                                    type="MPI",
                                    homogeneous=conf$homo,
                                    outfile=conf$outfile,
                                    verbose=conf$verbose)
    } else {
      stop("Unrecognized cluster type\n")
    }
    
    clusterEvalQ(cl=clus, expr=library("parallel"))
    clusterEvalQ(cl=clus, expr=library("limma"))
    clusterEvalQ(cl=clus, expr=library("igraph"))
    clusterEvalQ(cl=clus, expr=library("CausalR"))
    clusterEvalQ(cl=clus, expr=library("testit"))
    clusterExport(cl=clus, varlist=ls(.GlobalEnv), envir=.GlobalEnv)
    clusterExport(cl=clus,
                  varlist=c("causal.boot", "deg", "deg.fdr", "is.empty"),
                  envir=.GlobalEnv)
    
    # Bootstrap Cross-Validation of Prediction Error
    clusterSetRNGStream(cl=clus, iseed=seed)
    obj.clus <- clusterApply(cl=clus,
                             fun=causal.boot,
                             x=nboot,
                             graph=graph,
                             express=express,
                             factor=factor1,
                             fdr=fdr,
                             parallel=parallel,
                             seed=NULL)
    causal.list <- list("nodes"=vector(mode="list", length=R),
                        "edges"=vector(mode="list", length=R))
    for (r in 1:R) {
      causal.list$nodes[[r]] <- obj.clus[[r]]$nodes
      causal.list$edges[[r]] <- obj.clus[[r]]$edges
    }
    
    # Stopping the cluster and cleaning all MPI states
    parallel::stopCluster(cl=clus)
  }
  
  # Preparations - Initializations
  nodes.names <- character(0)
  nodes.rep <- vector(mode="list", length=R)
  for (r in 1:R) {
    nodes.rep[[r]] <- causal.list$nodes[[r]]
    nodes.names <- unique(c(nodes.names, nodes.rep[[r]]))
  }
  nodes.names <- nodes.names[!is.na(nodes.names)]
  n <- ncol(express)
  p <- nrow(express)
  p.adj <- length(nodes.names)
  
  # Aggregating bootstrap replications results
  nodes.matrix <- matrix(data=NA, nrow=p.adj, ncol=R, dimnames=list(nodes.names, 1:R))
  edges.matrix <- matrix(data=NA, nrow=0, ncol=3, dimnames=list(NULL,c("from", "sign", "to")))
  for (r in 1:R) {
    w <- pmatch(x=causal.list$nodes[[r]], table=nodes.names)
    if (!is.null(w)) {
      nodes.matrix[w,r] <- causal.list$nodes[[r]]
      edges.matrix <- rbind(edges.matrix, causal.list$edges[[r]])
    }
  }
  edges.matrix <- unique(edges.matrix)
  edges.matrix <- edges.matrix[,c("from", "to", "type")]
  levels(edges.matrix) <- c(1,-1)
  
  # Filtering
  nodes.freq <- apply(nodes.matrix, 1, function(x){R - sum(is.na(x))})  
  causal.nodes <- names(which(nodes.freq >= FR*R))
  w <- intersect(pmatch(x=causal.nodes, table=edges.matrix[,"from"]), 
                 pmatch(x=causal.nodes, table=edges.matrix[,"to"]))
  w <- w[!is.na(w)]
  causal.edges <- edges.matrix[w,]
  
  # Returning the final results
  return(list("nodes"=causal.nodes,
              "edges"=causal.edges,
              "nodes.rep"=nodes.rep))
}
#===============================================================================================================================#




#===============================================================================================================================#
#==============#
# Usage         :
#==============#
#                   causal.boot(graph,
#                               express,
#                               factor,
#                               nboot,
#                               fdr,
#                               parallel,
#                               seed)
#
#==============#
# Description   :
#==============#
#                   Bootstrap replication function of CausalR at fixed FDR.
#
#==============#
# Arguments     :
#==============#
#
#==============#
# Values        :
#==============#
#
#===============================================================================================================================#

causal.boot <- function (graph,
                         express,
                         factor,
                         nboot,
                         fdr,
                         parallel,
                         seed) {
  
  #==============================================================================#
  # Initializations - Constants - Parameters
  #==============================================================================#
  n <- ncol(express)
  p <- nrow(express)
  
  #==============================================================================#
  # Definition of the treatment factor
  #==============================================================================#
  factor.lev <- levels(factor)
  factor.ng <- nlevels(factor)
  factor.def <- vector(mode="list", length=factor.ng)
  factor.tab <- numeric(factor.ng)
  for (g in 1:factor.ng) {
    factor.def[[g]] <- which(factor == factor.lev[g])
    factor.tab[g] <- length(factor.def[[g]])
  }
  
  #==============================================================================#
  # Retrieving sample names
  #==============================================================================#
  sample.names <- colnames(express)
  
  #==============================================================================#
  # Simplify the graph
  #==============================================================================#
  graph <- graph %>% simplify(edge.attr.comb = "sum")
  E(graph)$sign <- sign(E(graph)$sign)
  
  #==============================================================================#
  # Intersection of expression matrix variables and graph nodes
  # Common graph size and dimensionality of the data
  #==============================================================================#
  nodes.names <- igraph::V(graph)$name
  inter_names <- intersect(x=rownames(express), y=nodes.names)
  express <- express[inter_names,]
  graph <- delete_vertices(graph=graph, v=setdiff(x=nodes.names, y=inter_names))
  
  #==============================================================================#
  # Bootstrap replication loop for a fixed FDR value
  #==============================================================================#
  if (!parallel) {
    
    R <- length(nboot)
    nodes.output <- vector(mode="list", length=R)
    edges.output <- vector(mode="list", length=R)
    r <- 1
    while (r <= R) {
      
      if (R > 1) {
        cat("Replication: ", r, "\n")
      }
      cat("seed: ", seed[r], "\n", sep="")
      if (!is.null(seed[r])) {
        set.seed(seed[r])
      }
      
      #==============================================================================#
      # Define the Bootstrap quantities
      if (R == 1) {
        # No Bootstrap resampling
        samples.boot <- c(factor.def[[1]], factor.def[[2]])
      } else {
        # Bootstrap resampling
        samples.boot <- c(sample(x=factor.def[[1]], size=factor.tab[[1]], replace=TRUE),
                          sample(x=factor.def[[2]], size=factor.tab[[2]], replace=TRUE))
        # Make sure all groups qualify
        while (length(intersect(samples.boot, factor.def[[2]])) < factor.tab[[2]]/2 |
               length(intersect(samples.boot, factor.def[[1]])) < factor.tab[[1]]/2 ) {
          samples.boot <- c(sample(x=factor.def[[1]], size=factor.tab[[1]], replace=TRUE),
                            sample(x=factor.def[[2]], size=factor.tab[[2]], replace=TRUE))
        }
      }
      # Bootstrap quantities
      express.boot <- express[,samples.boot,drop=FALSE]
      factor.boot <- factor[samples.boot]
      factor.boot.lev <- levels(factor.boot)
      factor.boot.ng <- nlevels(factor.boot)
      factor.boot.def <- vector(mode="list", length=factor.boot.ng)
      factor.boot.tab <- numeric(factor.boot.ng)
      for (g in 1:factor.boot.ng) {
        factor.boot.def[[g]] <- which(factor.boot == factor.boot.lev[g])
        factor.boot.tab[g] <- length(factor.boot.def[[g]])
      }
      sample.names.boot <- sample.names[samples.boot]
      
      #==============================================================================#
      # Differential expression analysis to get observed train set DEGs
      fit.boot <- deg(express=express.boot, factor=factor.boot)
      
      # Table of FDR-controlled top ranked train set DEGs
      top.boot <- limma::topTable(fit=fit.boot, coef=1, adjust="BH", number=p, p.value=fdr, sort.by="p", resort.by="B")
      
      #==============================================================================#
      # CausalR
      
      # Creating result folder for each bootstrap
      dir.create(path=file.path(HOME.path, RESULT.path, paste0("CausalR", r), fsep=.Platform$file.sep),
                 showWarnings=TRUE,
                 recursive=TRUE)
      
      # Generating CausalR inputs
      causalR_edges <- igraph::as_data_frame(graph)[,c("from", "sign", "to")]
      causalR_edges[which(causalR_edges$sign == 1),2] <- "Activates"
      causalR_edges[which(causalR_edges$sign == -1),2] <- "Inhibits"
      write.table(causalR_edges, file.path(HOME.path, RESULT.path, paste0("CausalR", r, "/CausalR_knowledgebase.sif"), fsep=.Platform$file.sep), quote=FALSE, sep=" ")
      causal_exprs <- as.data.frame(cbind(as.character(rownames(top.boot)), sign(as.numeric(top.boot[,"logFC"]))))
      write.table(causal_exprs, file.path(HOME.path, RESULT.path, paste0("CausalR", r, "/DEG_list.txt"), fsep=.Platform$file.sep), quote=FALSE, sep=" ")
      
      # Loading CausalR network
      simulation_network <- CausalR::CreateCCG(filename=file.path(HOME.path, RESULT.path, paste0("CausalR", r, "/CausalR_knowledgebase.sif"), fsep=.Platform$file.sep))
      
      # Loading CausalR DEGs
      simulation_DEGs <- CausalR::ReadExperimentalData(fileName=file.path(HOME.path, RESULT.path, paste0("CausalR", r, "/DEG_list.txt"), fsep=.Platform$file.sep),
                                                       network=simulation_network)
      # Running CausalR
      invisible(capture.output(CausalR::runSCANR(network=simulation_network,
                                                 experimentalData=simulation_DEGs,
                                                 writeNetworkFiles="all",
                                                 outputDir=file.path(HOME.path, RESULT.path, paste0("CausalR", r), fsep=.Platform$file.sep))))
      
      # Combining all CausalR results into a single file
      obj <- as.numeric(system(paste0("ls ", eval(file.path(HOME.path, RESULT.path, paste0("CausalR", r), fsep=.Platform$file.sep)), " | grep ^cor* | wc -l"), intern = TRUE))
      if (obj != 0) {
        # Recording result
        tryCatch({
          system(paste0("cat ", file.path(HOME.path, RESULT.path, paste0("CausalR", r, "/cor*.sif"), fsep=.Platform$file.sep), " > ", file.path(HOME.path, RESULT.path, paste0("CausalR", r, "/CausalR_cor_combined_network.sif"), fsep=.Platform$file.sep)))
          if (testit::has_error(read.table(file.path(HOME.path, RESULT.path, paste0("CausalR", r, "/CausalR_cor_combined_network.sif"), fsep=.Platform$file.sep)))) {
            # If returned CausalR files are empty, collect all file names and grep the gene names from file names
            system(paste0("ls ",
                          eval(file.path(HOME.path, RESULT.path, paste0("CausalR", r, "/corExplainedNodes-CausalR_knowledgebase-DEG_list-delta5-top150-*.sif"), fsep=.Platform$file.sep)),
                          "| sed 's/corExplainedNodes-CausalR_knowledgebase-DEG_list-delta5-top150-//g'",
                          "| sed 's/+.sif//g' "," > ",
                          eval(file.path(HOME.path, RESULT.path, paste0("CausalR", r, "/list_nodes"), fsep=.Platform$file.sep))
            )
            )
            # Read in gene name and record
            causalR_significant_nodes_result <- read.table(file.path(HOME.path, RESULT.path, paste0("CausalR", r, "/list_nodes"), fsep=.Platform$file.sep))
            colnames(causalR_significant_nodes_result) <- c("from", "type", "to")
            tmp <- sub(pattern=file.path(HOME.path, RESULT.path, paste0("CausalR", r, "/"), fsep=.Platform$file.sep),
                       replacement="",
                       causalR_significant_nodes_result[,"from"])
            if (is.empty(tmp)) {
              nodes.output[[r]] <- NA
            } else {
              nodes.output[[r]] <- tmp
            }
            edges.output[[r]] <- matrix(data=NA, nrow=0, ncol=3, dimnames=list(NULL,c("from", "type", "to")))
          } else {
            # If returned CausalR files are not empty, use files content
            causalR_significant_nodes_result <- read.table(file.path(HOME.path, RESULT.path, paste0("CausalR", r, "/CausalR_cor_combined_network.sif"), fsep=.Platform$file.sep))
            colnames(causalR_significant_nodes_result) <- c("from", "type", "to")
            tmp <- unique(c(as.character(causalR_significant_nodes_result[,"from"]), as.character(causalR_significant_nodes_result[,"to"])))
            nodes.output[[r]] <- tmp
            edges.output[[r]] <- causalR_significant_nodes_result
          }
        })
      } else {
        # If returned CausalR files are nonexistent, collect all file names and grep the gene names from file names
        nodes.output[[r]] <- NA
        edges.output[[r]] <- matrix(data=NA, nrow=0, ncol=3, dimnames=list(NULL,c("from", "type", "to")))
      }
      
      # Cleanup
      system(paste0("rm -rf ", file.path(HOME.path, RESULT.path, paste0("CausalR", r), fsep=.Platform$file.sep)))
      
      r <- r + 1
      
    }
    
  } else {
    
    #==============================================================================#
    # Get the process ID of each slave R session
    pid <- Sys.getpid()
    
    #==============================================================================#
    # Define the Bootstrap quantities
    if (R == 1) {
      # No Bootstrap resampling
      samples.boot <- c(factor.def[[1]], factor.def[[2]])
    } else {
      # Bootstrap resampling
      samples.boot <- c(sample(x=factor.def[[1]], size=factor.tab[[1]], replace=TRUE),
                        sample(x=factor.def[[2]], size=factor.tab[[2]], replace=TRUE))
      # Make sure all groups qualify
      while (length(intersect(samples.boot, factor.def[[2]])) < factor.tab[[2]]/2 |
             length(intersect(samples.boot, factor.def[[1]])) < factor.tab[[1]]/2 ) {
        samples.boot <- c(sample(x=factor.def[[1]], size=factor.tab[[1]], replace=TRUE),
                          sample(x=factor.def[[2]], size=factor.tab[[2]], replace=TRUE))
      }
    }
    # Bootstrap quantities
    express.boot <- express[,samples.boot,drop=FALSE]
    factor.boot <- factor[samples.boot]
    factor.boot.lev <- levels(factor.boot)
    factor.boot.ng <- nlevels(factor.boot)
    factor.boot.def <- vector(mode="list", length=factor.boot.ng)
    factor.boot.tab <- numeric(factor.boot.ng)
    for (g in 1:factor.boot.ng) {
      factor.boot.def[[g]] <- which(factor.boot == factor.boot.lev[g])
      factor.boot.tab[g] <- length(factor.boot.def[[g]])
    }
    sample.names.boot <- sample.names[samples.boot]
    
    #==============================================================================#
    # Differential expression analysis to get observed train set DEGs
    fit.boot <- deg(express=express.boot, factor=factor.boot)
    
    # Table of FDR-controlled top ranked train set DEGs
    top.boot <- limma::topTable(fit=fit.boot, coef=1, adjust="BH", number=p, p.value=fdr, sort.by="p", resort.by="B")
    
    #==============================================================================#
    # CausalR
    
    # Creating result folder for each bootstrap
    dir.create(path=file.path(HOME.path, RESULT.path, paste0("CausalR", pid), fsep=.Platform$file.sep),
               showWarnings=TRUE,
               recursive=TRUE)
    
    # Generating CausalR inputs
    causalR_edges <- igraph::as_data_frame(graph)[,c("from", "sign", "to")]
    causalR_edges[which(causalR_edges$sign == 1),2] <- "Activates"
    causalR_edges[which(causalR_edges$sign == -1),2] <- "Inhibits"
    write.table(causalR_edges, file.path(HOME.path, RESULT.path, paste0("CausalR", pid, "/CausalR_knowledgebase.sif"), fsep=.Platform$file.sep), quote=FALSE, sep=" ")
    causal_exprs <- as.data.frame(cbind(as.character(rownames(top.boot)), sign(as.numeric(top.boot[,"logFC"]))))
    write.table(causal_exprs, file.path(HOME.path, RESULT.path, paste0("CausalR", pid, "/DEG_list.txt"), fsep=.Platform$file.sep), quote=FALSE, sep=" ")
    
    # Loading CausalR network
    simulation_network <- CausalR::CreateCCG(filename=file.path(HOME.path, RESULT.path, paste0("CausalR", pid, "/CausalR_knowledgebase.sif"), fsep=.Platform$file.sep))
    
    # Loading CausalR DEGs
    simulation_DEGs <- CausalR::ReadExperimentalData(fileName=file.path(HOME.path, RESULT.path, paste0("CausalR", pid, "/DEG_list.txt"), fsep=.Platform$file.sep),
                                                     network=simulation_network)
    # Running CausalR
    invisible(capture.output(CausalR::runSCANR(network=simulation_network,
                                               experimentalData=simulation_DEGs,
                                               writeNetworkFiles="all",
                                               outputDir=file.path(HOME.path, RESULT.path, paste0("CausalR", pid), fsep=.Platform$file.sep))))
    
    # Combining all CausalR results into a single file
    obj <- as.numeric(system(paste0("ls ", eval(file.path(HOME.path, RESULT.path, paste0("CausalR", pid), fsep=.Platform$file.sep)), " | grep ^cor* | wc -l"), intern = TRUE))
    if (obj != 0) {
      # Recording result
      tryCatch({
        system(paste0("cat ", file.path(HOME.path, RESULT.path, paste0("CausalR", pid, "/cor*.sif"), fsep=.Platform$file.sep), " > ", file.path(HOME.path, RESULT.path, paste0("CausalR", pid, "/CausalR_cor_combined_network.sif"), fsep=.Platform$file.sep)))
        if (testit::has_error(read.table(file.path(HOME.path, RESULT.path, paste0("CausalR", pid, "/CausalR_cor_combined_network.sif"), fsep=.Platform$file.sep)))) {
          # If returned CausalR files are empty, collect all file names and grep the gene names from file names
          system(paste0("ls ",
                        eval(file.path(HOME.path, RESULT.path, paste0("CausalR", pid, "/corExplainedNodes-CausalR_knowledgebase-DEG_list-delta5-top150-*.sif"), fsep=.Platform$file.sep)),
                        "| sed 's/corExplainedNodes-CausalR_knowledgebase-DEG_list-delta5-top150-//g'",
                        "| sed 's/+.sif//g' "," > ",
                        eval(file.path(HOME.path, RESULT.path, paste0("CausalR", pid, "/list_nodes"), fsep=.Platform$file.sep))))
          # Read in gene name and record
          causalR_significant_nodes_result <- read.table(file.path(HOME.path, RESULT.path, paste0("CausalR", pid, "/list_nodes"), fsep=.Platform$file.sep))
          colnames(causalR_significant_nodes_result) <- c("from", "type", "to")
          tmp <- sub(pattern=file.path(HOME.path, RESULT.path, paste0("CausalR", pid, "/"), fsep=.Platform$file.sep),
                     replacement="",
                     causalR_significant_nodes_result[,"from"])
          if (is.empty(tmp)) {
            nodes.output <- NA
          } else {
            nodes.output <- tmp
          }
          edges.output <- matrix(data=NA, nrow=0, ncol=3, dimnames=list(NULL,c("from", "type", "to")))
        } else {
          # If returned CausalR files are not empty, use files content
          causalR_significant_nodes_result <- read.table(file.path(HOME.path, RESULT.path, paste0("CausalR", pid, "/CausalR_cor_combined_network.sif"), fsep=.Platform$file.sep))
          colnames(causalR_significant_nodes_result) <- c("from", "type", "to")
          tmp <- unique(c(as.character(causalR_significant_nodes_result[,"from"]), as.character(causalR_significant_nodes_result[,"to"])))
          nodes.output <- tmp
          edges.output <- causalR_significant_nodes_result
        }
      })
    } else {
      # If returned CausalR files are nonexistent, collect all file names and grep the gene names from file names
      nodes.output <- NA
      edges.output <- matrix(data=NA, nrow=0, ncol=3, dimnames=list(NULL,c("from", "type", "to")))
    }
    
    # Cleanup
    system(paste0("rm -rf ", file.path(HOME.path, RESULT.path, paste0("CausalR", pid), fsep=.Platform$file.sep)))
    
  }
  
  #==============================================================================#
  # Returning replication results
  #==============================================================================#
  return(list("nodes"=nodes.output,
              "edges"=edges.output))
  
}
#===============================================================================================================================#




#===============================================================================================================================#
#==============#
# Usage         :
#==============#
#                   ipa(graph,
#                       express,
#                       factors,
#                       FR=0.1,
#                       R=100,
#                       fdr=10^(-2),
#                       conf=NULL,
#                       parallel=FALSE,
#                       seed=NULL)
#
#==============#
# Description   :
#==============#
#                   IPA function at fixed FDR.
#
#==============#
# Arguments     :
#==============#
#
#==============#
# Values        :
#==============#
#
#===============================================================================================================================#

ipa <- function (graph,
                 express,
                 factors,
                 FR=0.1,
                 R=100,
                 fdr=10^(-2),
                 conf=NULL,
                 parallel=FALSE,
                 seed=NULL) {
  
  seed <- seed[1]
  nboot <- 1:R
  
  factor1 <- factors[[1]]
  factor1.lev <- levels(factor1)
  factor1.ng <- nlevels(factor1)
  factor1.def <- vector(mode="list", length=factor1.ng)
  factor1.tab <- numeric(factor1.ng)
  for (g in 1:factor1.ng) {
    factor1.def[[g]] <- which(factor1 == factor1.lev[g])
    factor1.tab[g] <- length(factor1.def[[g]])
  }
  
  if (!parallel) {
    
    if (!is.null(seed)) {
      seed <- (0:(R-1)) + seed
    }
    
    ipa.list <- ipa.boot(graph=graph,
                         express=express,
                         factor=factor1,
                         nboot=nboot,
                         fdr=fdr,
                         parallel=parallel,
                         seed=seed)
    
  } else {
    
    # Parallel backend registration
    # To be used with the 'parallel' package:
    if (conf$type == "SOCKET") {
      clus <- parallel::makeCluster(spec=conf$spec,
                                    type="PSOCK",
                                    homogeneous=conf$homo,
                                    outfile=conf$outfile,
                                    verbose=conf$verbose)
    } else if (conf$type == "MPI") {
      clus <- parallel::makeCluster(spec=conf$spec,
                                    type="MPI",
                                    homogeneous=conf$homo,
                                    outfile=conf$outfile,
                                    verbose=conf$verbose)
    } else {
      stop("Unrecognized cluster type\n")
    }
    
    clusterEvalQ(cl=clus, expr=library("parallel"))
    clusterEvalQ(cl=clus, expr=library("limma"))
    clusterEvalQ(cl=clus, expr=library("igraph"))
    clusterExport(cl=clus, varlist=ls(.GlobalEnv), envir=.GlobalEnv)
    clusterExport(cl=clus,
                  varlist=c("ipa.boot", "deg", "deg.fdr", "ipa.pvalue", "is.empty"),
                  envir=.GlobalEnv)
    
    # Bootstrap Cross-Validation of Prediction Error
    clusterSetRNGStream(cl=clus, iseed=seed)
    ipa.list <- clusterApply(cl=clus,
                             fun=ipa.boot,
                             x=nboot,
                             graph=graph,
                             express=express,
                             factor=factor1,
                             fdr=fdr,
                             parallel=parallel,
                             seed=NULL)
    
    # Stopping the cluster and cleaning all MPI states
    parallel::stopCluster(cl=clus)
  }
  
  # Preparations - Initializations
  nodes.names <- character(0)
  for (r in 1:R) {
    nodes.names <- unique(c(nodes.names, ipa.list[[r]]))
  }
  nodes.names <- nodes.names[!is.na(nodes.names)]
  n <- ncol(express)
  p <- nrow(express)
  p.adj <- length(nodes.names)
  
  # Aggregating bootstrap replications results
  nodes.matrix <- matrix(data=NA, nrow=p.adj, ncol=R, dimnames=list(nodes.names, 1:R))
  for (r in 1:R) {
    w <- pmatch(x=ipa.list[[r]], table=nodes.names)
    if (!is.null(w)) {
      nodes.matrix[w,r] <- ipa.list[[r]]
    }
  }
  
  # Filtering
  nodes.freq <- apply(nodes.matrix, 1, function(x){R - sum(is.na(x))})  
  ipa.nodes <- names(which(nodes.freq >= FR*R))
  ipa.edges <- unique(igraph::as_data_frame(induced_subgraph(graph=graph, vids=ipa.nodes)))
  
  # Returning the final results
  return(list("nodes"=ipa.nodes,
              "edges"=ipa.edges,
              "nodes.rep"=ipa.list))
}
#===============================================================================================================================#




#===============================================================================================================================#
#==============#
# Usage         :
#==============#
#                   ipa.boot(graph,
#                            express,
#                            factor,
#                            nboot,
#                            fdr,
#                            parallel,
#                            seed)
#
#==============#
# Description   :
#==============#
#                   Bootstrap replication function of IPA at fixed FDR.
#
#==============#
# Arguments     :
#==============#
#
#==============#
# Values        :
#==============#
#
#===============================================================================================================================#

ipa.boot <- function(graph,
                     express,
                     factor,
                     nboot,
                     fdr,
                     conf,
                     parallel,
                     seed) {
  
  #==============================================================================#
  # Initializations - Constants - Parameters
  #==============================================================================#
  n <- ncol(express)
  p <- nrow(express)
  
  #==============================================================================#
  # Definition of the treatment factor
  #==============================================================================#
  factor.lev <- levels(factor)
  factor.ng <- nlevels(factor)
  factor.def <- vector(mode="list", length=factor.ng)
  factor.tab <- numeric(factor.ng)
  for (g in 1:factor.ng) {
    factor.def[[g]] <- which(factor == factor.lev[g])
    factor.tab[g] <- length(factor.def[[g]])
  }
  
  #==============================================================================#
  # Retrieving sample names
  #==============================================================================#
  sample.names <- colnames(express)
  
  #==============================================================================#
  # Simplify the graph
  #==============================================================================#
  graph <- graph %>% simplify(edge.attr.comb = "sum")
  E(graph)$sign <- sign(E(graph)$sign)
  
  #==============================================================================#
  # Intersection of expression matrix variables and graph nodes
  # Common graph size and dimensionality of the data
  #==============================================================================#
  nodes.names <- igraph::V(graph)$name
  inter_names <- intersect(x=rownames(express), y=nodes.names)
  express <- express[inter_names,]
  graph <- delete_vertices(graph=graph, v=setdiff(x=nodes.names, y=inter_names))
  
  #==============================================================================#
  # Bootstrap replication loop for a fixed FDR value
  #==============================================================================#
  if (!parallel) {
    
    R <- length(nboot)
    nodes.output <- vector(mode="list", length=R)
    r <- 1
    while (r <= R) {
      
      if (R > 1) {
        cat("Replication: ", r, "\n")
      }
      cat("seed: ", seed[r], "\n", sep="")
      if (!is.null(seed[r])) {
        set.seed(seed[r])
      }
      
      #==============================================================================#
      # Define the Bootstrap quantities
      if (R == 1) {
        # No Bootstrap resampling
        samples.boot <- c(factor.def[[1]], factor.def[[2]])
      } else {
        # Bootstrap resampling
        samples.boot <- c(sample(x=factor.def[[1]], size=factor.tab[[1]], replace=TRUE),
                          sample(x=factor.def[[2]], size=factor.tab[[2]], replace=TRUE))
        # Make sure all groups qualify
        while (length(intersect(samples.boot, factor.def[[2]])) < factor.tab[[2]]/2 |
               length(intersect(samples.boot, factor.def[[1]])) < factor.tab[[1]]/2 ) {
          samples.boot <- c(sample(x=factor.def[[1]], size=factor.tab[[1]], replace=TRUE),
                            sample(x=factor.def[[2]], size=factor.tab[[2]], replace=TRUE))
        }
      }
      # Bootstrap quantities
      express.boot <- express[,samples.boot,drop=FALSE]
      factor.boot <- factor[samples.boot]
      factor.boot.lev <- levels(factor.boot)
      factor.boot.ng <- nlevels(factor.boot)
      factor.boot.def <- vector(mode="list", length=factor.boot.ng)
      factor.boot.tab <- numeric(factor.boot.ng)
      for (g in 1:factor.boot.ng) {
        factor.boot.def[[g]] <- which(factor.boot == factor.boot.lev[g])
        factor.boot.tab[g] <- length(factor.boot.def[[g]])
      }
      sample.names.boot <- sample.names[samples.boot]
      
      #==============================================================================#
      # Differential expression for a fixed (user-defined) FDR value
      DEGs <- deg.fdr(fit=deg(express=express.boot, factor=factor.boot), fdr=fdr)$degs[[1]]
      
      #==============================================================================#
      # IPA
      
      # Initializations
      nodeNames <- igraph::V(graph)$name
      
      # Selecting nodes that have significant p-values
      IPAresult <- as.data.frame(matrix(NA,ncol=2, nrow=length(nodeNames)))
      IPAresult[,1] <- nodeNames
      colnames(IPAresult) <- c("node", "p-value")
      
      # Calculating p-values
      for (i in 1:length(nodeNames)) {
        IPAresult[i,2] <-  ipa.pvalue(r=nodeNames[i], graph=graph, degs=DEGs)
      }
      
      # Recording result
      tmp <- IPAresult[which(IPAresult$`p-value` <= 0.05),]$node
      #      tmp <- IPAresult[which(IPAresult$`p-value` <= fdr),]$node
      if (is.empty(tmp)) {
        nodes.output[[r]] <- NA
      } else {
        nodes.output[[r]] <- tmp
      }
      
      r <- r + 1
      
    }
    
  } else {
    
    #==============================================================================#
    # Define the Bootstrap quantities
    if (R == 1) {
      # No Bootstrap resampling
      samples.boot <- c(factor.def[[1]], factor.def[[2]])
    } else {
      # Bootstrap resampling
      samples.boot <- c(sample(x=factor.def[[1]], size=factor.tab[[1]], replace=TRUE),
                        sample(x=factor.def[[2]], size=factor.tab[[2]], replace=TRUE))
      # Make sure all groups qualify
      while (length(intersect(samples.boot, factor.def[[2]])) < factor.tab[[2]]/2 |
             length(intersect(samples.boot, factor.def[[1]])) < factor.tab[[1]]/2 ) {
        samples.boot <- c(sample(x=factor.def[[1]], size=factor.tab[[1]], replace=TRUE),
                          sample(x=factor.def[[2]], size=factor.tab[[2]], replace=TRUE))
      }
    }
    # Bootstrap quantities
    express.boot <- express[,samples.boot,drop=FALSE]
    factor.boot <- factor[samples.boot]
    factor.boot.lev <- levels(factor.boot)
    factor.boot.ng <- nlevels(factor.boot)
    factor.boot.def <- vector(mode="list", length=factor.boot.ng)
    factor.boot.tab <- numeric(factor.boot.ng)
    for (g in 1:factor.boot.ng) {
      factor.boot.def[[g]] <- which(factor.boot == factor.boot.lev[g])
      factor.boot.tab[g] <- length(factor.boot.def[[g]])
    }
    sample.names.boot <- sample.names[samples.boot]
    
    #==============================================================================#
    # Differential expression for a fixed (user-defined) FDR value
    DEGs <- deg.fdr(fit=deg(express=express.boot, factor=factor.boot), fdr=fdr)$degs[[1]]
    
    #==============================================================================#
    # IPA
    
    # Initializations
    nodeNames <- igraph::V(graph)$name
    
    # Selecting nodes that have significant p-values
    IPAresult <- as.data.frame(matrix(NA,ncol=2, nrow=length(nodeNames)))
    IPAresult[,1] <- nodeNames
    colnames(IPAresult) <- c("node", "p-value")
    
    # Calculating p-values
    for (i in 1:length(nodeNames)) {
      IPAresult[i,2] <-  ipa.pvalue(r=nodeNames[i], graph=graph, degs=DEGs)
    }
    
    # Recording result
    tmp <- IPAresult[which(IPAresult$`p-value` <= 0.05),]$node
    #    tmp <- IPAresult[which(IPAresult$`p-value` <= fdr),]$node
    if (is.empty(tmp)) {
      nodes.output <- NA
    } else {
      nodes.output <- tmp
    }
    
  }
  
  #==============================================================================#
  # Returning replication results
  #==============================================================================#
  return(nodes.output)
  
}
#===============================================================================================================================#




#===============================================================================================================================#
#==============#
# Usage         :
#==============#
#                   ipa.pvalue(r, graph, degs)
#
#==============#
# Description   :
#==============#
#                   IPA Subroutine for Fisher's exact test p-value
#
#==============#
# Arguments     :
#==============#
#
#==============#
# Values        :
#==============#
#
#===============================================================================================================================#

ipa.pvalue <- function(r, graph, degs) {
  
  if (r %in% V(graph)$name) {
    Vrg <- unique(igraph::as_data_frame(graph)[,2])
    D <- unique(degs)
    R <- unique(neighbors(graph,r,"out")$name)
    n <- length(Vrg)
    O <- unique(R[which(R %in% D)])
    a <- length(O)
    b <- length(unique(D[which(D %in% Vrg)])) - a
    c <- length(R) - a
    d <- n - a - b - c
    p.value <- fisher.test(as.table(matrix(base::c((a),(c),(b),(d)), nrow = 2)), alternative = "greater")$p.value
    
    return(p.value)
    
  } else {
    
    return(NA)
    
  }
  
}
#===============================================================================================================================#




#===============================================================================================================================#
#===============#
# Usage         :
#===============#
#                    list2array (list, rowtrunc=NULL, coltrunc=NULL, sub=NULL, fill=NA)
#
#===============#
# Description   :
#===============#
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#===============================================================================================================================#

list2array <- function (list, rowtrunc=NULL, coltrunc=NULL, sub=NULL, fill=NA) {
  
  if (!is.empty(list)) {
    if (is.null(sub)) {
      my.list <- list
    } else {
      L <- length(list)
      my.list <- vector(mode="list", length=L)
      for (i in 1:L) {
        my.list[[i]] <- list[[i]][[sub]]
      }
    }
    min.row <- min(vapply(my.list, function(x) nrow(x), numeric(1), na.rm = TRUE))
    max.row <- max(vapply(my.list, function(x) nrow(x), numeric(1), na.rm = TRUE))
    min.col <- min(vapply(my.list, function(x) ncol(x), numeric(1), na.rm = TRUE))
    max.col <- max(vapply(my.list, function(x) ncol(x), numeric(1), na.rm = TRUE))
    if (!is.null(coltrunc)) {
      if (coltrunc == "min") {
        adjusted.list <- lapply(my.list, function(x) {x[,1:min.col,drop=FALSE]})
      } else if (coltrunc == "max") {
        adjusted.list <- lapply(my.list, function(x) {cbind(x, matrix(data=fill, nrow=nrow(x), ncol=max.col - ncol(x)))})
      } else {
        adjusted.list <- lapply(my.list, function(x, coltrunc) {if (coltrunc <= ncol(x)) {
          x[,1:coltrunc,drop=FALSE]
        } else if (coltrunc > ncol(x)) {
          cbind(x, matrix(data=fill, nrow=nrow(x), ncol=coltrunc - ncol(x)))
        }
        }, coltrunc)
      }
    } else {
      adjusted.list <- lapply(my.list, function(x) {cbind(x, matrix(data=fill, nrow=nrow(x), ncol=max.col - ncol(x)))})
    }
    if (!is.null(rowtrunc)) {
      if (rowtrunc == "min") {
        adjusted.list <- lapply(adjusted.list, function(x) {x[1:min.row,,drop=FALSE]})
      } else if (rowtrunc == "max") {
        adjusted.list <- lapply(adjusted.list, function(x) {rbind(x, matrix(data=fill, nrow=max.row - nrow(x), ncol=ncol(x)))})
      } else {
        
        adjusted.list <- lapply(my.list, function(x, rowtrunc) {if (rowtrunc <= nrow(x)) {
          x[1:rowtrunc,,drop=FALSE]
        } else if (rowtrunc > nrow(x)) {
          rbind(x, matrix(data=fill, nrow=rowtrunc - nrow(x), ncol=ncol(x)))
        }
        }, rowtrunc)
      }
    } else {
      adjusted.list <- lapply(adjusted.list, function(x) {rbind(x, matrix(data=fill, nrow=max.row - nrow(x), ncol=ncol(x)))})
    }
    my.array <- array(data=fill, dim=c(nrow(adjusted.list[[1]]), ncol(adjusted.list[[1]]), length(adjusted.list)))
    for(i in 1:length(adjusted.list)) {
      my.array[,,i] <- adjusted.list[[i]]
    }
  } else {
    my.array <- array(data=fill, dim=c(0,0,0))
  }
  
  return(my.array)
  
}
#===============================================================================================================================#