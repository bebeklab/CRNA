# ===============================================================================================================================#
# ==============#
# Usage         :
# ==============#
#                   crna.sl(graph,
#                           express,
#                           factors,
#                           FC=1,
#                           FR=0.1,
#                           R=100,
#                           alpha=0.80,
#                           size=100,
#                           fdr=10^(-2),
#                           conf=NULL,
#                           parallel=FALSE,
#                           seed=NULL)
#
# ==============#
# Description   :
# ==============#
#                   CRNA function at fixed FDR.
#
# ==============#
# Arguments     :
# ==============#
#
# ==============#
# Values        :
# ==============#
#
# ===============================================================================================================================#
#' @export
crna.sl <- function(graph,
                    express,
                    factors,
                    FC = 1,
                    FR = 0.1,
                    R = 100,
                    alpha = 0.80,
                    size = 100,
                    fdr = 10^(-2),
                    conf = NULL,
                    parallel = FALSE,
                    seed = NULL) {
  seed <- seed[1]
  nboot <- 1:R
  n <- ncol(express)
  p <- nrow(express)

  factor1 <- factors[[1]]
  factor1.lev <- levels(factor1)
  factor1.ng <- nlevels(factor1)
  factor1.def <- vector(mode = "list", length = factor1.ng)
  factor1.tab <- numeric(factor1.ng)
  for (g in 1:factor1.ng) {
    factor1.def[[g]] <- which(factor1 == factor1.lev[g])
    factor1.tab[g] <- length(factor1.def[[g]])
  }

  if (!parallel) {
    if (!is.null(seed)) {
      seed <- (0:(R - 1)) + seed
    }

    crna.list <- crna.boot.sl(
      graph = graph,
      express = express,
      factor = factor1,
      nboot = nboot,
      alpha = alpha,
      size = size,
      fdr = fdr,
      parallel = parallel,
      seed = seed
    )
  } else {
    # Parallel backend registration
    # To be used with the 'parallel' package:
    if (conf$type == "SOCKET") {
      clus <- parallel::makeCluster(
        spec = conf$spec,
        type = "PSOCK",
        #rscript_args = c("-e", shQuote("getRversion()")),
        homogeneous = conf$homo,
        outfile = conf$outfile,
        verbose = conf$verbose
      )
    } else if (conf$type == "MPI") {
      clus <- parallel::makeCluster(
        spec = conf$spec,
        type = "MPI",
        #rscript_args = c("-e", shQuote("getRversion()")),
        homogeneous = conf$homo,
        outfile = conf$outfile,
        verbose = conf$verbose
      )
    } else {
      stop("Unrecognized cluster type\n")
    }
    clusterEvalQ(cl = clus, expr = library("parallel"))
    clusterEvalQ(cl = clus, expr = library("limma"))
    clusterEvalQ(cl = clus, expr = library("Matrix"))
    clusterEvalQ(cl = clus, expr = library("igraph"))
    clusterExport(cl = clus, varlist = ls(.GlobalEnv), envir = .GlobalEnv)
    clusterExport(
      cl = clus,
      varlist = c("crna.boot.sl", "deg", "deg.fdr", "rwr", "dfs", "is.empty"),
      envir = .GlobalEnv
    )
    clusterSetRNGStream(cl = clus, iseed = seed)
    obj.clus <- clusterApply(
      cl = clus,
      fun = crna.boot.sl,
      x = nboot,
      graph = graph,
      express = express,
      factor = factor1,
      alpha = alpha,
      size = size,
      fdr = fdr,
      parallel = parallel,
      seed = NULL
    )
    crna.list <- list(
      "edgelist" = vector(mode = "list", length = R),
      "scores.samples" = vector(mode = "list", length = R)
    )
    for (r in 1:R) {
      crna.list$edgelist[[r]] <- obj.clus[[r]]$edgelist
      crna.list$scores.samples[[r]] <- obj.clus[[r]]$scores.samples
    }
    # Stopping the cluster and cleaning all MPI states
    parallel::stopCluster(cl = clus)
  }

  # Aggregating replications results of activity scores by sample
  for (r in 1:R) {
    tmp <- colnames(crna.list$scores.samples[[r]])
    colnames(crna.list$scores.samples[[r]]) <- 1:ncol(crna.list$scores.samples[[r]])
    scores.samples <- t(crna.list$scores.samples[[r]])
    scores.samples <- cbind.data.frame("symbol" = tmp, scores.samples, stringsAsFactors = FALSE)
    scores.samples <- aggregate(. ~ symbol, data = scores.samples, mean, na.rm = TRUE)
    crna.list$scores.samples[[r]] <- as.matrix(t(scores.samples[, -1]))
    colnames(crna.list$scores.samples[[r]]) <- scores.samples[, 1]
  }
  crna.scores.samples <- matrix(data = 0, nrow = p, ncol = n, dimnames = list(rownames(express), colnames(express)))
  for (r in 1:R) {
    scores.samples <- crna.list$scores.samples[[r]]
    crna.scores.samples[rownames(scores.samples),colnames(scores.samples)] <- crna.scores.samples[rownames(scores.samples),colnames(scores.samples)] + scores.samples
  }

  # Calculating activity scores by treatment group
  crna.scores.treatments <- cbind(
    rowMeans(crna.scores.samples[, factor1.def[[1]]], na.rm = TRUE),
    rowMeans(crna.scores.samples[, factor1.def[[2]]], na.rm = TRUE)
  )
  rownames(crna.scores.treatments) <- rownames(crna.scores.samples)
  colnames(crna.scores.treatments) <- factor1.lev

  # Aggregating replications results of edges
  crna.edgelist <- matrix(data = NA, nrow = 0, ncol = 4, dimnames = list(NULL, c("from", "to", "edgesign", "weight")))
  for (r in 1:R) {
    crna.edgelist <- unique(rbind(crna.list$edgelist[[r]], crna.edgelist))
  }

  # Aggregating replications results of nodes
  nodes.names <- character(0)
  nodes.rep <- vector(mode = "list", length = R)
  for (r in 1:R) {
    nodes.rep[[r]] <- rownames(crna.list$scores.samples[[r]])
    nodes.names <- unique(c(nodes.names, nodes.rep[[r]]))
  }
  nodes.names <- nodes.names[!is.na(nodes.names)]
  p.adj <- length(nodes.names)

  # Filtering by frequency of occurrence (in replications) and by activity scores
  scores.samples.array <- array(data = NA, dim = c(p.adj, n, R), dimnames = list(nodes.names, colnames(express), 1:R))
  for (r in 1:R) {
    q <- ncol(crna.list$scores.samples[[r]])
    w <- pmatch(x = rownames(crna.list$scores.samples[[r]]), table = nodes.names)
    if (!is.null(w)) {
      scores.samples.array[w, 1:q, r] <- crna.list$scores.samples[[r]]
    }
  }
  scores.freq.samples <- apply(scores.samples.array, 1:2, function(x) {
    R - sum(is.na(x))
  })
  scores.freq.treatments <- cbind(
    rowMeans(scores.freq.samples[, factor1.def[[1]]], na.rm = TRUE),
    rowMeans(scores.freq.samples[, factor1.def[[2]]], na.rm = TRUE)
  )
  if (FC != 1) {
    selected.nodes.names <- names(which(((scores.freq.treatments[, 1] >= FR * R) | (scores.freq.treatments[, 2] >= FR * R)) &
      (abs(log(x = crna.scores.treatments[nodes.names, 2] / crna.scores.treatments[nodes.names, 1], base = FC)) >= 1)))
  } else {
    selected.nodes.names <- names(which((scores.freq.treatments[, 1] >= FR * R) | (scores.freq.treatments[, 2] >= FR * R)))
  }

  # Updating crna.scores.samples, crna.scores.treatments, and crna.states
  crna.scores.treatments <- crna.scores.treatments[selected.nodes.names, ]
  crna.scores.samples <- crna.scores.samples[selected.nodes.names, ]
  crna.states <- cbind.data.frame(
    "symbol" = rownames(crna.scores.treatments),
    "nodesign" = sign(rowSums(crna.scores.treatments, na.rm = TRUE))
  )

  # Updating graph
  RWR_JT_graph <- igraph::graph_from_data_frame(crna.edgelist, directed = TRUE)
  RWR_JT_graph <- igraph::induced.subgraph(RWR_JT_graph, rownames(crna.states))
  RWR_JT_graph <- igraph::simplify(
    graph = RWR_JT_graph,
    remove.multiple = TRUE,
    remove.loops = TRUE,
    edge.attr.comb = "max"
  )
  crna.edges <- igraph::as_data_frame(RWR_JT_graph)
  tmp <- igraph::V(RWR_JT_graph)$name
  if (is.empty(tmp)) {
    crna.nodes <- NA
  } else {
    crna.nodes <- tmp
  }

  # Updating by removing null activity score nodes
  crna.states <- crna.states[which(rownames(crna.states) %in% crna.nodes), ]
  crna.scores.samples <- crna.scores.samples[which(rownames(crna.scores.samples) %in% crna.nodes), ]
  crna.express.samples <- express[which(rownames(express) %in% crna.nodes), ]

  # Returning the final 'crna' object
  return(structure(
    list(
      "nodes" = crna.nodes,
      "edges" = crna.edges,
      "states" = crna.states,
      "scores.samples" = crna.scores.samples,
      "express.samples" = crna.express.samples,
      "nodes.rep" = nodes.rep
    ),
    class = "crna"
  ))
}
# ===============================================================================================================================#




# ===============================================================================================================================#
# ==============#
# Usage         :
# ==============#
#                   crna.as(graph,
#                           express,
#                           factors,
#                           FC=1,
#                           FR=0.1,
#                           R=100,
#                           alpha=0.80,
#                           quant=0.90,
#                           fdr=10^(-2),
#                           conf=NULL,
#                           parallel=FALSE,
#                           seed=NULL)
#
# ==============#
# Description   :
# ==============#
#                   CRNA function at fixed FDR.

#
# ==============#
# Arguments     :
# ==============#
#
# ==============#
# Values        :
# ==============#
#
# ===============================================================================================================================#

crna.as <- function(graph,
                    express,
                    factors,
                    FC = 1,
                    FR = 0.1,
                    R = 100,
                    alpha = 0.80,
                    quant = 0.90,
                    fdr = 10^(-2),
                    conf = NULL,
                    parallel = FALSE,
                    seed = NULL) {
  seed <- seed[1]
  nboot <- 1:R
  n <- ncol(express)
  p <- nrow(express)

  factor1 <- factors[[1]]
  factor1.lev <- levels(factor1)
  factor1.ng <- nlevels(factor1)
  factor1.def <- vector(mode = "list", length = factor1.ng)
  factor1.tab <- numeric(factor1.ng)
  for (g in 1:factor1.ng) {
    factor1.def[[g]] <- which(factor1 == factor1.lev[g])
    factor1.tab[g] <- length(factor1.def[[g]])
  }

  if (!parallel) {
    if (!is.null(seed)) {
      seed <- (0:(R - 1)) + seed
    }
    # debug(crna.bocrna.asot.as)
    crna.list <- crna.boot.as(
      graph = graph,
      express = express,
      factor = factor1,
      nboot = nboot,
      alpha = alpha,
      quant = quant,
      fdr = fdr,
      parallel = parallel,
      seed = seed
    )
  } else {
    # Parallel backend registration
    # To be used with the 'parallel' package:
    if (conf$type == "SOCKET") {
      clus <- parallel::makeCluster(
        spec = conf$spec,
        type = "PSOCK",
        #rscript_args = c("-e", shQuote("getRversion()")),
        homogeneous = conf$homo,
        outfile = conf$outfile,
        verbose = conf$verbose
      )
    } else if (conf$type == "MPI") {
      clus <- parallel::makeCluster(
        spec = conf$spec,
        type = "MPI",
        #rscript_args = c("-e", shQuote("getRversion()")),
        homogeneous = conf$homo,
        outfile = conf$outfile,
        verbose = conf$verbose
      )
    } else {
      stop("Unrecognized cluster type\n")
    }
    clusterEvalQ(cl = clus, expr = library("parallel"))
    clusterEvalQ(cl = clus, expr = library("limma"))
    clusterEvalQ(cl = clus, expr = library("Matrix"))
    clusterEvalQ(cl = clus, expr = library("igraph"))
    clusterExport(cl = clus, varlist = ls(.GlobalEnv), envir = .GlobalEnv)
    clusterExport(
      cl = clus,
      varlist = c("crna.boot.as", "deg", "deg.fdr", "rwr", "dfs", "is.empty"),
      envir = .GlobalEnv
    )
    clusterSetRNGStream(cl = clus, iseed = seed)
    obj.clus <- clusterApply(
      cl = clus,
      fun = crna.boot.as,
      x = nboot,
      graph = graph,
      express = express,
      factor = factor1,
      alpha = alpha,
      quant = quant,
      fdr = fdr,
      parallel = parallel,
      seed = NULL
    )
    crna.list <- list(
      "edgelist" = vector(mode = "list", length = R),
      "scores.samples" = vector(mode = "list", length = R)
    )
    for (r in 1:R) {
      crna.list$edgelist[[r]] <- obj.clus[[r]]$edgelist
      crna.list$scores.samples[[r]] <- obj.clus[[r]]$scores.samples
    }
    # Stopping the cluster and cleaning all MPI states
    parallel::stopCluster(cl = clus)
  }

  # Aggregating replications results of activity scores by sample
  for (r in 1:R) {
    tmp <- colnames(crna.list$scores.samples[[r]])
    colnames(crna.list$scores.samples[[r]]) <- 1:ncol(crna.list$scores.samples[[r]])
    scores.samples <- t(crna.list$scores.samples[[r]])
    scores.samples <- cbind.data.frame("symbol" = tmp, scores.samples, stringsAsFactors = FALSE)
    scores.samples <- aggregate(. ~ symbol, data = scores.samples, mean, na.rm = TRUE)
    crna.list$scores.samples[[r]] <- as.matrix(t(scores.samples[, -1]))
    colnames(crna.list$scores.samples[[r]]) <- scores.samples[, 1]
  }
  crna.scores.samples <- matrix(data = 0, nrow = p, ncol = n, dimnames = list(rownames(express), colnames(express)))
  for (r in 1:R) {
    scores.samples <- crna.list$scores.samples[[r]]
    crna.scores.samples[rownames(scores.samples),colnames(scores.samples)] <- crna.scores.samples[rownames(scores.samples),colnames(scores.samples)] + scores.samples
  }

  # Calculating activity scores by treatment group
  crna.scores.treatments <- cbind(
    rowMeans(crna.scores.samples[, factor1.def[[1]]], na.rm = TRUE),
    rowMeans(crna.scores.samples[, factor1.def[[2]]], na.rm = TRUE)
  )
  rownames(crna.scores.treatments) <- rownames(crna.scores.samples)
  colnames(crna.scores.treatments) <- factor1.lev

  # Aggregating replications results of edges
  crna.edgelist <- matrix(data = NA, nrow = 0, ncol = 4, dimnames = list(NULL, c("from", "to", "edgesign", "weight")))
  for (r in 1:R) {
    crna.edgelist <- unique(rbind(crna.list$edgelist[[r]], crna.edgelist))
  }

  # Aggregating replications results of nodes
  nodes.names <- character(0)
  nodes.rep <- vector(mode = "list", length = R)
  for (r in 1:R) {
    nodes.rep[[r]] <- rownames(crna.list$scores.samples[[r]])
    nodes.names <- unique(c(nodes.names, nodes.rep[[r]]))
  }
  nodes.names <- nodes.names[!is.na(nodes.names)]
  p.adj <- length(nodes.names)

  # Filtering by frequency of occurrence (in replications) and by activity scores
  scores.samples.array <- array(data = NA, dim = c(p.adj, n, R), dimnames = list(nodes.names, colnames(express), 1:R))
  for (r in 1:R) {
    q <- ncol(crna.list$scores.samples[[r]])
    w <- pmatch(x = rownames(crna.list$scores.samples[[r]]), table = nodes.names)
    if (!is.null(w)) {
      scores.samples.array[w, 1:q, r] <- crna.list$scores.samples[[r]]
    }
  }
  scores.freq.samples <- apply(scores.samples.array, 1:2, function(x) {
    R - sum(is.na(x))
  })
  scores.freq.treatments <- cbind(
    rowMeans(scores.freq.samples[, factor1.def[[1]]], na.rm = TRUE),
    rowMeans(scores.freq.samples[, factor1.def[[2]]], na.rm = TRUE)
  )
  if (FC != 1) {
    selected.nodes.names <- names(which(((scores.freq.treatments[, 1] >= FR * R) | (scores.freq.treatments[, 2] >= FR * R)) &
      (abs(log(x = crna.scores.treatments[nodes.names, 2] / crna.scores.treatments[nodes.names, 1], base = FC)) >= 1)))
  } else {
    selected.nodes.names <- names(which((scores.freq.treatments[, 1] >= FR * R) | (scores.freq.treatments[, 2] >= FR * R)))
  }

  # Updating crna.scores.samples, crna.scores.treatments, and crna.states
  crna.scores.treatments <- crna.scores.treatments[selected.nodes.names, ]
  crna.scores.samples <- crna.scores.samples[selected.nodes.names, ]
  crna.states <- cbind.data.frame(
    "symbol" = rownames(crna.scores.treatments),
    "nodesign" = sign(rowSums(crna.scores.treatments, na.rm = TRUE))
  )

  # Updating graph
  RWR_JT_graph <- igraph::graph_from_data_frame(crna.edgelist, directed = TRUE)
  RWR_JT_graph <- igraph::induced.subgraph(RWR_JT_graph, rownames(crna.states))
  RWR_JT_graph <- igraph::simplify(
    graph = RWR_JT_graph,
    remove.multiple = TRUE,
    remove.loops = TRUE,
    edge.attr.comb = "max"
  )
  crna.edges <- igraph::as_data_frame(RWR_JT_graph)
  tmp <- igraph::V(RWR_JT_graph)$name
  if (is.empty(tmp)) {
    crna.nodes <- NA
  } else {
    crna.nodes <- tmp
  }

  # Updating by removing null activity score nodes
  crna.states <- crna.states[which(rownames(crna.states) %in% crna.nodes), ]
  crna.scores.samples <- crna.scores.samples[which(rownames(crna.scores.samples) %in% crna.nodes), ]
  crna.express.samples <- express[which(rownames(express) %in% crna.nodes), ]

  # Returning the final 'crna' object
  return(structure(
    list(
      "nodes" = crna.nodes,
      "edges" = crna.edges,
      "states" = crna.states,
      "scores.samples" = crna.scores.samples,
      "express.samples" = crna.express.samples,
      "nodes.rep" = nodes.rep
    ),
    class = "crna"
  ))
}
# ===============================================================================================================================#




# ===============================================================================================================================#
# ==============#
# Usage         :
# ==============#
#                   crna.boot.sl(graph,
#                                express,
#                                factor,
#                                nboot,
#                                alpha,
#                                size,
#                                fdr,
#                                conf,
#                                parallel,
#                                seed)
#
# ==============#
# Description   :
# ==============#
#                   Bootstrap replication function of CRNA at fixed FDR.
#
# ==============#
# Arguments     :
# ==============#
#
# ==============#
# Values        :
# ==============#
#
# ===============================================================================================================================#

crna.boot.sl <- function(graph,
                         express,
                         factor,
                         nboot,
                         alpha,
                         size,
                         fdr,
                         conf,
                         parallel,
                         seed) {
  # ==============================================================================#
  # Definition of the treatment factor
  # ==============================================================================#
  factor.lev <- levels(factor)
  factor.ng <- nlevels(factor)
  factor.def <- vector(mode = "list", length = factor.ng)
  factor.tab <- numeric(factor.ng)
  for (g in 1:factor.ng) {
    factor.def[[g]] <- which(factor == factor.lev[g])
    factor.tab[g] <- length(factor.def[[g]])
  }

  # ==============================================================================#
  # Simplify the graph
  # ==============================================================================#
  graph <- graph %>% simplify(edge.attr.comb = "sum")
  E(graph)$sign <- sign(E(graph)$sign)
  # RWR_edges <- igraph::as_data_frame(graph)
  # RWR_edges <- RWR_edges[RWR_edges$sign != 0, ]
  # graph <- graph_from_data_frame(RWR_edges)
  # ==============================================================================#
  # Intersection of expression matrix variables and graph nodes
  # Common graph size and dimensionality of the data
  # ==============================================================================#
  nodes.names <- igraph::V(graph)$name
  inter_names <- intersect(x = rownames(express), y = nodes.names)
  express <- express[inter_names, ]
  graph <- delete_vertices(graph = graph, v = setdiff(x = nodes.names, y = inter_names))

  # ==============================================================================#
  # Convert the graph to an undirected graph, and then to an adjacency matrix
  # ==============================================================================#
  graph.ud <- igraph::as_adjacency_matrix(igraph::as.undirected(graph), attr = "weight")

  # ==============================================================================#
  # Bootstrap replication loop for a fixed FDR value
  # ==============================================================================#
  if (!parallel) {
    R <- length(nboot)
    edgelist.output <- vector(mode = "list", length = R)
    scores.samples.output <- vector(mode = "list", length = R)
    r <- 1
    while (r <= R) {
      if (R > 1) {
        cat("Replication: ", r, "/", R, "\n", sep = "")
      }
      cat("seed: ", seed[r], "\n", sep = "")
      if (!is.null(seed[r])) {
        set.seed(seed[r])
      }

      # ==============================================================================#
      # Define the Bootstrap quantities
      if (R == 1) {
        # No Bootstrap resampling
        samples.boot <- c(factor.def[[1]], factor.def[[2]])
      } else {
        # Bootstrap resampling
        samples.boot <- c(
          sample(x = factor.def[[1]], size = factor.tab[[1]], replace = TRUE),
          sample(x = factor.def[[2]], size = factor.tab[[2]], replace = TRUE)
        )
        # Make sure all groups qualify
        while (length(intersect(samples.boot, factor.def[[2]])) < factor.tab[[2]] / 2 |
          length(intersect(samples.boot, factor.def[[1]])) < factor.tab[[1]] / 2) {
          samples.boot <- c(
            sample(x = factor.def[[1]], size = factor.tab[[1]], replace = TRUE),
            sample(x = factor.def[[2]], size = factor.tab[[2]], replace = TRUE)
          )
        }
      }
      # Bootstrap quantities
      express.boot <- express[, samples.boot, drop = FALSE]
      factor.boot <- factor[samples.boot]
      factor.boot.lev <- levels(factor.boot)
      factor.boot.ng <- nlevels(factor.boot)
      factor.boot.def <- vector(mode = "list", length = factor.boot.ng)
      factor.boot.tab <- numeric(factor.boot.ng)
      for (g in 1:factor.boot.ng) {
        factor.boot.def[[g]] <- which(factor.boot == factor.boot.lev[g])
        factor.boot.tab[g] <- length(factor.boot.def[[g]])
      }
      n <- ncol(express.boot)
      p <- nrow(express.boot)
      sample.names.boot <- sample.names[samples.boot]

      # ==============================================================================#
      # Differential expression for a fixed (user-defined) FDR value
      DEGs <- deg.fdr(fit = deg(express = express.boot, factor = factor.boot), fdr = fdr)$degs[[1]]

      # ==============================================================================#
      # CRNA

      #---------------------------------- RWR ---------------------------------#
      affinity.scores <- rwr(
        graph = graph.ud,
        seedSet = DEGs,
        r = alpha,
        random.seed = FALSE,
        n.thread = NULL
      )
      affinity.scores <- sort(affinity.scores, decreasing = TRUE)
      sig_nodes_names <- names(affinity.scores[1:size])

      # Generating network with significant nodes only and clean up
      RWR_graph <- induced.subgraph(graph, sig_nodes_names)
      RWR_graph <- igraph::simplify(
        graph = RWR_graph,
        remove.multiple = TRUE,
        remove.loops = TRUE,
        edge.attr.comb = "max"
      )
      RWR_edges <- igraph::as_data_frame(RWR_graph)
      RWR_edges <- RWR_edges[RWR_edges$sign != 0, ]
      RWR_graph <- graph_from_data_frame(RWR_edges)

      #---------------------------------- DFS ---------------------------------#
      # Initializing DFS result receptor
      edgelist <- matrix(data = NA, nrow = 0, ncol = 4)
      colnames(edgelist) <- c("from", "to", "edgesign", "weight")
      activity.scores.samples <- matrix(data = 0, nrow = vcount(RWR_graph), ncol = length(sample.names.boot))
      rownames(activity.scores.samples) <- V(RWR_graph)$name
      colnames(activity.scores.samples) <- sample.names.boot

      # Scaling expression data
      express.scaled <- t(scale(t(express.boot), center = TRUE, scale = TRUE))

      # DFS graph for every candidate
      sig_nodes_names <- intersect(sig_nodes_names, igraph::V(RWR_graph)$name)

      for (i in 1:length(sig_nodes_names)) {
        RWR_DFS_graph <- dfs(sig_nodes_names[i], RWR_graph, DEGs)
        RWR_DFS_edges <- igraph::as_data_frame(RWR_DFS_graph)

        if (nrow(RWR_DFS_edges) > 0) {
          #-------------------------- Inputs for Junction Tree ------------------------#
          # Nodes and edges description
          RWR_DFS_nodes <- as.data.frame(as.character(V(RWR_DFS_graph)$name))
          RWR_DFS_nodes[, 2] <- as.character(NA)
          colnames(RWR_DFS_nodes) <- c("id", "type")
          RWR_DFS_nodes[, 2] <- "protein"
          RWR_DFS_nodes <- RWR_DFS_nodes[, c(2, 1)]
          RWR_DFS_edges <- RWR_DFS_edges[, c("from", "to", "sign")]
          RWR_DFS_edges[which(RWR_DFS_edges$sign == 1), 3] <- "-a>"
          RWR_DFS_edges[which(RWR_DFS_edges$sign == -1), 3] <- "-a|"

          # Subsetting expression data
          express.scaled.sub <- express.scaled[as.character(RWR_DFS_nodes$id), ]
          rownames(express.scaled.sub) <- as.character(RWR_DFS_nodes$id)
          colnames(express.scaled.sub) <- 1:ncol(express.scaled.sub)
          RWR_DFS_mRNA <- t(express.scaled.sub)

          #-------------------------- Junction Tree Results ---------------------------#
          RWR_DFS_genes <- colnames(RWR_DFS_mRNA)
          RWR_DFS_mRNA <- cbind(sample.names.boot, RWR_DFS_mRNA)
          colnames(RWR_DFS_mRNA) <- paste0("V", seq_len(ncol(RWR_DFS_mRNA)))
          colnames(RWR_DFS_nodes) <- paste0("V", seq_len(ncol(RWR_DFS_nodes)))
          colnames(RWR_DFS_edges) <- paste0("V", seq_len(ncol(RWR_DFS_edges)))
          JT.mat <- CRNA::JTinf(RWR_DFS_nodes, RWR_DFS_edges, RWR_DFS_mRNA, RWR_DFS_genes)
          JT.mat <- JT.mat[-1, ]

          simulated_as <- matrix(
            data = NA,
            nrow = nrow(JT.mat) / n,
            ncol = n,
            dimnames = list(JT.mat$V1[1:(nrow(JT.mat) / n)], sample.names.boot)
          )

          # Triming up all data
          for (j in 1:n) {
            tails <- j * nrow(JT.mat) / n
            heads <- tails - (nrow(JT.mat) / n - 1)
            simulated_as[, j] <- as.numeric(JT.mat[heads:tails, 2])
          }

          # Selecting only nodes which are not all zero
          simulated_as_sig <- simulated_as[which(rowSums(simulated_as) != 0), ]

          # Building activity score matrix and majority vote matrix
          activity.scores.samples[rownames(simulated_as_sig), ] <- activity.scores.samples[rownames(simulated_as_sig), ] + simulated_as_sig / degree(RWR_DFS_graph)[rownames(simulated_as_sig)]

          # Recording the subgraph for each DFS
          RWR_JT_graph <- igraph::induced.subgraph(RWR_DFS_graph, rownames(simulated_as_sig))
          edgelist <- rbind(edgelist, igraph::as_data_frame(RWR_JT_graph))
        }
      }

      if (nrow(edgelist) > 0) {
        activity.scores.samples <- activity.scores.samples[unique(c(edgelist[, "from"], edgelist[, "to"])), ]

        # Outputting results
        message("CRNA successfully run. Graph results are returned.\n")
        edgelist.output[[r]] <- unique(edgelist)
        scores.samples.output[[r]] <- activity.scores.samples
      } else {
        # Outputting results
        message("CRNA could not run (Unconnected DFS graph). NA results are returned\n")
        edgelist.output[[r]] <- NA
        scores.samples.output[[r]] <- NA
      }

      r <- r + 1
    }
  } else {
    # ==============================================================================#
    # Get the process ID of each slave R session
    pid <- Sys.getpid()

    # ==============================================================================#
    # Define the Bootstrap quantities
    if (R == 1) {
      # No Bootstrap resampling
      samples.boot <- c(factor.def[[1]], factor.def[[2]])
    } else {
      # Bootstrap resampling
      samples.boot <- c(
        sample(x = factor.def[[1]], size = factor.tab[[1]], replace = TRUE),
        sample(x = factor.def[[2]], size = factor.tab[[2]], replace = TRUE)
      )
      # Make sure all groups qualify
      while (length(intersect(samples.boot, factor.def[[2]])) < factor.tab[[2]] / 2 |
        length(intersect(samples.boot, factor.def[[1]])) < factor.tab[[1]] / 2) {
        samples.boot <- c(
          sample(x = factor.def[[1]], size = factor.tab[[1]], replace = TRUE),
          sample(x = factor.def[[2]], size = factor.tab[[2]], replace = TRUE)
        )
      }
    }
    # Bootstrap quantities
    express.boot <- express[, samples.boot, drop = FALSE]
    factor.boot <- factor[samples.boot]
    factor.boot.lev <- levels(factor.boot)
    factor.boot.ng <- nlevels(factor.boot)
    factor.boot.def <- vector(mode = "list", length = factor.boot.ng)
    factor.boot.tab <- numeric(factor.boot.ng)
    for (g in 1:factor.boot.ng) {
      factor.boot.def[[g]] <- which(factor.boot == factor.boot.lev[g])
      factor.boot.tab[g] <- length(factor.boot.def[[g]])
    }
    n <- ncol(express.boot)
    p <- nrow(express.boot)
    sample.names.boot <- sample.names[samples.boot]

    # ==============================================================================#
    # Differential expression for a fixed (user-defined) FDR value
    DEGs <- deg.fdr(fit = deg(express = express.boot, factor = factor.boot), fdr = fdr)$degs[[1]]

    # ==============================================================================#
    # CRNA

    #---------------------------------- RWR ---------------------------------#
    affinity.scores <- rwr(
      graph = graph.ud,
      seedSet = DEGs,
      r = alpha,
      random.seed = FALSE,
      n.thread = NULL
    )
    affinity.scores <- sort(affinity.scores, decreasing = TRUE)
    sig_nodes_names <- names(affinity.scores[1:size])

    # Generating network with significant nodes only and clean up
    RWR_graph <- induced.subgraph(graph, sig_nodes_names)
    RWR_graph <- igraph::simplify(
      graph = RWR_graph,
      remove.multiple = TRUE,
      remove.loops = TRUE,
      edge.attr.comb = "max"
    )
    RWR_edges <- igraph::as_data_frame(RWR_graph)
    RWR_edges <- RWR_edges[RWR_edges$sign != 0, ]
    RWR_graph <- graph_from_data_frame(RWR_edges)

    #---------------------------------- DFS ---------------------------------#
    # Initializing DFS result receptor
    edgelist <- matrix(data = NA, nrow = 0, ncol = 4)
    colnames(edgelist) <- c("from", "to", "edgesign", "weight")
    activity.scores.samples <- matrix(data = 0, nrow = vcount(RWR_graph), ncol = length(sample.names.boot))
    rownames(activity.scores.samples) <- V(RWR_graph)$name
    colnames(activity.scores.samples) <- sample.names.boot

    # Scaling expression data
    express.scaled <- t(scale(t(express.boot), center = TRUE, scale = TRUE))

    # DFS graph for every candidate
    sig_nodes_names <- intersect(sig_nodes_names, igraph::V(RWR_graph)$name)

    for (i in 1:length(sig_nodes_names)) {
      RWR_DFS_graph <- dfs(sig_nodes_names[i], RWR_graph, DEGs)
      RWR_DFS_edges <- igraph::as_data_frame(RWR_DFS_graph)

      if (nrow(RWR_DFS_edges) > 0) {
        #-------------------------- Inputs for Junction Tree ------------------------#
        # Nodes and edges description
        RWR_DFS_nodes <- as.data.frame(as.character(V(RWR_DFS_graph)$name))
        RWR_DFS_nodes[, 2] <- as.character(NA)
        colnames(RWR_DFS_nodes) <- c("id", "type")
        RWR_DFS_nodes[, 2] <- "protein"
        RWR_DFS_nodes <- RWR_DFS_nodes[, c(2, 1)]
        RWR_DFS_edges <- RWR_DFS_edges[, c("from", "to", "sign")]
        RWR_DFS_edges[which(RWR_DFS_edges$sign == 1), 3] <- "-a>"
        RWR_DFS_edges[which(RWR_DFS_edges$sign == -1), 3] <- "-a|"

        # Subsetting expression data
        express.scaled.sub <- express.scaled[as.character(RWR_DFS_nodes$id), ]
        rownames(express.scaled.sub) <- as.character(RWR_DFS_nodes$id)
        colnames(express.scaled.sub) <- 1:ncol(express.scaled.sub)
        RWR_DFS_mRNA <- t(express.scaled.sub)

        #-------------------------- Junction Tree Results ---------------------------#
        RWR_DFS_genes <- colnames(RWR_DFS_mRNA)
        RWR_DFS_mRNA <- cbind(sample.names.boot, RWR_DFS_mRNA)
        colnames(RWR_DFS_mRNA) <- paste0("V", seq_len(ncol(RWR_DFS_mRNA)))
        colnames(RWR_DFS_nodes) <- paste0("V", seq_len(ncol(RWR_DFS_nodes)))
        colnames(RWR_DFS_edges) <- paste0("V", seq_len(ncol(RWR_DFS_edges)))
        JT.mat <- CRNA::JTinf(RWR_DFS_nodes, RWR_DFS_edges, RWR_DFS_mRNA, RWR_DFS_genes)
        JT.mat <- JT.mat[-1, ]

        simulated_as <- matrix(
          data = NA,
          nrow = nrow(JT.mat) / n,
          ncol = n,
          dimnames = list(JT.mat$V1[1:(nrow(JT.mat) / n)], sample.names.boot)
        )

        # Triming up all data
        for (j in 1:n) {
          tails <- j * nrow(JT.mat) / n
          heads <- tails - (nrow(JT.mat) / n - 1)
          simulated_as[, j] <- as.numeric(JT.mat[heads:tails, 2])
        }

        # Selecting only nodes which are not all zero
        simulated_as_sig <- simulated_as[which(rowSums(simulated_as) != 0), ]

        # Building activity score matrix and majority vote matrix
        activity.scores.samples[rownames(simulated_as_sig), ] <- activity.scores.samples[rownames(simulated_as_sig), ] + simulated_as_sig / degree(RWR_DFS_graph)[rownames(simulated_as_sig)]

        # Recording the subgraph for each DFS
        RWR_JT_graph <- igraph::induced.subgraph(RWR_DFS_graph, rownames(simulated_as_sig))
        edgelist <- rbind(edgelist, igraph::as_data_frame(RWR_JT_graph))
      }
    }

    if (nrow(edgelist) > 0) {
      activity.scores.samples <- activity.scores.samples[unique(c(edgelist[, "from"], edgelist[, "to"])), ]

      # Outputting results
      message("CRNA successfully run. Graph results are returned.\n")
      edgelist.output <- unique(edgelist)
      scores.samples.output <- activity.scores.samples
    } else {
      # Outputting results
      message("CRNA could not run (Unconnected DFS graph). NA results are returned\n")
      edgelist.output <- NA
      scores.samples.output <- NA
    }
  }

  # ==============================================================================#
  # Returning replication results
  # ==============================================================================#
  return(list(
    "edgelist" = edgelist.output,
    "scores.samples" = scores.samples.output
  ))
}
# ===============================================================================================================================#




# ===============================================================================================================================#
# ==============#
# Usage         :
# ==============#
#                   crna.boot.as(graph,
#                                express,
#                                factor,
#                                nboot,
#                                alpha,
#                                quant,
#                                fdr,
#                                conf,
#                                parallel,
#                                seed)
#
# ==============#
# Description   :
# ==============#
#                   Bootstrap replication function of CRNA at fixed FDR.
#
# ==============#
# Arguments     :
# ==============#
#
# ==============#
# Values        :
# ==============#
#
# ===============================================================================================================================#

crna.boot.as <- function(graph,
                         express,
                         factor,
                         nboot,
                         alpha,
                         quant,
                         fdr,
                         conf,
                         parallel,
                         seed) {
  # ==============================================================================#
  # Definition of the treatment factor
  # ==============================================================================#

  factor.lev <- levels(factor)
  factor.ng <- nlevels(factor)
  factor.def <- vector(mode = "list", length = factor.ng)
  factor.tab <- numeric(factor.ng)
  for (g in 1:factor.ng) {
    factor.def[[g]] <- which(factor == factor.lev[g])
    factor.tab[g] <- length(factor.def[[g]])
  }

  # ==============================================================================#
  # Simplify the graph
  # ==============================================================================#
  graph <- graph %>% simplify(edge.attr.comb = "sum")
  E(graph)$sign <- sign(E(graph)$sign)
  # RWR_edges <- igraph::as_data_frame(graph)
  # RWR_edges <- RWR_edges[RWR_edges$sign != 0, ]
  # graph <- graph_from_data_frame(RWR_edges)
  # ==============================================================================#
  # Intersection of expression matrix variables and graph nodes
  # Common graph size and dimensionality of the data
  # ==============================================================================#
  nodes.names <- igraph::V(graph)$name
  inter_names <- intersect(x = rownames(express), y = nodes.names)
  express <- express[inter_names, ]
  graph <- delete_vertices(graph = graph, v = setdiff(x = nodes.names, y = inter_names))

  # ==============================================================================#
  # Convert the graph to an undirected graph, and then to an adjacency matrix
  # ==============================================================================#
  graph.ud <- igraph::as_adjacency_matrix(igraph::as.undirected(graph), attr = "weight")

  # ==============================================================================#
  # Bootstrap replication loop for a fixed FDR value
  # ==============================================================================#
  if (!parallel) {
    R <- length(nboot)
    edgelist.output <- vector(mode = "list", length = R)
    scores.samples.output <- vector(mode = "list", length = R)
    r <- 1
    while (r <= R) {
      if (R > 1) {
        cat("Replication: ", r, "/", R, "\n", sep = "")
      }
      cat("seed: ", seed[r], "\n", sep = "")
      if (!is.null(seed[r])) {
        set.seed(seed[r])
      }

      # ==============================================================================#
      # Define the Bootstrap quantities
      if (R == 1) {
        # No Bootstrap resampling
        samples.boot <- c(factor.def[[1]], factor.def[[2]])
      } else {
        # Bootstrap resampling
        samples.boot <- c(
          sample(x = factor.def[[1]], size = factor.tab[[1]], replace = TRUE),
          sample(x = factor.def[[2]], size = factor.tab[[2]], replace = TRUE)
        )
        # Make sure all groups qualify
        while (length(intersect(samples.boot, factor.def[[2]])) < factor.tab[[2]] / 2 |
          length(intersect(samples.boot, factor.def[[1]])) < factor.tab[[1]] / 2) {
          samples.boot <- c(
            sample(x = factor.def[[1]], size = factor.tab[[1]], replace = TRUE),
            sample(x = factor.def[[2]], size = factor.tab[[2]], replace = TRUE)
          )
        }
      }
      # Bootstrap quantities
      express.boot <- express[, samples.boot, drop = FALSE]
      factor.boot <- factor[samples.boot]
      factor.boot.lev <- levels(factor.boot)
      factor.boot.ng <- nlevels(factor.boot)
      factor.boot.def <- vector(mode = "list", length = factor.boot.ng)
      factor.boot.tab <- numeric(factor.boot.ng)
      for (g in 1:factor.boot.ng) {
        factor.boot.def[[g]] <- which(factor.boot == factor.boot.lev[g])
        factor.boot.tab[g] <- length(factor.boot.def[[g]])
      }
      n <- ncol(express.boot)
      p <- nrow(express.boot)
      sample.names.boot <- sample.names[samples.boot]

      # ==============================================================================#
      # Differential expression for a fixed (user-defined) FDR value
      DEGs <- deg.fdr(fit = deg(express = express.boot, factor = factor.boot), fdr = fdr)$degs[[1]]

      # ==============================================================================#
      # CRNA

      #---------------------------------- RWR ---------------------------------#
      affinity.scores <- rwr(
        graph = graph.ud,
        seedSet = DEGs,
        r = alpha,
        random.seed = FALSE,
        n.thread = NULL
      )
      w <- which(affinity.scores >= quantile(x = affinity.scores, probs = quant))
      sig_nodes_names <- names(affinity.scores[w])

      # Generating network with significant nodes only and clean up
      RWR_graph <- induced.subgraph(graph, sig_nodes_names)
      RWR_graph <- igraph::simplify(
        graph = RWR_graph,
        remove.multiple = TRUE,
        remove.loops = TRUE,
        edge.attr.comb = "max"
      )
      RWR_edges <- igraph::as_data_frame(RWR_graph)
      RWR_edges <- RWR_edges[RWR_edges$sign != 0, ]
      RWR_graph <- graph_from_data_frame(RWR_edges)

      #---------------------------------- DFS ---------------------------------#
      # Initializing DFS result receptor
      edgelist <- matrix(data = NA, nrow = 0, ncol = 4)
      colnames(edgelist) <- c("from", "to", "edgesign", "weight")
      activity.scores.samples <- matrix(data = 0, nrow = vcount(RWR_graph), ncol = length(sample.names.boot))
      rownames(activity.scores.samples) <- V(RWR_graph)$name
      colnames(activity.scores.samples) <- sample.names.boot

      # Scaling expression data
      express.scaled <- t(scale(t(express.boot), center = TRUE, scale = TRUE))

      # DFS graph for every candidate
      sig_nodes_names <- intersect(sig_nodes_names, igraph::V(RWR_graph)$name)

      for (i in 1:length(sig_nodes_names)) {
        RWR_DFS_graph <- dfs(sig_nodes_names[i], RWR_graph, DEGs)
        RWR_DFS_edges <- igraph::as_data_frame(RWR_DFS_graph)

        if (nrow(RWR_DFS_edges) > 0) {
          #-------------------------- Inputs for Junction Tree ------------------------#
          # Nodes and edges description
          RWR_DFS_nodes <- as.data.frame(as.character(V(RWR_DFS_graph)$name))
          RWR_DFS_nodes[, 2] <- as.character(NA)
          colnames(RWR_DFS_nodes) <- c("id", "type")
          RWR_DFS_nodes[, 2] <- "protein"
          RWR_DFS_nodes <- RWR_DFS_nodes[, c(2, 1)]
          RWR_DFS_edges <- RWR_DFS_edges[, c("from", "to", "sign")]
          RWR_DFS_edges[which(RWR_DFS_edges$sign == 1), 3] <- "-a>"
          RWR_DFS_edges[which(RWR_DFS_edges$sign == -1), 3] <- "-a|"

          # Subsetting expression data
          express.scaled.sub <- express.scaled[as.character(RWR_DFS_nodes$id), ]
          rownames(express.scaled.sub) <- as.character(RWR_DFS_nodes$id)
          colnames(express.scaled.sub) <- 1:ncol(express.scaled.sub)
          RWR_DFS_mRNA <- t(express.scaled.sub)

          #-------------------------- Junction Tree Results ---------------------------#
          RWR_DFS_genes <- colnames(RWR_DFS_mRNA)
          RWR_DFS_mRNA <- cbind(sample.names.boot, RWR_DFS_mRNA)
          colnames(RWR_DFS_mRNA) <- paste0("V", seq_len(ncol(RWR_DFS_mRNA)))
          colnames(RWR_DFS_nodes) <- paste0("V", seq_len(ncol(RWR_DFS_nodes)))
          colnames(RWR_DFS_edges) <- paste0("V", seq_len(ncol(RWR_DFS_edges)))
          JT.mat <- CRNA::JTinf(RWR_DFS_nodes, RWR_DFS_edges, RWR_DFS_mRNA, RWR_DFS_genes)
          JT.mat <- JT.mat[-1, ]

          simulated_as <- matrix(
            data = NA,
            nrow = nrow(JT.mat) / n,
            ncol = n,
            dimnames = list(JT.mat$V1[1:(nrow(JT.mat) / n)], sample.names.boot)
          )

          # Triming up all data
          for (j in 1:n) {
            tails <- j * nrow(JT.mat) / n
            heads <- tails - (nrow(JT.mat) / n - 1)
            simulated_as[, j] <- as.numeric(JT.mat[heads:tails, 2])
          }

          # Selecting only nodes which are not all zero
          simulated_as_sig <- simulated_as[which(rowSums(simulated_as) != 0), ]

          # Building activity score matrix and majority vote matrix
          activity.scores.samples[rownames(simulated_as_sig), ] <- activity.scores.samples[rownames(simulated_as_sig), ] + simulated_as_sig / degree(RWR_DFS_graph)[rownames(simulated_as_sig)]

          # Recording the subgraph for each DFS
          RWR_JT_graph <- igraph::induced.subgraph(RWR_DFS_graph, rownames(simulated_as_sig))
          edgelist <- rbind(edgelist, igraph::as_data_frame(RWR_JT_graph))
        }
      }

      if (nrow(edgelist) > 0) {
        activity.scores.samples <- activity.scores.samples[unique(c(edgelist[, "from"], edgelist[, "to"])), ]

        # Outputting results
        message("CRNA successfully run. Graph results are returned.\n")
        edgelist.output[[r]] <- unique(edgelist)
        scores.samples.output[[r]] <- activity.scores.samples
      } else {
        # Outputting results
        message("CRNA could not run (Unconnected DFS graph). NA results are returned\n")
        edgelist.output[[r]] <- NA
        scores.samples.output[[r]] <- NA
      }

      r <- r + 1
    }
  } else {
    # ==============================================================================#
    # Get the process ID of each slave R session
    pid <- Sys.getpid()

    # ==============================================================================#
    # Define the Bootstrap quantities
    if (B == 1) {
      # No Bootstrap resampling
      samples.boot <- c(factor.def[[1]], factor.def[[2]])
    } else {
      # Bootstrap resampling
      samples.boot <- c(
        sample(x = factor.def[[1]], size = factor.tab[[1]], replace = TRUE),
        sample(x = factor.def[[2]], size = factor.tab[[2]], replace = TRUE)
      )
      # Make sure all groups qualify
      while (length(intersect(samples.boot, factor.def[[2]])) < factor.tab[[2]] / 2 |
        length(intersect(samples.boot, factor.def[[1]])) < factor.tab[[1]] / 2) {
        samples.boot <- c(
          sample(x = factor.def[[1]], size = factor.tab[[1]], replace = TRUE),
          sample(x = factor.def[[2]], size = factor.tab[[2]], replace = TRUE)
        )
      }
    }
    # Bootstrap quantities
    express.boot <- express[, samples.boot, drop = FALSE]
    factor.boot <- factor[samples.boot]
    factor.boot.lev <- levels(factor.boot)
    factor.boot.ng <- nlevels(factor.boot)
    factor.boot.def <- vector(mode = "list", length = factor.boot.ng)
    factor.boot.tab <- numeric(factor.boot.ng)
    for (g in 1:factor.boot.ng) {
      factor.boot.def[[g]] <- which(factor.boot == factor.boot.lev[g])
      factor.boot.tab[g] <- length(factor.boot.def[[g]])
    }
    n <- ncol(express.boot)
    p <- nrow(express.boot)
    sample.names.boot <- sample.names[samples.boot]

    # ==============================================================================#
    # Differential expression for a fixed (user-defined) FDR value
    DEGs <- deg.fdr(fit = deg(express = express.boot, factor = factor.boot), fdr = fdr)$degs[[1]]

    # ==============================================================================#
    # CRNA

    #---------------------------------- RWR ---------------------------------#
    affinity.scores <- rwr(
      graph = graph.ud,
      seedSet = DEGs,
      r = alpha,
      random.seed = FALSE,
      n.thread = NULL
    )
    w <- which(affinity.scores >= quantile(x = affinity.scores, probs = quant))
    sig_nodes_names <- names(affinity.scores[w])

    # Generating network with significant nodes only and clean up
    RWR_graph <- induced.subgraph(graph, sig_nodes_names)
    RWR_graph <- igraph::simplify(
      graph = RWR_graph,
      remove.multiple = TRUE,
      remove.loops = TRUE,
      edge.attr.comb = "max"
    )
    RWR_edges <- igraph::as_data_frame(RWR_graph)
    RWR_edges <- RWR_edges[RWR_edges$sign != 0, ]
    RWR_graph <- graph_from_data_frame(RWR_edges)

    #---------------------------------- DFS ---------------------------------#
    # Initializing DFS result receptor
    edgelist <- matrix(data = NA, nrow = 0, ncol = 4)
    colnames(edgelist) <- c("from", "to", "sign", "weight")
    activity.scores.samples <- matrix(data = 0, nrow = vcount(RWR_graph), ncol = length(sample.names.boot))
    rownames(activity.scores.samples) <- V(RWR_graph)$name
    colnames(activity.scores.samples) <- sample.names.boot

    # Scaling expression data
    express.scaled <- t(scale(t(express.boot), center = TRUE, scale = TRUE))

    # DFS graph for every candidate
    sig_nodes_names <- intersect(sig_nodes_names, igraph::V(RWR_graph)$name)

    for (i in 1:length(sig_nodes_names)) {
      RWR_DFS_graph <- dfs(sig_nodes_names[i], RWR_graph, DEGs)
      RWR_DFS_edges <- igraph::as_data_frame(RWR_DFS_graph)

      if (nrow(RWR_DFS_edges) > 0) {
        #-------------------------- Inputs for Junction Tree ------------------------#
        # Nodes and edges description
        RWR_DFS_nodes <- as.data.frame(as.character(V(RWR_DFS_graph)$name))
        RWR_DFS_nodes[, 2] <- as.character(NA)
        colnames(RWR_DFS_nodes) <- c("id", "type")
        RWR_DFS_nodes[, 2] <- "protein"
        RWR_DFS_nodes <- RWR_DFS_nodes[, c(2, 1)]
        RWR_DFS_edges <- RWR_DFS_edges[, c("from", "to", "sign")]
        RWR_DFS_edges[which(RWR_DFS_edges$sign == 1), 3] <- "-a>"
        RWR_DFS_edges[which(RWR_DFS_edges$sign == -1), 3] <- "-a|"

        # Subsetting expression data
        express.scaled.sub <- express.scaled[as.character(RWR_DFS_nodes$id), ]
        rownames(express.scaled.sub) <- as.character(RWR_DFS_nodes$id)
        colnames(express.scaled.sub) <- 1:ncol(express.scaled.sub)
        RWR_DFS_mRNA <- t(express.scaled.sub)

        #-------------------------- Junction Tree Results ---------------------------#
        RWR_DFS_genes <- colnames(RWR_DFS_mRNA)
        RWR_DFS_mRNA <- cbind(sample.names.boot, RWR_DFS_mRNA)
        colnames(RWR_DFS_mRNA) <- paste0("V", seq_len(ncol(RWR_DFS_mRNA)))
        colnames(RWR_DFS_nodes) <- paste0("V", seq_len(ncol(RWR_DFS_nodes)))
        colnames(RWR_DFS_edges) <- paste0("V", seq_len(ncol(RWR_DFS_edges)))
        JT.mat <- CRNA::JTinf(RWR_DFS_nodes, RWR_DFS_edges, RWR_DFS_mRNA, RWR_DFS_genes)
        JT.mat <- JT.mat[-1, ]

        simulated_as <- matrix(
          data = NA,
          nrow = nrow(JT.mat) / n,
          ncol = n,
          dimnames = list(JT.mat$V1[1:(nrow(JT.mat) / n)], sample.names.boot)
        )

        # Triming up all data
        for (j in 1:n) {
          tails <- j * nrow(JT.mat) / n
          heads <- tails - (nrow(JT.mat) / n - 1)
          simulated_as[, j] <- as.numeric(JT.mat[heads:tails, 2])
        }

        # Selecting only nodes which are not all zero
        simulated_as_sig <- simulated_as[which(rowSums(simulated_as) != 0), ]

        # Building activity score matrix and majority vote matrix
        activity.scores.samples[rownames(simulated_as_sig), ] <- activity.scores.samples[rownames(simulated_as_sig), ] + simulated_as_sig / degree(RWR_DFS_graph)[rownames(simulated_as_sig)]

        # Recording the subgraph for each DFS
        RWR_JT_graph <- igraph::induced.subgraph(RWR_DFS_graph, rownames(simulated_as_sig))
        edgelist <- rbind(edgelist, igraph::as_data_frame(RWR_JT_graph))
      }
    }

    if (nrow(edgelist) > 0) {
      activity.scores.samples <- activity.scores.samples[unique(c(edgelist[, "from"], edgelist[, "to"])), ]

      # Outputting results
      message("CRNA successfully run. Graph results are returned.\n")
      edgelist.output <- unique(edgelist)
      scores.samples.output <- activity.scores.samples
    } else {
      # Outputting results
      message("CRNA could not run (Unconnected DFS graph). NA results are returned\n")
      edgelist.output <- NA
      scores.samples.output <- NA
    }
  }

  # ==============================================================================#
  # Returning replication results
  # ==============================================================================#
  return(list(
    "edgelist" = edgelist.output,
    "scores.samples" = scores.samples.output
  ))
}
# ===============================================================================================================================#




# ===============================================================================================================================#
# ==============#
# Usage         :
# ==============#
#                   crna.tuning.sl(graph,
#                                 express,
#                                 factors,
#                                 B=50,
#                                 alpha=seq(0.10, 0.90, by=0.10),
#                                 size=seq(from=1, to=100, by=1),
#                                 fdr=10^(-2),
#                                 lag=2, span=0.40, degree=2, family="gaussian",
#                                 conf=NULL,
#                                 parallel=FALSE,
#                                 seed=NULL)
#
# ==============#
# Description   :
# ==============#
#                   CRNA parameter tuning function at fixed FDR.

#
# ==============#
# Arguments     :
# ==============#
#
# ==============#
# Values        :
# ==============#
#
# ===============================================================================================================================#

crna.tuning.sl <- function(graph,
                           express,
                           factors,
                           B = 50,
                           alpha = seq(0.10, 0.90, by = 0.10),
                           size = seq(from = 1, to = 100, by = 1),
                           fdr = 10^(-2),
                           lag = 2, span = 0.40, degree = 2, family = "gaussian",
                           conf = NULL,
                           parallel = FALSE,
                           seed = NULL) {
  # ==============================================================================#
  # Initializations - Constants - Parameters
  # ==============================================================================#
  n <- ncol(express)
  p <- nrow(express)

  # ==============================================================================#
  # Definitions of the factors
  # ==============================================================================#
  factor1 <- factors[[1]]
  factor1.lev <- levels(factor1)
  factor1.ng <- nlevels(factor1)
  factor1.def <- vector(mode = "list", length = factor1.ng)
  factor1.tab <- numeric(factor1.ng)
  for (g in 1:factor1.ng) {
    factor1.def[[g]] <- which(factor1 == factor1.lev[g])
    factor1.tab[g] <- length(factor1.def[[g]])
  }

  # ==============================================================================#
  # Simplify the graph
  # ==============================================================================#
  graph <- graph %>% simplify(edge.attr.comb = "sum")
  E(graph)$sign <- sign(E(graph)$sign)
  # RWR_edges <- igraph::as_data_frame(graph)
  # RWR_edges <- RWR_edges[RWR_edges$sign != 0, ]
  # graph <- graph_from_data_frame(RWR_edges)
  # ==============================================================================#
  # Intersection of expression matrix variables and graph nodes
  # Common graph size and dimensionality of the data
  # ==============================================================================#
  nodes.names <- igraph::V(graph)$name
  inter_names <- intersect(x = rownames(express), y = nodes.names)
  express <- express[inter_names, ]
  graph <- delete_vertices(graph = graph, v = setdiff(x = nodes.names, y = inter_names))

  # ==============================================================================#
  # Convert the graph to an undirected graph, and then to an adjacency matrix
  # ==============================================================================#
  graph.ud <- igraph::as_adjacency_matrix(igraph::as.undirected(graph), attr = "weight")

  # ==============================================================================#
  # Tuning of RWR parameters size (graph model size) and alpha (restart probability)
  # ==============================================================================#
  # Computation of OOB PE as a function of RWR parameters for fixed FDR
  mce.array <- rwr.tuning.sl(
    graph = graph.ud,
    express = express,
    factor = factor1,
    B = B,
    alpha = alpha,
    size = size,
    fdr = fdr,
    conf = conf,
    parallel = parallel,
    seed = seed
  )
  PE.mu <- round(apply(X = mce.array, MARGIN = c(1, 2, 3), FUN = mean, na.rm = TRUE), 2)
  PE.se <- round(apply(X = mce.array, MARGIN = c(1, 2, 3), FUN = sd, na.rm = TRUE), 2)

  # ==============================================================================#
  # Selection of optimal parameters minimizing the PE surface
  # ==============================================================================#
  alpha.cv <- rwr.id.sl(pe.mu = PE.mu, pe.se = PE.se, size = size, fdr = fdr, lag = lag, span = span, degree = degree, family = family)
  size.cv <- rwr.id.sl(pe.mu = PE.mu, pe.se = PE.se, alpha = alpha, fdr = fdr, lag = lag, span = span, degree = degree, family = family)

  return(list("PE.mu" = PE.mu, "PE.se" = PE.se, "alpha.cv" = alpha.cv, "size.cv" = size.cv))
}
# ===============================================================================================================================#




# ===============================================================================================================================#
# ==============#
# Usage         :
# ==============#
#                   crna.tuning.as(graph,
#                                 express,
#                                 factors,
#                                 B=50,
#                                 P=100,
#                                 alpha=seq(0.10, 0.90, by=0.10),
#                                 quant=seq(0.01, 0.99, by=0.01),
#                                 fdr=10^(-2),
#                                 lag=2, span=0.30, degree=2, family="gaussian",
#                                 conf=NULL,
#                                 parallel=FALSE,
#                                 seed=NULL)
#
# ==============#
# Description   :
# ==============#
#                   CRNA parameter tuning function at fixed FDR.

#
# ==============#
# Arguments     :
# ==============#
#
# ==============#
# Values        :
# ==============#
#
# ===============================================================================================================================#

crna.tuning.as <- function(graph,
                           express,
                           factors,
                           B = 50,
                           P = 100,
                           alpha = seq(0.10, 0.90, by = 0.10),
                           quant = seq(0.01, 0.99, by = 0.01),
                           fdr = 10^(-2),
                           lag = 2, span = 0.30, degree = 2, family = "gaussian",
                           conf = NULL,
                           parallel = FALSE,
                           seed = NULL) {
  # ==============================================================================#
  # Initializations - Constants - Parameters
  # ==============================================================================#
  n <- ncol(express)
  p <- nrow(express)

  # ==============================================================================#
  # Definitions of the factors
  # ==============================================================================#
  factor1 <- factors[[1]]
  factor1.lev <- levels(factor1)
  factor1.ng <- nlevels(factor1)
  factor1.def <- vector(mode = "list", length = factor1.ng)
  factor1.tab <- numeric(factor1.ng)
  for (g in 1:factor1.ng) {
    factor1.def[[g]] <- which(factor1 == factor1.lev[g])
    factor1.tab[g] <- length(factor1.def[[g]])
  }

  # ==============================================================================#
  # Simplify the graph
  # ==============================================================================#
  graph <- graph %>% simplify(edge.attr.comb = "sum")
  E(graph)$sign <- sign(E(graph)$sign)
  # RWR_edges <- igraph::as_data_frame(graph)
  # RWR_edges <- RWR_edges[RWR_edges$sign != 0, ]
  # graph <- graph_from_data_frame(RWR_edges)
  # ==============================================================================#
  # Intersection of expression matrix variables and graph nodes
  # Common graph size and dimensionality of the data
  # ==============================================================================#
  nodes.names <- igraph::V(graph)$name
  inter_names <- intersect(x = rownames(express), y = nodes.names)
  express <- express[inter_names, ]
  graph <- delete_vertices(graph = graph, v = setdiff(x = nodes.names, y = inter_names))

  # ==============================================================================#
  # Convert the graph to an undirected graph, and then to an adjacency matrix
  # ==============================================================================#
  graph.ud <- igraph::as_adjacency_matrix(igraph::as.undirected(graph), attr = "weight")

  # ==============================================================================#
  # Tuning of RWR parameters size (graph model size) and alpha (restart probability)
  # ==============================================================================#
  # Computation of OOB PE as a function of RWR parameters for fixed FDR
  mse.array <- rwr.tuning.as(
    graph = graph.ud,
    express = express,
    factor = factor1,
    B = B,
    P = P,
    alpha = alpha,
    quant = quant,
    fdr = fdr,
    conf = conf,
    parallel = parallel,
    seed = seed
  )
  PE.mu <- round(apply(X = mse.array, MARGIN = c(1, 2, 3), FUN = mean, na.rm = TRUE), 2)
  PE.se <- round(apply(X = mse.array, MARGIN = c(1, 2, 3), FUN = sd, na.rm = TRUE), 2)

  # ==============================================================================#
  # Selection of optimal parameters minimizing the PE surface
  # ==============================================================================#
  alpha.cv <- rwr.id.as(pe.mu = PE.mu, pe.se = PE.se, quant = quant, fdr = fdr, lag = lag, span = span, degree = degree, family = family)
  quant.cv <- rwr.id.as(pe.mu = PE.mu, pe.se = PE.se, alpha = alpha, fdr = fdr, lag = lag, span = span, degree = degree, family = family)

  return(list("PE.mu" = PE.mu, "PE.se" = PE.se, "alpha.cv" = alpha.cv, "quant.cv" = quant.cv))
}
# ===============================================================================================================================#




# ===============================================================================================================================#
# ==============#
# Usage         :
# ==============#
#                   rwr.tuning.sl(graph,
#                                 express,
#                                 factor,
#                                 B,
#                                 alpha,
#                                 size,
#                                 fdr,
#                                 conf,
#                                 parallel,
#                                 seed)
#
# ==============#
# Description   :
# ==============#
#                   RWR tuning function used for the estimation of the
#                   tuning parameter 'alpha' (restart probability) and 'size' (graph size)
#
# ==============#
# Arguments     :
# ==============#
#
# ==============#
# Values        :
# ==============#
#                   The function returns an array of Prediction Error for each parameter value.
#
# ===============================================================================================================================#

#' RWR tuning function used for the estimation of the tuning parameter 'alpha' (restart probability) and 'size' (graph size)
#'
#' @param graph Protein-protein interaction graph.
#' @param express Expression data.
#' @param factor Treatment factors.
#' @param B The number of bootstrap samples.
#' @param alpha The restart probability.
#' @param size The graph size.
#' @param fdr False Discovery Rate.
#' @param conf Cluster configuration.
#' @param parallel Binary value that is either TRUE or FALSE, enabling parallelization.
#' @param seed Seed for reproducibility.
#'
#' @return An array of Prediction Error for each parameter value.
#' @export
#'
#' @examples
rwr.tuning.sl <- function(graph,
                          express,
                          factor,
                          B,
                          alpha,
                          size,
                          fdr,
                          conf,
                          parallel,
                          seed) {
  seed <- seed[1]
  nboot <- 1:B

  if (!parallel) {
    if (!is.null(seed)) {
      seed <- (0:(B - 1)) + seed
    }

    mce.list <- cv.rwr.tuning.sl(
      graph = graph,
      express = express,
      factor = factor,
      nboot = nboot,
      alpha = alpha,
      size = size,
      fdr = fdr,
      parallel = parallel,
      seed = seed
    )
  } else {
    # Parallel backend registration
    # To be used with the 'foreach' and 'doParallel' packages:
    # require(doParallel)
    # doParallel::registerDoParallel(cores=cpus)

    # Parallel backend registration
    # To be used with the 'parallel' package:
    if (conf$type == "SOCKET") {
      clus <- parallel::makeCluster(
        spec = conf$spec,
        type = "PSOCK",
        homogeneous = conf$homo,
        outfile = conf$outfile,
        verbose = conf$verbose
      )
    } else if (conf$type == "MPI") {
      clus <- parallel::makeCluster(
        spec = conf$spec,
        type = "MPI",
        homogeneous = conf$homo,
        outfile = conf$outfile,
        verbose = conf$verbose
      )
    } else {
      stop("Unrecognized cluster type\n")
    }

    clusterEvalQ(cl = clus, expr = library("parallel"))
    clusterEvalQ(cl = clus, expr = library("limma"))
    clusterEvalQ(cl = clus, expr = library("Matrix"))
    clusterEvalQ(cl = clus, expr = library("igraph"))
    clusterEvalQ(cl = clus, expr = library("glmnet"))
    clusterExport(cl = clus, varlist = ls(.GlobalEnv), envir = .GlobalEnv)
    clusterExport(
      cl = clus,
      varlist = c("cv.rwr.tuning.sl", "deg", "deg.fdr", "rwr", "dfs", "cv.class", "mce", "is.empty"),
      envir = .GlobalEnv
    )

    # Bootstrap Cross-Validation of Prediction Error
    # To be used with the 'foreach' and 'doParallel' packages:
    # set.seed(seed)
    # mce.list <- foreach::foreach(b=1:B,
    #                              .inorder=TRUE,
    #                              .packages=c(""parallel","limma","Matrix","igraph"")) %dopar% {cv.rwr.tuning.sl(graph=graph,
    #                                                                                                             express=express,
    #                                                                                                             factor=factor,
    #                                                                                                             alpha=alpha,
    #                                                                                                             size=size,
    #                                                                                                             fdr=fdr,
    #                                                                                                             parallel=parallel,
    #                                                                                                             seed=NULL))}

    # Bootstrap Cross-Validation of Prediction Error
    clusterSetRNGStream(cl = clus, iseed = seed)
    mce.list <- clusterApply(
      cl = clus,
      fun = cv.rwr.tuning.sl,
      x = nboot,
      graph = graph,
      express = express,
      factor = factor,
      alpha = alpha,
      size = size,
      fdr = fdr,
      parallel = parallel,
      seed = NULL
    )

    # Stopping the cluster and cleaning all MPI states
    parallel::stopCluster(cl = clus)
  }

  mce.boot <- array(data = NA, dim = c(dim(mce.list[[1]])[1], dim(mce.list[[1]])[2], dim(mce.list[[1]])[3], 0))
  for (b in 1:B) {
    mce.boot <- abind::abind(mce.boot, mce.list[[b]])
  }

  return(mce.boot)
}
# ===============================================================================================================================#




# ===============================================================================================================================#
# ==============#
# Usage         :
# ==============#
#                   rwr.tuning.as(graph,
#                                express,
#                                factor,
#                                B,
#                                P,
#                                quant,
#                                alpha,
#                                fdr,
#                                conf,
#                                parallel,
#                                seed)
#
# ==============#
# Description   :
# ==============#
#                   RWR tuning function used for the estimation of the
#                   tuning parameter 'alpha' (restart probability) and 'quant' (quantile)
#
# ==============#
# Arguments     :
# ==============#
#
# ==============#
# Values        :
# ==============#
#                   The function returns an array of Prediction Error for each parameter value.
#
# ===============================================================================================================================#

rwr.tuning.as <- function(graph,
                          express,
                          factor,
                          B,
                          P,
                          quant,
                          alpha,
                          fdr,
                          conf,
                          parallel,
                          seed) {
  seed <- seed[1]
  nboot <- 1:B

  if (!parallel) {
    if (!is.null(seed)) {
      seed <- (0:(B - 1)) + seed
    }

    mse.list <- cv.rwr.tuning.as(
      graph = graph,
      express = express,
      factor = factor,
      nboot = nboot,
      P = P,
      alpha = alpha,
      quant = quant,
      fdr = fdr,
      parallel = parallel,
      seed = seed
    )
  } else {
    # Parallel backend registration
    # To be used with the 'foreach' and 'doParallel' packages:
    # require(doParallel)
    # doParallel::registerDoParallel(cores=cpus)

    # Parallel backend registration
    # To be used with the 'parallel' package:
    if (conf$type == "SOCKET") {
      clus <- parallel::makeCluster(
        spec = conf$spec,
        type = "PSOCK",
        homogeneous = conf$homo,
        outfile = conf$outfile,
        verbose = conf$verbose
      )
    } else if (conf$type == "MPI") {
      clus <- parallel::makeCluster(
        spec = conf$spec,
        type = "MPI",
        homogeneous = conf$homo,
        outfile = conf$outfile,
        verbose = conf$verbose
      )
    } else {
      stop("Unrecognized cluster type\n")
    }

    clusterEvalQ(cl = clus, expr = library("parallel"))
    clusterEvalQ(cl = clus, expr = library("limma"))
    clusterEvalQ(cl = clus, expr = library("Matrix"))
    clusterEvalQ(cl = clus, expr = library("igraph"))
    clusterExport(cl = clus, varlist = ls(.GlobalEnv), envir = .GlobalEnv)
    clusterExport(
      cl = clus,
      varlist = c("cv.rwr.tuning.as", "deg", "deg.fdr", "rwr", "dfs", "mse", "is.empty"),
      envir = .GlobalEnv
    )

    # Bootstrap Cross-Validation of Prediction Error
    # To be used with the 'foreach' and 'doParallel' packages:
    # set.seed(seed)
    # mse.list <- foreach::foreach(b=1:B,
    #                              .inorder=TRUE,
    #                              .packages=c(""parallel","limma","Matrix","igraph"")) %dopar% {cv.rwr.tuning.as(graph=graph,
    #                                                                                                             express=express,
    #                                                                                                             factor=factor,
    #                                                                                                             alpha=alpha,
    #                                                                                                             quant=quant,
    #                                                                                                             P=P,
    #                                                                                                             fdr=fdr,
    #                                                                                                             parallel=parallel,
    #                                                                                                             seed=NULL))}

    # Bootstrap Cross-Validation of Prediction Error
    clusterSetRNGStream(cl = clus, iseed = seed)
    mse.list <- clusterApply(
      cl = clus,
      fun = cv.rwr.tuning.as,
      x = nboot,
      graph = graph,
      express = express,
      factor = factor,
      P = P,
      alpha = alpha,
      quant = quant,
      fdr = fdr,
      parallel = parallel,
      seed = NULL
    )

    # Stopping the cluster and cleaning all MPI states
    parallel::stopCluster(cl = clus)
  }

  mse.boot <- array(data = NA, dim = c(dim(mse.list[[1]])[1], dim(mse.list[[1]])[2], dim(mse.list[[1]])[3], 0))
  for (b in 1:B) {
    mse.boot <- abind::abind(mse.boot, mse.list[[b]])
  }

  return(mse.boot)
}
# ===============================================================================================================================#




# ===============================================================================================================================#
# ==============#
# Usage         :
# ==============#
#                   cv.rwr.tuning.sl(graph,
#                                    express,
#                                    factor,
#                                    nboot,
#                                    alpha,
#                                    size,
#                                    fdr,
#                                    parallel,
#                                    seed)
#
# ==============#
# Description   :
# ==============#
#                   Computes the Out-Of-Bootstrap Cross-Validation Missclassification Error (MCE)
#                   between test set predicted and observed sample labels for each parameter value.
#
# ==============#
# Arguments     :
# ==============#
#
# ==============#
# Values        :
# ==============#
#                   Returns an array of MCE.
#
# ===============================================================================================================================#

cv.rwr.tuning.sl <- function(graph,
                             express,
                             factor,
                             nboot,
                             alpha,
                             size,
                             fdr,
                             parallel,
                             seed) {
  # ==============================================================================#
  # Initializations - Constants - Parameters
  # ==============================================================================#
  n <- ncol(express)
  p <- nrow(express)

  # ==============================================================================#
  # Definition of the treatment factor
  # ==============================================================================#
  factor.lev <- levels(factor)
  factor.ng <- nlevels(factor)
  factor.def <- vector(mode = "list", length = factor.ng)
  factor.tab <- numeric(factor.ng)
  for (g in 1:factor.ng) {
    factor.def[[g]] <- which(factor == factor.lev[g])
    factor.tab[g] <- length(factor.def[[g]])
  }

  # ==============================================================================#
  # OOB loop for a fixed FDR value
  # ==============================================================================#
  if (!parallel) {
    B <- length(nboot)
    MCE <- vector(mode = "list", length = B)
    b <- 1
    while (b <= B) {
      cat("Bootstrap: ", b, "\n")
      cat("seed: ", seed[b], "\n", sep = "")
      if (!is.null(seed[b])) {
        set.seed(seed[b])
      }

      # ==============================================================================#
      # Define the bootstrap quantities
      # Bootstrap samples (training and test sets)
      samples.train <- c(
        sample(x = factor.def[[1]], size = factor.tab[[1]], replace = TRUE),
        sample(x = factor.def[[2]], size = factor.tab[[2]], replace = TRUE)
      )
      samples.test <- setdiff(x = 1:n, samples.train)

      # Make sure all groups qualify
      while (length(intersect(samples.train, factor.def[[2]])) < factor.tab[[2]] / 2 |
        length(intersect(samples.train, factor.def[[1]])) < factor.tab[[1]] / 2 |
        length(intersect(samples.test, factor.def[[2]])) < factor.tab[[2]] / 4 |
        length(intersect(samples.test, factor.def[[1]])) < factor.tab[[1]] / 4) {
        samples.train <- c(
          sample(x = factor.def[[1]], size = factor.tab[[1]], replace = TRUE),
          sample(x = factor.def[[2]], size = factor.tab[[2]], replace = TRUE)
        )
        samples.test <- setdiff(x = 1:n, samples.train)
      }

      # Training set
      express.train <- express[, samples.train, drop = FALSE]
      factor.train <- factor[samples.train]

      # Test set
      express.test <- express[, samples.test, drop = FALSE]
      factor.test <- factor[samples.test]

      # ==============================================================================#
      # Initialize the MCE array
      MCE[[b]] <- array(
        data = NA,
        dim = c(length(size), length(alpha), length(fdr)),
        dimnames = list(as.character(size), as.character(alpha), as.character(fdr))
      )

      # Observed test set sample labels
      observed_test_labels <- as.numeric(factor[samples.test])

      # ==============================================================================#
      # Outer loop: for each FDR value
      for (f in 1:length(fdr)) {
        cat("fdr: ", fdr[f], "\n")

        # ==============================================================================#
        # Top d observed all samples DEGs at fixed FDR
        all_DEGs <- deg.fdr(fit = deg(express = express, factor = factor), fdr = fdr[f])$degs[[1]]
        d <- length(all_DEGs)

        # ==============================================================================#
        # Differential expression analysis to get observed train set DEGs
        observed_train_DEGs <- deg.fdr(fit = deg(express = express.train, factor = factor.train), fdr = 1)$degs[[1]]

        # Top d observed train set DEGs names that are in the observed graph
        observed_train_DEGs <- observed_train_DEGs[1:d]

        # ==============================================================================#
        # Inner loop: for each restart probability value
        for (a in 1:length(alpha)) {
          cat("alpha: ", alpha[a], "\n")

          # ==============================================================================#
          # RWR on train set seeds (DEGS) to get train observed affinity scores
          observed_train_as <- rwr(
            graph = graph,
            r = alpha[a],
            seedSet = observed_train_DEGs,
            random.seed = FALSE,
            n.thread = NULL
          )
          observed_sorted_train_as <- sort(observed_train_as, decreasing = TRUE)

          # ==============================================================================#
          # Second loop: for each size of the model
          for (q in 1:length(size)) {
            cat("size: ", size[q], "\n")

            # ==============================================================================#
            # Get the names of observed train set nodes whose affinity scores are in the top selected nodes
            observed_train_nodes <- names(observed_sorted_train_as[1:size[q]])

            # ==============================================================================#
            # Build a classifier on the training set
            # Get the classifier predicted test set labels
            object <- tryCatch(
              {
                cv.class(
                  x = t(express[observed_train_nodes, samples.train, drop = FALSE]),
                  y = as.numeric(factor.train),
                  K = 3,
                  onese = FALSE,
                  nlambda = 100
                )
              },
              error = function(w) {
                NULL
              }
            )

            if (!is.null(object)) {
              class.fit <- cv.class(
                x = t(express[observed_train_nodes, samples.train, drop = FALSE]),
                y = as.numeric(factor.train),
                K = 3,
                onese = FALSE,
                nlambda = 100
              )
              predicted_test_labels <- as.numeric(predict(
                object = class.fit$fit,
                newx = t(express[observed_train_nodes, samples.test]),
                s = class.fit$lambda.min,
                type = "class"
              ))
            } else {
              class.fit <- NA
              predicted_test_labels <- rep(NA, length(samples.test))
            }

            # ==============================================================================#
            # Misclassification Error
            MCE[[b]][q, a, f] <- mce(
              observed = observed_test_labels,
              predicted = predicted_test_labels
            )
          }
        }
      }

      b <- b + 1
    }
  } else {
    # ==============================================================================#
    # Define the bootstrap quantities
    # Bootstrap samples (training and test sets)
    samples.train <- c(
      sample(x = factor.def[[1]], size = factor.tab[[1]], replace = TRUE),
      sample(x = factor.def[[2]], size = factor.tab[[2]], replace = TRUE)
    )
    samples.test <- setdiff(x = 1:n, samples.train)

    # Make sure all groups qualify
    while (length(intersect(samples.train, factor.def[[2]])) < factor.tab[[2]] / 2 |
      length(intersect(samples.train, factor.def[[1]])) < factor.tab[[1]] / 2 |
      length(intersect(samples.test, factor.def[[2]])) < factor.tab[[2]] / 4 |
      length(intersect(samples.test, factor.def[[1]])) < factor.tab[[1]] / 4) {
      samples.train <- c(
        sample(x = factor.def[[1]], size = factor.tab[[1]], replace = TRUE),
        sample(x = factor.def[[2]], size = factor.tab[[2]], replace = TRUE)
      )
      samples.test <- setdiff(x = 1:n, samples.train)
    }

    # Training set
    express.train <- express[, samples.train, drop = FALSE]
    factor.train <- factor[samples.train]

    # Test set
    express.test <- express[, samples.test, drop = FALSE]
    factor.test <- factor[samples.test]

    # ==============================================================================#
    # Initialize the MCE array
    MCE <- array(
      data = NA,
      dim = c(length(size), length(alpha), length(fdr)),
      dimnames = list(as.character(size), as.character(alpha), as.character(fdr))
    )

    # Observed test set sample labels
    observed_test_labels <- as.numeric(factor[samples.test])

    # ==============================================================================#
    # Outer loop: for each FDR value
    for (f in 1:length(fdr)) {
      cat("fdr: ", fdr[f], "\n")

      # ==============================================================================#
      # Top d observed all samples DEGs at fixed FDR
      all_DEGs <- deg.fdr(fit = deg(express = express, factor = factor), fdr = fdr[f])$degs[[1]]
      d <- length(all_DEGs)

      # ==============================================================================#
      # Differential expression analysis to get observed train set DEGs
      observed_train_DEGs <- deg.fdr(fit = deg(express = express.train, factor = factor.train), fdr = 1)$degs[[1]]

      # Top d observed train set DEGs names that are in the observed graph
      observed_train_DEGs <- observed_train_DEGs[1:d]

      # ==============================================================================#
      # Inner loop: for each restart probability value
      for (a in 1:length(alpha)) {
        cat("alpha: ", alpha[a], "\n")

        # ==============================================================================#
        # RWR on train set seeds (DEGS) to get train observed affinity scores
        observed_train_as <- rwr(
          graph = graph,
          r = alpha[a],
          seedSet = observed_train_DEGs,
          random.seed = FALSE,
          n.thread = NULL
        )
        observed_sorted_train_as <- sort(observed_train_as, decreasing = TRUE)

        # ==============================================================================#
        # Second loop: for each size of the model
        for (q in 1:length(size)) {
          cat("size: ", size[q], "\n")

          # ==============================================================================#
          # Get the names of observed train set nodes whose affinity scores are in the top selected nodes
          observed_train_nodes <- names(observed_sorted_train_as[1:size[q]])

          # ==============================================================================#
          # Build a classifier on the training set
          # Get the classifier predicted test set labels
          object <- tryCatch(
            {
              cv.class(
                x = t(express[observed_train_nodes, samples.train, drop = FALSE]),
                y = as.numeric(factor.train),
                K = 3,
                onese = FALSE,
                nlambda = 100
              )
            },
            error = function(w) {
              NULL
            }
          )

          if (!is.null(object)) {
            class.fit <- cv.class(
              x = t(express[observed_train_nodes, samples.train, drop = FALSE]),
              y = as.numeric(factor.train),
              K = 3,
              onese = FALSE,
              nlambda = 100
            )
            predicted_test_labels <- as.numeric(predict(
              object = class.fit$fit,
              newx = t(express[observed_train_nodes, samples.test]),
              s = class.fit$lambda.min,
              type = "class"
            ))
          } else {
            class.fit <- NA
            predicted_test_labels <- rep(NA, length(samples.test))
          }

          # ==============================================================================#
          # Misclassification Error
          MCE[q, a, f] <- mce(
            observed = observed_test_labels,
            predicted = predicted_test_labels
          )
        }
      }
    }
  }

  # ==============================================================================#
  # Returning OOB results
  # ==============================================================================#
  return(MCE)
}
# ===============================================================================================================================#




# ===============================================================================================================================#
# ==============#
# Usage         :
# ==============#
#                   cv.rwr.tuning.as(graph,
#                                    express,
#                                    factor,
#                                    nboot,
#                                    P,
#                                    alpha,
#                                    quant,
#                                    fdr,
#                                    parallel,
#                                    seed)
#
# ==============#
# Description   :
# ==============#
#                   Computes the Out-Of-Bootstrap Cross-Validation Mean Squared Error (MSE)
#                   between test set predicted and observed affinity scores (as) for each parameter value.
#
# ==============#
# Arguments     :
# ==============#
#
# ==============#
# Values        :
# ==============#
#                   Returns an array of MSE.
#
# ===============================================================================================================================#

cv.rwr.tuning.as <- function(graph,
                             express,
                             factor,
                             nboot,
                             P,
                             alpha,
                             quant,
                             fdr,
                             parallel,
                             seed) {
  # ==============================================================================#
  # Initializations - Constants - Parameters
  # ==============================================================================#
  n <- ncol(express)
  p <- nrow(express)

  # ==============================================================================#
  # Definition of the treatment factor
  # ==============================================================================#
  factor.lev <- levels(factor)
  factor.ng <- nlevels(factor)
  factor.def <- vector(mode = "list", length = factor.ng)
  factor.tab <- numeric(factor.ng)
  for (g in 1:factor.ng) {
    factor.def[[g]] <- which(factor == factor.lev[g])
    factor.tab[g] <- length(factor.def[[g]])
  }

  # ==============================================================================#
  # OOB loop for a fixed FDR value
  # ==============================================================================#
  if (!parallel) {
    B <- length(nboot)
    MSE <- vector(mode = "list", length = B)
    b <- 1
    while (b <= B) {
      cat("Bootstrap: ", b, "\n")
      cat("seed: ", seed[b], "\n", sep = "")
      if (!is.null(seed[b])) {
        set.seed(seed[b])
      }

      # ==============================================================================#
      # Define the bootstrap quantities
      # Bootstrap samples (training and test sets)
      samples.train <- c(
        sample(x = factor.def[[1]], size = factor.tab[[1]], replace = TRUE),
        sample(x = factor.def[[2]], size = factor.tab[[2]], replace = TRUE)
      )
      samples.test <- setdiff(x = 1:n, samples.train)

      # make sure all groups qualify
      while (length(intersect(samples.train, factor.def[[2]])) < factor.tab[[2]] / 2 |
        length(intersect(samples.train, factor.def[[1]])) < factor.tab[[1]] / 2 |
        length(intersect(samples.test, factor.def[[2]])) < factor.tab[[2]] / 4 |
        length(intersect(samples.test, factor.def[[1]])) < factor.tab[[1]] / 4) {
        samples.train <- c(
          sample(x = factor.def[[1]], size = factor.tab[[1]], replace = TRUE),
          sample(x = factor.def[[2]], size = factor.tab[[2]], replace = TRUE)
        )
        samples.test <- setdiff(x = 1:n, samples.train)
      }

      # Training set
      express.train <- express[, samples.train]
      factor.train <- factor[samples.train]

      # Test set
      express.test <- express[, samples.test]
      factor.test <- factor[samples.test]

      # ==============================================================================#
      # Initialize the MCE array
      MSE[[b]] <- array(
        data = NA,
        dim = c(length(quant), length(alpha), length(fdr)),
        dimnames = list(as.character(quant), as.character(alpha), as.character(fdr))
      )

      # ==============================================================================#
      # Outer loop: for each FDR value
      for (f in 1:length(fdr)) {
        cat("fdr: ", fdr[f], "\n")

        # ==============================================================================#
        # Top d observed all samples DEGs at fixed FDR
        all_DEGs <- deg.fdr(fit = deg(express = express, factor = factor), fdr = fdr[f])$degs[[1]]
        d <- length(all_DEGs)

        # ==============================================================================#
        # Differential expression analysis to get observed test set DEGs
        observed_test_DEGs <- deg.fdr(fit = deg(express = express.test, factor = factor.test), fdr = 1)$degs[[1]]

        # Top d observed test set DEGs names that are in the observed graph
        observed_test_DEGs <- observed_test_DEGs[1:d]

        # ==============================================================================#
        # Inner loop: for each restart probability value
        for (a in 1:length(alpha)) {
          cat("alpha: ", alpha[a], "\n")

          # ==============================================================================#
          # RWR on test set seeds (DEGS) to get test observed affinity scores
          observed_test_affinity_scores <- rwr(
            graph = graph,
            r = alpha[a],
            seedSet = observed_test_DEGs,
            random.seed = FALSE,
            n.thread = NULL
          )

          # ==============================================================================#
          # Permutation loop:
          # Permute training samples and labels P times
          # Get the training null distribution of affinity scores
          affinity_scores_null_distributions <- matrix(0, nrow = p, ncol = P)
          rownames(affinity_scores_null_distributions) <- rownames(graph)

          for (k in 1:P) {
            # Permute the training labels
            shuffled_train <- sample(x = samples.train, size = length(samples.train), replace = FALSE)

            # Differential expression analysis to get shuffled train set DEGs
            fit.shuffled.train <- deg(express = express.train[, shuffled_train], factor = factor.train)

            # Table of top ranked shuffled trained set DEGs
            top.shuffled.train <- limma::topTable(fit = fit.shuffled.train, coef = 1, adjust = "BH", number = p, p.value = 1, sort.by = "p", resort.by = "B")

            # Only keep the nodes that are in the graph
            shuffled_train_DEGs <- intersect(rownames(top.shuffled.train), rownames(graph))
            shuffled_train_DEGs <- shuffled_train_DEGs[1:d]

            # Null training affinity scores
            affinity_scores_null_distributions[, k] <- rwr(
              graph = graph,
              r = alpha[a],
              seedSet = shuffled_train_DEGs,
              random.seed = FALSE,
              n.thread = NULL
            )
          }
          affinity_scores_null_distributions <- as.vector(affinity_scores_null_distributions)

          # ==============================================================================#
          # Second loop: for each probability distribution value
          for (q in 1:length(quant)) {
            cat("quant: ", quant[q], "\n")

            # ==============================================================================#
            # Observed threshold corresponding to the quantile of combined training null affinity scores
            observed_training_threshold <- quantile(x = affinity_scores_null_distributions, probs = quant[q])

            # ==============================================================================#
            # Get the names of predicted test set nodes whose affinity scores are above the threshold
            test_predicted_newseeds <- names(observed_test_affinity_scores)[observed_test_affinity_scores > observed_training_threshold]

            # ==============================================================================#
            # Get predicted RWR graph with new test set nodes as new seeds to get predicted test set affinity scores
            if (length(test_predicted_newseeds) == 0) {
              predicted_test_affinity_scores <- numeric(p)
            } else {
              predicted_test_affinity_scores <- rwr(
                graph = graph,
                r = alpha[a],
                seedSet = test_predicted_newseeds,
                random.seed = FALSE,
                n.thread = NULL
              )
            }

            # ==============================================================================#
            # Mean Squared Error
            MSE[[b]][q, a, f] <- mse(
              observed = observed_test_affinity_scores[test_predicted_newseeds],
              predicted = predicted_test_affinity_scores[test_predicted_newseeds]
            )
          }
        }
      }

      b <- b + 1
    }
  } else {
    # ==============================================================================#
    # Define the bootstrap quantities
    # Bootstrap samples (training and test sets)
    samples.train <- c(
      sample(x = factor.def[[1]], size = factor.tab[[1]], replace = TRUE),
      sample(x = factor.def[[2]], size = factor.tab[[2]], replace = TRUE)
    )
    samples.test <- setdiff(x = 1:n, samples.train)

    # make sure all groups qualify
    while (length(intersect(samples.train, factor.def[[2]])) < factor.tab[[2]] / 2 |
      length(intersect(samples.train, factor.def[[1]])) < factor.tab[[1]] / 2 |
      length(intersect(samples.test, factor.def[[2]])) < factor.tab[[2]] / 4 |
      length(intersect(samples.test, factor.def[[1]])) < factor.tab[[1]] / 4) {
      samples.train <- c(
        sample(x = factor.def[[1]], size = factor.tab[[1]], replace = TRUE),
        sample(x = factor.def[[2]], size = factor.tab[[2]], replace = TRUE)
      )
      samples.test <- setdiff(x = 1:n, samples.train)
    }

    # Training set
    express.train <- express[, samples.train]
    factor.train <- factor[samples.train]

    # Test set
    express.test <- express[, samples.test]
    factor.test <- factor[samples.test]

    # ==============================================================================#
    # Initialize the MCE array
    MSE <- array(
      data = NA,
      dim = c(length(quant), length(alpha), length(fdr)),
      dimnames = list(as.character(quant), as.character(alpha), as.character(fdr))
    )

    # ==============================================================================#
    # Outer loop: for each FDR value
    for (f in 1:length(fdr)) {
      cat("fdr: ", fdr[f], "\n")

      # ==============================================================================#
      # Top d observed all samples DEGs at fixed FDR
      all_DEGs <- deg.fdr(fit = deg(express = express, factor = factor), fdr = fdr[f])$degs[[1]]
      d <- length(all_DEGs)

      # ==============================================================================#
      # Differential expression analysis to get observed test set DEGs
      observed_test_DEGs <- deg.fdr(fit = deg(express = express.test, factor = factor.test), fdr = 1)$degs[[1]]

      # Top d observed test set DEGs names that are in the observed graph
      observed_test_DEGs <- observed_test_DEGs[1:d]

      # ==============================================================================#
      # Inner loop: for each restart probability value
      for (a in 1:length(alpha)) {
        cat("alpha: ", alpha[a], "\n")

        # ==============================================================================#
        # RWR on test set seeds (DEGS) to get test observed affinity scores
        observed_test_affinity_scores <- rwr(
          graph = graph,
          r = alpha[a],
          seedSet = observed_test_DEGs,
          random.seed = FALSE,
          n.thread = NULL
        )

        # ==============================================================================#
        # Permutation loop:
        # Permute training samples and labels P times
        # Get the training null distribution of affinity scores
        affinity_scores_null_distributions <- matrix(0, nrow = p, ncol = P)
        rownames(affinity_scores_null_distributions) <- rownames(graph)

        for (k in 1:P) {
          # Permute the training labels
          shuffled_train <- sample(x = samples.train, size = length(samples.train), replace = FALSE)

          # Differential expression analysis to get shuffled train set DEGs
          fit.shuffled.train <- deg(express = express.train[, shuffled_train], factor = factor.train)

          # Table of top ranked shuffled trained set DEGs
          top.shuffled.train <- limma::topTable(fit = fit.shuffled.train, coef = 1, adjust = "BH", number = p, p.value = 1, sort.by = "p", resort.by = "B")

          # Only keep the nodes that are in the graph
          shuffled_train_DEGs <- intersect(rownames(top.shuffled.train), rownames(graph))
          shuffled_train_DEGs <- shuffled_train_DEGs[1:d]

          # Null training affinity scores
          affinity_scores_null_distributions[, k] <- rwr(
            graph = graph,
            r = alpha[a],
            seedSet = shuffled_train_DEGs,
            random.seed = FALSE,
            n.thread = NULL
          )
        }
        affinity_scores_null_distributions <- as.vector(affinity_scores_null_distributions)

        # ==============================================================================#
        # Second loop: for each probability distribution value
        for (q in 1:length(quant)) {
          cat("quant: ", quant[q], "\n")

          # ==============================================================================#
          # Observed threshold corresponding to the quantile of combined training null affinity scores
          observed_training_threshold <- quantile(x = affinity_scores_null_distributions, probs = quant[q])

          # ==============================================================================#
          # Get the names of predicted test set nodes whose affinity scores are above the threshold
          test_predicted_newseeds <- names(observed_test_affinity_scores)[observed_test_affinity_scores > observed_training_threshold]

          # ==============================================================================#
          # Get predicted RWR graph with new test set nodes as new seeds to get predicted test set affinity scores
          if (length(test_predicted_newseeds) == 0) {
            predicted_test_affinity_scores <- numeric(p)
          } else {
            predicted_test_affinity_scores <- rwr(
              graph = graph,
              r = alpha[a],
              seedSet = test_predicted_newseeds,
              random.seed = FALSE,
              n.thread = NULL
            )
          }

          # ==============================================================================#
          # Mean Squared Error
          MSE[q, a, f] <- mse(
            observed = observed_test_affinity_scores[test_predicted_newseeds],
            predicted = predicted_test_affinity_scores[test_predicted_newseeds]
          )
        }
      }
    }
  }

  # ==============================================================================#
  # Returning OOB results
  # ==============================================================================#
  return(MSE)
}
# ===============================================================================================================================#




# ===============================================================================================================================#
# ==============#
# Usage         :
# ==============#
#                   rwr.id.sl(pe.mu,
#                             pe.se,
#                             alpha=NULL,
#                             size=NULL,
#                             fdr,
#                             lag=2, span=0.40, degree=2, family="gaussian")
#
# ==============#
# Description   :
# ==============#
#                   Return the minimum of Cross-Validated Prediction Error curve as a function of RWR parameters.
#
# ==============#
# Arguments     :
# ==============#
#
# ==============#
# Values        :
# ==============#
#
# ===============================================================================================================================#

rwr.id.sl <- function(pe.mu,
                      pe.se,
                      alpha = NULL,
                      size = NULL,
                      fdr,
                      lag = 2, span = 0.40, degree = 2, family = "gaussian") {
  q.range <- as.numeric(dimnames(pe.mu)[[1]])
  a.range <- as.numeric(dimnames(pe.mu)[[2]])
  f.range <- as.numeric(dimnames(pe.mu)[[3]])

  if ((is.null(size)) && !(is.null(alpha))) {
    q.min <- matrix(
      data = NA, nrow = length(a.range), ncol = length(f.range),
      dimnames = list(a.range, f.range)
    )
    q.1se <- matrix(
      data = NA, nrow = length(a.range), ncol = length(f.range),
      dimnames = list(a.range, f.range)
    )
    i.min <- matrix(
      data = NA, nrow = length(a.range), ncol = length(f.range),
      dimnames = list(a.range, f.range)
    )
    i.1se <- matrix(
      data = NA, nrow = length(a.range), ncol = length(f.range),
      dimnames = list(a.range, f.range)
    )
    for (j in 1:length(fdr)) {
      if (length(f.range) == 1) {
        f <- 1
      } else {
        f <- pmatch(x = fdr[j], table = f.range)
        f <- f[!is.na(f)]
      }
      for (i in 1:length(alpha)) {
        if (length(a.range) == 1) {
          a <- 1
        } else {
          a <- pmatch(x = alpha[i], table = a.range)
          a <- a[!is.na(a)]
        }
        z <- zeroslope(y = pe.mu[, a, f], x = q.range, lag = lag, span = span, degree = degree, family = family, minimum = TRUE)$min
        if (is.empty(z)) {
          i.min[a, f] <- which.min(pe.mu[, a, f])
        } else {
          i.min[a, f] <- z
        }
        w <- which(pe.mu[, a, f] <= pe.mu[i.min[a, f], a, f] + pe.se[i.min[a, f], a, f])
        if (i.min[a, f] == length(q.range)) {
          i.1se[a, f] <- max(w, na.rm = TRUE)
        } else {
          i.1se[a, f] <- min(w, na.rm = TRUE)
        }
        q.min[a, f] <- q.range[i.min[a, f]]
        q.1se[a, f] <- q.range[i.1se[a, f]]
      }
    }
    return(list("qmin" = q.min, "qmin.id" = i.min, "q1se" = q.1se, "q1se.id" = i.1se))
  } else if ((is.null(alpha)) && !(is.null(size))) {
    a.min <- matrix(
      data = NA, nrow = length(q.range), ncol = length(f.range),
      dimnames = list(q.range, f.range)
    )
    a.1se <- matrix(
      data = NA, nrow = length(q.range), ncol = length(f.range),
      dimnames = list(q.range, f.range)
    )
    i.min <- matrix(
      data = NA, nrow = length(q.range), ncol = length(f.range),
      dimnames = list(q.range, f.range)
    )
    i.1se <- matrix(
      data = NA, nrow = length(q.range), ncol = length(f.range),
      dimnames = list(q.range, f.range)
    )
    for (j in 1:length(fdr)) {
      if (length(f.range) == 1) {
        f <- 1
      } else {
        f <- pmatch(x = fdr[j], table = f.range)
        f <- f[!is.na(f)]
      }
      for (i in 1:length(size)) {
        if (length(q.range) == 1) {
          q <- 1
        } else {
          q <- pmatch(x = size[i], table = q.range)
          q <- q[!is.na(q)]
        }
        i.min[q, f] <- which.min(pe.mu[q, , f])
        w <- which(pe.mu[q, , f] <= pe.mu[q, i.min[q, f], f] + pe.se[q, i.min[q, f], f])
        if (i.min[q, f] == length(q.range)) {
          i.1se[q, f] <- max(w, na.rm = TRUE)
        } else {
          i.1se[q, f] <- min(w, na.rm = TRUE)
        }
        a.min[q, f] <- a.range[i.min[q, f]]
        a.1se[q, f] <- a.range[i.1se[q, f]]
      }
    }
    return(list("amin" = a.min, "amin.id" = i.min, "a1se" = a.1se, "a1se.id" = i.1se))
  } else if ((is.null(size)) && (is.null(alpha))) {
    stop("Arguments 'alpha' and 'size' cannot be both NULL) \n")
  } else {
    stop("Something is wrong in the inputs of argument 'alpha' or 'size') \n")
  }
}
# ===============================================================================================================================#




# ===============================================================================================================================#
# ==============#
# Usage         :
# ==============#
#                   rwr.id.as(pe.mu,
#                             pe.se,
#                             alpha=NULL,
#                             quant=NULL,
#                             fdr,
#                             lag=2, span=0.40, degree=2, family="gaussian")
#
# ==============#
# Description   :
# ==============#
#                   Return the minimum of Cross-Validated Prediction Error curve as a function of RWR parameters.
#
# ==============#
# Arguments     :
# ==============#
#
# ==============#
# Values        :
# ==============#
#
# ===============================================================================================================================#

rwr.id.as <- function(pe.mu,
                      pe.se,
                      alpha = NULL,
                      quant = NULL,
                      fdr,
                      lag = 2, span = 0.40, degree = 2, family = "gaussian") {
  q.range <- as.numeric(dimnames(pe.mu)[[1]])
  a.range <- as.numeric(dimnames(pe.mu)[[2]])
  f.range <- as.numeric(dimnames(pe.mu)[[3]])

  if ((is.null(quant)) && !(is.null(alpha))) {
    q.min <- matrix(
      data = NA, nrow = length(a.range), ncol = length(f.range),
      dimnames = list(a.range, f.range)
    )
    q.1se <- matrix(
      data = NA, nrow = length(a.range), ncol = length(f.range),
      dimnames = list(a.range, f.range)
    )
    i.min <- matrix(
      data = NA, nrow = length(a.range), ncol = length(f.range),
      dimnames = list(a.range, f.range)
    )
    i.1se <- matrix(
      data = NA, nrow = length(a.range), ncol = length(f.range),
      dimnames = list(a.range, f.range)
    )
    for (j in 1:length(fdr)) {
      if (length(f.range) == 1) {
        f <- 1
      } else {
        f <- pmatch(x = fdr[j], table = f.range)
        f <- f[!is.na(f)]
      }
      for (i in 1:length(alpha)) {
        if (length(a.range) == 1) {
          a <- 1
        } else {
          a <- pmatch(x = alpha[i], table = a.range)
          a <- a[!is.na(a)]
        }
        z <- zeroslope(y = pe.mu[, a, f], x = q.range, lag = lag, span = span, degree = degree, family = family, minimum = TRUE)$min
        if (is.empty(z)) {
          i.min[a, f] <- which.min(pe.mu[, a, f])
        } else {
          i.min[a, f] <- z
        }
        w <- which(pe.mu[, a, f] <= pe.mu[i.min[a, f], a, f] + pe.se[i.min[a, f], a, f])
        if (i.min[a, f] == length(q.range)) {
          i.1se[a, f] <- max(w, na.rm = TRUE)
        } else {
          i.1se[a, f] <- min(w, na.rm = TRUE)
        }
        q.min[a, f] <- q.range[i.min[a, f]]
        q.1se[a, f] <- q.range[i.1se[a, f]]
      }
    }
    return(list("qmin" = q.min, "qmin.id" = i.min, "q1se" = q.1se, "q1se.id" = i.1se))
  } else if ((is.null(alpha)) && !(is.null(quant))) {
    a.min <- matrix(
      data = NA, nrow = length(q.range), ncol = length(f.range),
      dimnames = list(q.range, f.range)
    )
    a.1se <- matrix(
      data = NA, nrow = length(q.range), ncol = length(f.range),
      dimnames = list(q.range, f.range)
    )
    i.min <- matrix(
      data = NA, nrow = length(q.range), ncol = length(f.range),
      dimnames = list(q.range, f.range)
    )
    i.1se <- matrix(
      data = NA, nrow = length(q.range), ncol = length(f.range),
      dimnames = list(q.range, f.range)
    )
    for (j in 1:length(fdr)) {
      if (length(f.range) == 1) {
        f <- 1
      } else {
        f <- pmatch(x = fdr[j], table = f.range)
        f <- f[!is.na(f)]
      }
      for (i in 1:length(quant)) {
        if (length(q.range) == 1) {
          q <- 1
        } else {
          q <- pmatch(x = quant[i], table = q.range)
          q <- q[!is.na(q)]
        }
        i.min[q, f] <- which.min(pe.mu[q, , f])
        w <- which(pe.mu[q, , f] <= pe.mu[q, i.min[q, f], f] + pe.se[q, i.min[q, f], f])
        if (i.min[q, f] == length(q.range)) {
          i.1se[q, f] <- max(w, na.rm = TRUE)
        } else {
          i.1se[q, f] <- min(w, na.rm = TRUE)
        }
        a.min[q, f] <- a.range[i.min[q, f]]
        a.1se[q, f] <- a.range[i.1se[q, f]]
      }
    }
    return(list("amin" = a.min, "amin.id" = i.min, "a1se" = a.1se, "a1se.id" = i.1se))
  } else if ((is.null(quant)) && (is.null(alpha))) {
    stop("Arguments 'alpha' and 'quant' cannot be both NULL) \n")
  } else {
    stop("Something is wrong in the inputs of argument 'alpha' or 'quant') \n")
  }

  invisible()
}
# ===============================================================================================================================#



# ===============================================================================================================================#
# ===============================================================================================================================#
# ==============#
# Usage         :
# ==============#
#                    deg(express,
#                        factor,
#                        contrasts="contr.treatment")
#
# ==============#
# Description   :
# ==============#
#                   Differential Expression Analysis for a Complete Randomized Design (CRD).
#                   Analysis by moderated one-way ANOVA and Emmpirical Bayes.
#                   Typically, use the default options for "contrasts":
#                      'treatment' for unordered factors
#                      'orthogonal polynomials' for ordered factors.
#                   options(contrasts=c("contr.treatment", "contr.poly"))
#
# ==============#
# Arguments     :
# ==============#
#
# ==============#
# Values        :
# ==============#
#                   Return a limma fit of differentially expressed nodes.
#
# ===============================================================================================================================#

deg <- function(express,
                factor,
                contrasts = "contr.treatment") {
  factor.lev <- levels(factor)
  # Design matrix
  design <- model.matrix(
    object = ~ 1 + factor,
    data = model.frame(~ 1 + factor),
    contrasts.arg = list(factor = contrasts)
  )
  rownames(design) <- factor
  colnames(design) <- c("Intercept", paste(factor.lev[2], "-", factor.lev[1], sep = ""))
  # Contrast matrix
  cont.matrix.factor <- limma::makeContrasts("T-C" = c(0, 1), levels = c("Intercept", "Factor"))
  colnames(cont.matrix.factor) <- paste(factor.lev[2], "-", factor.lev[1], sep = "")
  rownames(cont.matrix.factor) <- colnames(design)
  # Fit the EBLM
  fit <- limma::lmFit(object = express, design = design, block = NULL, correlation = NULL)
  fit <- limma::contrasts.fit(fit, cont.matrix.factor)
  fit <- limma::eBayes(fit, proportion = 0.01, robust = FALSE)

  return(fit)
}
# ===============================================================================================================================#




# ===============================================================================================================================#
# ==============#
# Usage         :
# ==============#
#                    deg.fdr(fit,
#                            fdr)
#
# ==============#
# Description   :
# ==============#
#                   Differential Expression Analysis for a Complete Randomized Design (CRD).
#                   Analysis by moderated one-way ANOVA and Emmpirical Bayes.
#                   Checks for different number of differentially expressed nodes for corresponding fdr values.
#
# ==============#
# Arguments     :
# ==============#
#
# ==============#
# Values        :
# ==============#
#                   Return a list of differentially expressed nodes with corresponding fdr.
#
# ===============================================================================================================================#

deg.fdr <- function(fit,
                    fdr) {
  degs <- vector(mode = "list", length = length(fdr))
  degs.toptable <- vector(mode = "list", length = length(fdr))
  p <- nrow(fit$coefficients)
  d <- numeric(length(fdr))
  for (f in 1:length(fdr)) {
    degs.toptable[[f]] <- limma::topTable(fit = fit, coef = 1, adjust = "BH", number = p, p.value = fdr[f], lfc = 0, sort.by = "p", resort.by = "B")
    degs[[f]] <- rownames(degs.toptable[[f]])
    d[f] <- length(degs[[f]])
  }

  return(list("degstoptable" = degs.toptable, "degs" = degs, "d" = d, "fdr" = fdr))
}
# ===============================================================================================================================#




# ===============================================================================================================================#
# ==============#
# Usage         :
# ==============#
#                   rwr(graph,
#                       seedSet,
#                       r=0.75,
#                       random.seed=FALSE,
#                       n.thread=NULL)
#
# ==============#
# Description   :
# ==============#
#                   Undirected Random Walk w/ Restart graph function, that uses
#                   an original graph of node interaction and a set of seeds,
#                   parametrized by an initial vector of probabilities,
#                   and a restart probability.
#
# ==============#
# Arguments     :
# ==============#
#
# ==============#
# Values        :
# ==============#
#                   Return a vector of normalized (by node degree) affinity scores (steady-state probabilities).
#
# ===============================================================================================================================#

rwr <- function(graph,
                seedSet,
                r = 0.75,
                random.seed = FALSE,
                n.thread = NULL) {
  seedSet <- intersect(seedSet, colnames(graph))
  graph.norm <- as.matrix(graph %*% Matrix::Diagonal(x = Matrix::colSums(graph)^-1))
  nd <- apply(as.matrix(graph), 1, function(x) length(which(x != 0)))

  if (random.seed) { # Random seed scheme

    nSeed <- length(seedSet)

    cl <- parallel::makeCluster(n.thread)
    parallel::clusterEvalQ(cl = cl, expr = library("parallel"))
    parallel::clusterEvalQ(cl = cl, expr = library("Matrix"))
    parallel::clusterExport(
      cl = cl,
      varlist = c("graph.norm", "r", "nSeed"),
      envir = environment()
    )

    res.random <- parallel::parSapply(cl = cl, X = 1:1000, FUN = function() {
      p0 <- Matrix(as.numeric(ifelse(rownames(graph.norm) %in% sample(rownames(graph.norm), nSeed, replace = FALSE), 1, 0)),
        ncol = 1,
        nrow = nrow(graph.norm)
      )
      p0 <- p0 / sum(p0)
      pt <- Matrix(1 / nrow(graph.norm),
        ncol = 1,
        nrow = nrow(graph.norm)
      )
      delta <- 1
      count <- 1
      while (delta > 1e-16 && count < 100) {
        px <- (1 - r) * graph.norm %*% pt + r * p0
        delta <- sum(abs(px - pt))
        count <- count + 1
        pt <- px
      }
      pt <- as.numeric(log2(pt) / (log2(nd + 1)))
      return(as.vector(pt))
    })
    parallel::stopCluster(cl = cl)
    rownames(res.random) <- colnames(as.matrix(graph))

    return(res.random)
  } else { # Original scheme

    p0 <- Matrix(ifelse(rownames(graph.norm) %in% seedSet, 1, 0),
      ncol = 1,
      nrow = nrow(graph.norm)
    )
    p0 <- p0 / sum(p0)
    pt <- Matrix(1 / nrow(graph.norm),
      ncol = 1,
      nrow = nrow(graph.norm)
    )
    delta <- 1
    count <- 1
    while (delta > 1e-16 && count < 100) {
      px <- (1 - r) * graph.norm %*% pt + r * p0
      delta <- sum(abs(px - pt))
      count <- count + 1
      pt <- px
    }
    pt <- as.numeric(log2(pt) / (log2(nd + 1)))
    names(pt) <- colnames(as.matrix(graph))

    return(pt)
  }
}
# ===============================================================================================================================#




# ===============================================================================================================================#
# ==============#
# Usage         :
# ==============#
#                   dfs(Root, Network, degs)
#
# ==============#
# Description   :
# ==============#
#                   Return the Depth First Search graph
#
# ==============#
# Arguments     :
# ==============#
#
# ==============#
# Values        :
# ==============#
#
# ===============================================================================================================================#

dfs <- function(Root,
                Network,
                degs) {
  record_epaths <- c()
  record_npaths <- list()
  npaths_counter <- 1
  passedNodes <- c(Root)
  passedEdges <- c()
  regulatedDEG_count <- 0
  regulatedDEG <- c()

  # stop_point <- rownames(DEG_List)[-which(rownames(DEG_List) %in% V(Network)$name)]
  ### Find r1
  possible_r1 <- setdiff(unique(neighbors(Network, Root, "out")$name), passedNodes) ### That's r1

  ### if Root has childs(r1)
  if (length(possible_r1) > 0) { ### Root has r1
    for (r1_iterator in 1:length(possible_r1)) {
      r1 <- possible_r1[r1_iterator]
      passedNodes <- c(passedNodes, r1)
      possible_r2 <- setdiff(unique(neighbors(Network, r1, "out")$name), passedNodes) ### That's r2
      ### if r1 has childs(r2)
      if (length(possible_r2) > 0) { ### r1 has r2
        for (r2_iterator in 1:length(possible_r2)) {
          r2 <- possible_r2[r2_iterator]
          passedNodes <- c(passedNodes, r2)
          possible_r3 <- setdiff(unique(neighbors(Network, r2, "out")$name), passedNodes) ### That's r2
          ### if r2 has childs(r3)
          if (length(possible_r3) > 0) { ### r2 has r3
            for (r3_iterator in 1:length(possible_r3)) {
              r3 <- possible_r3[r3_iterator]
              passedNodes <- c(passedNodes, r3)
              possible_r4 <- setdiff(unique(neighbors(Network, r3, "out")$name), passedNodes) ### That's r3
              ### if r3 has childs(r4)
              if (length(possible_r4) > 0) { ### r3 has r4
                for (r4_iterator in 1:length(possible_r4)) {
                  r4 <- possible_r4[r4_iterator]
                  passedNodes <- c(passedNodes, r4)
                  possible_r5 <- setdiff(unique(neighbors(Network, r4, "out")$name), passedNodes) ### That's r5

                  if (length(possible_r5) > 0) { ### r3 has r5
                    for (r5_iterator in 1:length(possible_r5)) {
                      r5 <- possible_r5[r5_iterator]
                      passedNodes <- c(passedNodes, r5)
                      record_npaths[[npaths_counter]] <- c(Root, r1, r2, r3, r4, r5)
                      npaths_counter <- npaths_counter + 1
                    }
                  } else {
                    record_npaths[[npaths_counter]] <- c(Root, r1, r2, r3, r4)
                    npaths_counter <- npaths_counter + 1
                  }
                }
              } else {
                record_npaths[[npaths_counter]] <- c(Root, r1, r2, r3)
                npaths_counter <- npaths_counter + 1
              }
            }
          } else {
            record_npaths[[npaths_counter]] <- c(Root, r1, r2)
            npaths_counter <- npaths_counter + 1
          }
        }
      } else {
        record_npaths[[npaths_counter]] <- c(Root, r1)
        npaths_counter <- npaths_counter + 1
      }
    }
  }

  if (length(record_npaths) > 0) {
    for (each_path in 1:length(record_npaths)) {
      while (!tail(record_npaths[[each_path]], n = 1) %in% degs & length(record_npaths[[each_path]]) > 1) {
        record_npaths[[each_path]] <- head(record_npaths[[each_path]], -1)
      }
    }
    record_npaths <- unique(record_npaths[lapply(record_npaths, length) > 1])
  }

  if (length(record_npaths) > 0) {
    for (each_path in 1:length(record_npaths)) {
      vp <- rep(record_npaths[[each_path]], each = 2)
      vp <- vp[-c(1, length(vp))]
      passedEdges <- c(passedEdges, get.edge.ids(Network, vp))
    }
  }

  return(subgraph.edges(Network, eids = unique(passedEdges)))
}
# ===============================================================================================================================#




# ===============================================================================================================================#
# ==============#
# Usage         :
# ==============#
#                   cv.class(x, y, K, onese, nlambda)
#
# ==============#
# Description   :
# ==============#
#                   Cross-Validated Penalized Logistic Regression model.
#                   Logistic regression model is fit with Elasticnet Penalty for p > n situations and variable selection.
#
# ==============#
# Arguments     :
# ==============#
#
# ==============#
# Values        :
# ==============#
#                   The function returns the cross-validated fit, vectors of tunig parameters, and minimizers.
#
# ===============================================================================================================================#

cv.class <- function(x,
                     y,
                     K,
                     onese,
                     nlambda) {
  cv.fit <- glmnet::cv.glmnet(
    x = x,
    y = y,
    lambda = NULL,
    alpha = 0,
    nlambda = nlambda,
    nfolds = pmax(3, K),
    type.measure = "class",
    family = "binomial",
    maxit = 1e5
  )
  cv.errmu <- cv.fit$cvm
  cv.errse <- cv.fit$cvsd
  enlambda <- cv.fit$lambda

  if (all(is.na(cv.errmu)) || any(is.nan(cv.errmu)) || is.empty(cv.errmu)) {
    index.min <- NA
    index.1se <- NA
  } else {
    index.min <- min(which(cv.errmu == min(cv.errmu, na.rm = TRUE)))
    w <- (cv.errmu <= cv.errmu[index.min] + cv.errse[index.min])
    if (all(is.na(w)) || is.empty(w)) {
      index.1se <- NA
    } else {
      index.1se <- min(which(w))
    }
  }
  if (is.empty(index.min) || is.empty(index.1se)) {
    varsel <- NA
  } else {
    enlambda.min <- enlambda[index.min]
    enlambda.1se <- enlambda[index.1se]
    fit <- glmnet::glmnet(
      x = x,
      y = y,
      alpha = 0,
      family = "binomial",
      maxit = 1e5
    )
    w <- apply(X = fit$beta, MARGIN = 2, FUN = function(x) {
      sum(!(is.na(x)) & (x != 0))
    })
    if (all(w == 0)) {
      varsel <- NA
    } else {
      if (onese) {
        cv.coef <- coef(object = fit, s = enlambda.1se)
        lambda.min <- enlambda.1se
        i.min <- index.1se
      } else {
        cv.coef <- coef(object = fit, s = enlambda.min)
        lambda.min <- enlambda.min
        i.min <- index.min
      }
      varsel <- rownames(cv.coef)[vapply(cv.coef, FUN = function(x) {
        !(is.na(x)) & (x != 0)
      }, FUN.VALUE = logical(1))]
      if (is.empty(varsel)) {
        varsel <- NA
      }
    }
  }
  return(list("fit" = fit, "lambda" = enlambda, "i.min" = i.min, "lambda.min" = lambda.min, "varsel" = varsel))
}
# ===============================================================================================================================#




# ===============================================================================================================================#
# ===============#
# Usage         :
# ===============#
#                    zeroslope(y, x, lag=1, span=0.10, degree=2, family="gaussian", minimum=TRUE)
#
# ===============#
# Description   :
# ===============#
#
# ===============#
# Arguments     :
# ===============#
#
# ===============#
# Values        :
# ===============#
#
# ===============================================================================================================================#

zeroslope <- function(y, x, lag, span, degree, family, minimum) {
  if (anyNA(x)) {
    stop("'x' cannot contain NA values. Exiting ... \n\n")
  } else {
    y <- y[order(x)] # reorder the data in ascending values of x
    x <- x[order(x)] # do the same for x
    na <- is.na(y)
    wa <- which(na)
    if (!is.empty(wa)) {
      xc <- x[-wa]
      yc <- y[-wa]
    } else {
      xc <- x
      yc <- y
    }
    fitc <- loess(yc ~ as.numeric(xc), na.action = "na.omit", span = span, degree = degree, family = family)$fitted
    loe <- rep(NA, length(x))
    loe[!na] <- fitc
    d <- diff(x = loe, lag = lag) / diff(x = x, lag = lag)
    d.sign <- diff(x = sign(d), lag = lag)
    zs <- which(d.sign != 0) + lag
    if (minimum) {
      w <- which.min(loe[zs])
    } else {
      w <- which.max(loe[zs])
    }
    return(list("loess" = loe, "min" = zs[w]))
  }
}
# ===============================================================================================================================#




# ===============================================================================================================================#
# ==============#
# Usage         :
# ==============#
#                    mce(observed, predicted)
#
# ==============#
# Description   :
# ==============#
#                   Computes the Misclassification Error of two input vectors: observed vs. predicted
#
# ==============#
# Arguments     :
# ==============#
#
# ==============#
# Values        :
# ==============#
#                   Return a scalar
#
# ===============================================================================================================================#

mce <- function(observed, predicted) {
  return(sum(abs(observed - predicted)))
}
# ===============================================================================================================================#




# ===============================================================================================================================#
# ==============#
# Usage         :
# ==============#
#                    mse(observed, predicted)
#
# ==============#
# Description   :
# ==============#
#                   Computes the Mean Squarred Error of two input vectors: observed vs. predicted
#
# ==============#
# Arguments     :
# ==============#
#
# ==============#
# Values        :
# ==============#
#                   Return a scalar
#
# ===============================================================================================================================#

mse <- function(observed, predicted) {
  return(mean(observed - predicted)^2)
}
# ===============================================================================================================================#




# ===============================================================================================================================#
# ===============#
# Usage         :
# ===============#
#                    is.empty(x)
#
# ===============#
# Description   :
# ===============#
#
# ===============#
# Arguments     :
# ===============#
#
# ===============#
# Values        :
# ===============#
#
# ===============================================================================================================================#

is.empty <- function(x) {
  if (is.null(x)) {
    return(TRUE)
  } else if (is.vector(x, mode = "any")) {
    na <- is.na(x)
    if (length(which(na)) != 0) {
      return(FALSE) # NA is a non-empty value
    } else {
      x <- x[!na]
      if (is.character(x)) {
        if (length(x) == 0) {
          return(TRUE)
        } else if (length(x) > 0) {
          return(all(vapply(as.list(x), function(x) {
            x == ""
          }, logical(1))))
        }
      } else {
        if (length(x) == 0) {
          return(TRUE)
        } else {
          return(FALSE)
        }
      }
    }
  } else if (is.matrix(x) || is.data.frame(x)) {
    return(((nrow(x) == 0) || (ncol(x) == 0)))
  }
}
# ===============================================================================================================================#




#==========================================================================================#
#===============#
# Usage         :
#===============#
#                    CRNA.news()
#
#===============#
# Description   :
#===============#
#                    Function to display the log file NEWS of updates of the CRNA package.
#
#===============#
# Arguments     :
#===============#
#
#===============#
# Values        :
#===============#
#
#==========================================================================================#

CRNA.news <- function(...) {
   
   newsfile <- file.path(system.file(package="CRNA"), "NEWS")
   file.show(newsfile)
   
}
#==========================================================================================#
