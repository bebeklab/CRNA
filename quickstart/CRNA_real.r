knitr::opts_chunk$set(echo = TRUE)
library(CRNA)
load("/home/fxe35/your_workspace.RData") # includes objects like ncbi_mapping
.libPaths( c( "/home/fxe35/R/x86_64-pc-linux-gnu-library/4.2" , .libPaths() ) )
.libPaths( c( .libPaths(), "/usr/local/easybuild_allnodes/software/R/4.2.1-foss-2021b/lib64/R/library", "/usr/local/easybuild_allnodes/software/R-bundle-Bioconductor/3.15-foss-2022a-R-4.2.2", "/usr/local/easybuild_allnodes/software/arrow-R/8.0.0-foss-2022a-R-4.2.2") )
library(edgeR)
if ((.Platform$OS.type == "windows") || (.Platform$OS.type == "unix")) {
  HOME.path <- Sys.getenv("HOME")
  setwd(dir=file.path(HOME.path, "RESULTS/METHODS/CRNA/REAL", fsep = .Platform$file.sep))
} else {
  warning("OS not recognized \n")
}

CODE.path <- "CODES"
RESULT.path <- "RESULTS/METHODS/CRNA/REAL"
DATA.path <- "DATA"


library(GEOquery)
library(DESeq2)

##  START OF CUSTOMIZATION ##

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE149428", "file=GSE149428_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tblt <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")


GEOquery::getGEOSuppFiles('GSE149428',fetch_files = T)

sample_info<-getGEO(filename="~/GSE149428_series_matrix.txt.gz")

smlt <-  sample_info@phenoData@data$geo_accession

Group_Name1="Mefloquine hydrochloride"
Group_Name2="DMSO"
Group_Name3="Tamoxifen citrate"
Group_Name4="Withaferin A"

group_names <- c(Group_Name1, Group_Name2, Group_Name3, Group_Name4)
group_labels <- c("Mef", "WT", "Tamoxifen", "Withaferin")
combinations <- combn(group_names, 2, simplify = FALSE)

## END OF CUSTOMIZATION ##

options(max.print=1e6)
options(digits=12)
options(echo=TRUE)

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
library("OmnipathR")
library("Biobase")
library("ROCR")
library("CausalR")

for (i in 1:length(combinations)) {
    

    # LOOP THE FOLLOWING FOR EVERY POSSIBLE COMBINATION:
    group_1 <- combinations[[i]][1]
    group_2 <- combinations[[i]][2]

    sel2<-c( grep(group_2,sample_info@phenoData@data$`treatment:ch1`) ,
            grep(group_1,sample_info@phenoData@data$`treatment:ch1`)
    )
    groups<-NULL
    groups<-c(rep(group_labels[which(group_names == group_2)], length(grep(group_2,sample_info@phenoData@data$`treatment:ch1`))) , 
            rep(group_labels[which(group_names == group_1)],length(grep(group_1,sample_info@phenoData@data$`treatment:ch1`))))

    gs<-factor(groups)

    sml <- smlt[sel2]

    tbl <- tblt[ ,sel2]

    #dim(tbl)
    #length(sel2)

    # group membership for samples
    sample_info_subset <- data.frame(Group = gs, row.names = colnames(tbl))

    # pre-filter low count genes
    # keep genes with at least N counts > 10, where N = size of smallest group
    keep <- rowSums( tbl >= 10 ) >= min(table(gs))
    tbl <- tbl[keep, ]


    y <- DGEList(counts = tbl, group = sample_info_subset$Group,genes =row.names(tbl) )

    keep <- filterByExpr(y)
    y<- y[keep,]
    y <- calcNormFactors(y)



    if (require("parallel")) {
    print("'parallel' is attached correctly \n")
    } else {
    stop("'parallel' must be attached first \n")
    }

    if (.Platform$OS.type == "windows") {
    # Windows OS
    cpus <- detectCores(logical = TRUE)
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
        cpus <- as.numeric(Sys.getenv("SLURM_NTASKS"))
        if (type == "SOCKET") {
        conf <- list("spec"=rep("localhost", 40),
                    "type"="SOCKET",
                    "homo"=TRUE,
                    "verbose"=TRUE,
                    "outfile"=file.path(getwd(), "output_SOCK.txt", fsep=.Platform$file.sep))
        } else if (type == "MPI") {
        if (require("Rmpi")) {
            print("Rmpi is loaded correctly \n")
        } else {
            stop("Rmpi must be installed first \n")
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
    E <- cpm(y, log=FALSE, prior.count=2)

    gene_data <- strsplit(row.names(E), "\\|")

    # Create a data frame with two columns
    gene_df <- data.frame(
    GeneName = sapply(gene_data, function(x) x[1]),
    EnsemblID = sapply(gene_data, function(x) x[2])
    )
    y$genes<-gene_df$GeneName

    row.names(E)<-y$genes
    # Import the graph data from the data frame and convert it to an 'igraph' object
    PPI <- readRDS("PPI_graph_2021.RDS")
    PPI.ig <- igraph::graph_from_data_frame(d=PPI, directed=TRUE) %>% simplify(edge.attr.comb = "sum")
    E(PPI.ig)$sign <- sign(E(PPI.ig)$sign)
    igraph::vcount(PPI.ig)  # |V|=3743 initial nodes
    igraph::ecount(PPI.ig)  # |E|=12086 initial edges
    
    rownames(E) <- as.character(rownames(E))
    ncbi_mapping$entrezgene_id <- as.character(ncbi_mapping$entrezgene_id)

    # Create a named vector for mapping NCBI Gene IDs to gene names
    id_to_name <- setNames(ncbi_mapping$hgnc_symbol, ncbi_mapping$entrezgene_id)

    # Replace NCBI Gene IDs with gene names where available, otherwise keep the original ID
    new_row_names <- sapply(rownames(E), function(id) {
    if (!is.na(id_to_name[id])) {
        id_to_name[id]
    } else {
        id
    }
    })

    rownames(E) <- new_row_names
    #==========================================================================================#
    # Intersection of expression matrix genes and graph nodes
    # Final Expression Data matrix
    # Final Knowledge Base graph
    #==========================================================================================#
    #----------------------- Updating final Matrix dimensions ---------------------#
    nodes_names <- igraph::V(PPI.ig)$name
    inter_names <- intersect(x=rownames(E), y=nodes_names)
    E <- E[inter_names,]
    n <- ncol(E)  # n=27
    p <- nrow(E)  # p=3053
    #--------------------- Updating final Graph characteristics -------------------#
    PPI.ig <- delete_vertices(graph=PPI.ig, v=setdiff(x=nodes_names, y=inter_names))
    e <- igraph::ecount(PPI.ig)  # e=9968 final edges
    p <- igraph::vcount(PPI.ig)  # p=3053 final nodes
    dist <- igraph::distances(
      graph = PPI.ig,
      v = igraph::V(PPI.ig)$name,
      to = igraph::V(PPI.ig)$name,
      mode = "out",
      weights = NULL)[1,]
    dist <- dist[!is.infinite(dist)]
    #for (f in 1:length(fdr.set)) {
    #  cat("fdr: ", fdr.set[f], "\n", sep="")
    #  cat("d/p: ", d[[f]]/p, "\n", sep="")
    #}
    cat("Median degree (\\nu): ", median(igraph::degree(PPI.ig)), "\n", sep="") #Median degree (\nu): 3
    cat("Median path length (\\lambda): ", median(dist), "\n", sep="") #Median path length (\lambda): 3.658
    cat("e: ", e, "\n", sep="") #e: 11741
    cat("p: ", p, "\n", sep="") #p: 3635
    cat("e/p: ", e/p, "\n", sep="") #e/p: 3.22998624484
    rm( inter_names, nodes_names)

    X <- E 
    cat("col: ", length(colnames(X)), "\n", sep="")
    cat("row: ", length(rownames(X)), "\n", sep="")
    colnames_df <- colnames(X)
    cur_trt <- sample_info_subset$Group[1]
    ct <- 1
    for (i in seq_along(colnames_df)) {
    # Find corresponding Sample ID in sample_info

        sample_id <- colnames_df[i]
        if (sample_info_subset$Group[match(sample_id, rownames(sample_info_subset))] != cur_trt) {
            cur_trt <- sample_info_subset$Group[match(sample_id, rownames(sample_info_subset))]
            ct <- 1
        }
        
        group_name <- sample_info_subset$Group[match(sample_id, rownames(sample_info_subset))]

        # Create new column name
        new_col_name <- paste("Sample", group_name, ct, sep = "_")

        # Rename the column
        colnames(X)[i] <- new_col_name
        ct <- ct + 1
    }

    sample.names <- colnames(X)
    myFactor<- gsub("[^a-zA-Z]", "", sub(".*\\.", "", colnames(X) ))

    TF<-as.factor(myFactor)
    TF


    a <- 1
    f <- 1
    qmin <- quant.cv.as$qmin[a,f]
    qmin.id <- quant.cv.as$qmin.id[a,f]
    q1se <- quant.cv.as$q1se[a,f]
    q1se.id <- quant.cv.as$q1se.id[a,f]
    seed <- 777
    CRNA.output <- CRNA::crna.as(graph=PPI.ig,
                        express=X,
                        factors=list(TF),
                        alpha=0.9,
                        quant=qmin[[f]],
                        FC=FC,
                        FR=FR,
                        R=R,
                        fdr=fdr.set[f],
                        conf=conf,
                        parallel=FALSE,
                        seed=seed)
    save.image(file=file.path(HOME.path, RESULT.path, "CRNA_REAL_RNASEQ.RData", fsep = .Platform$file.sep))
    # Visualization of inferred network graph - Exporting lists to Cytoscape
    cat("TEST3\n")
    write.table(CRNA.output$edges, 
                file.path(HOME.path, RESULT.path, paste("CRNA_real_reduced_edgelist_", group_labels[which(group_names == group_1)], "_vs_", group_labels[which(group_names == group_2)], "_fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep),
                quote=FALSE, 
                row.names=FALSE, 
                sep="\t")
    write.table(CRNA.output$states, 
                file.path(HOME.path, RESULT.path, paste("CRNA_real_reduced_nodelist_", group_labels[which(group_names == group_1)], "_vs_", group_labels[which(group_names == group_2)], "_fdr=", format(x=fdr.set[f], nsmall=2, scientific=TRUE), ".txt", sep=""), fsep=.Platform$file.sep), 
                quote=FALSE, 
                row.names=FALSE, 
                sep="\t")
}