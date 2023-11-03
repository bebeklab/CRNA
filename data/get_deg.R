library(edgeR)
library(limma)

## GSE11352
gse11352 <- readRDS('~/Research/crna/GSE11352.RDS')
ppi<-readRDS('~/Research/crna/PPI_graph.RDS')

expr_mat <- exprs(gse11352)
#expr_mat <- expr_mat[rownames(expr_mat) %in% union(ppi$source_name, ppi$target_name),]
pheno <- as.data.frame(phenoData(gse11352)@data)
idx <- pheno$time_point %in% c('24hr timepoint', '48hr timepoint')
pheno <- pheno[idx,]
expr_mat <- expr_mat[,match(pheno$geo_accession, table=colnames(expr_mat))]
tx.status <- as.factor(pheno$treatment)
design <- model.matrix(~0+tx.status)
fit <- lmFit(expr_mat, design = design)
cont.matrix <- makeContrasts(Treatment=tx.statusTx-tx.statusnone, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
res <- topTable(fit2, adjust="BH", sort.by = 'logFC', number = Inf)
res$id <- rownames(res)

require(data.table)
res_ft <- sapply(1:nrow(res), simplify = FALSE, function(i){
  if(grepl(pattern='///', res$id[i])){
    uniq_names <- sapply(trimws(unlist(strsplit(res$id[i], split = '///', fixed=TRUE))), simplify = FALSE, function(k){
      return(data.table('logFC' = res$logFC[i],
                        'AveExpr' = res$AveExpr[i],
                        't' = res$t[i],
                        'P.Value' = res$P.Value[i],
                        'adj.P.Val' = res$adj.P.Val[i],
                        'B' = res$B[i],
                        'id' = k))
    })
    return(rbindlist(uniq_names))
  }else{
    return(data.table('logFC' = res$logFC[i],
                      'AveExpr' = res$AveExpr[i],
                      't' = res$t[i],
                      'P.Value' = res$P.Value[i],
                      'adj.P.Val' = res$adj.P.Val[i],
                      'B' = res$B[i],
                      'id' = res$id[i]))
  }
})
res_ft <- rbindlist(res_ft)
res_ft <- res_ft[res_ft$id %in% union(ppi$source_name, ppi$target_name),]
write.csv(res_ft, file='~/Research/crna/DEG_GSE11352.csv', col.names = TRUE, row.names = FALSE,
          append = FALSE, quote=FALSE)

## GSE60880
gse.ft <- readRDS('~/Research/crna/GSE60880.RDS')
ppi<-readRDS('~/Research/crna/PPI_graph.RDS')

expr_mat <- exprs(gse.ft)
pheno <- as.data.frame(phenoData(gse.ft)@data)
idx <- pheno$time %in% c('1', '2')
pheno <- pheno[idx,]
expr_mat <- expr_mat[,match(pheno$geo_accession, table=colnames(expr_mat))]
tx.status <- as.factor(pheno$treatment)
design <- model.matrix(~0+tx.status)
fit <- lmFit(expr_mat, design = design)
cont.matrix <- makeContrasts(Treatment=tx.statusEGF-tx.statusnone, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
res <- topTable(fit2, adjust="BH", sort.by = 'logFC', number = Inf)
res$id <- rownames(res)

require(data.table)
res_ft <- sapply(1:nrow(res), simplify = FALSE, function(i){
  if(grepl(pattern='///', res$id[i])){
    uniq_names <- sapply(trimws(unlist(strsplit(res$id[i], split = '///', fixed=TRUE))), simplify = FALSE, function(k){
      return(data.table('logFC' = res$logFC[i],
                        'AveExpr' = res$AveExpr[i],
                        't' = res$t[i],
                        'P.Value' = res$P.Value[i],
                        'adj.P.Val' = res$adj.P.Val[i],
                        'B' = res$B[i],
                        'id' = k))
    })
    return(rbindlist(uniq_names))
  }else{
    return(data.table('logFC' = res$logFC[i],
                      'AveExpr' = res$AveExpr[i],
                      't' = res$t[i],
                      'P.Value' = res$P.Value[i],
                      'adj.P.Val' = res$adj.P.Val[i],
                      'B' = res$B[i],
                      'id' = res$id[i]))
  }
})
res_ft <- rbindlist(res_ft)
res_ft <- res_ft[res_ft$id %in% union(ppi$source_name, ppi$target_name),]
write.csv(res_ft, file='~/Research/crna/DEG_GSE60880.csv', col.names = TRUE, row.names = FALSE,
          append = FALSE, quote=FALSE)

