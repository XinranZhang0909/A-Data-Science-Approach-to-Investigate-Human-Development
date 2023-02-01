doLimma <- function(exprsMat, cellTypes, exprs_pct = 0.05){
  cellTypes <- droplevels(as.factor(cellTypes))
  tt <- list()
  for (i in 1:nlevels(cellTypes)) {
    tmp_celltype <- (ifelse(cellTypes == levels(cellTypes)[i], 1, 0))
    design <- stats::model.matrix(~tmp_celltype)
    meanExprs <- do.call(cbind, lapply(c(0,1), function(i){
      Matrix::rowMeans(exprsMat[, tmp_celltype == i, drop = FALSE])
    }))
    meanPct <- do.call(cbind, lapply(c(0,1), function(i){
      Matrix::rowSums(exprsMat[, tmp_celltype == i, drop = FALSE] > 0)/sum(tmp_celltype == i)
    }))
    #keep <- meanPct[,2] > exprs_pct
    y <- methods::new("EList")
    #y$E <- exprsMat[keep, ]
    y$E <- exprsMat
    fit <- limma::lmFit(y, design = design)
    fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
    tt[[i]] <- limma::topTable(fit, n = Inf, adjust.method = "BH", coef = 2)
    if (!is.null(tt[[i]]$ID)) {
      tt[[i]] <- tt[[i]][!duplicated(tt[[i]]$ID),]
      rownames(tt[[i]]) <- tt[[i]]$ID
    }
    tt[[i]]$meanExprs.1 <- meanExprs[rownames(tt[[i]]), 1]
    tt[[i]]$meanExprs.2 <- meanExprs[rownames(tt[[i]]), 2]
    tt[[i]]$meanPct.1 <- meanPct[rownames(tt[[i]]), 1]
    tt[[i]]$meanPct.2 <- meanPct[rownames(tt[[i]]), 2]
  }
  names(tt) <- levels(cellTypes)
  return(tt)
}

###### adapted from https://github.com/csoneson/conquer_comparison/tree/master/scripts
doEdgeR <- function(countMatrix, cellTypes) {
  
  message("edgeRQLFDetRate")
  cty <- droplevels(as.factor(cellTypes))
  df <- list()
  
  for (i in 1:nlevels(cty)) {
    grp <- (ifelse(cty == levels(cty)[i], 1, 0))
    names(grp) <- names(cty)
    
    dge <- DGEList(countMatrix, group = cty)
    dge <- calcNormFactors(dge)
    cdr <- scale(colMeans(countMatrix > 0))
    design <- model.matrix(~ cdr + grp)
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design)
    qlf <- glmQLFTest(fit)
    tt <- topTags(qlf, n = Inf)
    
    #plotBCV(dge)
    #plotQLDisp(fit)
    #hist(tt$table$PValue, 50)
    #hist(tt$table$FDR, 50)
    #limma::plotMDS(dge, col = as.numeric(as.factor(L$condt)), pch = 19)
    #plotSmear(qlf)
    
    df[[i]] = tt$table
    
  }
  names(df) <- levels(cty)
  return(df)
}

###### Adapted from https://github.com/csoneson/conquer_comparison/tree/master/scripts
doMAST <- function(exprsMat, cellTypes) {
  # input must be normalised, log-transformed data
  message("MAST")
  cty <- droplevels(as.factor(cellTypes))
  names(cty) <- colnames(exprsMat)
  df <- list()
  for (i in 1:nlevels(cty)) {
    grp <- (ifelse(cty == levels(cty)[i], 1, 0))
    names(grp) <- names(cty)
    cdr <- scale(colMeans(exprsMat > 0))
    sca <- FromMatrix(exprsArray = exprsMat, 
                      cData = data.frame(wellKey = names(grp), 
                                         grp = grp, cdr = cdr))
    zlmdata <- zlm(~cdr + grp, sca)
    mast <- lrTest(zlmdata, "grp")
    
    df[[i]] = data.frame(pval = mast[, "hurdle", "Pr(>Chisq)"],
                         lambda = mast[, "cont", "lambda"],
                         row.names = names(mast[, "hurdle", "Pr(>Chisq)"]))
    df[[i]]$fdr <- stats::p.adjust(df[[i]]$pval, method="BH")
    b <- getLogFC(zlmdata)
    df[[i]] <- cbind(df[[i]], b[b$contrast=="grp",c("logFC", "varLogFC","z"), with=F])
    df[[i]] <- df[[i]][order(df[[i]]$pval),]
  }
  names(df) <- levels(cty)
  return(df)
  
}

###### Adapted from https://github.com/csoneson/conquer_comparison/tree/master/scripts
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(genefilter))

doTest <- function(exprsMat, cellTypes) {
  # input must be normalised, log-transformed data
  message("t-test")
  cty <- droplevels(as.factor(cellTypes))
  names(cty) <- colnames(exprsMat)
  
  tt <- list()
  for (i in 1:nlevels(cty)) {
    tmp_celltype <- (ifelse(cty == levels(cty)[i], 1, 0))
    
    tt[[i]] <- t(apply(exprsMat, 1, function(x) {
      x1 <- x[tmp_celltype == 0]
      x2 <- x[tmp_celltype == 1]
      
      res <- stats::t.test(x2, y=x1)
      return(c(stats=res$statistic,
               pvalue=res$p.value))
    }))
    tt[[i]] <- as.data.frame(tt[[i]])
    tt[[i]]$adj.pvalue <- stats::p.adjust(tt[[i]]$pvalue, method = "BH")
  }
  names(tt) <- levels(cty)
  return(tt)
}

###### Adapted from https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_Wilcoxon.R
suppressPackageStartupMessages(library(edgeR))

doWilcoxon <- function(exprsMat, cellTypes) {
  # input must be normalised, log-transformed data
  message("Wilcoxon")
  cty <- droplevels(as.factor(cellTypes))
  names(cty) <- colnames(exprsMat)
  
  tt <- list()
  for (i in 1:nlevels(cty)) {
    tmp_celltype <- (ifelse(cty == levels(cty)[i], 1, 0))
    
    tt[[i]] <- t(apply(exprsMat, 1, function(x) {
      
      res <- stats::wilcox.test(x ~ tmp_celltype)
      c(stats=res$statistic,
        pvalue=res$p.value)
    }))
    tt[[i]] <- as.data.frame(tt[[i]])
    tt[[i]]$adj.pvalue <- stats::p.adjust(tt[[i]]$pvalue, method = "BH")
  }
  names(tt) <- levels(cty)
  return(tt)
}


###### Adapted from https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_voomlimma.R
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))

doVoom <- function(exprsMat, cellTypes) {
  # input must be count data
  message("Voom Limma")
  cty <- droplevels(as.factor(cellTypes))
  names(cty) <- colnames(exprsMat)
  
  tt <- list()
  for (i in 1:nlevels(cty)) {
    tmp_celltype <- (ifelse(cty == levels(cty)[i], 1, 0))
    
    design <- stats::model.matrix(~tmp_celltype)
    
    y <- methods::new("EList")
    #y$E <- exprsMat[keep, ]
    y$E <- exprsMat
    vm <- limma::voom(y, design = design)
    fit <- limma::lmFit(vm, design = design)
    fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
    tt[[i]] <- limma::topTable(fit, n = Inf, adjust.method = "BH", coef = 2)
    if (!is.null(tt[[i]]$ID)) {
      tt[[i]] <- tt[[i]][!duplicated(tt[[i]]$ID),]
      rownames(tt[[i]]) <- tt[[i]]$ID
    }
    
  }
  names(tt) <- levels(cty)
  return(tt)
}
