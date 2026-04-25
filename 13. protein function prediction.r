suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(DOSE)
  library(igraph)
})
df <- read.csv("results/Common_Significant_Proteins.csv")

symbols <- sub("\\..*$", "", as.character(df$Protein))

sig_symbols <- unique(symbols)

test <- suppressMessages(
  bitr(
    sig_symbols,
    fromType = "SYMBOL",
    toType   = c("ENTREZID","SYMBOL"),
    OrgDb    = org.Hs.eg.db
  )
)

unmapped <- setdiff(sig_symbols, unique(test$SYMBOL))
if (length(unmapped) > 0) {
  alias_map <- tryCatch(
      suppressMessages(
        bitr(
          unmapped,
          fromType = "ALIAS",
          toType   = c("ENTREZID","SYMBOL"),
          OrgDb    = org.Hs.eg.db
        )
      ),
      error = function(e) return(NULL)
    )
  if (!is.null(alias_map) && nrow(alias_map) > 0) {
    alias_map <- alias_map[, c("SYMBOL","ENTREZID")]
    test <- unique(rbind(test[, c("SYMBOL","ENTREZID")], alias_map))
  }
}

if ("HGNChelper" %in% rownames(installed.packages())) {
  suppressWarnings(suppressPackageStartupMessages(library(HGNChelper)))
  unmapped2 <- setdiff(sig_symbols, unique(test$SYMBOL))
  if (length(unmapped2) > 0) {
    fix <- checkGeneSymbols(unmapped2, unmapped.as.na = TRUE)
    corrected <- unique(fix$Suggested.Symbol[!is.na(fix$Suggested.Symbol)])
    corrected <- setdiff(corrected, unique(test$SYMBOL))
    if (length(corrected) > 0) {
      corr_map <- tryCatch(
        suppressMessages(
          bitr(
            corrected,
            fromType = "SYMBOL",
            toType   = c("ENTREZID","SYMBOL"),
            OrgDb    = org.Hs.eg.db
          )
        ),
        error = function(e) return(NULL)
      )
      if (!is.null(corr_map) && nrow(corr_map) > 0) {
        test <- unique(rbind(test[, c("SYMBOL","ENTREZID")], corr_map[, c("SYMBOL","ENTREZID")]))
      }
    }
  }
}

test_entrez <- unique(test$ENTREZID)
ego <- suppressMessages(
  enrichGO(
    gene          = test_entrez,
    OrgDb         = org.Hs.eg.db,
    ont           = "ALL",
    pAdjustMethod = "BH",
    minGSSize     = 1,
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
)

print(
  dotplot(ego, showCategory = 20)
)
kk <- suppressMessages(
  enrichKEGG(
    gene         = test_entrez,
    organism     = "hsa",
    pvalueCutoff = 0.05
  )
)

print(
  dotplot(kk, showCategory = 20)
)

ego2 <- pairwise_termsim(ego)
print(cnetplot(ego, showCategory = 20))

dores <- suppressMessages(
  enrichDO(
    gene = as.character(test_entrez),
    pvalueCutoff = 0.05
  )
)

print(dotplot(dores, showCategory = 20))




