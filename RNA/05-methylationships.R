# Methylation and RNA relationships
# Plot all overalps and test for correlations between placent-brain overlaping DMRs
# Ben Laufer

.libPaths("/share/lasallelab/programs/DMRichR/R_4.1")

if (!requireNamespace("Vennerable", quietly = TRUE))
  BiocManager::install("js229/Vennerable") 

packages <- c("DMRichR", "Vennerable",
              "TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db",
              "UpSetR", "magrittr", "tidyverse")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

setwd("/share/lasallelab/Ben/PEBBLES/")
dir.create("methylationships")

# Overlaps ----------------------------------------------------------------

dir.create("methylationships/overlaps")

loadDMRs <- function(name){
  load(glue::glue("DNA/DMRs/{name}/RData/DMRs.RData"))
  assign(name, sigRegions, envir = .GlobalEnv)
  rm(sigRegions, regions)
}

contrasts <- c("placenta_male", "placenta_female", "brain_male", "brain_female")

purrr::walk(contrasts, loadDMRs)

DMRichR::annotationDatabases("mm10")

#' getSymbol
#' @description Get gene symbols as character vector from excel column 
#' @param file An excel spreadsheet from \code{DMRichR} (DMRs_annotated.xlsx)
#' @return Vector of gene symbols
#' @importFrom dplyr select filter
#' @importFrom purrr pluck
#' @importFrom magrittr %>% 
#' @export getSymbol
getSymbol <- function(regions = sigRegions,
                      TxDb = TxDb,
                      annoDb = annoDb){
  regions %>% 
    DMRichR::annotateRegions(TxDb = TxDb,
                             annoDb = annoDb) %>%
    #dplyr::filter(annotation != "Distal Intergenic") %>% 
    purrr::pluck("geneSymbol")
}

contrasts <- c("placenta_male",
               "placenta_female",
               "brain_male",
               "brain_female")

print(glue::glue("Annotating all DMRs"))

DNA <- contrasts %>%
  purrr::set_names() %>% 
  purrr::map(function(contrast){
    getSymbol(get(contrast),
              TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
              annoDb = "org.Mm.eg.db")}) %>%
  purrr::set_names(c("placenta_male_WGBS",
                     "placenta_female_WGBS",
                     "brain_male_WGBS",
                     "brain_female_WGBS"))

print(glue::glue("Loading DEGs"))

RNA <- contrasts %>% 
  purrr::set_names() %>% 
  purrr::map(function(contrast){readxl::read_xlsx(glue::glue("RNA/{contrast}/sig_DEGs.xlsx")) %>%
      dplyr::filter(!is.na(SYMBOL)) %>% 
      purrr::pluck("SYMBOL")}) %>% 
  purrr::set_names(c("placenta_male_RNA-seq",
                     "placenta_female_RNA-seq",
                     "brain_male_RNA-seq",
                     "brain_female_RNA-seq"))

purrr::walk(contrasts, function(contrast){
  
  print(glue::glue("Generating gene symbol overlaps for {contrast} samples"))
  
  v <- c(DNA,RNA)[grep(contrast, names(c(DNA,RNA)))] %>%
    purrr::set_names(c("WGBS", "RNA-seq")) %>%
    Vennerable::Venn()
  
  print(glue::glue("Saving gene symbol Venn of overlaps for {contrast} samples"))
  
  v@IntersectionSets[["11"]] %>% 
    tibble::as_tibble() %>% 
    purrr::set_names("SYMBOL") %>% 
    openxlsx::write.xlsx(glue::glue("methylationships/overlaps/{contrast}_symbol_overlaps.xlsx"))
  
  svg(glue::glue("methylationships/overlaps/{contrast}_symbols_Venn.svg"),
      height = 8.5,
      width = 12)
  
  plot(v)
  
  dev.off()
})

purrrIds <- function(DEGs){
  purrr::map_dfc(c("SYMBOL", "ENSEMBL", "GENENAME", "ENTREZID", "CHR"),
                 function(column){
                   DEGs$SYMBOL %>%
                     AnnotationDbi::mapIds(org.Mm.eg.db,
                                           keys = .,
                                           column = column,
                                           keytype = 'SYMBOL') %>%
                     as.data.frame() %>%
                     tibble::remove_rownames() %>%
                     purrr::set_names(column)
                 })
}

tidyr::crossing(tissue = c("placenta", "brain"),
                sex = c("male", "female")) %>% 
  purrr::pmap_dfr(function(tissue, sex){
    readxl::read_xlsx(glue::glue("methylationships/overlaps/{tissue}_{sex}_symbol_overlaps.xlsx")) %>%
      dplyr::mutate(Tissue = Hmisc::capitalize(tissue),
                    Sex = Hmisc::capitalize(sex))
  }) %>% 
  dplyr::left_join(purrrIds(.) %>%
                     dplyr::distinct()) %>% 
  dplyr::select(Tissue, Sex,
                Gene = SYMBOL,
                Description = GENENAME,
                ENSEMBL, ENTREZID, CHR) %>% 
  openxlsx::write.xlsx("methylationships/overlaps/Overlaps_supplementary_table.xlsx")

# Correlations ------------------------------------------------------------

dir.create("methylationships/correlations")

rm(list=ls())
DMRichR::annotationDatabases("mm10")

tidyr::crossing(tissue = c("placenta", "brain"),
                sex = c("male", "female")) %>% 
  purrr::pwalk(function(tissue, sex){
    
    # Get tissue overlaps
    dmrs <- readxl::read_xlsx(glue::glue("methylationships/{sex}_tissue_overlaps.xlsx")) %>%
      plyranges::as_granges()

    # Get matching genes that are expressed
    rna <- readxl::read_xlsx(glue::glue("RNA/{tissue}_{sex}/plotData/voomLogCPM.xlsx")) %>%
      dplyr::filter(geneSymbol %in% dmrs$geneSymbol) %>%
      magrittr::set_names(names(.) %>% stringr::word(1, sep = "-")) %>%
      tibble::column_to_rownames("geneSymbol") %>% 
      as.matrix()
    
    # Get smoothed values
    load(glue::glue("DNA/DMRs/{tissue}_{sex}/RData/bsseq.RData"))
    
    rownames(pData(bs.filtered.bsseq)) <- bs.filtered.bsseq %>%
      pData() %>%
      rownames() %>%
      stringr::word(1, sep = "-")
    
    smoothed <- bs.filtered.bsseq %>%
      bsseq::getMeth(BSseq = .,
                     regions = dmrs %>% plyranges::filter(geneSymbol %in% rownames(rna)),
                     type = "smooth",
                     what = "perRegion") %>% 
      as.data.frame(check.names = FALSE) %>%
      dplyr::bind_cols(dmrs %>%
                         plyranges::filter(geneSymbol %in% rownames(rna)) %>%
                         tibble::as_tibble() %>%
                         dplyr::select(geneSymbol),
                       .) %>%
      tibble::column_to_rownames("geneSymbol") %>% 
      as.matrix()
    
    # Match order
    genes.idx <- match(rownames(smoothed), rownames(rna))
    samples.idx <- match(colnames(smoothed), colnames(rna))
    
    rna <- rna[genes.idx, samples.idx]
    
    stopifnot(rownames(smoothed) == rownames(rna))
    stopifnot(colnames(smoothed) == colnames(rna))
    
    # Test correlation
    moduleTraitCor <- WGCNA::cor(t(smoothed), t(rna))
    moduleTraitPvalue <- WGCNA::corPvalueStudent(moduleTraitCor, ncol(smoothed)) 
    
    # Plot
    pdf(glue::glue("methylationships/correlations/{tissue}_{sex}_cor_tissue_overlaps.pdf"), height = 10, width = 12)
    
    textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                        signif(moduleTraitPvalue, 1), ")", sep = "")
    dim(textMatrix) <- dim(moduleTraitCor)
    par(mar = c(6, 8.5, 3, 3))
    
    WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                          xLabels = colnames(moduleTraitCor),
                          yLabels = rownames(moduleTraitCor),
                          ySymbols = rownames(moduleTraitCor),
                          colorLabels = FALSE,
                          colors = WGCNA::blueWhiteRed(50),
                          textMatrix = textMatrix,
                          setStdMargins = FALSE,
                          cex.text = 1,
                          zlim = c(-1,1),
                          main = paste("DMR-RNA Correlations"))
    
    dev.off()
  })
