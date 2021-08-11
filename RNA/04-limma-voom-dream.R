# RNA-seq pipeline using limma-voom-dream
# Ben Laufer

# Modifies and expands on these references:
#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
#https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
#https://www.bioconductor.org/packages/devel/bioc/vignettes/variancePartition/inst/doc/dream.html

# Load packages -----------------------------------------------------------

setwd("/Users/blaufer/Box Sync/PEBBLES/RNA")

#BiocManager::install("ben-laufer/DMRichR")
packages <- c("edgeR", "tidyverse", "RColorBrewer", "org.Mm.eg.db", "AnnotationDbi", "EnhancedVolcano",
              "enrichR", "openxlsx", "glue", "Glimma", "DMRichR", "magrittr", "variancePartition",
              "UpSetR", "ComplexUpset")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

# Needed or else "EnrichR website not responding"
enrichR:::.onAttach()

# Create log
sink("RNA-seq_log.txt", type = "output", append = FALSE, split = TRUE)

# Set up parallelization
param <- SnowParam(12, "SOCK", progressbar = TRUE)
register(param)

# Pipeline ----------------------------------------------------------------

# To test and develop, assign the variable tissue and then just run the main sections
tidyr::crossing(tissue = c("placenta", "brain"),
                sex = c("male", "female")) %>% 
  purrr::pwalk(function(tissue, sex){
    
    dir.create(glue::glue("{tissue}_{sex}"))
    dir.create(glue::glue("{tissue}_{sex}/plotData"))
    dir.create(glue::glue("{tissue}_{sex}/QC"))
  
  # Count Matrix ------------------------------------------------------------
  
  #name <- gsub( "(?:[^_]+_){4}([^_ ]+)*$","", files)
  
  # STAR quantMode geneCounts output:
  #column 1: gene ID
  #column 2: counts for unstranded RNA-seq
  #column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
  #column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)
  
  # KAPA mRNA HyperPrep Kit reads are reverse stranded, so select column 4
  # Confirm by looking at the N_noFeature line for the 3rd and 4th column and pick the column with the lowest count.
  
  sampleNames <- list.files(path = glue::glue(getwd(), "/GeneCounts"),
                            pattern = "*.ReadsPerGene.out.tab") %>%
    stringr::str_split_fixed("_", n = 3) %>%
    tibble::as_tibble() %>%
    tidyr::unite(Name, c(V1:V2), sep = "-") %>%
    dplyr::select(Name) %>% 
    purrr::flatten_chr()
  
  # Note: if using ensembl GTF file, use this section and swap in other ensembl comments
  # Could alternatively use edgeR::readDGE() but that calls to the slower read.delim()
  # ensemblIDs <- list.files(path = glue::glue(getwd(), "/GeneCounts"),
  #                          pattern = "*.ReadsPerGene.out.tab", full.names = TRUE)[1] %>% 
  #   data.table::fread(select = 1) %>%
  #   purrr::flatten_chr()
  # 
  # countMatrix <- list.files(path = glue::glue(getwd(), "/GeneCounts"),
  #                           pattern = "*.ReadsPerGene.out.tab", full.names = TRUE) %>%
  #   purrr::map_dfc(data.table::fread, select = 4, data.table = FALSE) %>%
  #   magrittr::set_colnames(sampleNames) %>% 
  #   magrittr::set_rownames(ensemblIDs)
  
  
  # Could alternatively use edgeR::readDGE() but that calls to the slower read.delim()
  geneSymbols <- list.files(path = glue::glue(getwd(), "/GeneCounts"),
                            pattern = "*.ReadsPerGene.out.tab", full.names = TRUE)[1] %>% 
    data.table::fread(select = 1) %>%
    purrr::flatten_chr()
  
  countMatrix <- list.files(path = glue::glue(getwd(), "/GeneCounts"),
                            pattern = "*.ReadsPerGene.out.tab", full.names = TRUE) %>%
    purrr::map_dfc(data.table::fread, select = 4, data.table = FALSE) %>%
    magrittr::set_colnames(sampleNames) %>% 
    magrittr::set_rownames(geneSymbols)
  
  # Remove meta info
  countMatrix <- countMatrix[-c(1:4),]
  
  # Design Matrix -----------------------------------------------------------
  
  designMatrix <- readxl::read_xlsx("sample_info.xlsx") %>%
    dplyr::rename(group = Treatment) %>% 
    dplyr::mutate_if(is.character, as.factor) %>%
    dplyr::mutate(Name = as.character(Name))
  
  # # Recode sex
  # designMatrix$Sex <- as.character(designMatrix$Sex)
  # designMatrix$Sex[designMatrix$Sex  == "F"] <- "0"
  # designMatrix$Sex[designMatrix$Sex  == "M"] <- "1"
  # designMatrix$Sex <- as.factor(designMatrix$Sex)
  
  samples.idx <- pmatch(designMatrix$Name, colnames(countMatrix))
  designMatrix <- designMatrix[order(samples.idx),]
  
  # Preprocessing -----------------------------------------------------------
  
  print(glue::glue("Preprocessing {sex} {tissue} samples"))
 
  # Select sample subset
  designMatrix <- designMatrix %>%
    dplyr::filter(Tissue == tissue & Sex == sex)
  
  countMatrix <- countMatrix %>%
    dplyr::select(contains(designMatrix$Name)) %>% 
    as.matrix()
  
  # Create DGE list and calculate normalization factors
  countMatrix <- countMatrix %>%
    DGEList() %>%
    calcNormFactors()
  
  # Reorder design matrix 
  samples.idx <- pmatch(designMatrix$Name, rownames(countMatrix$samples))
  designMatrix <- designMatrix[order(samples.idx),]
  stopifnot(rownames(countMatrix$samples) == designMatrix$Name)
  
  designMatrix %>%
    openxlsx::write.xlsx(glue::glue("{tissue}_{sex}/plotData/designMatrix.xlsx"))
  
  # Add sample info from design matrix to DGE list
  countMatrix$samples <- countMatrix$samples %>%
    tibble::add_column(designMatrix %>% dplyr::select(-Name))
  
  # ensembl 
  # # Add gene info
  # countMatrix$genes <- purrr::map_dfc(c("SYMBOL", "GENENAME", "ENTREZID", "CHR"), function(column){
  #   rownames(countMatrix$counts) %>% 
  #     AnnotationDbi::mapIds(org.Mm.eg.db,
  #                           keys = .,
  #                           column = column,
  #                           keytype = 'ENSEMBL') %>%
  #     as.data.frame() %>% 
  #     tibble::remove_rownames() %>% 
  #     purrr::set_names(column)
  # })
  
  # Raw density of log-CPM values
  
  L <- mean(countMatrix$samples$lib.size) * 1e-6
  M <- median(countMatrix$samples$lib.size) * 1e-6
  
  logCPM <- cpm(countMatrix, log = TRUE)
  logCPM.cutoff <- log2(10/M + 2/L)
  nsamples <- ncol(countMatrix)
  col <- brewer.pal(nsamples, "Paired")
  
  pdf(glue::glue("{tissue}_{sex}/QC/density_plot.pdf"), height = 8.5, width = 11)
  par(mfrow = c(1,2))
  
  plot(density(logCPM[,1]), col = col[1], lwd = 2, las = 2, main = "", xlab = "")
  title(main = "A. Raw data", xlab = "Log-cpm")
  abline(v = logCPM.cutoff, lty = 3)
  for (i in 2:nsamples){
    den <- density(logCPM[,i])
    lines(den$x, den$y, col = col[i], lwd = 2)
  }
  legend("topright", designMatrix$Name, text.col = col, bty = "n", cex = 0.5)
  
  # Filter genes with low expression
  
  rawCount <- dim(countMatrix)
  
  keep.exprs <- filterByExpr(countMatrix,
                             group = countMatrix$samples$group,
                             lib.size = countMatrix$samples$lib.size)
  
  countMatrix <- countMatrix[keep.exprs,, keep.lib.sizes = FALSE] %>%
    calcNormFactors() 
  
  filterCount <- dim(countMatrix)
  
  print(glue::glue("{100 - round((filterCount[1]/rawCount[1])*100)}% of genes were filtered from {rawCount[2]} samples, \\
             where there were {rawCount[1]} genes before filtering and {filterCount[1]} genes after filtering for {tissue}"))
  
  # Filtered density plot of log-CPM values 
  logCPM <- cpm(countMatrix, log = TRUE)
  plot(density(logCPM[,1]), col = col[1], lwd = 2, las =2 , main = "", xlab = "")
  title(main = "B. Filtered data", xlab = "Log-cpm")
  abline(v = logCPM.cutoff, lty = 3)
  for (i in 2:nsamples){
    den <- density(logCPM[,i])
    lines(den$x, den$y, col = col[i], lwd = 2)
  }
  legend("topright", designMatrix$Name, text.col = col, bty = "n", cex = 0.5)
  dev.off()
  
  # Interactive MDS plot
  Glimma::glMDSPlot(countMatrix,
                    groups = designMatrix,
                    path = getwd(),
                    folder = glue::glue("{tissue}_{sex}/QC/"),
                    html = "MDS-Plot",
                    launch = FALSE)

  # Surrogate variables analysis --------------------------------------------
  
  # # Create model matrices, with null model for svaseq, and don't force a zero intercept
  # mm <- model.matrix(~group + Litter,
  #                    data = designMatrix)
  # 
  # mm0 <- model.matrix(~1 + Litter,
  #                     data = designMatrix)
  # 
  # # svaseq requires normalized data that isn't log transformed
  # cpm <- cpm(countMatrix, log = FALSE)
  # 
  # # Calculate number of surrogate variables
  # nSv <- num.sv(cpm,
  #               mm,
  #               method = "leek")
  # 
  # # Estimate surrogate variables
  # svObj <- svaseq(cpm,
  #                 mm,
  #                 mm0,
  #                 n.sv = nSv)
  # 
  # # Update model to include surrogate variables
  # mm <- model.matrix(~Treatment + svObj$sv,
  #                    data = designMatrix)
  
  # Voom transformation and calculation of variance weights -----------------
  
  print(glue::glue("Normalizing {sex} {tissue} samples"))
  
  # Design
  model <- ~ group + (1|Litter)
  
  # Voom
  
  pdf(glue::glue("{tissue}_{sex}/QC/voom_mean-variance_trend.pdf"), height = 8.5, width = 11)
  voomLogCPM <- variancePartition::voomWithDreamWeights(countMatrix,
                                                        model,
                                                        designMatrix,
                                                        plot = TRUE)
  dev.off()
  
  # Boxplots of logCPM values before and after normalization
  pdf(glue::glue("{tissue}_{sex}/QC/normalization_boxplots.pdf"), height = 8.5, width = 11)
  par(mfrow=c(1,2))
  
  boxplot(logCPM, las = 2, col = col, main = "")
  title(main = "A. Unnormalised data", ylab = "Log-cpm")
  
  boxplot(voomLogCPM$E, las = 2, col = col, main = "")
  title(main = "B. Normalised data", ylab = "Log-cpm")
  
  dev.off()
  
  # Fitting linear models in limma ------------------------------------------
  
  print(glue::glue("Testing {sex} {tissue} samples for differential expression"))
  
  # Weight standard errors of log fold changes by within litter correlation 
  fit <- variancePartition::dream(voomLogCPM,
                                  model,
                                  designMatrix) #,
  #ddf = "Kenward-Roger") # small sample method
  
  head(coef(fit))
  
  # Save normalized expression values for correlation analyses
  voomLogCPM$E %>%
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "geneSymbol") %>%
    openxlsx::write.xlsx(glue::glue("{tissue}_{sex}/plotData/voomLogCPM.xlsx"))
  
  # Create DEG tibble -------------------------------------------------------
  
  print(glue::glue("Creating DEG list of {sex} {tissue} samples"))
  
  dfit <- fit %>%
    contrasts.fit(coef = "groupPCB")
  
  # Final model plot
  pdf(glue::glue("{tissue}_{sex}/QC/final_model_mean-variance_trend.pdf"),
      height = 8.5, width = 11)
  
  plotSA(dfit, main = "Final model: Mean-variance trend")
  
  dev.off()
  
  # Interactive MA plot
  # Glimma::glimmaMA(dfit,
  #                  dge = countMatrix,
  #                  path = getwd(),
  #                  folder = glue::glue("{tissue}_{sex}/QC/"),
  #                  html = "MA-Plot",
  #                  launch = FALSE)
  
  # Top differentially expressed genes
  DEGs <- dfit %>%
    topTable(sort.by = "P", n = Inf) %>% 
    rownames_to_column() %>% 
    tibble::as_tibble() %>%
    dplyr::rename(SYMBOL = rowname) %>% #ensgene
    dplyr::mutate(FC = dplyr::case_when(logFC > 0 ~ 2^logFC,
                                        logFC < 0 ~ -1/(2^logFC))) %>%
    dplyr::select(SYMBOL, FC, logFC, P.Value, adj.P.Val, AveExpr, t, z.std) %T>%  #ensgene
    openxlsx::write.xlsx(file = glue::glue("{tissue}_{sex}/plotData/DEGs.xlsx")) %>% 
    dplyr::filter(P.Value < 0.05) %T>%
    openxlsx::write.xlsx(file = glue::glue("{tissue}_{sex}/sig_DEGs.xlsx"))
  
  # Ontologies and Pathways -------------------------------------------------
  
  print(glue::glue("Performing GO and pathway analysis of {sex} {tissue} samples"))

  tryCatch({
  DEGs %>% 
    dplyr::select(SYMBOL) %>%
    purrr::flatten() %>%
    enrichR::enrichr(c("GO_Biological_Process_2018",
                       "GO_Cellular_Component_2018",
                       "GO_Molecular_Function_2018",
                       "KEGG_2019_Mouse",
                       "Panther_2016",
                       "Reactome_2016",
                       "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")) %>% 
    purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis = "")) %T>% # %>% 
    #purrr::map(~ dplyr::filter(., Adjusted.P.value < 0.05)) %>% 
    #purrr::map(~ dplyr::filter(., stringr::str_detect(Genes, ";"))) %>% 
    openxlsx::write.xlsx(file = glue::glue("{tissue}_{sex}/enrichr.xlsx")) %>%
    DMRichR::slimGO(tool = "enrichR",
                    annoDb = "org.Mm.eg.db",
                    plots = FALSE) %T>%
    openxlsx::write.xlsx(file = glue::glue("{tissue}_{sex}/rrvgo_enrichr.xlsx")) %>%
    DMRichR::GOplot() %>%
    ggplot2::ggsave(glue::glue("{tissue}_{sex}/enrichr_plot.pdf"),
                    plot = .,
                    device = NULL,
                    height = 8.5,
                    width = 10) },
  error = function(error_condition) {
    print(glue::glue("ERROR: Gene Ontology pipe didn't finish for {sex} {tissue}"))
  })
  print(glue::glue("The pipeline has finished for {sex} {tissue} samples"))
})

sink()

# Main Plot ---------------------------------------------------------------

## Volcano Plot -----------------------------------------------------------

volcanoPlot <- function(tissue, sex){
  readxl::read_xlsx(glue::glue("{tissue}_{sex}/plotData/DEGs.xlsx")) %>% 
    EnhancedVolcano::EnhancedVolcano(title = "",
                                     subtitle = "",
                                     labSize = 5,
                                     lab = .$SYMBOL, # 
                                     x = 'logFC',
                                     y = 'P.Value', # P.Value 'adj.P.Val'
                                     col = c("grey30", "royalblue", "royalblue", "red2"),
                                     pCutoff = 0.05,
                                     FCcutoff = 0.0,
                                     legendPosition = "none",
                                     raster = FALSE) +
    ggplot2::coord_cartesian(xlim = c(-2, 2),
                             ylim = c(0, 8))
}

## Heatmap ----------------------------------------------------------------

DEGheatmap <- function(tissue, sex){
  
  DEGs <- readxl::read_xlsx(glue::glue("{tissue}_{sex}/sig_DEGs.xlsx")) 
  
  voomLogCPM <- readxl::read_xlsx(glue::glue("{tissue}_{sex}/plotData/voomLogCPM.xlsx")) %>%
    tibble::column_to_rownames("geneSymbol") %>%
    as.matrix()
  
  designMatrix <- readxl::read_xlsx(glue::glue("{tissue}_{sex}/plotData/designMatrix.xlsx"))
  
  voomLogCPM[which(rownames(voomLogCPM) %in% DEGs$SYMBOL),] %>% #DEGs$ensgene
    as.matrix() %>% 
    ComplexHeatmap::pheatmap(.,
                             scale = "row",
                             annotation_col = designMatrix %>%
                               tibble::column_to_rownames(var = "Name") %>% 
                               dplyr::select(Treatment = group),
                             color = RColorBrewer::brewer.pal(11, name = "RdBu") %>%
                               rev(),
                             show_colnames = FALSE,
                             show_rownames = FALSE,
                             main = glue::glue("{nrow(DEGs)} Differentially Expressed Genes"),
                             border_color = NA,
                             use_raster = TRUE,
                             annotation_colors = list(Treatment = c("PCB" = "#F8766D",
                                                                    "Control" = "#619CFF"))) %>%
    grid::grid.grabExpr()
}

## GO ---------------------------------------------------------------------

GOplot <- function(tissue, sex){
  readxl::read_xlsx(glue::glue("{tissue}_{sex}/rrvgo_enrichr.xlsx")) %>% 
    dplyr::mutate("Database" = dplyr::recode_factor(`Gene Ontology`,
                                                    "Biological Process" = "GO Biological Process",
                                                    "Cellular Component" = "GO Cellular Component",
                                                    "Molecular Function" = "GO Molecular Function")) %>%
    dplyr::select(-`Gene Ontology`) %>% 
    dplyr::mutate(Term = stringr::str_trim(Term)) %>%
    dplyr::mutate(Term = Hmisc::capitalize(Term)) %>%
    dplyr::mutate(Term = stringr::str_trunc(Term, 50)) %>%
    dplyr::group_by(Database) %>%
    dplyr::slice(1:7) %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(Term = factor(.$Term, levels = unique(.$Term[order(forcats::fct_rev(.$Database), .$`-log10.p-value`)]))) %>%
    ggplot2::ggplot(ggplot2::aes(x = Term,
                                 y = `-log10.p-value`,
                                 fill = Database,
                                 group = Database)) +
    ggplot2::geom_bar(stat = "identity",
                      position = position_dodge(),
                      color = "Black",
                      size = 0.25) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = c(0, 0, 0.1, 0),
                                breaks = scales::pretty_breaks(4)) +
    ggsci::scale_fill_d3() +
    ggplot2::labs(y = expression("-log"[10](italic(p)))) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = element_text(size = 20),
                   axis.title.y = element_blank(),
                   axis.title.x = element_text(size = 12),
                   legend.position = "none")
}

## Run --------------------------------------------------------------------

cowplot::plot_grid(volcanoPlot("placenta", "female"),
                   volcanoPlot("placenta", "male"),
                   volcanoPlot("brain", "female"),
                   volcanoPlot("brain", "male"),
                   DEGheatmap("placenta", "female"),
                   DEGheatmap("placenta", "male"),
                   DEGheatmap("brain", "female"),
                   DEGheatmap("brain", "male"),
                   GOplot("placenta", "female"),
                   GOplot("placenta", "male"),
                   GOplot("brain", "female"),
                   GOplot("brain", "male"),
                   ncol = 4,
                   #align = "v",
                   labels = c("A)", "B)", "C)", "D)",
                              "E)", "F)", "G)", "H)",
                              "I)", "J)", "K)", "L)"),
                   label_size = 24,
                   rel_heights = c(2, 1.4, 1.5),
                   scale = 0.95) %>%
  ggplot2::ggsave("RNA-seq figure.pdf",
                  plot = ., 
                  width = 30,
                  height = 20)

# Tables ------------------------------------------------------------------

## DEGs -------------------------------------------------------------------

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
    readxl::read_xlsx(glue::glue("{tissue}_{sex}/sig_DEGs.xlsx")) %>%
      dplyr::mutate(Tissue = Hmisc::capitalize(tissue),
                    Sex = Hmisc::capitalize(sex))
  }) %>% 
  dplyr::left_join(purrrIds(.) %>%
                     dplyr::distinct()) %>% 
  dplyr::select(Tissue, Sex,
                Gene = SYMBOL,
                FC, logFC, AveExpr, t, z.std,
                P.Value, adj.P.Val,
                Description = GENENAME,
                ENSEMBL, ENTREZID, CHR) %>% 
  openxlsx::write.xlsx("DEGs_supplementary_table.xlsx")

## GO ---------------------------------------------------------------------

tidyr::crossing(tissue = c("placenta", "brain"),
                sex = c("male", "female")) %>% 
  purrr::pmap_dfr(function(tissue, sex){
    readxl::read_xlsx(glue::glue("{tissue}_{sex}/rrvgo_enrichr.xlsx")) %>% 
      dplyr::mutate(Tissue = Hmisc::capitalize(tissue),
                    Sex = Hmisc::capitalize(sex))
  }) %>% 
  dplyr::select(Tissue, Sex, "Gene Ontology", Term, "-log10.p-value") %>% 
  openxlsx::write.xlsx("GO_supplementary_table.xlsx")

# UpSet of overlaps -------------------------------------------------------

print(glue::glue("Creating UpSet plot of all gene symbol overlaps"))

geneUpsetPlot <- function(list = list){
  list %>%
    UpSetR::fromList() %>%
    ComplexUpset::upset(.,
                        names(.),
                        n_intersections = 40, # Default from UpSetR
                        width_ratio = 0.2,
                        height_ratio = 0.4,
                        base_annotations = list(
                          'Intersection size'= intersection_size(
                            counts = TRUE
                          ) +
                            theme(
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.line = element_line(colour = 'black'),
                              axis.text = element_text(size = 10),
                              axis.title = element_text(size = 10)
                            ) +
                            scale_y_continuous(expand = c(0, 0, 0.1, 0))
                        ),
                        sort_sets = FALSE,
                        queries = list(
                          # upset_query(
                          #   intersect = c("Placenta Female", "Placenta Male",
                          #                 "Brain Female", "Brain Male"),
                          #   color = '#E41A1C',
                          #   fill = '#E41A1C',
                          #   only_components = c('intersections_matrix', 'Intersection size')
                          # ),
                          upset_query(
                            intersect = c("Placenta Female", "Placenta Male"),
                            color = '#FF7F00',
                            fill = '#FF7F00',
                            only_components = c('intersections_matrix', 'Intersection size')
                          ),
                          # upset_query(
                          #   intersect = c("Brain Female", "Brain Male"),
                          #   color = '#4DAF4A',
                          #   fill = '#4DAF4A',
                          #   only_components = c('intersections_matrix', 'Intersection size')
                          # ),
                          upset_query(set = 'Placenta Female', fill = '#FF7F00'),
                          upset_query(set = 'Placenta Male', fill = '#FF7F00'),
                          upset_query(set = 'Brain Female', fill = '#4DAF4A'),
                          upset_query(set = 'Brain Male', fill = '#4DAF4A')
                        ),
                        matrix = intersection_matrix(
                          geom = geom_point(
                            shape = 'circle filled',
                            size = 3.5,
                            stroke = 0.45
                          )
                        ),
                        set_sizes = (
                          upset_set_size(geom = geom_bar(color = 'black')
                          )
                        ),
                        themes = upset_modify_themes( # names(upset_themes) 
                          list(
                            'intersections_matrix' = theme(
                              axis.text = element_text(size = 10),
                              axis.title = element_blank()
                            ),
                            'overall_sizes' = theme(
                              axis.title = element_text(size = 10),
                              axis.text = element_text(angle = 45,
                                                       vjust = 0.5,
                                                       hjust = 1,
                                                       size = 6)
                            )
                          )
                        )
    )
}

pdf(glue::glue("UpSet_RNA-seq.pdf"),
    height = 4, width = 6, onefile = FALSE)

print({
  c("placenta_male", "placenta_female", "brain_male", "brain_female") %>% 
    purrr::set_names() %>% 
    purrr::map(function(contrast){readxl::read_xlsx(glue::glue("{contrast}/sig_DEGs.xlsx")) %>%
        purrr::pluck("SYMBOL")}) %>% #ensgene
    purrr::set_names(c("Placenta Male", "Placenta Female", "Brain Male", "Brain Female")) %>% 
    rev() %>% 
    geneUpsetPlot()
})

dev.off()
