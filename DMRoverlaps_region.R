# DMR Overlaps
# Plots and enrichment analyses for overlaps between datasets from the genomic coordinate and gene symbol perspectives
# Ben Laufer

# Packages ----------------------------------------------------------------

.libPaths("/share/lasallelab/programs/DMRichR/R_4.0")

if (!requireNamespace("Vennerable", quietly = TRUE))
  BiocManager::install("js229/Vennerable") 

packages <- c("DMRichR", "ChIPpeakAnno", "Vennerable", "TxDb.Mmusculus.UCSC.mm10.knownGene", "ggplot2",
              "org.Mm.eg.db", "kableExtra", "regioneR", "BSgenome.Mmusculus.UCSC.mm10.masked", "UpSetR",
              "ComplexUpset", "metap", "pheatmap", "magrittr", "Hmisc", "tidyverse")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))
enrichR:::.onAttach() 

# Load DMRs ---------------------------------------------------------------

setwd("/share/lasallelab/Ben/PEBBLES/DNA/DMRs/")
dir.create("overlaps")

loadDMRs <- function(name){
  load(glue::glue("{name}/RData/DMRs.RData"))
  assign(name, sigRegions, envir = .GlobalEnv)
  rm(sigRegions, regions)
}

contrasts <- c("placenta_male", "placenta_female", "brain_male", "brain_female")
  
purrr::walk(contrasts, loadDMRs)

# Genomic coordinate ------------------------------------------------------

dir.create("overlaps/coordinate")

purrr::walk(c("male", "female"), function(sex){
  
  print(glue::glue("Obtaining genomic coordinate overlaps for {sex} samples"))
  
  res <- suppressMessages(ChIPpeakAnno::makeVennDiagram(Peaks = list(get(paste0("placenta_",sex)),
                                                                     get(paste0("brain_",sex))), 
                                                        NameOfPeaks = c(paste0("placenta_",sex),
                                                                        paste0("brain_",sex))))

  
  print(glue::glue("Plotting Venn of genomic coordinate overlaps for {sex} samples"))
  # ref: https://support.bioconductor.org/p/67429/
  venn_cnt2venn <- function(venn_cnt){
    n <- which(colnames(venn_cnt) == "Counts") - 1
    SetNames <- colnames(venn_cnt)[1:n]
    Weight <- venn_cnt[,"Counts"]
    names(Weight) <- apply(venn_cnt[,1:n], 1, paste, collapse="")
    Vennerable::Venn(SetNames = SetNames, Weight = Weight)
  }
  v <- venn_cnt2venn(res$vennCounts)
  
  svg(glue::glue("overlaps/coordinate/{sex}_venn_coordinate.svg"),
      height = 8.5,
      width = 12)
  
  plot(v)
  
  dev.off()
  
  print(glue::glue("Annotating genomic coordinate overlaps for {sex} samples"))
  
  shared <- c(get(paste0("placenta_",sex)),
              get(paste0("brain_",sex))) %>% 
    GenomicRanges::sort() %>%
    plyranges::reduce_ranges() %>%
    plyranges::filter_by_overlaps(get(paste0("placenta_",sex))) %>% 
    plyranges::filter_by_overlaps(get(paste0("brain_",sex)))
  
  sharedAnnotated <- shared %>%
    ChIPseeker::annotatePeak(TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                             annoDb = "org.Mm.eg.db",
                             overlap = "all",
                             verbose = FALSE) %>%
    dplyr::as_tibble() %>%
    dplyr::select(seqnames, start, end, width, annotation, geneSymbol = SYMBOL, gene = GENENAME)
  
  sharedOverlaps <- . %>%
    plyranges::filter_by_overlaps(shared) %>%
    DMRichR::annotateRegions(TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                             annoDb = "org.Mm.eg.db") %>%
    dplyr::select(geneSymbol, difference, `p-value`)
  
  brainShared <- get(paste0("brain_",sex)) %>%
    sharedOverlaps()

  placentaShared <- get(paste0("placenta_",sex)) %>% 
    sharedOverlaps
  
  sharedAnnotated %>% 
    dplyr::left_join(placentaShared) %>%
    dplyr::rename("placenta p-value" = `p-value`, "placenta Difference" = difference) %>% 
    dplyr::left_join(brainShared) %>% 
    dplyr::rename("brain p-value" = `p-value`, "brain Difference" = difference) %>% 
    dplyr::rename(chr = "seqnames", symbol = geneSymbol, Name = "gene") %>% 
    dplyr::rename_with(Hmisc::capitalize) %>% 
    dplyr::mutate(Name = Hmisc::capitalize(Name)) %>% 
    dplyr::mutate(Annotation = gsub(" \\(.*","", Annotation)) %>% 
    dplyr::arrange(`Placenta p-value`) %>% 
    openxlsx::write.xlsx(file = glue::glue("overlaps/coordinate/{sex}_coordinate_overlaps_annotated.xlsx"))
  
  print(glue::glue("Creating a kable of genomic coordinate overlaps for {sex} samples"))
  
  readxl::read_xlsx(glue::glue("overlaps/coordinate/{sex}_coordinate_overlaps_annotated.xlsx"),
                    col_types = c(rep("text", 7), rep(c("text","numeric"), 2))) %>%
    dplyr::mutate_if(is.numeric, function(x) {
      formatC(x, digits = 1, format = "e", drop0trailing = TRUE)
      }) %>% 
    dplyr::mutate("Placenta Difference" = paste0(`Placenta Difference`, "%"),
                  "Brain Difference" = paste0(`Brain Difference`, "%")) %>% 
    kbl(align = c(rep("l",7), rep(c("r", "l"),2)),
        col.names = c("Chr", "Start", "End", "Width",
                      "Region", "Symbol", "Name",
                      "Difference", "p-value", "Difference", "p-value")) %>%
    kable_classic(full_width = FALSE, html_font = "Cambria") %>%
    column_spec(6, italic = TRUE) %>% 
    add_header_above(c("Coordinates" = 4, "Annotation" = 3, "Placenta" = 2, "Brain" = 2)) %>%
    kableExtra::save_kable(file = glue::glue("overlaps/coordinate/{sex}_coordinate_overlaps.html"))
  
  print(glue::glue("{sex} genomic coordinate overlaps pipeline is complete"))
})

# regioneR ----------------------------------------------------------------

dir.create("overlaps/regioneR")

getConsensus <- . %>%
  plyranges::as_granges() %>% 
  GenomicRanges::sort() %>%
  plyranges::reduce_ranges()

placenta_both <- c(placenta_male, placenta_female) %>% 
  getConsensus

brain_both <- c(brain_male, brain_female) %>% 
  getConsensus

purrr::walk(c("male", "female", "both"), function(sex){
  
  print(glue::glue("Counting overlaps for {sex} samples"))
  
  print(regioneR::numOverlaps(A = get(paste0("brain_",sex)), 
                              B = get(paste0("placenta_",sex)),
                              count.once = TRUE))
  
  print(glue::glue("Running permutation test for {sex} samples"))
  
  pt <- regioneR::overlapPermTest(A = get(paste0("brain_",sex)), 
                                  B = get(paste0("placenta_",sex)), 
                                  alternative = "greater",
                                  genome = "mm10",
                                  ntimes = 10000,
                                  count.once = TRUE)
  
  pdf(glue::glue("overlaps/regioneR/regioneR_{sex}_Brain_Placenta_Overlap.pdf"),
      height = 7.50,
      width = 11.50)
  
  plot(pt)
  
  dev.off()
  
  print(glue::glue("Calculating local Z-score for {sex} samples"))
  
  lz <- regioneR::localZScore(A = get(paste0("brain_",sex)), 
                              B = get(paste0("placenta_",sex)),
                              pt = pt, 
                              window = 10*mean(width(get(paste0("brain_",sex)))), 
                              step = mean(width(get(paste0("brain_",sex))))/2,
                              count.once = TRUE)
  
  pdf(glue::glue("overlaps/regioneR/regioneR_{sex}_Brain_Placenta_Overlap_Zscore.pdf"),
      height = 7.50,
      width = 11.50)
  
  plot(lz)
  
  dev.off()
  
  print(glue::glue("Saving data for {sex} samples"))
  
  save(pt, lz, file = glue::glue("overlaps/regioneR/regioneR_{sex}_brain_placenta_overlap.RData"))
  
  print(glue::glue("regioneR Brain DMR within Placenta DMR overlap finished for {sex} samples"))
})

# Gene symbol -------------------------------------------------------------

getSymbol <- function(regions = sigRegions,
                      TxDb = TxDb,
                      annoDb = annoDb){
  regions %>% 
    DMRichR::annotateRegions(TxDb = TxDb,
                             annoDb = annoDb) %>%
    #dplyr::filter(annotation != "Distal Intergenic") %>% 
    dplyr::distinct() %>% 
    purrr::pluck("geneSymbol") %>%
    na.omit()
}

contrasts <- c("placenta_male",
               "placenta_female",
               "brain_male",
               "brain_female")

print(glue::glue("Annotating all DMRs"))

genes <- contrasts %>%
  purrr::set_names() %>% 
  purrr::map(function(contrast){
    getSymbol(get(contrast),
              TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
              annoDb = "org.Mm.eg.db")
  })

dir.create("overlaps/symbol")

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
                              axis.text = element_text(size = 20),
                              axis.title = element_text(size = 20)
                            ) +
                            scale_y_continuous(expand = c(0, 0, 0.1, 0))
                        ),
                        sort_sets = FALSE,
                        queries = list(
                          upset_query(
                            intersect = c("Placenta Female", "Placenta Male",
                                          "Brain Female", "Brain Male"),
                            color = '#E41A1C',
                            fill = '#E41A1C',
                            only_components = c('intersections_matrix', 'Intersection size')
                          ),
                          upset_query(
                            intersect = c("Placenta Female", "Placenta Male"),
                            color = '#FF7F00',
                            fill = '#FF7F00',
                            only_components = c('intersections_matrix', 'Intersection size')
                          ),
                          upset_query(
                            intersect = c("Brain Female", "Brain Male"),
                            color = '#4DAF4A',
                            fill = '#4DAF4A',
                            only_components = c('intersections_matrix', 'Intersection size')
                          ),
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
                              axis.text = element_text(size = 20),
                              axis.title = element_blank()
                            ),
                            'overall_sizes' = theme(
                              axis.title = element_text(size = 14)
                            )
                          )
                        )
    )
}

pdf(glue::glue("overlaps/symbol/UpSet_Genes_All_Symbol.pdf"),
    height = 6, width = 10, onefile = FALSE)

print({
  genes %>%
    magrittr::set_names(c("Placenta Male", "Placenta Female",
                          "Brain Male", "Brain Female")) %>% 
    rev() %>% 
    geneUpsetPlot()
})

dev.off()

purrr::walk(c("all", "male", "female"), function(contrast){
  
  print(glue::glue("Generating gene symbol overlaps for {contrast} samples"))
  
  if(contrast == "male"){
    contrasts <- c("placenta_male", "brain_male")
  }else if(contrast == "female"){
    contrasts <- c("placenta_female", "brain_female")
  }
  
  genes <- genes[contrasts]
  
  v <- Vennerable::Venn(genes,
                        SetNames = contrasts)
  
  print(glue::glue("Saving gene symbol Venn of overlaps for {contrast} samples"))
  
  svg(glue::glue("overlaps/symbol/{contrast}_symbols_Venn.svg"),
      height = 8.5,
      width = 12)
  
  plot(v)
  
  dev.off()
  
  print(glue::glue("Performing GO analysis for {contrast} samples"))
  
  intersects <- dplyr::case_when(contrast == "all" ~ "1111",
                                 contrast == "male" ~ "11",
                                 contrast == "female" ~ "11")
  
  v@IntersectionSets[[intersects]] %T>%
    openxlsx::write.xlsx(glue::glue("overlaps/symbol/{contrast}_symbol_overlaps_annotated.xlsx")) %>%
    enrichR::enrichr(c("GO_Biological_Process_2018",
                       "GO_Cellular_Component_2018",
                       "GO_Molecular_Function_2018",
                       "KEGG_2019_Mouse",
                       "Panther_2016",
                       "Reactome_2016",
                       "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")) %T>% # %>% 
    #purrr::map(~ dplyr::filter(., Adjusted.P.value < 0.05)) %T>%
    openxlsx::write.xlsx(file = glue::glue("overlaps/symbol/{contrast}_symbol_overlaps_enrichr.xlsx")) %>%
    DMRichR::slimGO(tool = "enrichR",
                    annoDb = "org.Mm.eg.db",
                    plots = TRUE) %T>%
    openxlsx::write.xlsx(file =glue::glue("overlaps/symbol/{contrast}_symbol_overlaps_enrichr_rrvgo_results.xlsx")) %>% 
    DMRichR::GOplot() %>% 
    ggplot2::ggsave(glue::glue("overlaps/symbol/{contrast}_symbol_overlaps_enrichr_plot.pdf"),
                    plot = .,
                    device = NULL,
                    height = 8.5,
                    width = 10)
  
  print(glue::glue("Finished gene symbol overlap pipe for {contrast} samples"))
  
})

# Meta p-value ------------------------------------------------------------

print(glue::glue("Running meta p-value analysis for female and male samples"))

contrasts <- c("female", "male") 

contrasts %>% 
  purrr::set_names() %>% 
  purrr::map(function(sex){
    DMRichR::read_excel_all(glue::glue("overlaps/symbol/{sex}_symbol_overlaps_enrichr.xlsx")) %>%
      data.table::rbindlist(idcol = "Gene Ontology") %>%
      dplyr::as_tibble()}) %>%
  purrr::map(~ dplyr::select(., Term, P.value, "Gene Ontology")) %>% 
  purrr::map2(contrasts, ~ data.table::setnames(.x, "P.value", .y)) %>% 
  purrr::reduce(dplyr::inner_join, by = c("Term", "Gene Ontology")) %>%
  dplyr::rowwise(Term, "Gene Ontology") %>% 
  dplyr::mutate(meta_p = metap::sumlog(dplyr::c_across(where(is.numeric)))$p) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(meta_p) %>%
  dplyr::mutate(P.value = stats::p.adjust(meta_p, method = 'fdr')) %>% 
  dplyr::filter(P.value < 0.05, female < 0.05, male < 0.05) %>% 
  split(.$`Gene Ontology`) %>%
  purrr::map(~ dplyr::select(., -`Gene Ontology`)) %T>%
  openxlsx::write.xlsx(file = "overlaps/symbol/meta_p_all_symbol_overlaps_enrichr.xlsx") %>%
  DMRichR::slimGO(tool = "enrichR",
                  annoDb = "org.Mm.eg.db",
                  plots = FALSE) %T>%
  openxlsx::write.xlsx(file = "overlaps/symbol/meta_p_all_symbol_overlaps_enrichr_rrvgo_results.xlsx")

## Production plot --------------------------------------------------------

plotData <- DMRichR::read_excel_all("overlaps/symbol/meta_p_all_symbol_overlaps_enrichr.xlsx") %>%
  data.table::rbindlist(idcol = "Gene Ontology") %>%
  dplyr::as_tibble() %>%
  dplyr::filter(`Gene Ontology` %in% c( "Panther_2016", "RNA-Seq_Disease_Gene_and_Drug_S")) %>% # "KEGG_2019_Mouse",
  dplyr::filter(P.value < 0.05, female < 0.05, male < 0.05) %>% 
  dplyr::mutate("-log10.p-value" = -log10(P.value)) %>% 
  dplyr::select("Gene Ontology", Term, "-log10.p-value") %>%
  dplyr::bind_rows(readxl::read_xlsx("overlaps/symbol/meta_p_all_symbol_overlaps_enrichr_rrvgo_results.xlsx"), .) %>%
  dplyr::mutate("Database" = dplyr::recode_factor(`Gene Ontology`,
                                                  "Biological Process" = "GO Biological Process",
                                                  "Cellular Component" = "GO Cellular Component",
                                                  "Molecular Function" = "GO Molecular Function",
                                                  #"KEGG_2019_Mouse" = "KEGG Pathways",
                                                  "Panther_2016" = "Panther Pathways",
                                                  "RNA-Seq_Disease_Gene_and_Drug_S" = "GEO RNA-seq Disease and Drug")) %>%
  dplyr::select(-`Gene Ontology`) %>% 
  dplyr::mutate(Term = stringr::str_trim(Term)) %>%
  dplyr::mutate(Term = Hmisc::capitalize(Term)) %>%
  dplyr::mutate(Term = stringr::str_remove(Term, "Homo sapiens.*$")) %T>%
  openxlsx::write.xlsx("overlaps/symbol/meta_p_all_symbol_overlaps_enrichr_plot_table.xlsx") %>%
  dplyr::mutate(Term = stringr::str_trunc(Term, 50)) %>% 
  dplyr::group_by(Database) %>%
  dplyr::slice(1:7) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(Term = stringr::str_replace(Term, "up", "Up"),
                Term = stringr::str_replace(Term, "down", "Down")) %>% 
  dplyr::mutate(Term = factor(.$Term, levels = unique(.$Term[order(forcats::fct_rev(.$Database), .$`-log10.p-value`)]))) 

metaPlot <- function(data){
  data %>% 
    ggplot2::ggplot(aes(x = Term,
                        y = `-log10.p-value`,
                        fill = Database,
                        group = Database)) +
    ggplot2::geom_bar(stat = "identity",
                      position = position_dodge(),
                      color = "Black") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggsci::scale_fill_d3() +
    ggplot2::labs(y = expression("-log"[10](italic(q)))) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = element_text(size = 40),
                   axis.title.y = element_blank(),
                   axis.title.x = element_text(size = 25),
                   legend.position = "none") 
} 

cowplot::plot_grid(plotData %>%
                     dplyr::filter(Database != "GEO RNA-seq Disease and Drug") %>% 
                     metaPlot() +
                     ggplot2::scale_y_continuous(breaks = c(0,5,10),
                                                 expand = c(0, 0)) +
                     ggplot2::theme(axis.title.x = element_blank()),
                   plotData %>%
                     dplyr::filter(Database == "GEO RNA-seq Disease and Drug") %>% 
                     metaPlot() + 
                     ggplot2::scale_fill_manual(values = "#9467BDFF"), # ggsci::pal_d3()(5)
                   align = "v",
                   ncol = 1,
                   rel_heights = c(11,4)) %>%
  ggplot2::ggsave(glue::glue("overlaps/symbol/meta_p_all_production_symbol_overlaps_enrichr_plot.pdf"),
                  plot = ., 
                  width = 14,
                  height = 16)

# Heatmap -----------------------------------------------------------------

# Modified from UpSetR
# https://github.com/hms-dbmi/UpSetR/blob/master/R/fromList.R
fromList2 <- function(input){
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x){x <- as.vector(match(elements, x))}))
  data[is.na(data)] <- as.integer(0); data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) !=0), ]
  names(data) <- names(input)
  rownames(data) <- elements %>%
    stringr::str_to_lower() %>%
    Hmisc::capitalize()
  return(data)
}

# Modified from: https://stackoverflow.com/a/60177668
make_italics <- function(x) {
  as.expression(lapply(rownames(x), function(y) bquote(italic(.(y)))))
}

## Make list --------------------------------------------------------------

geneList <- c("female", "male") %>% 
  purrr::set_names() %>% 
  purrr::map_dfr(function(sex){
    DMRichR::read_excel_all(glue::glue("overlaps/symbol/{sex}_symbol_overlaps_enrichr.xlsx")) %>%
      data.table::rbindlist(idcol = "Gene Ontology") %>%
      dplyr::as_tibble() %>% 
      dplyr::mutate(Genes = Genes %>% 
                      purrr::map(~ stringr::str_split(., pattern = ";")) %>%
                      purrr::flatten()) %>% 
      dplyr::mutate("Database" = dplyr::recode_factor(`Gene Ontology`,
                                                      "GO_Biological_Process_2018" = "GO Biological Process",
                                                      "GO_Cellular_Component_2018" = "GO Cellular Component",
                                                      "GO_Molecular_Function_2018" = "GO Molecular Function",
                                                      "Panther_2016" = "Panther Pathways",
                                                      "RNA-Seq_Disease_Gene_and_Drug_S" = "GEO RNA-seq Disease and Drug")) %>% 
      dplyr::mutate(Term = stringr::str_remove(.$Term, "\\(GO.*")) %>%
      dplyr::mutate(Term = stringr::str_remove(Term, "Homo sapiens.*$")) %>%
      dplyr::mutate(Term = stringr::str_trim(Term)) %>%
      dplyr::mutate(Term = Hmisc::capitalize(Term)) %>%
      dplyr::select(Term, Database, Genes)
  }, .id = "Sex") %>% 
  dplyr::group_by(Term) %>%
  dplyr::filter(dplyr::n() > 1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Sex = Hmisc::capitalize(Sex)) %>% 
  tidyr::pivot_wider(names_from = Sex, values_from = Genes) %>%
  dplyr::inner_join(readxl::read_xlsx("overlaps/symbol/meta_p_all_symbol_overlaps_enrichr_plot_table.xlsx"),
                    .,
                    by = c("Term", "Database")) %>% 
  dplyr::mutate(Shared = purrr::map2(Female, Male, intersect)) %T>%
  openxlsx::write.xlsx("overlaps/symbol/enrichr_meta_p.xlsx") %>% 
  dplyr::filter(Database == "GEO RNA-seq Disease and Drug") %>% 
  dplyr::select(Term, `-log10.p-value`, Database, Shared) %>% 
  dplyr::slice(1:7) %>%
  purrr::pluck("Shared") %>% 
  purrr::set_names(c("MeCP2 Hypothalamus Knockout Up",
                     "Topotecan Cortical neurons 300 nM Down",
                     "MeCP2 Hypothalamus Transgenic Down",
                     "LPS Neuron Down",
                     "TAF15 Striatum Knockdown Down",
                     "Bicuculin Hippocampus 20 uM Down",
                     "MeCP2 Visual Cortex Knockout Up"))

## Main plot --------------------------------------------------------------

geneList %>% 
  magrittr::extract(c("MeCP2 Hypothalamus Knockout Up", "MeCP2 Hypothalamus Transgenic Down")) %>% 
  magrittr::set_names(c("MeCP2 Knockout Up", "MeCP2 Transgenic Down")) %>% 
  fromList2() %>% 
  pheatmap::pheatmap(angle_col = 45,
                     legend = FALSE,
                     labels_row = make_italics(.),
                     border_color = "black",
                     treeheight_col = 10, 
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     filename = "overlaps/symbol/Rett_heatmap.pdf",
                     width = 1.1,
                     height = 8)

## Supplementary Plot -----------------------------------------------------

geneList %>% 
  fromList2() %>% 
  pheatmap::pheatmap(angle_col = 45,
                     legend = FALSE,
                     labels_row = make_italics(.),
                     border_color = "black",
                     treeheight_col = 10, 
                     cluster_rows = FALSE,
                     filename = "overlaps/symbol/full_NDD_heatmap.pdf",
                     width = 2,
                     height = 12)
