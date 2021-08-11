# Functional DMR enrichment spreadsheets and plots
# GO, HOMER, PANTHER pathways, MEME AME, CpG and gene annotations, and chromHMM
# Recreates Figure 2, Supplementary Figure 1, and Figure 3
# Ben Laufer

BiocManager::install('robertamezquita/marge', ref = 'master')
packages <- c( "magrittr", "marge", "tidyverse", "ggsci", "cowplot")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

setwd("/Users/blaufer/Box Sync/PEBBLES/DNA/DMRs")

# Figure 2 ----------------------------------------------------------------

## GO ---------------------------------------------------------------------

GO <- tidyr::crossing(tissue = c("placenta", "brain"),
                      sex = c("female", "male")) %>% 
  purrr::pmap_dfr(function(tissue, sex){
    readxl::read_xlsx(glue::glue("{tissue}_{sex}/Ontologies/GOfuncR_slimmed_results.xlsx")) %>%
      dplyr::mutate(tissue = Hmisc::capitalize(tissue),
                    sex = Hmisc::capitalize(sex)) %>%
      dplyr::select(tissue, sex, Term, `Gene Ontology`, `-log10.p-value`) %>% 
      dplyr::left_join(readxl::read_xlsx(glue::glue("{tissue}_{sex}/Ontologies/GOfuncR.xlsx")) %>%
                         dplyr::select(Term = node_name, FWER = FWER_overrep))}) %T>%
  openxlsx::write.xlsx("slimmed_GOfuncR_enrichments.xlsx", overwrite = TRUE) %>%
  tidyr::unite("Contrast", c(tissue, sex), sep = " ") %>%
  dplyr::mutate("Database" = dplyr::recode_factor(`Gene Ontology`,
                                                  "Biological Process" = "GO Biological Process",
                                                  "Cellular Component" = "GO Cellular Component",
                                                  "Molecular Function" = "GO Molecular Function")) %>%
  dplyr::select(-`Gene Ontology`) %>% 
  dplyr::mutate(Term = stringr::str_trim(Term)) %>%
  dplyr::mutate(Term = Hmisc::capitalize(Term)) %>%
  dplyr::mutate(Term = stringr::str_trunc(Term, 40)) 

GOplot <- function(contrast = contrast){
  GO %>%
    dplyr::filter(Contrast == !!contrast) %>% 
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
    ggplot2::theme(text = element_text(size = 14),
                   axis.title.y = element_blank(),
                   axis.title.x = element_text(size = 12),
                   legend.position = "none") %>%
    return()
}

## HOMER ------------------------------------------------------------------

homer <- tidyr::crossing(tissue = c("placenta", "brain"),
                         sex = c("female", "male")) %>% 
  purrr::pmap_dfr(function(tissue, sex){
    marge::read_known_results(glue::glue("{tissue}_{sex}/Extra/HOMER/both")) %>%
      dplyr::mutate(tissue = Hmisc::capitalize(tissue),
                    sex = Hmisc::capitalize(sex)) %>%
      dplyr::select(tissue, sex,
                    motif_name, motif_family, experiment, accession, consensus,
                    "-log10.p-value" = log_p_value, fdr, tgt_num, tgt_pct, bgd_num, bgd_pct)
  }) %T>%
  openxlsx::write.xlsx("homer_enrichments.xlsx", overwrite = TRUE) %>%
  tidyr::unite("Contrast", c(tissue, sex), sep = " ") %>% 
  dplyr::mutate(Motif = paste0(motif_name," (",motif_family,")"),
                Contrast = as.factor(Contrast)) %>%
  dplyr::select(Contrast, Motif, `-log10.p-value`) %>%
  dplyr::group_by(Contrast) %>%
  dplyr::slice(1:10) %>%
  dplyr::arrange(dplyr::desc(`-log10.p-value`)) %>% 
  dplyr::ungroup()

homerPlot <- function(contrast = contrast){
  homer %>%
    dplyr::filter(Contrast == !!contrast) %>% 
    ggplot2::ggplot(ggplot2::aes(x = reorder(Motif,`-log10.p-value`),
                                 y = `-log10.p-value`,
                                 fill = Contrast)) +
    ggplot2::geom_bar(stat = "identity",
                      position = ggplot2::position_dodge(),
                      color = "Black",
                      size = 0.25) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = c(0, 0, 0.1, 0),
                                breaks = scales::pretty_breaks(3)) +
    ggsci::scale_fill_simpsons() +
    ggplot2::labs(y = expression("-log"[10](italic(p)))) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = element_text(size = 14),
                   axis.title.y = element_blank(),
                   axis.title.x = element_text(size = 12),
                   legend.position = "none") %>% 
    return()
}

## Plot --------------------------------------------------------------------

cowplot::plot_grid(GOplot("Placenta Female"),
                   GOplot("Placenta Male"),
                   GOplot("Brain Female"),
                   GOplot("Brain Male"),
                   homerPlot("Placenta Female"),
                   homerPlot("Placenta Male"),
                   homerPlot("Brain Female"),
                   homerPlot("Brain Male"),
                   ncol = 4,
                   #align = "v",
                   labels = c("A)", "B)", "C)", "D)",
                              "E)", "F)", "G)", "H)"),
                   label_size = 16,
                   rel_heights = c(1.8, 1.3),
                   scale = c(rep(0.95,4), rep(0.75,4))) %>%
  ggplot2::ggsave("DMR functional enrichments.pdf",
                  plot = ., 
                  width = 16,
                  height = 8)

# Supplementary Figure 2 --------------------------------------------------

## PANTHER Pathways -------------------------------------------------------

PANTHER <- tidyr::crossing(tissue = c("placenta", "brain"),
                           sex = c("female", "male")) %>% 
  purrr::pmap_dfr(function(tissue, sex){
    readxl::read_xlsx(glue::glue("{tissue}_{sex}/Ontologies/enrichr.xlsx"), sheet = "Panther_2016") %>%
      dplyr::mutate("-log10.p-value" = -log10(Adjusted.P.value),
                    Tissue = Hmisc::capitalize(tissue),
                    Sex = Hmisc::capitalize(sex))}) %>%
  dplyr::filter(Adjusted.P.value < 0.1) %>% 
  dplyr::select(Tissue, Sex, Term, "-log10.p-value", Combined.Score, Overlap, Genes) %T>%
  openxlsx::write.xlsx("PANTHER_enrichments.xlsx", overwrite = TRUE) %>%
  tidyr::unite("Contrast", c(Tissue, Sex), sep = " ") %>%
  dplyr::mutate(Term = stringr::str_trim(Term)) %>%
  dplyr::mutate(Term = Hmisc::capitalize(Term)) %>%
  dplyr::mutate(Term = stringr::str_remove(Term, "Homo sapiens.*$"))

pathwayPlot <- function(contrast = contrast){
  PANTHER %>%
    dplyr::filter(Contrast == !!contrast) %>% 
    dplyr::mutate(Database = "Panther Pathways") %>% 
    dplyr::mutate(Term = stringr::str_trunc(Term, 50)) %>% 
    dplyr::slice(1:10) %>%
    dplyr::mutate(Term = factor(.$Term, levels = rev(.$Term))) %>%
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
                                breaks = scales::pretty_breaks(2)) +
    ggplot2::scale_fill_manual(values = ggsci::pal_d3()(4)[4]) + 
    ggplot2::labs(y = expression("-log"[10](italic(q)))) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = element_text(size = 14),
                   axis.title.y = element_blank(),
                   axis.title.x = element_text(size = 12),
                   legend.position = "none") %>%
    return()
}

## MEME AME ---------------------------------------------------------------

.libPaths("/share/lasallelab/programs/DMRichR/R_4.1")
AnnotationHub::setAnnotationHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_4.1")
BiocManager::install('robertamezquita/marge', ref = 'master')
packages <- c("memes", "BSgenome.Mmusculus.UCSC.mm10", "tidyverse")

# Run this part on the cluster
setwd("/share/lasallelab/Ben/PEBBLES/DNA/DMRs/")
dir.create("memes")

#options(meme_bin = "/Users/blaufer/opt/miniconda3/pkgs/meme-5.3.0-py37pl5262he48d6b8_2/bin") 
options(meme_bin = "/software/meme/5.3.3/lssc0-linux/bin") 
memes::check_meme_install()
options(meme_db = "motif_databases/METHYLCYTOSINE/yin2017.meme")
goi <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
contrasts <- c("placenta_male", "placenta_female", "brain_male", "brain_female")

methylMeme <- function(sigRegions = sigRegions,
                       regions = regions,
                       goi = goi){
  sigRegions %>%
    memes::get_sequence(goi) %>%
    memes::runAme(control = regions %>% 
                    memes::get_sequence(goi))
}

parallel::mclapply(contrasts, 
                   function(contrast){
                     load(glue::glue("../{contrast}/RData/DMRs.RData"))
                     
                     sigRegions %>% 
                       methylMeme(regions = regions,
                                  goi = goi)  %>%
                       openxlsx::write.xlsx(file = glue::glue("memes/{contrast}_meme_methylcytosine_results.xlsx"))
                   },
                   mc.cores = length(contrasts),
                   mc.silent = TRUE)

# Run on desktop
setwd("/Users/blaufer/Box Sync/PEBBLES/DNA/DMRs")

ame <- tidyr::crossing(tissue = c("placenta", "brain"),
                       sex = c("female", "male")) %>% 
  purrr::pmap_dfr(function(tissue, sex){
    readxl::read_xlsx(glue::glue("memes/{tissue}_{sex}_meme_methylcytosine_results.xlsx")) %>%
      dplyr::rename("Motif ID" = motif_id) %>% 
      dplyr::mutate("Transcription Factor" = stringr::word(`Motif ID`, sep = "-"),
                    Domain = dplyr::case_when(stringr::str_detect(`Motif ID`, "eDBD") ~ "Extended DNA-binding Domain",
                                              stringr::str_detect(`Motif ID`, "FL") ~ "Full Length"),
                    Methylated = dplyr::case_when(stringr::str_detect(`Motif ID`, "methyl") ~ "Methylated",
                                                  TRUE ~ "Not Methylated"),
                    "-log10.E-value" = -log10(evalue),
                    Tissue = Hmisc::capitalize(tissue),
                    Sex = Hmisc::capitalize(sex)) %>%
      dplyr::select(Tissue, Sex, `Motif ID`, `Transcription Factor`, Domain, Methylated, `-log10.E-value`, adj.pvalue)}) %T>%
  openxlsx::write.xlsx("meme_methylcytosine_enrichments.xlsx", overwrite = TRUE) %>%
  tidyr::unite("Contrast", c(Tissue, Sex), sep = " ")

amePlot <- function(contrast = contrast){
  ame %>% 
    dplyr::filter(Contrast == !!contrast) %>% 
    dplyr::slice(1:10) %>% 
    ggplot2::ggplot(ggplot2::aes(x = reorder(`Motif ID`,`-log10.E-value`),
                                 y = `-log10.E-value`,
                                 fill = Methylated)) +
    ggplot2::geom_bar(stat = "identity",
                      position = ggplot2::position_dodge(),
                      color = "Black",
                      size = 0.25) +
    ggplot2::scale_fill_manual(values = c("Methylated" = "black",
                                          "Not Methylated" = "white")) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = c(0, 0, 0.1, 0),
                                breaks = scales::pretty_breaks(2)) +
    ggplot2::labs(y = expression("-log"[10]("E-value"))) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = element_text(size = 14),
                   axis.title.y = element_blank(),
                   axis.title.x = element_text(size = 12),
                   legend.position = "none") %>% 
    return()
}

## Plot --------------------------------------------------------------------

cowplot::plot_grid(pathwayPlot("Placenta Female"),
                   pathwayPlot("Placenta Male"),
                   pathwayPlot("Brain Female"),
                   pathwayPlot("Brain Male"),
                   amePlot("Placenta Female"),
                   amePlot("Placenta Male"),
                   amePlot("Brain Female"),
                   amePlot("Brain Male"),
                   ncol = 4,
                   align = "v",
                   labels = c("A)", "B)", "C)", "D)",
                              "E)", "F)", "G)", "H)"),
                   label_size = 16,
                   rel_heights = c(10, 10),
                   scale = 0.95) %>%
  ggplot2::ggsave("DMR functional enrichments supplement.pdf",
                  plot = ., 
                  width = 18,
                  height = 5)

# Figure 3 ----------------------------------------------------------------

## chromHMM ---------------------------------------------------------------

# Run chromHMM.R first
dir.create("chromHMM")

data <- tidyr::crossing(source = c("placenta", "brain"),
                        sex = c("female", "male")) %>% 
  purrr::pmap_dfr(function(source, sex){
    readr::read_tsv(glue::glue("{source}_{sex}/Extra/LOLA/All DMRs/Development_ChromHMM/allEnrichments.tsv")) %>%
      dplyr::mutate(source = source,
                    sex = sex,
                    contrast = paste0(sex,"_",source),
                    tissue = tissue %>% 
                      Hmisc::capitalize() %>% 
                      stringr::str_replace_all("_", " "),
                    timepoint = stringr::str_split(description," ") %>%
                      sapply(magrittr::extract2, 3)) %>%
      dplyr::as_tibble() %>%
      dplyr::select(contrast, source, sex, qValue, oddsRatio, tissue, timepoint, antibody) %>%
      dplyr::mutate(antibody = as.factor(antibody)) %>% 
      dplyr::mutate(timepoint = dplyr::case_when(timepoint == "GD0" ~ "PD0",
                                                 TRUE ~ timepoint)) %>% 
      dplyr::mutate(antibody = dplyr::recode_factor(antibody,
                                                    "Tss" = "Active TSS",
                                                    "TssFlnk" = "Flanking Active TSS",
                                                    "Tx" = "Strong Transcription",
                                                    "TxWk" = "Weak Transcription",
                                                    "EnhG"= "Genic Enhancer",
                                                    "Enh" = "Enhancers",
                                                    "EnhLo" = "Weak Enhancer",
                                                    "EnhPois" = "Poised Enhancers",
                                                    "EnhPr" = "Primed Enhancer",
                                                    "TssBiv" = "Bivalent TSS",
                                                    "ReprPC" = "Repressed PolyComb",
                                                    "ReprPCWk" = "Weak Repressed PolyComb",
                                                    "QuiesG" = "Quiescent Gene",
                                                    "Quies" = "Quiescent 1",
                                                    "Quies2" = "Quiescent 2",
                                                    "Quies3" = "Quiescent 3",
                                                    "Quies4" = "Quiescent 4",
                                                    "Het" = "Heterochromatin"),
                    itemRgb = dplyr::case_when(antibody == "Active TSS" ~ "#FF0000",
                                               antibody == "Flanking Active TSS" ~ "#FF4500",
                                               antibody == "Strong Transcription" ~ "#008000",
                                               antibody == "Weak Transcription" ~ "#3F9A50",
                                               antibody == "Genic Enhancer" ~ "#AADF07",
                                               antibody == "Enhancers" ~ "#FFDF00",
                                               antibody == "Weak Enhancer" ~ "#FFFF80",
                                               antibody == "Poised Enhancers" ~ "#BDB76B",
                                               antibody == "Primed Enhancer" ~ "#BDB76B",
                                               antibody == "Bivalent TSS" ~ "#CD5C5C",
                                               antibody == "Repressed PolyComb" ~ "#8937DF",
                                               antibody == "Weak Repressed PolyComb" ~ "#9750E3",
                                               antibody == "Quiescent Gene" ~ "#808080",
                                               antibody == "Quiescent 1" ~ "#DCDCDC",
                                               antibody == "Quiescent 2" ~ "#DCDCDC",
                                               antibody == "Quiescent 3" ~ "#DCDCDC",
                                               antibody == "Quiescent 4" ~ "#DCDCDC",
                                               antibody == "Heterochromatin" ~ "#4B0082")) 
  }) %>%
  dplyr::mutate(contrast = as.factor(contrast)) %>% 
  dplyr::mutate(contrast = dplyr::recode_factor(contrast,
                                                "male_brain" = "Brain Male",
                                                "female_brain" = "Brain Female",
                                                "male_placenta" = "Placenta Male",
                                                "female_placenta" = "Placenta Female")) %>%
  dplyr::filter(tissue %in% c("Embryonic facial prominence", "Neural tube",
                              "Forebrain","Midbrain", "Hindbrain")) %>% 
  dplyr::mutate(tissue = dplyr::recode_factor(tissue,
                                              "Embryonic facial prominence" = "Facial Prominence",
                                              "Neural tube" = "Neural Tube",
                                              "Forebrain" = "Forebrain",
                                              "Midbrain" = "Midbrain",
                                              "Hindbrain" = "Hindbrain")) %>% 
  dplyr::filter(tissue == "Forebrain") %>% 
  dplyr::arrange(qValue)

data %>%
  dplyr::select(-contrast) %>%
  openxlsx::write.xlsx("chromHMM/chromHMM_forebrain_enrichments.xlsx")

### Bar plot ---------------------------------------------------------------

barData <- data %>%
  dplyr::group_by(contrast, antibody, itemRgb) %>%
  dplyr::slice_max(oddsRatio) %>%
  ungroup() %>%
  dplyr::select(contrast, qValue, oddsRatio, antibody, itemRgb, timepoint) %>%
  dplyr::mutate(fold = dplyr::case_when(oddsRatio < 1 ~ -1/oddsRatio,
                                        oddsRatio >= 1 ~ oddsRatio)) %>%
  dplyr::mutate(signif = dplyr::case_when(qValue <= 0.05 ~ 1,
                                          qValue > 0.05 ~ 0)) %>%
  dplyr::mutate(contrast = forcats::fct_rev(contrast),
                antibody = forcats::fct_rev(antibody))

chromHMM <- (barData %>% 
               ggplot(aes(x = antibody,
                          y = fold,
                          fill = antibody)) +
               geom_bar(stat = "identity",
                        color = "Black") +
               coord_flip() +
               labs(y = "Fold Enrichment",
                    x = element_blank()) +
               theme_classic() + 
               theme(axis.text = element_text(size = 20),
                     axis.title = element_text(size = 20),
                     strip.text = element_text(size = 20),
                     legend.position = "none") + 
               scale_fill_manual(values = c("#4B0082","#DCDCDC","#DCDCDC","#DCDCDC","#DCDCDC","#808080","#9750E3",
                                            "#8937DF","#CD5C5C","#BDB76B","#BDB76B","#FFFF80","#FFDF00","#AADF07",
                                            "#3F9A50","#008000","#FF4500","#FF0000")) +
               facet_grid(~contrast) +
               # scale_y_continuous(limits = c(-7,7),
               #                    breaks = c(-6,-4,-2,0,2,4,6)) +
               geom_hline(yintercept = 0) +
               geom_text(data = barData[(barData$signif == 1 & barData$fold > 0), ],
                         label = "*",
                         size = 8,
                         show.legend = FALSE,
                         nudge_y = 0.5,
                         nudge_x = -0.25) +
               geom_text(data = barData[(barData$signif == 1 & barData$fold < 0), ],
                         label = "*",
                         size = 8,
                         show.legend = FALSE,
                         nudge_y = -0.5,
                         nudge_x = -0.25)) %T>% 
  ggsave("chromHMM/Mouse development chromHMM enrichments forebrain summary.pdf",
         plot = .,
         width = 14,
         height = 5.5)

## DMRichments -------------------------------------------------------------

DMRichPlot <- function(data = data,
                       type = c("CpG", "genic")){
  
  stopifnot(type %in% c("CpG", "genic"))
  print(glue::glue("Plotting {type} annotation results"), "\n")
  
  data <- data %>%
    dplyr::filter(Type == !!type) %>% 
    dplyr::mutate(OR = dplyr::case_when(OR < 1 ~ -1/OR,
                                        OR >= 1 ~ OR)) %>%
    dplyr::mutate(signif = dplyr::case_when(q_value <= 0.05 ~ 1,
                                            q_value > 0.05 ~ 0)) 
  
  p <- ggplot(data = data,
              aes(x = Annotation,
                  y = OR,
                  fill = Annotation)) +
    geom_bar(stat = "identity", 
             color = "Black") +
    coord_flip() +
    labs(y = "Fold Enrichment",
         x = element_blank()) +
    theme_classic() + 
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          strip.text = element_text(size = 20),
          legend.position = "none") +
    scale_y_continuous(expand = c(0.1, 0.1)) + 
    scale_x_discrete(limits = data$Annotation %>%
                       forcats::as_factor() %>% 
                       levels() %>%
                       rev()) + 
    geom_hline(yintercept = 0) +
    geom_text(data = data[(data$signif == 1 & data$OR > 0), ],
              label = "*",
              size = 8,
              show.legend = FALSE,
              nudge_y = 0.5,
              nudge_x = -0.25) +
    geom_text(data = data[(data$signif == 1 & data$OR < 0), ],
              label = "*",
              size = 8,
              show.legend = FALSE,
              nudge_y = -0.5,
              nudge_x = -0.25)
  
  if(type == "CpG"){
    p <- p +
      scale_fill_manual(values = c("forestgreen", "goldenrod2", "dodgerblue", "blue3"),
                        breaks = data$Annotation %>%
                          forcats::as_factor() %>% 
                          levels(),
                        name = "Annotation")
  }else if(type == "genic"){
    p <- p +
      scale_fill_manual(values = data$Annotation %>%
                          forcats::as_factor() %>% 
                          nlevels() %>% 
                          wesanderson::wes_palette("Zissou1", n = ., type = "continuous") %>%
                          rev(),
                        breaks = data$Annotation %>%
                          forcats::as_factor(),
                        name = "Annotation")
  }
  
  p +
    facet_grid(~Contrast) %>% 
    return()
}

dir.create("DMRichments")

data <- tidyr::crossing(tissue = c("placenta", "brain"),
                        sex = c("female", "male"),
                        type = c("CpG", "genic")) %>% 
  purrr::pmap_dfr(function(tissue, sex, type){
    readxl::read_xlsx(glue::glue("{tissue}_{sex}/DMRichments/All DMRs_{type}_enrichments.xlsx")) %>%
      dplyr::mutate(Tissue = Hmisc::capitalize(tissue),
                    Sex = Hmisc::capitalize(sex),
                    Type = type) %>%
      dplyr::select(Tissue, Sex, Type, Annotation, OR, q_value = fdr)}) %T>%
  openxlsx::write.xlsx("DMRichments/DMRichments.xlsx") %>% 
  tidyr::unite("Contrast", c(Tissue, Sex), sep = " ") %>%
  dplyr::mutate(Contrast = factor(Contrast, levels = c("Placenta Female", "Placenta Male", "Brain Female", "Brain Male")))

CpG <- data %>%
  DMRichPlot(type = "CpG") 

Genic <- data %>%
  DMRichPlot(type = "genic") 

## Plot -------------------------------------------------------------------

cowplot::plot_grid(CpG +
                     scale_y_continuous(limits = c(-3.5,5.5)) +
                     theme(axis.title.x = element_blank(), 
                           axis.text.x = element_blank(), 
                           axis.ticks.x= element_blank()),
                   Genic +
                     scale_y_continuous(limits = c(-3.5,5.5)) +
                     theme(axis.title.x = element_blank(), 
                           axis.text.x = element_blank(), 
                           axis.ticks.x= element_blank(),
                           strip.background = element_blank(),
                           strip.text.x = element_blank()),
                   chromHMM +
                     scale_y_continuous(limits = c(-3.5,5.5)) +
                     theme(strip.text.x = element_blank()), 
                   align = "v",
                   nrow = 3,
                   rel_heights = c((4+1)/28, 6/28, 18/28)) %>%
  ggplot2::ggsave("DMRichments.pdf",
                  plot = ., 
                  width = 14,
                  height = 10)
