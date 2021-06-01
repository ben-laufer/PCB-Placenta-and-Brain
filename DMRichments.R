# DMRichments plot
# Plot CpG, gene region, and chromHMM results together for all pairwise contrasts
# Ben Laufer

packages <-  c("tidyverse", "magrittr", "cowplot")

stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

setwd("/Users/blaufer/Box Sync/PEBBLES/DNA/DMRs")

# chromHMM ----------------------------------------------------------------

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

## Bar plot ---------------------------------------------------------------

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

# DMRichments -------------------------------------------------------------

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

## Run --------------------------------------------------------------------

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
  DMRichPlot(type = "CpG") %T>% 
  ggplot2::ggsave("DMRichments/CpG_enrichments.pdf",
                  plot = ., 
                  width = 12,
                  height = 2.125)

Genic <- data %>%
  DMRichPlot(type = "genic") %T>%
  ggplot2::ggsave(glue::glue("DMRichments/genic_enrichments.pdf"),
                  plot = ., 
                  width = 12,
                  height = 3)

# Combined plot -----------------------------------------------------------

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
  ggplot2::ggsave(glue::glue("DMRichments.pdf"),
                  plot = ., 
                  width = 14,
                  height = 10)
