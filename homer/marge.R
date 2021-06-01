# Tidy and plot results from multiple HOMER analyses using MARGE
# Ben Laufer

BiocManager::install('robertamezquita/marge', ref = 'master')
packages <-  c("marge", "tidyverse", "magrittr")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

setwd("/Users/blaufer/Box Sync/PEBBLES/DNA/DMRs")

dir.create("HOMER")

# Tidy and table ----------------------------------------------------------

terms <- tidyr::crossing(tissue = c("placenta", "brain"),
                sex = c("female", "male")) %>% 
  purrr::pmap_dfr(function(tissue, sex){
    marge::read_known_results(glue::glue("{tissue}_{sex}/Extra/HOMER/both")) %>%
      dplyr::mutate(tissue = Hmisc::capitalize(tissue),
                    sex = Hmisc::capitalize(sex)) %>%
      dplyr::select(tissue, sex,
                    motif_name, motif_family, experiment, accession, consensus,
                    "-log10.p-value" = log_p_value, fdr, tgt_num, tgt_pct, bgd_num, bgd_pct)
    }) %T>%
  openxlsx::write.xlsx("HOMER/homer_enrichments.xlsx") %>%
  tidyr::unite("Contrast", c(tissue, sex), sep = " ") %>% 
  dplyr::mutate(Motif = paste0(motif_name," (",motif_family,")"),
                Contrast = as.factor(Contrast)) %>%
  dplyr::select(Contrast, Motif, `-log10.p-value`) %>%
  dplyr::group_by(Contrast) %>%
  dplyr::slice(1:10) %>%
  dplyr::arrange(dplyr::desc(`-log10.p-value`)) %>% 
  dplyr::ungroup()

terms$Motif %>% 
  table() %>%
  dplyr::as_tibble() %>%
  dplyr::arrange(dplyr::desc(n))

# Plot --------------------------------------------------------------------

purrr::walk(levels(terms$Contrast),
            function(contrast){
              (terms %>%
                dplyr::filter(Contrast == !!contrast) %>% 
                ggplot2::ggplot(ggplot2::aes(x = reorder(Motif,`-log10.p-value`),
                                             y = `-log10.p-value`,
                                             fill = Contrast)) +
                ggplot2::geom_bar(stat = "identity",
                                  position = ggplot2::position_dodge(),
                                  color = "Black") +
                ggplot2::coord_flip() +
                ggplot2::scale_y_continuous(expand = c(0, 0)) +
                ggsci::scale_fill_simpsons() +
                ggplot2::labs(y = expression("-log"[10](p))) +
                ggplot2::theme_classic() +
                 ggplot2::theme(text = element_text(size = 40),
                                axis.title.y = element_blank(),
                                axis.title.x = element_text(size = 25),
                                legend.position = "none")) %>% 
                ggplot2::ggsave(glue::glue("HOMER/{contrast}_homer_plot.pdf"),
                                plot = .,
                                device = NULL,
                                height = 5,
                                width = 9)
            })
 