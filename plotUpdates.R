# Update GO plots and heatmaps from DMRichR for figures
# Ben Laufer

GOplot <- function(slimmedGO = slimmedGO){
  
  slimmedGO %>% 
    dplyr::group_by(`Gene Ontology`) %>%
    dplyr::slice(1:7) %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(Term = stringr::str_trim(.$Term)) %>%
    dplyr::mutate(Term = Hmisc::capitalize(.$Term)) %>%
    dplyr::mutate(Term = stringr::str_trunc(.$Term, 40, side = "right")) %>% 
    dplyr::mutate(Term = factor(.$Term, levels = unique(.$Term[order(forcats::fct_rev(.$`Gene Ontology`), .$`-log10.p-value`)]))) %>% 
    ggplot2::ggplot(aes(x = Term,
                        y = `-log10.p-value`,
                        fill = `Gene Ontology`,
                        group = `Gene Ontology`)) +
    ggplot2::geom_bar(stat = "identity",
                      position = position_dodge(),
                      color = "Black") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggsci::scale_fill_d3() +
    ggplot2::labs(y = expression("-log"[10](p))) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = element_text(size = 40),
                   axis.title.y = element_blank(),
                   axis.title.x = element_text(size = 25)) %>% 
    return()
}

library(org.Mm.eg.db)
annoDb <- "org.Mm.eg.db"
library(tidyverse)
setwd("/Users/benlaufer/Box/PEBBLES/DNA/DMRs")
dir.create("GOplots")

purrr::walk(c("brain_female", "brain_male", "placenta_female", "placenta_male"),
            function(x){
              
              (readxl::read_xlsx(glue::glue("{x}/Ontologies/GOfuncR_slimmed_results.xlsx")) %>% 
                 GOplot() + 
                 theme(legend.position = "none")) %>%
                ggplot2::ggsave(glue::glue("GOplots/{x}_GOfuncR_plot.pdf"),
                                plot = .,
                                device = NULL,
                                height = 12,
                                width = 12)
              
            })

# Heatmaps ----------------------------------------------------------------

#' smoothPheatmap2
#' @description Plot a heatmap of normalized individual smoothed methylation value z scores for selected regions (i.e. significant DMRs)
#' @param bs.filtered.bsseq Smoothed \code{bsseq} object
#' @param sigRegions \code{GRanges} object of regions to plot a heatmap for
#' @param testCovariate The factor tested for differences between groups
#' @param filename Character specifying the name of the heatmap image file
#' @param ... Additional arguments passed onto \code{pheatmap()}
#' @return Saves a pdf image of the heatmap in the DMR folder
#' @import pheatmap pheatmap
#' @importFrom dplyr select_if
#' @importFrom RColorBrewer brewer.pal
#' @importFrom bsseq getMeth
#' @importClassesFrom bsseq BSseq 
#' @importMethodsFrom bsseq pData
#' @importFrom glue glue
#' @importFrom magrittr %>% set_colnames
#' @references \url{https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/}
#' @export smoothPheatmap2
#' 
smoothPheatmap2 <- function(bs.filtered.bsseq = bs.filtered.bsseq,
                            sigRegions = sigRegions,
                            testCovariate = testCovariate,
                            annotation_colors = annotation_colors,
                            filename = "DMRs.pdf",
                            ...){
  cat("\n[DMRichR] DMR heatmap \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  pdf(filename)
  plot <- bsseq::getMeth(BSseq = bs.filtered.bsseq,
                 regions = sigRegions,
                 type = "smooth",
                 what = "perRegion") %>% 
    as.matrix() %>%
    ComplexHeatmap::pheatmap(.,
                       scale = "row",
                       annotation_col = pData(bs.filtered.bsseq) %>%
                         as.data.frame() %>%
                         dplyr::select(Treatment),
                       color = RColorBrewer::brewer.pal(11, name = "RdBu") %>%
                         rev(),
                       show_colnames = F,
                       border_color = NA,
                       main = glue::glue("{scales::comma(length(sigRegions))} DMRs"),
                       fontsize = 16,
                       width = 11,
                       height = 8.5,
                       annotation_colors = annotation_colors,
                       row_km = 2,
                       column_km = 2,
                       border = TRUE,
                       use_raster = FALSE,
                       ...
    ) 
  
  draw(plot)
  
  dev.off()
}

.libPaths("/share/lasallelab/programs/DMRichR/R_4.0")

packages <- c("DMRichR", "ComplexHeatmap")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

setwd("/share/lasallelab/Ben/PEBBLES/DNA/DMRs/")
dir.create("heatmaps")
testCovariate <- "Treatment"

contrasts <- c("brain_female", "brain_male", "placenta_female", "placenta_male")

mclapply(contrasts, function(contrast){
  
  load(glue::glue("{contrast}/RData/bsseq.RData"))
  load(glue::glue("{contrast}/RData/DMRs.RData"))
  filename <- glue::glue("./heatmaps/{contrast}_heatmap.pdf")
  
  annotation_colors <- list(Treatment = c("PCB" = "#F8766D",
                                          "Control" = "#619CFF"))

  sigRegions %>%
    smoothPheatmap2(bs.filtered.bsseq = bs.filtered.bsseq,
                    testCovariate = testCovariate,
                    annotation_colors = annotation_colors,
                    filename = filename)
  
},
mc.cores = length(contrasts),
mc.silent = TRUE)
