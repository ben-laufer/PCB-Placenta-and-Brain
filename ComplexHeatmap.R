# ComplexHeatmaps
# Ben Laufer

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
