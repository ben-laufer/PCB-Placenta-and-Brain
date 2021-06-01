# Global methylation statistics
# Ben Laufer

# Packages ----------------------------------------------------------------

.libPaths("/share/lasallelab/programs/DMRichR/R_3.6")
suppressPackageStartupMessages(require("DMRichR", quietly = TRUE))

# Load --------------------------------------------------------------------

dir.create("/share/lasallelab/Ben/PEBBLES/DNA/DMRs/global")
setwd("/share/lasallelab/Ben/PEBBLES/DNA/DMRs/global")
system("cp ../cytosine_reports/*.gz ./")

# Global variables --------------------------------------------------------

testCovariate <- "Treatment"
adjustCovariate <- c("Sex","Tissue")
matchCovariate <- NULL
coverage <- 1
perGroup <- 0.75
genome <- "mm10"
cores <- 10

# Load and process samples ------------------------------------------------

cat("\n[DMRichR] Loading Bismark genome-wide cytosine reports \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
start_time <- Sys.time()

bs.filtered <- DMRichR::processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                                       meta = openxlsx::read.xlsx("../sample_info_master.xlsx", colNames = TRUE) %>%
                                         dplyr::mutate_if(is.character, as.factor),
                                       testCovariate = testCovariate,
                                       adjustCovariate = adjustCovariate,
                                       matchCovariate = NULL,
                                       coverage = 1,
                                       cores = cores,
                                       perGroup = perGroup
) 

glue::glue("Saving Rdata...")
save(bs.filtered, file = "PEBBLES_all_bismark.RData")

glue::glue("Individual smoothing timing...")
end_time <- Sys.time()
end_time - start_time

# Individual smoothed values ----------------------------------------------

cat("\n[DMRichR] Smoothing individual methylation values \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
start_time <- Sys.time()

bs.filtered.bsseq <- bsseq::BSmooth(bs.filtered,
                                    BPPARAM = MulticoreParam(workers = cores,
                                                             progressbar = TRUE))

bs.filtered.bsseq

glue::glue("Saving Rdata...")
save(bs.filtered.bsseq, file = "PEBBLES_all_bsseq.RData")

glue::glue("Individual smoothing timing...")
end_time <- Sys.time()
end_time - start_time

system("rm *.txt.gz")

# Global stats ------------------------------------------------------------

tidyr::crossing(tissue = c("placenta", "brain"),
                sex = c("male", "female")) %>% 
  purrr::pwalk(function(tissue, sex){
    
    GenomeInfoDb::seqlevelsStyle(bs.filtered.bsseq) <- "UCSC"
    # Remove sex chromosomes for analyses including both sexes and add Sex to the model
    # bs.filtered.bsseq <- GenomeInfoDb::dropSeqlevels(bs.filtered.bsseq,
    #                                                  c("chrX", "chrY", "chrM"),
    #                                                  pruning.mode = "coarse")
    
    clean <- . %>%
      dplyr::filter(Tissue == tissue & Sex == sex) %>%
      dplyr::mutate("Treatment" = dplyr::recode_factor(Treatment,
                                                       "PCB" = "PCB",
                                                       "Control" = "Control"))
    
    # Global methylation ------------------------------------------------------
    
    print(glue::glue("Testing for global methylation differences in {sex} {tissue}"))
    
    print(glue::glue("Obtaining smoothed global methylation levels from {sex} {tissue}"))
    
    global <- data.frame(DelayedMatrixStats::colMeans2(getMeth(BSseq = bs.filtered.bsseq,
                                                               type = "smooth",
                                                               what = "perBase"),
                                                       na.rm = TRUE))
    global$sample <- sampleNames(bs.filtered.bsseq)
    names(global) <- c("CpG_Avg", "sample")
    global <- dplyr::as_tibble(cbind(global, data.frame(pData(bs.filtered.bsseq))), rownames = NULL) %>%
      clean()
    
    print(glue::glue("Testing for differences in {sex} {tissue} using a linear model"))
    
    ANOVA <- global %>%
      aov(CpG_Avg ~ Treatment + Litter, data = .)
    
    list("globalInput" = global,
         "ANOVA" = broom::tidy(ANOVA))%>%
      openxlsx::write.xlsx(glue::glue("{tissue}_{sex}_globalStats_lm.xlsx"))
  })
