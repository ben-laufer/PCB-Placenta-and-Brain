# Test PEBBLES DMRs for overlaps with human array PCB studies and ASD WGBS studies by genomic coordinate 
# Load liftOver module too if on cluster (module load ucsc-liftover/377)
# Ben Laufer

.libPaths("/share/lasallelab/programs/DMRichR/R_3.6")
packages <- c("DMRichR", "magrittr", "Biostrings", "BSgenome.Hsapiens.UCSC.hg38.masked")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

setwd("/share/lasallelab/Ben/PEBBLES/DNA/DMRs/")

# Create and liftOver consensus DMRs --------------------------------------

loadDMRs <- function(name){
  load(glue::glue("{name}/RData/DMRs.RData"))
  assign(name, sigRegions, envir = .GlobalEnv)
  rm(sigRegions, regions)
}

contrasts <- c("placenta_male", "placenta_female", "brain_male", "brain_female")

purrr::walk(contrasts, loadDMRs)

getConsensus <- . %>%
  plyranges::as_granges() %>% 
  GenomicRanges::sort() %>%
  plyranges::reduce_ranges()

c(placenta_male, placenta_female) %>% 
  getConsensus %>%
  DMRichR::gr2bed("placenta_consensus_mm10.bed")

c(brain_male, brain_female) %>% 
  getConsensus %>%
  DMRichR::gr2bed("brain_consensus_mm10.bed")

# R version of liftOver isn't the same algorithm and gives worse results for regions

system("rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToHg38.over.chain.gz .")
system("gunzip mm10ToHg38.over.chain.gz")

regions <- c("placenta_consensus", "brain_consensus")

genomes <- "hg38"

tidyr::crossing(regions,
                genomes) %>% 
  purrr::pwalk(function(regions, genomes){
    
    print(glue::glue("Lifting over {regions} to {genomes}"))
    
    #/Users/blaufer/opt/miniconda3/bin/liftOver is path for desktop conda install
    
    system(glue::glue("liftOver \\
                       -minMatch=0.1 \\
                       {regions}_mm10.bed \\
                       mm10To{Hmisc::capitalize(genomes)}.over.chain \\
                       {regions}_{genomes}.bed \\
                       {regions}_unmapped_{genomes}.txt"))
  })

system("rm mm10ToHg38.over.chain.gz")

loadBED <- function(name){
  sigRegions <- readr::read_tsv(paste0(name,"_hg38.bed"),
                                col_names = c("chr", "start", "end"))  %>% 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  assign(name, sigRegions, envir = .GlobalEnv)
  rm(sigRegions)
}

contrasts <- c("placenta_consensus", "brain_consensus")
purrr::walk(contrasts, loadBED)

# Test function -----------------------------------------------------------

overlaps <- function(contrast = contrast,
                     dataset = dataset){
  
  print(glue::glue("Testing for {contrast} within {dataset}"))
  
  pt <- regioneR::overlapPermTest(A = get(contrast), 
                                  B = get(dataset), 
                                  alternative = "greater",
                                  genome = "hg38",
                                  ntimes = 10000,
                                  count.once = TRUE,
                                  mc.set.seed = FALSE, 
                                  force.parallel = TRUE)
  
  tibble::tibble(contrast = !!contrast,
                 dataset = !!dataset,
                 overlaps = pt$numOverlaps$observed,
                 zscore = pt$numOverlaps$zscore,
                 pvalue = pt$numOverlaps$pval)
}

# PCB ---------------------------------------------------------------------

setwd("../PCB_overlaps")

Michigan <- readr::read_csv("Curtis_Michigan_TableS1.csv") %>%
  dplyr::select(CpGs = CPG.Labels) %>% 
  dplyr::distinct() %>% 
  dplyr::filter(!is.na(CpGs)) %>% 
  DMRichR::arrayLift("EPIC")

Anniston <- readxl::read_xlsx("Pittman_Anniston_Tables.xlsx", skip = 2) %>%
  dplyr::select(CpGs = ProbeID) %>%
  dplyr::filter(stringr::str_detect(CpGs, "cg")) %>% 
  dplyr::distinct() %>% 
  dplyr::filter(!is.na(CpGs)) %>% 
  DMRichR::arrayLift("EPIC")

PCBresults <- tidyr::crossing(contrast = c("placenta_consensus", "brain_consensus"),
                           dataset = c("Michigan", "Anniston")) %>% 
  purrr::pmap_dfr(overlaps) %>% 
  dplyr::mutate(p.adjust = p.adjust(pvalue)) 

save(PCBresults, file = "PCBresults.RData")

# ASD ---------------------------------------------------------------------

setwd("../ASD_overlaps")

loadDMRs <- function(name){
  load(glue::glue("{name}.RData"))
  assign(name, sigRegions, envir = .GlobalEnv)
  rm(sigRegions, regions)
}

contrasts <- c("Rett_brain_female", "Dup15_brain_male")
purrr::walk(contrasts, loadDMRs)

ASDresults <- tidyr::crossing(contrast = c("placenta_consensus", "brain_consensus"),
                           dataset = c("Rett_brain_female", "Dup15_brain_male")) %>% 
  purrr::pmap_dfr(overlaps) %>% 
  dplyr::mutate(p.adjust = p.adjust(pvalue)) 

save(ASDresults, file = "ASDresults.RData")
