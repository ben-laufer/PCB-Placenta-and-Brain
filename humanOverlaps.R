# Test PEBBLES DMRs for overlaps with human array PCB studies and ASD WGBS studies by genomic coordinate 
# Load liftOver module too if on cluster (module load ucsc-liftover/377)
# Ben Laufer

.libPaths("/share/lasallelab/programs/DMRichR/R_3.6")
packages <- c("DMRichR", "magrittr", "Biostrings", "BSgenome.Hsapiens.UCSC.hg38.masked")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

setwd("/share/lasallelab/Ben/PEBBLES/DNA/DMRs/")
dir.create("humanOverlaps")

# Create and liftOver consensus DMRs --------------------------------------

loadDMRs <- function(name){
  load(glue::glue("{name}/RData/DMRs.RData"))
  assign(glue::glue("{name}_sigRegions"), sigRegions, envir = .GlobalEnv)
  assign(glue::glue("{name}_regions"), regions, envir = .GlobalEnv)
  rm(sigRegions, regions)
}

contrasts <- c("placenta_male", "placenta_female", "brain_male", "brain_female")

purrr::walk(contrasts, loadDMRs)

dir.create("liftOver")

getConsensus <- . %>%
  plyranges::as_granges() %>% 
  GenomicRanges::sort() %>%
  plyranges::reduce_ranges()

c(placenta_male_sigRegions, placenta_female_sigRegions) %>% 
  getConsensus %>%
  DMRichR::gr2bed("liftOver/placenta_consensus_sigRegions_mm10.bed")

c(brain_male_sigRegions, brain_female_sigRegions) %>% 
  getConsensus %>%
  DMRichR::gr2bed("liftOver/brain_consensus_sigRegions_mm10.bed")

c(placenta_male_regions, placenta_female_regions) %>% 
  getConsensus %>%
  DMRichR::gr2bed("liftOver/placenta_consensus_regions_mm10.bed")

c(brain_male_regions, brain_female_regions) %>% 
  getConsensus %>%
  DMRichR::gr2bed("liftOver/brain_consensus_regions_mm10.bed")

tidyr::crossing(contrasts = c("placenta_male", "placenta_female", "brain_male", "brain_female"),
                region = c("sigRegions", "regions")) %>%
  purrr::pwalk(function(contrasts, region){
    get(glue::glue("{contrasts}_{region}")) %>% 
      DMRichR::gr2bed(glue::glue("liftOver/{contrasts}_{region}_mm10.bed"))
  })

# R version of liftOver isn't the same algorithm and gives worse results for regions

system("rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToHg38.over.chain.gz ./liftOver")
system("gunzip liftOver/mm10ToHg38.over.chain.gz")

tidyr::crossing(tissue = c("placenta", "brain"),
                sex = c("female", "male", "consensus"),
                region = c("sigRegions", "regions"),
                genome = c("hg38")) %>% 
  purrr::pwalk(function(tissue, sex, region, genome){
    
    print(glue::glue("Lifting over {tissue}_{sex}_{region} to {genome}"))
    
    #/Users/blaufer/opt/miniconda3/bin/liftOver is path for desktop conda install
    
    system(glue::glue("liftOver \\
                       -minMatch=0.1 \\
                       liftOver/{tissue}_{sex}_{region}_mm10.bed \\
                       liftOver/mm10To{Hmisc::capitalize(genome)}.over.chain \\
                       liftOver/{tissue}_{sex}_{region}_{genome}.bed \\
                       liftOver/{tissue}_{sex}_{region}_unmapped_{genome}.txt"))
  })

system("rm liftOver/mm10ToHg38.over.chain")

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
                                  count.once = TRUE)
  
  tibble::tibble(contrast = !!contrast,
                 dataset = !!dataset,
                 overlaps = pt$numOverlaps$observed,
                 zscore = pt$numOverlaps$zscore,
                 pvalue = pt$numOverlaps$pval)
}

# PCB ---------------------------------------------------------------------

Michigan <- readr::read_csv("humanOverlaps/Curtis_Michigan_TableS1.csv") %>%
  dplyr::select(CpGs = CPG.Labels) %>% 
  dplyr::distinct() %>% 
  dplyr::filter(!is.na(CpGs)) %>% 
  DMRichR::arrayLift("EPIC")

Anniston <- readxl::read_xlsx("humanOverlaps/Pittman_Anniston_Tables.xlsx", skip = 2) %>%
  dplyr::select(CpGs = ProbeID) %>%
  dplyr::filter(stringr::str_detect(CpGs, "cg")) %>% 
  dplyr::distinct() %>% 
  dplyr::filter(!is.na(CpGs)) %>% 
  DMRichR::arrayLift("EPIC")

PCBresults <- tidyr::crossing(contrast = c("placenta_consensus", "brain_consensus"),
                           dataset = c("Michigan", "Anniston")) %>% 
  purrr::pmap_dfr(overlaps) %>% 
  dplyr::mutate(p.adjust = p.adjust(pvalue)) 

save(PCBresults, file = "humanOverlaps/PCBresults.RData")

# ASD ---------------------------------------------------------------------

loadDMRs <- function(name){
  load(glue::glue("{name}.RData"))
  assign(name, sigRegions, envir = .GlobalEnv)
  rm(sigRegions, regions)
}

contrasts <- c("Rett_brain_female", "Dup15_brain_male")
purrr::walk(contrasts, loadDMRs)

ASDresults <- tidyr::crossing(contrast = c("placenta_consensus", "brain_consensus"),
                           dataset = c("Rett_brain_female", "Dup15_brain_male",
                                       "ASD_placenta_male")) %>% 
  purrr::pmap_dfr(overlaps) %>% 
  dplyr::mutate(p.adjust = p.adjust(pvalue)) 

save(ASDresults, file = "humanOverlaps/ASDresults.RData")
