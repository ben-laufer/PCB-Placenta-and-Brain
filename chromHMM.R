# Create LOLA collection of 18 state, 10 mark chromHMM model of mouse development
# Ben Laufer

# References:
# http://databio.org/regiondb
# https://www.nature.com/articles/s42003-021-01756-4#data-availability
# https://publications.wenglab.org/mouse_epigenomes/18state_10marks/

.libPaths("/share/lasallelab/programs/DMRichR/R_4.0")

# LOLA version with "two.sided" option for Fisher's exact test
BiocManager::install("ben-laufer/LOLA")

packages <- c("DMRichR", "LOLA", "simpleCache", "qvalue", "parallel", "tidyverse", "magrittr")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

cores <- 60

# Create database ---------------------------------------------------------

dir.create("/share/lasallelab/programs/LOLA/mm10/Development_ChromHMM/", recursive = TRUE)
setwd("/share/lasallelab/programs/LOLA/mm10/Development_ChromHMM/")

message("Downloading bed files")
system("wget -m https://publications.wenglab.org/mouse_epigenomes/18state_10marks/")
system("rm publications.wenglab.org/mouse_epigenomes/18state_10marks/index.html")
system("gunzip publications.wenglab.org/mouse_epigenomes/18state_10marks/*")

message("Splitting bed files")

parallel::mclapply(list.files("publications.wenglab.org/mouse_epigenomes/18state_10marks/",
                              full.names = TRUE),
                   function(file){
                     LOLA::splitFileIntoCollection(filename = file,
                                                   splitCol = 4,
                                                   filenamePrepend = file %>%
                                                     basename() %>%
                                                     gsub(pattern = ".bed", replacement = "") %>% 
                                                     paste0("_"),
                                                   collectionFolder = "regions")},
                   mc.cores = cores,
                   mc.silent = TRUE)

system("rm -r publications.wenglab.org/")

message("Creating indices")

tibble::tibble("collector" = "Ben Laufer",
               "date" = Sys.Date(),
               "source" = "van der Velde et al. (2021)",
               "description" = "18 state, 10 mark chromHMM model of mouse development") %>%
  readr::write_tsv("collection.txt")

list.files("regions") %>% 
  tibble::as_tibble_col(column_name = "filename") %>%
  dplyr::mutate(tissue = stringr::str_extract(filename, "^[[:alpha:]]+"),
                tissue = dplyr::case_when(tissue == "embryonic" ~ "embryonic_facial_prominence",
                                          tissue == "neural" ~ "neural_tube",
                                          TRUE ~ tissue),
                antibody = stringr::str_extract(filename, "_[^_]+$") %>%
                  stringr::str_remove("_") %>%
                  stringr::str_remove(".bed"),
                day = stringr::str_extract(filename, "[0-9]+"),
                day = dplyr::case_when(day %in% c(11:16) ~ paste0("GD",day,".5"),
                                       day == 0 ~ "PD0"),
                description = paste(antibody, tissue, day)) %>% 
  readr::write_tsv("index.txt")

message("Caching collection")

LOLA::loadRegionDB(dbLocation = "/share/lasallelab/programs/LOLA/mm10",
                   useCache = TRUE,
                   limit = NULL,
                   collections = "Development_ChromHMM")

message("Done creating collection")

# Run LOLA ----------------------------------------------------------------

regionDb <- LOLA::loadRegionDB(dbLocation = "/share/lasallelab/programs/LOLA/mm10",
                               useCache = TRUE,
                               limit = NULL,
                               collections = "Development_ChromHMM")

purrr::walk(contrasts <- c("placenta_male", "placenta_female", "brain_male", "brain_female"),
            function(contrast){
              
              setwd(glue::glue("/share/lasallelab/Ben/PEBBLES/DNA/DMRs/{contrast}"))
              
              load("RData/DMRs.RData")
              
              dmrList <- sigRegions %>% 
                DMRichR::dmrList()
              
              dir.create("Extra/LOLA")
              setwd("Extra/LOLA")
              
              purrr::walk(seq_along(dmrList),
                          function(x){
                            
                            dir.create(names(dmrList)[x])
                            setwd(names(dmrList)[x])
                            
                            dmrList[x] %>%
                              LOLA::runLOLA(userSets = .,
                                            userUniverse = regions,
                                            regionDB = regionDb,
                                            minOverlap = 1,
                                            cores = cores,
                                            redefineUserSets = FALSE,
                                            direction = "both") %T>%
                              LOLA::writeCombinedEnrichment(combinedResults = .,
                                                            outFolder = "Development_ChromHMM",
                                                            includeSplits = FALSE) 
                            
                            setwd("../")
                            
                          })
            })

rm(regionDb)
