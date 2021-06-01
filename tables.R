# Tables for top DMRs and all GO terms
# Ben Laufer

packages <- c("flextable", "officer", "tidyverse", "magrittr")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

setwd("/Users/blaufer/Box Sync/PEBBLES/DNA/DMRs")
dir.create("tables")

# Top 10 DMRs -------------------------------------------------------------

c("placenta_female", "placenta_male", "brain_female", "brain_male") %>% 
  purrr::set_names() %>%
  purrr::map_dfr(function(contrast){
    readxl::read_xlsx(glue::glue("{contrast}/DMRs/DMRs_annotated.xlsx"))},
    .id = "Contrast") %>% 
  tidyr::separate(Contrast, into = c("Tissue", "Sex")) %>% 
  dplyr::select(Tissue, Sex,
                Chr = seqnames,
                Region = annotation, Symbol = geneSymbol, Name = gene,
                Difference = difference, `p-value`) %T>%
  openxlsx::write.xlsx("tables/DMRs.xlsx") %>%
  #dplyr::mutate(Name = stringr::str_trunc(Name, 45)) %>% 
  dplyr::mutate(Tissue = Hmisc::capitalize(Tissue),
                Sex = Hmisc::capitalize(Sex)) %>% 
  dplyr::mutate(Tissue = forcats::as_factor(Tissue),
                Sex = forcats::as_factor(Sex)) %>% 
  dplyr::group_by(Tissue, Sex) %>%
  dplyr::slice(1:10) %>%
  dplyr::ungroup() %>% 
  #dplyr::select(-Tissue, -Sex) %>% 
  dplyr::mutate(Region = gsub(" \\(.*","", Region)) %>% 
  dplyr::mutate(Name = Hmisc::capitalize(.$Name)) %>% 
  dplyr::mutate(Difference = paste0(Difference, "%")) %>% 
  dplyr::mutate_at(vars(`p-value`),
                   function(x) {
                     formatC(x, digits = 1, format = "e", drop0trailing = TRUE)
                   }) %>% 
  dplyr::mutate(Chr = stringr::str_remove(Chr, "chr")) %>% 
  flextable::flextable() %>%
  flextable::theme_booktabs() %>%
  flextable::autofit() %>%
  flextable::align() %>%
  flextable::italic(j = c("Symbol", "Name"), italic = TRUE, part = "body") %>%
  flextable::hline(i = c(10,20,30), border = officer::fp_border(color = "black"), part = "body") %>%
  flextable::vline(border = officer::fp_border(color = "black")) %>% 
  flextable::padding(padding = 1, part = "all" ) %>% 
  flextable::fit_to_width(7.5) %>%
  flextable::save_as_docx(path = "tables/DMRs.docx")

# Gene Symbol Overlaps ----------------------------------------------------

c("placenta_female", "placenta_male", "brain_female", "brain_male") %>% 
  purrr::set_names() %>%
  purrr::map(function(contrast){
    readxl::read_xlsx(glue::glue("{contrast}/DMRs/DMRs_annotated.xlsx")) %>%
      dplyr::select(Chromosome = seqnames, Symbol = geneSymbol, Name = gene) %>%
      dplyr::distinct() %>%
      na.omit()}) %>% 
  purrr:::reduce(dplyr::semi_join) %>% 
  openxlsx::write.xlsx("tables/geneSymbolOverlapsAnnotated.xlsx")

# Genomic Coordinate Overlaps ---------------------------------------------

setwd("/Users/benlaufer/Box/PEBBLES/DNA/DMRs/overlaps/coordinate/")

readxl::read_xlsx(glue::glue("overlaps/coordinate/male_coordinate_overlaps_annotated.xlsx"),
                  col_types = c(rep("text", 7), rep(c("text","numeric"), 2))) %>%
  dplyr::mutate_if(is.numeric, function(x) {
    formatC(x, digits = 1, format = "e", drop0trailing = TRUE)
  }) %>% 
  dplyr::mutate("Placenta Difference" = paste0(`Placenta Difference`, "%"),
                "Brain Difference" = paste0(`Brain Difference`, "%")) %>%
  openxlsx::write.xlsx("overlaps/coordinate/male_coordinate_overlaps_annotated_reduced.xlsx")

# Manual fix in excel first
readxl::read_xlsx(glue::glue("overlaps/coordinate/female_coordinate_overlaps_annotated.xlsx"),
                  col_types = c(rep("text", 7), rep(c("text","numeric"), 2))) %>%
  dplyr::mutate_if(is.numeric, function(x) {
    formatC(x, digits = 1, format = "e", drop0trailing = TRUE)
  }) %>% 
  dplyr::mutate("Placenta Difference" = paste0(`Placenta Difference`, "%"),
                "Brain Difference" = paste0(`Brain Difference`, "%")) %>% 
  tibble::add_column("Sex" = "Female", .before = 1) %>%
  dplyr::bind_rows(readxl::read_xlsx(glue::glue("overlaps/coordinate/male_coordinate_overlaps_annotated_reduced.xlsx"),
                                     col_types = "text") %>%
                     tibble::add_column("Sex" = "Male", .before = 1)) %>% 
  dplyr::mutate(Name = stringr::str_trunc(Name, 35)) %>% 
  dplyr::mutate(Chr = stringr::str_remove(Chr, "chr")) %>% 
  dplyr::select(-Width) %>% 
  flextable::flextable() %>%
  flextable::add_header_row(top = TRUE, values = c("blank", "Placenta","Brain"), colwidths = c(7,2,2)) %>% 
  flextable::theme_booktabs() %>%
  flextable::autofit() %>%
  flextable::align() %>%
  flextable::italic(j = c("Symbol", "Name"), italic = TRUE, part = "body") %>%
  flextable::hline(i = 20, border = officer::fp_border(color = "black"), part = "body") %>%
  flextable::vline(border = officer::fp_border(color = "black")) %>% 
  flextable::padding(padding = 1, part = "all" ) %>% 
  flextable::fit_to_width(12) %>%
  flextable::save_as_docx(path = "tables/overlaps.docx")

# GO ----------------------------------------------------------------------

tidyr::crossing(source = c("placenta", "brain"),
                sex = c("female", "male")) %>% 
  purrr::pmap_dfr(function(source, sex){
    readxl::read_xlsx(glue::glue("{source}_{sex}/Ontologies/GOfuncR_slimmed_results.xlsx")) %>%
      dplyr::mutate(Source = source,
                    Sex = sex) %>%
      dplyr::select(Source, Sex, Term, `Gene Ontology`, `-log10.p-value`) %>% 
      dplyr::left_join(readxl::read_xlsx(glue::glue("{source}_{sex}/Ontologies/GOfuncR.xlsx")) %>%
                         dplyr::select(Term = node_name, FWER = FWER_overrep))}) %>%
  openxlsx::write.xlsx("slimmed_GOfuncR_results.xlsx")
