library(tidyverse)

combine_enrichment_files <- function(..., n){
  # input file you want to combine
  files <- c(...)
  # check csv file
  if(all(str_detect(string = files, pattern = ".csv"))){
    cat("Input file is csv file")
  }else{
    stop("Please Input csv file!")
  }
  
  # check files numbers
  number_of_file <- length(files)
  
  # read csv files
  out <- do.call(rbind, purrr::map(files, function(x){
    read_delim(x, col_names = T, delim = ",") %>%
      dplyr::select(ID, ONTOLOGY,Count, pvalue) %>%
      dplyr::mutate(Group = x)
  }))
  
  max_length = max(out$Count) + 50
  
  out %>%
    tidyr::pivot_wider(names_from = Group, values_from = c(Count, pvalue)) %>%
    na.omit() %>%
    purrr::set_names(c("ID", "ONTOLOGY", "Up_Counts", "Down_Counts", "Up_pvalue", "Down_pvalue")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(total_number = sum(c(Up_Counts, Down_Counts))) %>%
    dplyr::mutate(total_log = sum(c(-log10(Up_pvalue) + -log10(Down_pvalue)))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(total_number)) %>%
    dplyr::slice_head(n = n) %>%
    dplyr::arrange(ONTOLOGY, desc(total_number)) %>%
    dplyr::mutate(start = 0, 
                  end = max_length,
                  Up_percent = Up_Counts / total_number,
                  Down_percent = Down_Counts / total_number) %>% 
    as.data.frame() -> out_file
  
  rownames(out_file) <- out_file$ID
  
  return(out_file)
}

  

