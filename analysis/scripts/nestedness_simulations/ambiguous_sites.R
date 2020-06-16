library(rjson)

# _______Filtering out ambiguous sites________ #

  # Produce a mapping from a child site to its level1 parent site
get_child_to_par <- function(nested_json){
  intermediate_list <- lapply(nested_json$Child_Map, unlist)
  child_to_par = c()
  its_names = c()
  for (parent in names(intermediate_list)){
    children = unlist(intermediate_list[parent])
    its_names <- c(its_names, children)
    child_to_par <- c(child_to_par,rep(parent, length(children)))
  }
  names(child_to_par) <- its_names
  return(child_to_par)
}

# Add an ambiguity classification
get_ambiguity <- function(qsimu_path, qsite_num, qerr_rate, classif, nesting){
  if (nesting != "nested" || classif != "FP") {return(FALSE)}
  parent_site_num = as.integer(child_to_par[as.character(qsite_num)])
  if (is.na(parent_site_num)) {return(FALSE)}
  up_one = parent_site_num
  while (! is.na(up_one)){
    up_one = as.integer(child_to_par[as.character(up_one)])
    if (! is.na(up_one)) parent_site_num = up_one
  }
  parent_classif = filter(nested_data, simu_path == qsimu_path & site_num == parent_site_num & err_rate == qerr_rate)$classif 
  if (parent_classif == "TP") {return(TRUE)}
  return(FALSE)
}

#____Run____#
nested_json <- fromJSON(file=argv$nested_json)

data <- data %>% mutate(classif = pmap_chr(list(res_has_call, res_is_correct), get_classif))
nested_data <- filter(data, nesting == "nested")
# Add ambiguity
child_to_par <- get_child_to_par(nested_json)
data <- data %>% mutate(ambiguous = pmap_chr(list(simu_path, site_num, err_rate, classif, nesting), get_ambiguity))
data_noambi = filter(data, ambiguous == FALSE)