library(dplyr)
library(tidyr)
library(stringr)

parse_gff_annotations <- function(gff_df, prefix = "q") {
  gff_df <- gff_df %>%
    mutate(
      ID = sub("^ID=", "", str_extract(attributes, "ID=[^;\\n]+")),
      parent = sub("^Parent=", "", str_extract(attributes, "Parent=[^;\\n]+")),
      IPR_ID = str_extract(attributes, "IPR[0-9]{6}")
    )

  # Subset for annotation-bearing rows
  subset <- gff_df %>%
    select(seqname = 1, ID, parent, start = 4, end = 5, IPR_ID) %>%
    filter(!is.na(IPR_ID))

  parents <- gff_df %>%
    select(ID, parent)

  list(subset = subset, parents = parents)
}

merge_dnds_with_annotations <- function(dnds_df, gff_info, id_col_prefix) {
  subset <- gff_info$subset
  parents <- gff_info$parents

  id_col <- paste0(id_col_prefix, "_id")

  dnds_df <- dnds_df %>%
    left_join(parents, by = setNames("ID", id_col)) %>%
    left_join(subset, by = setNames("parent", id_col)) %>%
    left_join(subset, by = setNames("ID", "parent"), suffix = c(".x", ".y")) %>%
    distinct(!!sym(id_col), .keep_all = TRUE)

  # Count NA values and drop the worse column
  na_x <- sum(is.na(dnds_df$IPR_ID.x))
  na_y <- sum(is.na(dnds_df$IPR_ID.y))

  if (na_x > na_y) {
    dnds_df <- dnds_df %>% select(-ends_with(".x"))
  } else if (na_y > na_x) {
    dnds_df <- dnds_df %>% select(-ends_with(".y"))
  }

  return(dnds_df)
}
