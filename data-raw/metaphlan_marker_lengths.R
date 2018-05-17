## Metaphlan marker gene lengths
# Parse the gene ids and lengths from the Metaphlan marker info file
li <- readLines(file.path("~/downloads", "mpa_v20_m200_marker_info.txt"))
marker_lengths <- li %>%
    str_match("^(\\S+)\\t.*'len': ([0-9]+),.*$") %>%
    as_tibble %>%
    select(Gene = V2, Length = V3) %>%
    mutate(Length = as.integer(Length))
# Some gene identifiers have more than one location span, separated by a comma:
marker_lengths %>%
    filter(str_detect(Gene, ","))
# The gene identifiers in the clade_profiles don't have this form; they only
# have the first span (up to and excluding the comma). Ao we'll want to adjust
# accordingly by matching up to the first comma (or end of string if no commas)
marker_lengths <- marker_lengths %>%
    mutate(Gene_org = Gene,
        Gene = str_extract(Gene, "^.+?(?=,|$)"))
devtools::use_data(marker_lengths, internal = TRUE)
