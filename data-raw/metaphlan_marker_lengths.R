## Metaphlan marker gene lengths

# Example line we need to parse:
# gi|483970126|ref|NZ_KB891629.1|:c6456-5752	{'ext': set(['GCF_000226995', 'GCF_000355695', 'GCF_000373585']), 'score': 3.0, 'clade': 's__Streptomyces_sp_KhCrAH_244', 'len': 705, 'taxon': 'k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Streptomycetaceae|g__Streptomyces|s__Streptomyces_sp_KhCrAH_244'}

# Parse the clade and gene ids and lengths from the Metaphlan marker info file
li <- readLines(file.path("~/downloads", "mpa_v20_m200_marker_info.txt"))
marker_lengths_full <- li %>%
    str_match("^(\\S+)\\t.*'len': ([0-9]+),.*'taxon': '(\\S+)'.*$") %>%
    as_tibble %>%
    select(Gene = V2, Length = V3, Clade = V4) %>%
    mutate(Length = as.integer(Length))
# A small fraction of gene identifiers have more than one location span,
# separated by a comma:
marker_lengths_full %>%
    filter(str_detect(Gene, ","))
# The gene identifiers in the clade_profiles may have this form, but sometimes
# get truncated at the first comma. To deal with this, we can add copies with
# the truncated gene ids for these markers.
marker_lengths_trunc <- marker_lengths_full %>%
    filter(str_detect(Gene, ",")) %>%
    mutate(Gene = str_extract(Gene, "^.+?(?=,|$)"))
marker_lengths_full <- bind_rows(marker_lengths_full, marker_lengths_trunc)
# Subset to just the Gene ids and lengths for use in parsing clade profiles
marker_lengths <- marker_lengths_full %>%
    select(Gene, Length)
devtools::use_data(marker_lengths, internal = TRUE, overwrite = TRUE)
