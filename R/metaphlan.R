
# TODO: allow passing a list of file names and sample name (or perhaps a named
# list of file names)
parse_clade_profiles <- function (fn, filter = TRUE, combine_markers = TRUE) {
    # Read the clade profiles
    li <- readLines(fn)
    # Drop the header line
    li <- li[-1]
    # Split into clades (pairs of rows)
    clades <- split(li, (seq(li)+1) %/% 2)
    # Parse each clade into a tibble
    tb <- clades %>%
        map(parse_clade, combine_markers = combine_markers) %>%
        bind_rows
    # Return the tb with Clade as first column
    tb %>%
        select(Clade, Rank, everything())
}

# li is a list of two lines from the metaphlan2 clade_profiles output
parse_clade <- function (li, combine_markers = TRUE) {
    ## Parse the two lines for this taxon.
    # The first line is a tab-separated list of the gene identifiers, with the
    # first position empty. The second line is a tab-separated list of the
    # clade name followed by the normalized gene abundance (number of reads per
    # bp * 1000)
    li.split <- li %>% 
        map(str_split, "\t", simplify = TRUE)
    clade <- li.split[[2]][1]
    genes <- li.split[[1]][-1]
    tax_rank <- clade %>%
        str_extract("[a-z]__[\\w]+$") %>%
        str_sub(1, 1) %>%
        switchv(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", 
            f = "Family", g = "Genus", s = "Species", t = "Strain")
    abundances = as.double(li.split[[2]][-1])
    ## Build a data frame of marker gene info for the clade
    # The gene lengths are stored in the marker_lengths df; the number of reads
    # mapped for the gene can be back-calculated from the normalized abundance.
    tb <- tibble(Clade = clade, Rank = tax_rank, 
        Gene = genes, Abundance = abundances) %>%
        left_join(marker_lengths, by = "Gene") %>%
        mutate(Reads = Abundance * Length / 1000)
    # Combine markers
    if (combine_markers) {
        tb <- tb %>%
            group_by(Clade, Rank) %>%
            summarize(Length = sum(Length), Reads = sum(Reads)) %>%
            mutate(Weight = (Reads / Length) * 1000) %>%
            ungroup()
    }
    tb
}
