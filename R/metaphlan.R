
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
    # Filter viruses and eukaryotes, since don't have marker length info
    if (filter) {
        tb <- tb %>% 
            filter(!str_detect(Clade, "k__Viruses|k__Eukaryota"))
    }
    # Return the tb with Clade as first column
    tb %>%
        select(Clade, Rank, everything())
}

# li is a list of two lines from the metaphlan2 clade_profiles output
parse_clade <- function (li, combine_markers = TRUE) {
    # Parse the two lines for this taxon.
    # Split each line at tabs; the first line starts with a blank and then gives
    # the gene identifiers, the second line gives the Taxon and the normalized gene
    # abundance (number of reads per bp * 1000)
    li.ss <- li %>% 
        map(str_split, "\t", simplify = TRUE)
    clade <- li.ss[[2]][1]
    genes <- li.ss[[1]][-1]
    tax_rank <- clade %>%
        str_extract("[a-z]__[\\w]+$") %>%
        str_sub(1, 1) %>%
        switchv(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", 
            f = "Family", g = "Genus", s = "Species", t = "Strain")
    abundances = as.double(li.ss[[2]][-1])
    # Parse the gene id to get the gene length in bp
    tb <- genes %>%
        str_match(".+\\|:c?([0-9]+)-([0-9]+)") %>%
        as_tibble %>%
        rename(Gene = V1, Start = V2, Stop = V3) %>%
        mutate(Start = as.integer(Start), 
            Stop = as.integer(Stop),
            Length = abs(Start - Stop) + 1)
    # Back-calculate the number of mapped reads from the normalized abundance. Add
    # taxonomy for combining data across taxa
    tb <- tb %>%
        mutate(Abundance = abundances,
            Reads = Abundance * Length / 1000,
            Clade = clade,
            Rank = tax_rank)
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
