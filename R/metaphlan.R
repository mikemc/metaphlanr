# TODO: allow passing a list of file names and sample name (or perhaps a named
# list of file names)
parse_clade_profiles <- function (fn, combine_markers = TRUE) {
    ## Helper functions
    parse_clade <- function (li) {
        # li is a list of two lines from the metaphlan2 clade_profiles output
        ## Parse the two lines for this taxon.
        # The first line is a tab-separated list of the gene identifiers, with the
        # first position empty. The second line is a tab-separated list of the
        # clade name followed by the normalized gene abundance (number of reads per
        # bp * 1000)
        li.split <- li %>% 
            map(str_split, "\t", simplify = TRUE)
        clade <- li.split[[2]][1]
        genes <- li.split[[1]][-1]
        abundances = as.double(li.split[[2]][-1])
        tibble(Clade = clade,
            Gene = genes, Abundance = abundances)
    }
    tax_rank <- function (clade) {
        # param `clade` is a string of the form
        # 'k__...|p__...|c__...|o__...|f__...|g__...|s__...|t__...', up to the
        # taxonomic rank
        clade %>%
            str_extract("[k,p,c,o,f,g,s,t](?=__[\\w]+$)") %>%
            switchv(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", 
                f = "Family", g = "Genus", s = "Species", t = "Strain")
    }
    ## Main
    # Read the clade profiles
    li <- readLines(fn)
    # Drop the header line
    li <- li[-1]
    # Split into clades (pairs of rows)
    clade_lists <- split(li, (seq(li)+1) %/% 2)
    # Parse each clade into a tibble
    tb <- clade_lists %>%
        map(parse_clade) %>%
        bind_rows
    # Add the smallest taxonomic rank of the clade
    tb <- tb %>%
        mutate(Rank = tax_rank(Clade))
    # Get gene lengths
    tb <- tb %>%
        left_join(marker_lengths, 
            by = "Gene",
            suffix = c("", ".db")) %>%
        mutate(Reads = Abundance * Length / 1000)
    # Reads should all be integers unless something went wrong
    reads_are_ints <- tb %$% all(abs(Reads - round(Reads)) < 1e-4)
    if (reads_are_ints) {
        tb <- tb %>%
            mutate(Reads = as.integer(Reads))
    } else {
        warning("Reads may not all be integers")
    }
    # Combine markers.
    if (combine_markers) {
        tb <- tb %>%
            group_by(Clade, Rank) %>%
            summarize(Length = sum(Length), Reads = sum(Reads)) %>%
            # Weight = reads per kilobase
            mutate(Weight = (Reads / Length) * 1000) %>%
            ungroup()
    }
    # Return the tb with Clade as first column
    tb %>%
        select(Clade, Rank, everything())
}


#' Parse metaphlan clade strings to a taxonomy table
#' 
#' @param clade. Vector of clade strings
#' @param derep. Whether should only have one row per unique clade string
#' 
parse_taxonomy <- function (clade, derep = TRUE) {
    # The clade string has the form
    # "k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanosphaera|s__Methanosphaera_stadtmanae|t__GCF_000012545"
    # truncated at whatever the smallest rank is. We want to be able to parse
    # the clade string for any possible smallest rank, from Kingdom to Strain.
    rank_letters <- c("k", "p", "c", "o", "f", "g", "s", "t")
    tax_pattern  <- paste0("(?:", rank_letters, "__(\\w+))?") %>%
        paste(collapse = "\\|?")
    if (derep) {
        tax <- clade %>%
            unique %>%
            str_match(tax_pattern)
    } else {
        tax <- clade %>%
            str_match(tax_pattern)
    }
    colnames(tax) <- c("Clade", "Kingdom", "Phylum", "Class", "Order", "Family",
        "Genus", "Species", "Strain")
    tax %>% as_tibble
}

