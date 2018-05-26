# TODO: allow passing a list of file names and sample name (or perhaps a named
# list of file names)

#' Parse clade profiles produced by metaphlan2.
#'
#' @param fn. Path to metaphlan `clade_profiles` output
#' @param combine_markers. Whether to combine markers within each clade.
#' @param stat_q. Equivalent to metaphlan `stat_q` option.
#'
#' @import dplyr
#' @import tidyr
#' @export
parse_clade_profiles <- function (fn, combine_markers = TRUE, stat_q = 0) {
    ## Helper functions
    parse_clade <- function (li, stat_q = 0) {
        # li is a list of two lines from the metaphlan2 clade_profiles output
        ## Parse the two lines for this taxon.
        # The first line is a tab-separated list of the gene identifiers, with the
        # first position empty. The second line is a tab-separated list of the
        # clade name followed by the normalized gene abundance (number of reads per
        # 1000 bp)
        li.split <- li %>% 
            purrr::map(stringr::str_split, "\t", simplify = TRUE)
        clade <- li.split[[2]][1]
        genes <- li.split[[1]][-1]
        abundances = as.double(li.split[[2]][-1])
        tb <- tibble(Clade = clade,
            Gene = genes, Abundance = abundances)
        # If stat_q is set, filter out markers in the tails of Abundance,
        # (hopefully) following the "robust average" method of metaphlan2
        if (stat_q > 0) {
            # tb <- tb %>%
            #     mutate(Percentile = ecdf(Abundance)(Abundance)) %>%
            #     filter(Percentile >= stat_q, Percentile <= 1 - stat_q) %>%
            #     select(-Percentile)
            quantile_range <- quantile(tb$Abundance, c(stat_q, 1 - stat_q))
            tb <- tb %>%
                filter(Abundance > quantile_range[1], 
                    Abundance < quantile_range[2])
        }
        tb
    }
    tax_rank <- function (clade) {
        # param `clade` is a string of the form
        # 'k__...|p__...|c__...|o__...|f__...|g__...|s__...|t__...', up to the
        # taxonomic rank
        clade %>%
            stringr::str_extract("[k,p,c,o,f,g,s,t](?=__[\\w]+$)") %>%
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
        purrr::map(parse_clade, stat_q = stat_q) %>%
        bind_rows
    # Get gene lengths and read counts
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
    # Combine markers
    if (combine_markers) {
        tb <- tb %>%
            group_by(Clade) %>%
            summarize(Length = sum(Length), Reads = sum(Reads)) %>%
            mutate(Abundance = (Reads / Length) * 1000)
    }
    # Add the smallest taxonomic rank of the clade
    tb <- tb %>%
        mutate(Rank = tax_rank(Clade))
    # Return the tb with Clade as first column; rename Abundance to Weight to
    # avoid any possible confusion with absolute or relative abundance.
    tb %>%
        select(Clade, Rank, everything()) %>%
        rename(Weight = Abundance)
}


#' Parse metaphlan clade strings to a taxonomy table
#' 
#' @param clade. Vector of clade strings
#' @param derep. Whether should only have one row per unique clade string
#'
#' @import dplyr
#' @import tidyr
#' @export
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
            stringr::str_match(tax_pattern)
    } else {
        tax <- clade %>%
            stringr::str_match(tax_pattern)
    }
    colnames(tax) <- c("Clade", "Kingdom", "Phylum", "Class", "Order", "Family",
        "Genus", "Species", "Strain")
    tax %>% as_tibble
}

