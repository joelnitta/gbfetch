#' Extract a gene sequence from a genbank file
#'
#' If duplicates are detected via `detect_dups` and `rename` is TRUE, duplicates
#' will be differentiated by appending a numeral to the name, e.g.
#' "accession-gene-1", "accession-gene-2", etc.
#'
#' See https://github.com/gschofl/biofiles
#'
#' @param gene_name Name of gene to extract. Must match genbank data exactly,
#' case-sensitive.
#' @param rec An instance of the gbRecord class.
#' Typically a (sub)genome containing the gene of interest.
#' @param rename Logical; should the sequence be renamed as "accession-gene"?
#' @param detect_dups Logical; should duplicates be detected? If TRUE, a text
#' warning will be returned in the case that this gene has duplicate copies
#' in the genome; if false, all copies will be returned.
#'
#' @return object of class DNAbin
#' @examples
#' \dontrun{
#' gb_file <- reutils::efetch(uid = "KY427346", db = "nuccore",
#'   rettype = "gbwithparts", retmode = "text")
#' # Parse GenBank record text file
#' this_rec <- biofiles::gbRecord(gb_file)
#' fetch_gene("rbcL", this_rec)
#' }
#'
fetch_gene <- function (gene_name, rec, rename = TRUE, detect_dups = TRUE) {

  # Check for reutils
  if (!requireNamespace("reutils", quietly = TRUE)) {
    stop("Package \"reutils\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # Check for biofiles
  if (!requireNamespace("biofiles", quietly = TRUE)) {
    stop("Package \"biofiles\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # Check input
  assertthat::assert_that(assertthat::is.string(gene_name))
  assertthat::assert_that(inherits(rec, what = "gbRecord"))

  pattern <- paste0("^", gene_name, "$")

  # Filter for the gene of interest
  filtered_rec <- biofiles::filter(rec, key = "gene", gene = pattern)

  # Extract the sequence
  seq <- biofiles::getSequence(biofiles::ft(filtered_rec))

  # If no gene detected, return NULL.
  if(is.null(seq)) return (NULL)

  # Currently we are only interested in one gene per query, so return
  # a text warning if this is not the case. May want to make more
  # general in the future (return multiple hits).
  if(length(seq) > 1 & isTRUE(detect_dups)) return ("Multiple genes detected")

  # Convert to ape::DNAbin
  seq <- ape::as.DNAbin(seq)

  # Optionally rename
  if (isTRUE(rename)) {
    if(length(seq) > 1) {
      names(seq) <- paste(biofiles::getAccession(rec), gene_name, 1:length(seq), sep = "-")
    } else {
      names(seq) <- paste(biofiles::getAccession(rec), gene_name, sep = "-")
    }
  }

  seq
}

#' Extract a list of gene sequences from a genbank file
#'
#' @param gene Characte vector; gene names to extract.
#' @param accession String; GenBank accession number of (full or partial) genome
#' containing the gene(s) of interest.
#' @param rename Logical; should the sequence be renamed as "accession-gene"?
#' @param detect_dups Logical; should duplicates be detected? If TRUE, a text
#' warning will be returned in the case that this gene has duplicate copies
#' in the genome instead of the sequence; if false, all copies will be returned,
#' differentiated by appending numbers to the gene name.
#'
#' @return A list
#'
#' @examples
#' \dontrun{
#' # KY427346 is the GenBank accession no. for the Diplazium striatum plastome
#' # https://www.ncbi.nlm.nih.gov/nuccore/KY427346
#'
#' genes_to_get <- c("rbcL", "matK", "psbA")
#'
#' fetch_gene_from_genome(genes_to_get, "KY427346")
#' }
#'
#' @export
fetch_gene_from_genome <- function (gene, accession, rename = TRUE, detect_dups = TRUE) {

  # Check for reutils
    if (!requireNamespace("reutils", quietly = TRUE)) {
      stop("Package \"reutils\" needed for this function to work. Please install it.",
           call. = FALSE)
    }

  # Check for biofiles
    if (!requireNamespace("biofiles", quietly = TRUE)) {
      stop("Package \"biofiles\" needed for this function to work. Please install it.",
           call. = FALSE)
    }

  # Check input
  assertthat::assert_that(is.character(gene))
  assertthat::assert_that(assertthat::is.string(accession))
  assertthat::assert_that(is.logical(rename))
  assertthat::assert_that(is.logical(detect_dups))

  # Download GenBank record as text
  # (do this before the loop so it isn't downloaded
  # each time for every gene).
  gb_file <- reutils::efetch(
    uid = accession, db = "nuccore",
    rettype = "gbwithparts", retmode = "text")

  # Parse GenBank record text file
  rec <- biofiles::gbRecord(gb_file)

  # Use purrr to "vectorize" fetch_gene, and return output named by gene.
  if(length(gene) > 1) {
    gene %>% magrittr::set_names(., gene) %>%
      purrr::map(fetch_gene,
                 rec = rec, rename = rename, detect_dups = detect_dups)
  } else {
    fetch_gene(gene, rec = rec, rename = rename, detect_dups = detect_dups)
  }
}

#' Assemble a list of genes from annotated genbank files
#'
#' Sometimes genes may be repeated within a genome, such as chloroplast
#' genes in the inverted repeat. If `drop_dups` is set to TRUE, these will
#' be excluded from the results (with a warning). This is useful for
#' assembling chloroplast gene matrices of single-copy genes.
#'
#' When running in parallel (`parallel` option is set to TRUE),
#' it may be necessary to set the parallel backend first using
#' \code{\link[future]{plan}}, or the code will still run sequentially.
#'
#' @param accessions Vector of genbank accession numbers of (partial)
#' genomes including the genes of interest
#' @param genes Vector of gene names to assemble
#' @param parallel Logical; should \code{\link[furrr]{future_map}} be used to
#' fetch genes in parallel?
#' @param drop_dups Logical; should genes with duplicate copies be excluded
#' from the results?
#'
#' @return List. Each item in the list is a gene, which contains a list of
#' sequences of class DNAbin.
#'
#' @examples
#' \dontrun{
#' # KP136830 is the GenBank accession no. for the Cystopteris protrusa plastome
#' # https://www.ncbi.nlm.nih.gov/nuccore/KP136830
#'
#' # KP136830 is the GenBank accession no. for the Diplazium striatum plastome
#' # https://www.ncbi.nlm.nih.gov/nuccore/KY427346
#'
#' # Assemble a list of DNA sequences for three genes from these two species.
#' # Note that psbA is duplicated since it is in the Inverted Repeat.
#'
#' assemble_gene_set(
#'   c("KP136830", "KY427346"),
#'   c("accD", "atpA", "psbA", "not_a_proper_gene_name")
#'   )
#'
#' assemble_gene_set(
#'   c("KP136830", "KY427346"),
#'   c("accD", "atpA", "psbA", "not_a_proper_gene_name"),
#'   drop_dups = FALSE
#'   )
#' }
#'
#' @export
assemble_gene_set <- function (accessions, genes,
                               parallel = FALSE, drop_dups = TRUE) {

  assertthat::assert_that(is.character(accessions))
  assertthat::assert_that(is.character(genes))
  assertthat::assert_that(is.logical(parallel))

  # Extract target genes from each plastome.
  # Wei et al 2017 data (83 coding genes from 40 plastomes)
  # took 23 min using furrr with four CPUs.
  # future::plan("multiprocess")
  if (isTRUE(parallel)) {
    all_genes <- furrr::future_map(
      accessions,
      ~ fetch_gene_from_genome(
        accession = ., gene = genes, detect_dups = drop_dups, rename = TRUE)
    ) %>%
      purrr::transpose()
  } else {
    all_genes <- purrr::map(
      accessions,
      ~ fetch_gene_from_genome(
        accession = ., gene = genes, detect_dups = drop_dups, rename = TRUE)
    ) %>%
      purrr::transpose()
  }

  # Some genes are dublicated because they are in the IR (Inverted Repeat).
  # Check for these.
  repeats_detected <- purrr::map_lgl(all_genes,
                                     ~ purrr::map(., class) %>%
                                       unique %>%
                                       stringr::str_detect("character") %>%
                                       any)

  if(isTRUE(drop_dups) & sum(repeats_detected) > 0) {
    print(paste(
      "These genes were duplicated within at least one accession and dropped:",
      genes[repeats_detected]
    ))
    all_genes <- all_genes[!repeats_detected]
  }

  # Prepare gene list for output
  all_genes %>%
    # Remove NULL values (when a gene wasn't found)
    purrr::map(purrr::compact) %>%
    purrr::keep(~ length(.x) > 0) %>%
    # Remove redundant names
    purrr::map(~ purrr::set_names(., NULL)) %>%
    # Convert each list of sequence lists within a gene
    # into a single list per gene
    purrr::map(jntools::flatten_DNA_list)

}
