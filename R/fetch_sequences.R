# Helper functions ----

# Function to download genbanks sequences and read in as fasta.
# fetch_seqs_from_id(rep(c("383212727", "383212725"),10), FALSE)
fetch_seqs_from_id <- function (id, simple_names = TRUE, .pb = NULL) {

  # Check input
  assertthat::assert_that(is.character(id))
  assertthat::assert_that(is.logical(simple_names))

  # Run progress bar when this function is looped
  if ((!is.null(.pb)) && inherits(.pb, "Progress") && (.pb$i < .pb$n)) .pb$tick()$print()

  # Grab sequences by ID with rentrez
  seqs_text <- rentrez::entrez_fetch(
    db="nucleotide",
    id = id,
    rettype="fasta")

  # Write out result as plain text file
  temp_fasta_file <- tempfile()
  write(seqs_text, temp_fasta_file)

  # Read back in seqs as ape::DNAbin
  seqs <- ape::read.FASTA(temp_fasta_file)
  unlink(temp_fasta_file)

  # Rename by accession only
  # (don't keep ".1" part of name).
  if(isTRUE(simple_names))
    names(seqs) <- names(seqs) %>%
    stringr::str_split(" ") %>% purrr::map_chr(1) %>%
    stringr::str_split("\\.") %>% purrr::map_chr(1)

  seqs
}

# Main function ----

#' Fetch DNA sequences from GenBank
#'
#' Query GenBank using the same format as searches on the
#' \href{https://www.ncbi.nlm.nih.gov/nucleotide/}{NCBI nucleotide database} and
#' download the sequences directly into R.
#'
#'  \code{\link[rentrez]{entrez_search}} is used to obtain a vector of IDs from the
#' `query`, then downloads the corresponding DNA sequences
#' from the IDs. However, \code{\link[rentrez]{entrez_search}} will fail if too many
#' IDs are used as input (more than 200-300 or so). Therefore, \code{fetch_sequences}
#' splits the IDs into chunks (a list of vectors), and loops over the list.
#'
#' @param query String used to query NCBI GenBank. For more about the NCBI
#' query format see
#' https://www.ncbi.nlm.nih.gov/books/NBK3837/#EntrezHelp.Entrez_Searching_Options
#' @param simple_names Logical; should the sequence names be simplified to
#' the GenBank accession number only?
#' @param chunk_size Number of ids to use for each chunk. Changing this
#' doesn't tend to affect the results, but lower values have more accurate
#' progress bars.
#' @param .pb Internal agument used for setting the progress bar; don't change this.
#' @param ... Additional arguments, not used by this function but meant for enabling
#' tracking if this function is used as part of a \code{\link[drake]{drake_plan}}.
#'
#' @return List
#'
#' @examples
#' \dontrun{
#' fetch_sequences("Crepidomanes minutum[ORGN] AND rbcl[Gene]")
#' }
#' @export
fetch_sequences <- function(query, simple_names = TRUE, chunk_size = 100, .pb = NULL, ...) {

  # Check input
  assertthat::assert_that(assertthat::is.string(query))
  assertthat::assert_that(is.logical(simple_names))
  assertthat::assert_that(assertthat::is.number(chunk_size))

  # Do an initial search without downloading any IDs to see how many hits
  # we get.
  initial_genbank_results <- rentrez::entrez_search(
    db = "nucleotide",
    term = query,
    use_history=TRUE
  )

  assertthat::assert_that(
    initial_genbank_results$count > 0,
    msg = "Query resulted in no hits")

  # Download IDs with maximum set to 1 more than the total number of hits.
  genbank_results <- rentrez::entrez_search(
    db = "nucleotide",
    term = query,
    use_history=FALSE,
    retmax = initial_genbank_results$count + 1
  )

  # Extract vector of genbank IDs, split into chunks
  id <- genbank_results$id
  n <- length(id)
  r <- rep(1:ceiling(n/chunk_size), each=chunk_size)[1:n]
  id_list <- split(id, r) %>% magrittr::set_names(NULL)

  # Set progress bar
  pb <- dplyr::progress_estimated(length(id_list))

  # Loop over chunked vector
  results <- purrr::map(id_list, fetch_seqs_from_id, simple_names = simple_names, .pb = pb)

  # Combine results into single DNA object
  jntools::flatten_DNA_list(results)

}
