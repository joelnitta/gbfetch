# Helper functions ----

# Extract name of selected rank from taxonomic classification dataframe
name_from_classification <- function (taxonomy, rank_select) {
  tibble::as_tibble(taxonomy) %>%
    dplyr::filter(rank == rank_select) %>%
    dplyr::pull(name)
}

#' Fetch metadata for a GenBank sequence
#'
#' Metadata is obtained in two steps, first with rentrez::entrez_summary(),
#' then with taxize::classification() to get the species name.
#'
#' To see other possible values to use for `col_select`, run
#' rentrez::entrez_summary() with `db = nucleotide` for a valid ID,
#' e.g. rentrez::entrez_summary(db = "nucleotide", id = "383212727")
#'
#' @param id Character vector of IDs for GenBank records.
#' Output of rentrez::entrez_search().
#' @param col_select Character vector of metadata columns to include
#' in results. Must include at least `taxid` so that species
#' can be included. Default values include:
#' \itemize{
#'   \item{gi}{Genbank GI number}
#'   \item{caption}{Genbank accession number}
#'   \item{taxid}{Taxon ID (can use to query with taxize)}
#'   \item{title}{Sequence title}
#'   \item{slen}{Sequence length}
#'   \item{subname}{Misc. data (specimen, collection country, etc), separated
#'   by |}
#'   \item{subtype}{Column names of misc. data, separated by |}
#' }
#' @param .pb Internal agument used for setting the progress bar; don't change this.
#'
#' @return Tibble
#'
#' @examples
#' \dontrun{
#' fetch_metadata_from_id("383212727")
#' fetch_metadata_from_id(c("383212727", "383212725"))
#' }
fetch_metadata_from_id <- function (
  id,
  col_select = c("gi", "caption", "taxid", "title", "slen", "subtype", "subname"),
  .pb = NULL) {

  # Check input
  assertthat::assert_that(is.character(id))
  assertthat::assert_that(is.character(col_select))
  assertthat::assert_that("taxid" %in% col_select)

  # Set progress bar when this function is looped
  if ((!is.null(.pb)) && inherits(.pb, "Progress") && (.pb$i < .pb$n)) .pb$tick()$print()

  # Fetch metadata and extract selected columns
  if (length(id) == 1) {
    rentrez_results <- rentrez::entrez_summary(db = "nucleotide", id = id) %>%
      magrittr::extract(col_select) %>%
      tibble::as_tibble()
  } else {
    rentrez_results <- rentrez::entrez_summary(db = "nucleotide", id = id) %>%
      purrr::map_dfr(magrittr::extract, col_select)
  }

  # The metadata from rentrez::entrez_summary() don't include proper taxonomy.
  # There is a taxid field though, which we can use to lookup taxonomy with
  # taxize.
  dplyr::mutate(
    rentrez_results,
    taxonomy = taxize::classification(taxid, db = "ncbi"),
    species = purrr::map_chr(
      taxonomy,
      name_from_classification,
      rank_select = "species")
  ) %>%
    dplyr::select(-taxonomy)

}

#' Fetch metadata from GenBank sequences in chunks
#'
#' fetch_metadata() queries the NCBI API online, so we should do this
#' in smaller chunks to avoid errors.
#'
#' @param id Character vector of IDs for GenBank records.
#' Output of rentrez::entrez_search().
#' @param col_select Vector of column names to return.
#' @param chunk_size Number of rows to use for each chunk.
#'
#' @return A list of two items, `results` and `error`.
#'
#' @examples
#' \dontrun{
#' fetch_metadata_chunked(c("383212727", "383212725"))
#' }
fetch_metadata_chunked <- function(
  id,
  col_select = c("gi", "caption", "taxid", "title", "slen", "subtype", "subname"),
  chunk_size = 100) {

  # Check input
  assertthat::assert_that(is.character(id))
  assertthat::assert_that(is.character(col_select))
  assertthat::assert_that("taxid" %in% col_select)

  # Split input vector into chunks
  n <- length(id)
  r <- rep(1:ceiling(n/chunk_size), each=chunk_size)[1:n]
  id_list <- split(id, r) %>% magrittr::set_names(NULL)

  # Set progress bar
  pb <- dplyr::progress_estimated(length(id_list))

  fetch_metadata_from_id_safely <- purrr::safely(fetch_metadata_from_id)

  # Loop over vector chunks
  purrr::map(
    id_list,
    fetch_metadata_from_id_safely,
    .pb = pb)

}

# Main function ----

#' Fetch metadata from GenBank
#'
#' Sequences downloaded from GenBank with \code{\link{fetch_sequences}} only include
#' the title and/or accession number. Use \code{fetch_metadata} to obtain other
#' useful metadata associated with the sequences.
#'
#'  \code{\link[rentrez]{entrez_search}} is used to obtain a vector of IDs from the
#' `query`, then \code{\link[rentrez]{entrez_summary}} is used to download metadata
#' from the IDs. However, \code{\link[rentrez]{entrez_summary}} will fail if too many
#' IDs are used as input (more than 200-300 or so). Therefore, \code{fetch_metadata}
#' splits the IDs into chunks (a list of vectors), and loops over the list.
#'
#' Sometimes errors are encountered during the loop due to the API rejecting
#' the request, internet connectivity, etc. To avoid this, the loop repeats
#' until it finishes or the number of repeats reaches `max_tries`, upon
#' which it quits with an error.
#'
#' @param query String used to query NCBI GenBank. For more about the NCBI
#' query format see
#' https://www.ncbi.nlm.nih.gov/books/NBK3837/#EntrezHelp.Entrez_Searching_Options
#' @param chunk_size Number of ids to use for each chunk. Changing this
#' doesn't tend to affect the results, but lower values have more accurate
#' progress bars.
#' @param max_tries Maximum number of times to attempt the loop.
#'
#' @return Tibble of metadata resulting from Genbank query.
#' Columns include: \describe{
#'   \item{gi}{Genbank GI number}
#'   \item{accession}{Genbank accession number}
#'   \item{taxid}{Taxon ID (can use to query with taxize)}
#'   \item{title}{Sequence title}
#'   \item{slen}{Sequence length}
#'   \item{subname}{Misc. data (specimen, collection country, etc), separated
#'   by |}
#'   \item{subtype}{Column names of misc. data, separated by |}
#'   \item{species}{Species name}
#' }
#' @examples
#' fetch_metadata("rbcl[Gene] AND Crepidomanes[ORGN]")
#' @export
fetch_metadata <- function(
  query,
  chunk_size = 10,
  max_tries = 10) {

  # Check input
  assertthat::assert_that(assertthat::is.string(query))
  assertthat::assert_that(assertthat::is.number(chunk_size))
  assertthat::assert_that(assertthat::is.number(max_tries))
  assertthat::assert_that(chunk_size < 200)

  # Do an initial search without downloading any IDs to see how many hits
  # we get.
  initial_genbank_results <- rentrez::entrez_search(
    db = "nucleotide",
    term = query,
    use_history = TRUE
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

  # Extract vector of genbank IDs
  working_id <- genbank_results$id

  # Instantiate empty tibble to store results
  results <- tibble::tibble()

  # Instatiate counter to keep track of attempts
  # to fetch metadata
  this_try <- 0

  while(length(working_id) > 0) {

    this_try <- this_try + 1

    print(glue::glue("Fetching metadata, attempt {this_try} of {max_tries}"))

    if(this_try > max_tries) {
      print(glue::glue("Exceeded maximum number of tries ({max_tries}), quitting."))
      return(results)
    }

    this_result <-
      fetch_metadata_chunked(
        id = working_id,
        chunk_size = chunk_size) %>%
      purrr::transpose() %>%
      magrittr::extract("result") %>%
      purrr::flatten() %>%
      dplyr::bind_rows()

    results <- dplyr::bind_rows(results, this_result)

    working_id <- setdiff(working_id, results$gi)
  }

  # Confusingly, the entrez list actually uses "caption" to mean
  # accession, so swap this before selecting final columns.
  if ("caption" %in% colnames(results))
    results <- dplyr::rename(results, accession = caption)

  results

}
