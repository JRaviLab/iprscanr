.split_seqs <- function(fasta_path, outfolder) {

  cat("Splitting sequences\n", sep = "")

  tmp <- file.path(outfolder, "tmp")
  if (!(dir.exists(tmp))) {
    dir.create(tmp)
  }

  n_seqs <- nrow(Biostrings::fasta.index(fasta_path))
  i <- 0L
  batch_ct <- 0L
  # batch seqs into 30 to follow api reqs
  while (n_seqs > 0) {
    n_read <- ifelse(n_seqs < 30, n_seqs, 30L)
    batch <- Biostrings::readAAStringSet(fasta_path, nrec = n_read, skip = i)

    batch_name <- paste0("batch_", batch_ct, ".faa")
    out <- file.path(tmp, batch_name)
    Biostrings::writeXStringSet(batch, out)

    # substract seqs and add skipping index for next batch
    n_seqs <- n_seqs - n_read
    i <- i + n_read
    batch_ct <- batch_ct + 1L
  }
  # return the tmp outfolder path
  print('Exiting .split_seqs()')
  return(tmp)
}

.submit <- function(path2seq, outfile, user_email) {

  ### POST
  url_post <- 'https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run'
  fasta <- Biostrings::readAAStringSet(path2seq)
  fasta_header <- names(fasta)[[1]]
  sequence_str <- as.character(fasta[[1]])
  print(fasta)
  headers_post <- c("Accept: text/plain", "Content-Type: application/x-www-form-urlencoded")

  # job params
  data <- list(
    email=user_email,
    title="test",
    goterms="false",
    pathways="false",
    stype="p",
    appl="PfamA",
    sequence=sequence_str
  )

  res_post <- httr::POST(url_post, body = data, httr::add_headers(headers_post))
  cat("Request: \n", sep = "")
  print(res_post$request)
  if (as.integer(res_post$status) != 200L) {
    message <- paste0("POST error with sequence: ", fasta_header, "\n",
                      "Status code: ", as.integer(res_post$status), "\n")
    warning(message)
    return(NULL)
  }
  else {
    cat("POST status: ", as.integer(res_post$status), "\n", sep = "")
  }
  job_id <- rawToChar(res_post$content)
  cat("job_id: ", job_id, "\n", sep = "")

  ### GET
  url_get_base <- "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/"
  url_get <- paste0(url_get_base, job_id, "/tsv")
  header_get <- "Accept: text/tab-separated-values"

  status <- -1L
  n_polls <- 1L
  # Try poll 30 times OR until status is complete
  while (status != 200L && n_polls < 30L) {
    res_get <- httr::GET(url_get, httr::add_headers(header_get))
    status <- as.integer(res_get$status_code)
    cat("GET poll #", n_polls, "\t", "Status: ", status, "\n", sep = "")
    n_polls <- n_polls + 1L
    Sys.sleep(60L)
  }
  if (status != 200L) {
    message <- paste0("Failed to retrieve job results for <",
                      job_id,
                      ">",
                      "\n",
                      "Job exceeded 30 minutes runtime.")
    warning(message)
    return(NULL)
  }
  # parse and write tsv
  print(rawToChar(res_get$content))
  writeBin(res_get$content, paste0(outfile, ".tsv"))
  cat(job_id, " complete\n", "Output file located at ", outfile, "\n", sep = "")

  return(TRUE)
}

#' @export
ipr_submit <- function(path2seq, outfolder, user_email = "test@gmail.com") {

  outfolder <- file.path(outfolder)
  if (!(dir.exists(outfolder))) {
    dir.create(outfolder)
  }

  # test if multifasta
  n_seqs <- nrow(Biostrings::fasta.index(path2seq))

  if (n_seqs == 1L) {
    cat("Single sequence detected", "\n", sep = "")
    outfile <- file.path(outfolder, "iprout")
    .submit(path2seq, outfile, user_email)
    cat("Full submission complete. Results located at: ", outfolder, "\n", sep = "")
    return(outfolder)
  }

  cat("Multiple sequences detected", "\n", sep = "")
  split_seqs_folder <- .split_seqs(path2seq, outfolder)

  batch_i <- 1L
  batch_n <- length(list.files(split_seqs_folder))
  n_successes <- 0L
  for (batch in list.files(split_seqs_folder)) {
    sequence <- file.path(split_seqs_folder, batch)
    outfile <- file.path(outfolder, strsplit(batch, ".faa"[[1]]))
    success <- as.integer(.submit(sequence, outfile, user_email))
    cat("Batch # ", batch_i, "/", batch_n, " completed\n", sep = "")
    batch_i <- batch_i + 1L
    Sys.sleep(30L)
  }

  cat("Full submission complete. Results located at: ", outfolder, "\n", sep = "")
  return(outfolder)
}
