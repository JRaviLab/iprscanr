.split_seqs <- function(fasta_path, outfolder) {

  cat("Splitting sequences\n", sep = "")

  tmp <- file.path(outfolder, "tmp")
  if (!(dir.exists(tmp))) {
    dir.create(tmp)
  }

  n_seqs <- nrow(Biostrings::fasta.index(fasta_path))
  i <- 0L
  while (n_seqs > 0) {
    #n_read <- ifelse(n_seqs < 30, n_seqs, 30L)
    single_seq <- Biostrings::readAAStringSet(fasta_path, nrec = 1L, skip = i)
    seq_header <- names(single_seq)[[1]]
    # try to fix filenames
    if (grepl(" ", seq_header)) {
      seq_filename <- unlist(strsplit(seq_header, " "))[[1]]
    }
    else {
      seq_filename <- stringr::str_replace_all(seq_header, "[^[:alnum:]]", "")
    }
    out <- file.path(tmp, seq_filename)
    Biostrings::writeXStringSet(single_seq, out)

    # subtract seqs and add skipping index for next read
    n_seqs <- n_seqs - 1L
    i <- i + 1L
  }
  # return the tmp outfolder path
  print('Exiting .split_seqs()')
  return(tmp)
}

.construct_status_table <- function() {
  tb_status <<- tibble::new_tibble(
    list(
      job_id = c(NULL),
      status = c(NULL),
      polls = c(NULL),
      t_submit = c(NULL)
    )
  )
  return(tb_status)
}

.add_job_status_table <- function(job_code, submit_status) {
  new_job_row <- tibble::as_tibble_row(
    list(
      job_id = job_code,
      status = submit_status,
      polls = 0L,
      t_submit = as.POSIXct(Sys.time())
    )
  )
  tb_status <<- dplyr::bind_rows(tb_status, new_job_row)
  #return(job_code)
}

.submit <- function(path2seq, user_email) {

  fasta <- Biostrings::readAAStringSet(path2seq)
  sequence_str <- as.character(fasta[[1]])

  url_post <- 'https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run'
  headers_post <- c("Accept: text/plain", "Content-Type: application/x-www-form-urlencoded")
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
  job_id <- rawToChar(res_post$content)
  cat("Request info for ", job_id, ": \n", sep = "")
  print(res_post$request)

  status <- as.integer(res_post$status)
  if (as.integer(res_post$status) != 200L) {
    message <- paste0("POST error with sequence: ", fasta_header, "\n",
                      "Status code: ", as.integer(res_post$status), "\n")
    warning(message)
    .add_job_status_table(job_id, status)
    return("failed")
  }

  return(job_id)
}


.poll_job <- function(job_code) {

  url_get_base <- "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/"
  url_get <- paste0(url_get_base, job_code, "/tsv")
  header_get <- "Accept: text/tab-separated-values"

  res_get <- httr::GET(url_get, httr::add_headers(header_get))
  status <- as.integer(res_get$status_code)

  status <- ifelse(status == 200L, "complete", "running")

  return(status)

}

.update_status_table <- function(job_code, status_update, n_polls) {

  # mutate at
  tb_status <<- tb_status %>%
    mutate(status = ifelse(job_id == job_code, status_update, status)) %>%
    mutate(polls = ifelse(job_id == job_code, n_polls, polls))
}

.get_job_result <- function(job_code) {
  url_get_base <- "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/"
  url_get <- paste0(url_get_base, job_code, "/tsv")
  header_get <- "Accept: text/tab-separated-values"

  res_get <- httr::GET(url_get, httr::add_headers(header_get))
  # parse and write tsv
  plain_text <- rawToChar(res_get$content)
  tb_result <- tibble::as_tibble(utils::read.table(text = plain_text, header = TRUE, sep = "\t"))
  #writeBin(res_get$content, paste0(outfile, ".tsv"))
  #cat(job_id, " complete\n", "Output file located at ", outfile, "\n", sep = "")
  return(tb_result)

}

.check_status_table <- function() {
  print(tb_status)
  n_running <- tb_status %>%
    filter(status == "running") %>%
    count()
  # parse the count value int
  n_running <- as.integer(n_running[["n"]])
  return(n_running)

}

interproscan <- function(path2seq, outfolder, user_email = "test@gmail.com") {

  outfolder <- file.path(outfolder)
  if (!(dir.exists(outfolder))) {
    dir.create(outfolder)
  }

  tb_status <<- .construct_status_table()

  n_seqs <- nrow(Biostrings::fasta.index(path2seq))

  ### single seq submit
  if (n_seqs == 1L) {
    cat("Single sequence detected\n")
    outfile <- file.path(outfolder, "iprout")

    # handle POST fail
    submission_outcome <- .submit(path2seq, user_email)
    # .submit returns either "failed" or a job_id
    job_id <- ifelse(submission_outcome != "failed", submission_outcome, NULL)
    if (submission_outcome == "failed") {
      message <- cat("Single submission failed\n")
      warning(message)
      return(FALSE)
    }
    #cat("Full submission complete. Results located at: ", outfolder, "\n", sep = "")
    status <- "submitted"
    # continue while status does not equal 'complete' or 'failed'

    n_polls <- 1L
    while (!(status == "complete" | status == "failed")) {
      status <- .poll_job(job_id)
      .update_status_table(job_id, status, n_polls)
      # quit after 30 tries
      if (n_polls > 30) {
        status <- "failed"
      }
      n_polls <- n_polls + 1
    }

    if (status == "complete") {
      tb_result <- .get_job_result(job_id)
      write.table(tb_result,
                  paste0(out_folder,"iprout.tsv"),
                  sep="\t",
                  quote=FALSE,
                  row.names=FALSE)
    }
    else {
      return(FALSE)
    }
    return(outfolder)
  }

  ### multi-fasta submit
    cat("Multifasta detected\n")
    split_seqs_folder <- .split_seqs(path2seq, outfolder)
    fastas <- list.files(split_seqs_folder, pattern = ".faa")

    n_seqs <- nrow(Biostrings::fasta.index(path2seq))
    idx_start <- 1
    idx_stop <- ifelse(n_seqs >= 30, 30, n_seqs)
    while (n_seqs > 0) {

      while((.check_status_table()) < 30) {
        for (i in idx_start:idx_stop) {
         .submit(fastas[i], user_email)
        }
      }
    }

}


.old_submit <- function(path2seq, user_email) {
  #outfile

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
    Sys.sleep(30L)
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
  #print(rawToChar(res_get$content))
  #writeBin(res_get$content, paste0(outfile, ".tsv"))
  #cat(job_id, " complete\n", "Output file located at ", outfile, "\n", sep = "")

  return(TRUE)
}

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

  seq_i <- 1L
  seq_n <- length(list.files(split_seqs_folder))
  n_successes <- 0L
  for (seq_cur in list.files(split_seqs_folder)) {
    seq_cur_path <- file.path(split_seqs_folder, seq_cur)
    outfile <- file.path(outfolder, strsplit(seq_cur, ".faa"[[1]]))
    success <- as.integer(.submit(seq_cur_path, outfile, user_email))
    cat("Seq # ", seq_i, "/", seq_n, " completed\n", sep = "")
    seq_i <- seq_i + 1L
  }

  cat("Full submission complete. Results located at: ", outfolder, "\n", sep = "")
  return(outfolder)
}
