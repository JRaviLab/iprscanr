## Internal and user-facing functions to submit multifasta files to IPRscan

.split_seqs <- function(fasta_path, outfolder) {
  cat("Splitting sequences\n", sep = "")

  tmp <- file.path(outfolder, "tmp")
  if (!(dir.exists(tmp))) {
    dir.create(tmp)
  }

  # split the seqs
  n_seqs <- nrow(Biostrings::fasta.index(fasta_path))
  i <- 0L
  while (n_seqs > 0) {
    single_seq <- Biostrings::readAAStringSet(fasta_path, nrec = 1L, skip = i)
    seq_header <- names(single_seq)[[1]]

    # try to fix filenames
    if (grepl(" ", seq_header)) {
      seq_filename <- unlist(strsplit(seq_header, " "))[[1]]
    } else {
      # replace any alphanumeric characters otherwise
      seq_filename <- stringr::str_replace_all(seq_header, "[^[:alnum:]]", "")
    }
    # add a random suffix and '.faa' extension
    random_suffix <- paste0(sample(c(letters, LETTERS, 0:9),
                                   size = 7, replace = TRUE),
                            collapse = "")
    seq_filename <- paste0(seq_filename, "_", random_suffix, ".faa")
    cat("tmp: ", tmp, "seq_filename: ", seq_filename, "\n", sep = "")
    Biostrings::writeXStringSet(single_seq, file.path(tmp, seq_filename))

    # subtract seqs and add skipping index for next read
    n_seqs <- n_seqs - 1L
    i <- i + 1L
  }
  # return the tmp outfolder path
  print("Exiting .split_seqs()")
  return(tmp)
}

.construct_status_table <- function() {
  ### global
  tb_status <<- tibble::new_tibble(
    list(
      input_file = c(NULL),
      job_id = c(NULL),
      status = c(NULL),
      polls = c(NULL),
      t_submit = c(NULL),
      t_latest_poll = c(NULL)
    )
  )
  return(tb_status)
}

.add_job_status_table <- function(input_filename, job_code, submit_status) {
  new_job_row <- tibble::as_tibble_row(
    list(
      input_file = input_filename,
      job_id = job_code,
      status = submit_status,
      polls = 0L,
      t_submit = as.POSIXct(Sys.time()),
      t_latest_poll = as.POSIXct(Sys.time())
    )
  )
  # coerce tb_status to prevent type errors on bind_rows (not sure why, but sometimes it is converted to a double)
  tb_status[["t_latest_poll"]] <<- as.POSIXct(tb_status[["t_latest_poll"]])

  tb_status <<- dplyr::bind_rows(tb_status, new_job_row)
  return(NULL)
  # return(job_code)
}

.submit <- function(path2seq, user_email) {
  fasta <- Biostrings::readAAStringSet(path2seq)
  sequence_str <- as.character(fasta[[1]])

  url_post <- "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run"
  headers_post <- c("Accept: text/plain", "Content-Type: application/x-www-form-urlencoded")
  data <- list(
    email = user_email,
    title = "riprscan",
    goterms = "false",
    pathways = "false",
    stype = "p",
    appl = c("PfamA", "CDD"),
    sequence = sequence_str
  )
  ### global
  total_jobs <<- total_jobs + 1
  res_post <- httr::POST(url_post, body = data, httr::add_headers(headers_post))
  job_id <- rawToChar(res_post$content)
  cat("Request info for ", job_id, ": \n", sep = "")
  print(res_post$request)

  status <- ifelse(as.integer(res_post$status) == 200L, "running", "failed")
  if (status == "failed") {
    message <- paste0(
      "POST error with sequence: ", fasta_header, "\n",
      "Status code: ", as.integer(res_post$status), "\n"
    )
    warning(message)
    .add_job_status_table(basename(path2seq), job_id, status)
    return("failed")
  } else {
    .add_job_status_table(basename(path2seq), job_id, status)
    return(job_id)
  }
}


.poll_job <- function(job_code) {
  cat("Polling job: ", job_code, "\n", sep = "")
  url_get_base <- "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/"
  url_get <- paste0(url_get_base, job_code, "/tsv")
  header_get <- "Accept: text/tab-separated-values"
  res_get <- httr::GET(url_get, httr::add_headers(header_get))

  status <- as.integer(res_get$status_code)
  status <- ifelse(status == 200L, "completed", "running")

  # update N polls & latest poll time
  tb_status <<- tb_status %>%
    dplyr::mutate(polls = ifelse(job_id == job_code, polls + 1, polls)) %>%
    dplyr::mutate(t_latest_poll = ifelse(job_id == job_code,
                                         as.POSIXct(Sys.time()), t_latest_poll))

  return(status)
}

.update_status_table <- function(job_code, status_update) {
  # change status for a given job_id row
  tb_status <<- tb_status %>%
    dplyr::mutate(status = ifelse(job_id == job_code, status_update, status))
  cat("tb_status from .update_status_table()\n")
  print(tb_status)
}

.lookup_input_file <- function(job_code) {
  filename <- tb_status %>%
    dplyr::filter(job_id == job_code) %>%
    dplyr::pull(input_file)
  return(filename)
}

.lookup_polls <- function(job_code) {
  cat("Looking up poll count for: ", job_code, "\n", sep = "")
  poll_count <- tb_status %>%
    dplyr::filter(job_id == job_code) %>%
    dplyr::pull(polls)
  return(poll_count)
}

.get_running_jobs <- function() {
  print(tb_status)
  n_running <- tryCatch(
    expr = {
      n_running <- tb_status %>%
        dplyr::filter(status == "running") %>%
        dplyr::count()
    },
    error = function(e) {
      0
    }
  )
  # parse the count value int
  n_running <- tryCatch(
    expr = {
      as.integer(n_running[["n"]])
    },
    error = function(e) {
      0
    }
  )
  return(n_running)
}

.get_job_result <- function(job_code, outfolder, write_result = TRUE) {
  url_get_base <- "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/"
  url_get <- paste0(url_get_base, job_code, "/tsv")
  header_get <- "Accept: text/tab-separated-values"

  res_get <- httr::GET(url_get, httr::add_headers(header_get))
  # parse and write tsv
  plain_text <- rawToChar(res_get$content)
  print("Parsing job output text: ")
  print(plain_text)

  # return tibble or NULL if there was no annotation for sequence
  tb_result <- tryCatch(
    expr = {
      tb_result <- tibble::as_tibble(utils::read.table(text = plain_text, sep = "\t"))
    },
    error = function(e) {
      NULL
    }
  )
  if (is.null(tb_result)) {
    cat("No annotation for sequence, returning NULL\n")
    return(NULL)
  }
  else if (write_result == TRUE) {
    outfile <- unlist(strsplit((.lookup_input_file(job_code)), ".faa"))[[1]]
    outfile <- file.path(outfolder, paste0(outfile, ".tsv"))
    cat(job_code, " complete\n", "Output file located at ", outfolder, "\n", sep = "")
    writeBin(res_get$content, outfile)
    return(tb_result)
  }
  else {
    return(tb_result)
  }
}

.join_output <- function(outfolder) {
  files <- list.files(outfolder, pattern = ".tsv")
  # tb_init <- readr::read_tsv(file.path(outfolder, files[1]))
  # set colnames (adapted from from molevol_scipts/R/colnames_molevol.R)
  cols <- c(
    "InputFile", "SeqMD5Digest", "SLength", "Analysis",
    "DB.ID", "SignDesc", "StartLoc", "StopLoc", "Score",
    "Status", "RunDate", "IPRAcc", "IPRDesc"
  )
  col_classes <- c(
    "character", "character", "numeric", "character",
    "character", "character", "numeric", "numeric", "numeric",
    "character", "character", "character", "character"
  )

  # drop the first file used to build tb and iteratively bind the rest
  # files <- files[files != files[1]]
  tb_joined <- NULL
  idx_start <- 1L
  for (table in files) {
    input_file <- unlist(strsplit(table, ".tsv"))[[1]] # get input_file_name
    data <- read.table(file.path(outfolder, table),
                       sep = "\t", header = FALSE,
                       col.names = cols, colClasses = col_classes)
    tb_joined <- rbind(tb_joined, data)
    tb_joined[idx_start:nrow(tb_joined), "InputFile"] <- input_file
    idx_start <- nrow(tb_joined)
    cat("Joined", input_file, "to final table\n")
  }
  tb_joined <- tibble::as_tibble(tb_joined)

  cat("Writing joined table to", outfolder, "\n")
  write.table(tb_joined, file = file.path(outfolder, "ipr_joined.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  return(tb_joined)
}

submit_ipr <- function(path2seq, outfolder, user_email) {
  #' @export
  #' @title Submit IPRscan analysis for protein sequences (multifasta input)
  #' @param path2seq path to your sequence file
  #' @param outfolder location for IPRscan outputs
  #' @param email required email ID for job submission (to track and prevent accidental misuse of the API)
  #' @keywords domains, domain architectures, protein characterization

  total_jobs <<- 0
  outfolder <- file.path(outfolder)
  if (!(dir.exists(outfolder))) {
    dir.create(outfolder)
  }
  # copy input file to the output folder
  file.copy(from = path2seq, to = file.path(outfolder, "input.fa"))

  tb_status <<- .construct_status_table()

  n_seqs <- nrow(Biostrings::fasta.index(path2seq))

  ### single seq submit
  if (n_seqs == 1L) {
    cat("Single sequence detected\n")
    outfile <- file.path(outfolder, "iprout.tsv")

    # handle POST fail
    submission_outcome <- .submit(path2seq, user_email)
    # .submit returns either "failed" or a job_id
    job_id <- ifelse(submission_outcome != "failed", submission_outcome, NULL)
    if (submission_outcome == "failed") {
      message <- cat("Single submission failed\n")
      warning(message)
      return(NULL)
    }
    # cat("Full submission complete. Results located at: ", outfolder, "\n", sep = "")
    status <- "submitted"
    # continue while status does not equal 'complete' or 'failed'

    n_polls <- 1L
    while (!(status == "completed" | status == "failed")) {
      status <- .poll_job(job_id)
      .update_status_table(job_id, status)
      # quit after 180 tries; 15mins
      if (n_polls > 180) {
        status <- "failed"
      }
      n_polls <- n_polls + 1
      Sys.sleep(5)
    }

    if (status == "completed") {
      tb_result <- .get_job_result(job_id, outfolder, write_result = FALSE)
      cols <- c(
        "InputFile", "SeqMD5Digest", "SLength", "Analysis",
        "DB.ID", "SignDesc", "StartLoc", "StopLoc", "Score",
        "Status", "RunDate", "IPRAcc", "IPRDesc"
      )
      names(tb_result) <- cols
      .update_status_table(job_id, status)
      # write tb_results
      write.table(tb_result,
                  file.path(outfolder, "iprout.tsv"),
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE
      )
      # write tb_status
      write.table(tb_status,
                  file.path(outfolder, "api.log"),
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE
      )
      return(tb_result)
    }
    else { # if job failed
      .update_status_table(job_id, status)
      # write tb_status
      write.table(tb_status,
                  file.path(outfolder, "api.log"),
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE
      )
      return(NULL)
    }
  }
  else {
    ### multi-fasta submit
    cat("Multifasta detected\n")
    split_seqs_folder <- .split_seqs(path2seq, outfolder)
    fasta_files <- list.files(split_seqs_folder, pattern = ".faa")

    n_seqs <- nrow(Biostrings::fasta.index(path2seq))
    idx <- 1L
    batch_size_default <- 30L
    n_submitted <- 0L
    job_ids <- c() # job id vector
    while ((n_seqs > 0L) || ((.get_running_jobs()) > 0L)) {
      # only allow 30 job submissions at a time
      n_available <- 30L - .get_running_jobs() # .get_running_jobs() return N running jobs
      n_submissions <- min(n_seqs, batch_size_default, n_available) # minimum of available seqs for submission

      if (n_submissions > 0L) {
        # stop_idx has 2 conditions (1st batch submit && 2..N batch submit)
        stop_idx <- ifelse(n_submitted < 30, n_submissions, (idx + n_submissions - 1))

        cat("Submitting a batch of", (stop_idx - idx), "seqs\n")
        for (i in idx:stop_idx) {
          fasta <- file.path(split_seqs_folder, fasta_files[i])
          cat("Submitting the ", i, "th", "sequence: ", fasta, " for analysis\n", sep = "")
          job_id <- .submit(fasta, user_email)
          cat("Returned job id: ", job_id, "\n")
          if (job_id != "failed") { # .submit() will return either job_id or 'failed'
            job_ids <- append(job_ids, job_id)
          }
        }
        print("job_ids vec: ")
        print(job_ids)
        idx <- idx + n_submissions
        n_seqs <- n_seqs - n_submissions
        n_submitted <- n_submitted + n_submissions
      }

      # poll each job & update status table
      for (job in job_ids) {
        status <- .poll_job(job)
        .update_status_table(job, status)
        # write/get result (if completed)
        if (status == "completed") {
          .get_job_result(job, outfolder, write_result = TRUE)
          .update_status_table(job, status)
          job_ids <- job_ids[job_ids != job] # rm job id from vec when completed
        } else {
          n_polls <- .lookup_polls(job)
          cat("tb_status:\n")
          print(tb_status)
          cat("N polls for ", job, ": ", n_polls, "\n", sep = "")
          # forget jobs that are polled more than N times
          if (n_polls > 500) {
            job_ids <- job_ids[job_ids != job]
            .update_status_table(job_code = job, status_update = "failed")
          }
        }
      }
      cat("Number of sequences remaining: ", n_seqs, "\n", sep = "")
      cat("N running jobs: ", .get_running_jobs(), "\n", sep = "")
      cat("Total jobs: ", total_jobs, "\n", sep = "")
      Sys.sleep(5L)
    }
  }
  ### END job processing
  cat(
    "All jobs have been processed or forgotten\n",
    "Interproscan submission complete\n"
  )

  # coerce t_last_poll column to POSIX datetime format
  tb_status[["t_latest_poll"]] <<- as.POSIXct(tb_status[["t_latest_poll"]])

  # write api metadata table
  cat("Writing the api log\n")
  tb_status <<- tb_status %>%
    dplyr::mutate(job_duration = paste0(as.character(as.integer(t_latest_poll - t_submit)), "s"))
  write.table(tb_status, file = file.path(outfolder, "api.log"), sep = "\t", quote = FALSE, row.names = FALSE)

  # write and return the joined output
  cat("Joining results tables \n")
  tb_joined <- .join_output(outfolder)
  return(tb_joined)
}
