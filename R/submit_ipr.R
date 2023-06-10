## Internal and user-facing functions to submit FASTA files to InterProScan5 API

.set_env_ipr <- function() {
  env_ipr <<- new.env()
}

.validate_applications <- function(applications) {
  # use to verify applications input
  warning("not implemented yet")
  return()
}

.split_seqs <- function(fasta_path, outfolder) {
  cat("Splitting sequences\n", sep = "")

  split_folder <- file.path(outfolder, "split_seqs")
  if (!(dir.exists(split_folder))) {
    dir.create(split_folder)
  }

  # split the seqs
  n_seqs <- nrow(Biostrings::fasta.index(fasta_path))
  i <- 0L
  while (n_seqs > 0) {
    single_seq <- Biostrings::readAAStringSet(fasta_path, nrec = 1L, skip = i)
    seq_header <- names(single_seq)[[1]]

    # try to fix filenames
    # if there's a space, use as delimiter, and take first element
    if (grepl(" ", seq_header)) {
      seq_filename <- unlist(strsplit(seq_header, " "))[[1]]
    } else {
      # replace any alphanumeric characters otherwise
      seq_filename <- stringr::str_replace_all(seq_header, "[^[:alnum:]]", "")
    }
    # add a random suffix and '.faa' extension
    random_suffix <- paste0(
      sample(c(letters, LETTERS, 0:9), size = 7, replace = TRUE),
      collapse = ""
    )
    seq_filename <- paste0(seq_filename, "_", random_suffix, ".faa")
    cat("split_folder: ", split_folder, "seq_filename: ", seq_filename, "\n", sep = "")
    Biostrings::writeXStringSet(single_seq, file.path(split_folder, seq_filename))

    # subtract seqs and add skipping index for next read
    n_seqs <- n_seqs - 1L
    i <- i + 1L
  }
  print("Exiting .split_seqs()")
  return(split_folder)
}

.construct_status_table <- function() {
  ### global
  env_ipr$tb_status <<- tibble::new_tibble(
    list(
      input_file = c(NULL),
      job_id = c(NULL),
      status = c(NULL),
      polls = c(NULL),
      t_submit = c(NULL),
      t_latest_poll = c(NULL)
    )
  )
  return(env_ipr$tb_status)
}

.add_job_status_table <- function(input_filename, job_id, submit_status) {
  # get cur env for specifying job_id func param
  cur_env <- environment()
  new_job_row <- tibble::as_tibble_row(
    list(
      input_file = input_filename,
      job_id = get("job_id", cur_env),
      status = submit_status,
      polls = 0L,
      t_submit = as.POSIXct(Sys.time()),
      t_latest_poll = as.POSIXct(Sys.time())
    )
  )
  # coerce tb_status to prevent type errors on bind_rows (not sure why, but sometimes it is converted to a double)
  env_ipr$tb_status[["t_latest_poll"]] <<- as.POSIXct(env_ipr$tb_status[["t_latest_poll"]])

  env_ipr$tb_status <<- dplyr::bind_rows(env_ipr$tb_status, new_job_row)
  return(NULL)
  # return(job_code)
}

.submit <- function(path2seq, email) {
  fasta <- Biostrings::readAAStringSet(path2seq)
  sequence_str <- as.character(fasta[[1]])

  url_post <- "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run"
  headers_post <- c("Accept: text/plain", "Content-Type: application/x-www-form-urlencoded")
  data <- list(
    email = email,
    title = "riprscan",
    goterms = "false",
    pathways = "false",
    stype = "p",
    appl = paste0(env_ipr$.applications, collapse = ","),
    sequence = sequence_str
  )
  env_ipr$total_jobs <<- env_ipr$total_jobs + 1L
  res_post <- httr::POST(url_post, body = data, httr::add_headers(headers_post))
  job_id <- rawToChar(res_post$content)
  cat("Request info for ", job_id, ": \n", sep = "")
  print(res_post$request)

  status <- ifelse(as.integer(res_post$status) == 200L, "running", "failed")
  if (status == "failed") {
    message <- paste0(
      "POST error with sequence: ", names(fasta), "\n",
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


.poll_job <- function(job_id) {
  cat("Polling job: ", job_id, "\n", sep = "")
  url_get_base <- "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/"
  url_get <- paste0(url_get_base, job_id, "/tsv")
  header_get <- "Accept: text/tab-separated-values"
  res_get <- httr::GET(url_get, httr::add_headers(header_get))

  status <- as.integer(res_get$status_code)
  status <- ifelse(status == 200L, "completed", "running")

  # update N polls & latest poll time
  n_polls <- env_ipr$tb_status$polls[env_ipr$tb_status$job_id == job_id] 
  env_ipr$tb_status$polls[env_ipr$tb_status$job_id == job_id] <<- n_polls + 1
  env_ipr$tb_status$t_latest_poll[env_ipr$tb_status$job_id == job_id] <<- as.POSIXct(Sys.time())
  return(status)
}

.update_status_table <- function(job_id, status_update) {
  # change status for a given job_id row
  env_ipr$tb_status$status[env_ipr$tb_status$job_id == job_id] <<- status_update
  cat("env_ipr$tb_status from .update_status_table()\n")
  print(env_ipr$tb_status)
}

.lookup_input_file <- function(job_id) {
  # get cur env for specifying job_id func param
  cur_env <- environment()
  
  filename <- env_ipr$tb_status %>%
    dplyr::filter(job_id == get("job_id", cur_env)) %>%
    dplyr::pull(input_file)
  return(filename)
}

.lookup_polls <- function(job_id) {
  # get cur env for specifying job_id func param
  cur_env <- environment()
  
  cat("Looking up poll count for: ", job_id, "\n", sep = "")
  n_polls <- env_ipr$tb_status %>%
    dplyr::filter(job_id == get("job_id", cur_env)) %>%
    dplyr::pull(polls)
  cat("Polls for", job_id, ":", n_polls)
  return(n_polls)
}

.get_running_jobs <- function() {
  # handle case when tb_status is not yet init (or has NULL values)
  #   needed for iterative calls to .get_running_jobs
  n_running <- tryCatch(
    expr = {
      n_running <- env_ipr$tb_status %>%
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

.get_job_result <- function(job_id, outfolder, write_result = TRUE) {
  url_get_base <- "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/"
  url_get <- paste0(url_get_base, job_id, "/tsv")
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
  } else if (write_result == TRUE) {
    outfile <- unlist(strsplit((.lookup_input_file(job_id)), ".faa"))[[1]]
    outfile <- file.path(outfolder, "ipr_out", paste0(outfile, ".tsv"))
    cat(job_id, " complete\n", "Output file located at ", outfolder, "\n", sep = "")
    writeBin(res_get$content, outfile)
    return(tb_result)
  } else {
    return(tb_result)
  }
}

.join_output <- function(outfolder) {
  files <- list.files(file.path(outfolder, "ipr_out"), pattern = "\\.tsv$", full.names = TRUE)
  
  col_names <- c(
    "InputFile", "SeqMD5Digest", "SLength", "Analysis",
    "DB.ID", "SignDesc", "StartLoc", "StopLoc", "Score",
    "Status", "RunDate", "IPRAcc", "IPRDesc"
  )
  col_types <- c(
    "c", "c", "n", "c",
    "c", "c", "n", "n", "n",
    "c", "c", "c", "c"
  )
  # read all tables and bind the rows
  tb_joined <- dplyr::bind_rows(
    lapply(files, readr::read_tsv, col_names = col_names, col_types = col_types)
  )

  cat("Writing joined table to", outfolder, "\n")
  write.table(tb_joined, file = file.path(outfolder, "ipr_joined.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  return(tb_joined)
}

submit_ipr <- function(path2seq, outfolder, email, .applications = c("PfamA")) {
  #' @export
  #' @title Submit IPRscan analysis for protein sequences (multifasta input)
  #' @param path2seq path to your sequence file
  #' @param outfolder location for IPRscan outputs
  #' @param email required email ID for job submission (to track and prevent accidental misuse of the API)
  #' @param .applications a vector containing at least one of the annotation methods provided by InterProScan (https://www.ebi.ac.uk/Tools/services/rest/iprscan5/parameterdetails/appl)
  #' @keywords domains, domain architectures, protein characterization

  # init env_ipr for variables passed between multiple funcs
  .set_env_ipr()
  env_ipr$.applications <<- .applications
  env_ipr$total_jobs <<- 0L
  env_ipr$tb_status <<- .construct_status_table()
  
  outfolder <- file.path(outfolder)
  if (!(dir.exists(outfolder))) {
    dir.create(outfolder)
  }
  ipr_out <- file.path(outfolder, "ipr_out")
  if (!(dir.exists(ipr_out))) {
    dir.create(ipr_out)
  }
  file.copy(from = path2seq, to = file.path(outfolder, "input.fa"))
  
  n_seqs <- nrow(Biostrings::fasta.index(path2seq))

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

      cat("Submitting a batch of", n_submissions, "seqs\n")
      for (i in idx:stop_idx) {
        fasta <- file.path(split_seqs_folder, fasta_files[i])
        cat("Submitting sequence[", i, "]: ", fasta, " for analysis\n", sep = "")
        job_id <- .submit(fasta, email)
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
    for (job_id in job_ids) {
      status <- .poll_job(job_id)
      .update_status_table(job_id, status)
      # write/get result (if completed)
      if (status == "completed") {
        .get_job_result(job_id, outfolder, write_result = TRUE)
        .update_status_table(job_id, status)
        job_ids <- job_ids[job_ids != job_id] # rm job id from vec when completed
      } else {
        n_polls <- .lookup_polls(job_id)
        cat("env_ipr$tb_status:\n")
        print(env_ipr$tb_status)
        cat("N polls for ", job_id, ": ", n_polls, "\n", sep = "")
        # forget jobs that are polled more than N times
        if (n_polls > 500) {
          job_ids <- job_ids[job_ids != job_id]
          .update_status_table(job_id, status_update = "failed")
        }
      }
      Sys.sleep(3L)
    }
    cat("Number of sequences remaining: ", n_seqs, "\n", sep = "")
    cat("N running jobs: ", .get_running_jobs(), "\n", sep = "")
    cat("Total jobs: ", env_ipr$total_jobs, "\n", sep = "")
  }

  cat(
    "All jobs have been processed or forgotten\n",
    "Interproscan submission complete\n"
  )

  # coerce t_latest_poll column to POSIX datetime format
  env_ipr$tb_status[["t_latest_poll"]] <<- as.POSIXct(env_ipr$tb_status[["t_latest_poll"]])

  # write job polling metadata table
  cat("Writing the job_status table to", outfolder, "\n")
  env_ipr$tb_status <<- env_ipr$tb_status %>%
    dplyr::mutate(job_duration = paste0(as.character(as.integer(t_latest_poll - t_submit)), "s"))
  write.table(env_ipr$tb_status, file = file.path(outfolder, "job_status.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  # write and return the joined output
  cat("Joining results tables \n")
  tb_joined <- .join_output(outfolder)
  return(tb_joined)
}
