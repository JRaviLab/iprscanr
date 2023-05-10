library(Biostrings)
.split_seqs <- function(fasta_path, outfolder) {

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

#' @export
# To do: auto-detect multifasta param
ipr_submit <- function(path2seq, outfolder, outformat, email = "test@gmail.com", multifasta = FALSE, wait = TRUE) {

  if (!(dir.exists(outfolder))) {
    dir.create(outfolder)
  }

  pkg_ipr_py_path<- file.path("py_scripts", "iprscan5.py")
  script <- system.file(pkg_ipr_py_path, package = "rinterpro")
  #  if (!(exists(email))) {
  #    stop("Error: email is not set. Set email with rinterpro::set_email()")
  #  }
  pkg_py_path <- file.path(".venv", "bin", "python3")
  python <- system.file(pkg_py_path, package = "rinterpro")
  #script <- file.path(outfolder, "iprscan5.py")
  cmd <- paste(python, script,
               "--email", email,
               "--outformat", outformat,
               "--appl", "PfamA")

  ### single fasta
  if (multifasta == FALSE) {
    # single outfile
    outfile <- file.path(outfolder, "iprout")
    cmd <- paste(cmd, "--outfile", outfile, "--sequence", path2seq)

    system(cmd, wait = wait)

    # sleep to prevent recalls
    # ! memoise/cache and only call sleep if 30 have been submitted
    print("Job submitted. Sleeping for 30 seconds . . . ")
    Sys.sleep(30L)
  }
  ### multifasta
  else {
    print("multifasta mode")
    split_folder <- .split_seqs(path2seq, outfolder)
    print(split_folder)

    for (batch in list.files(split_folder)) {
      sequence <- file.path(split_folder, batch)
      outfile <- file.path(outfolder, strsplit(batch, ".faa"[[1]]))
      cmd <- paste(cmd, "--sequence", sequence, "--outfile", outfile, "--maxJobs", 30L)
      print(cmd)
      system(cmd, wait = wait)
      #print("Job submitted. Sleeping for 30 seconds . . . ")
      #Sys.sleep(30L)
    }
  }

}
