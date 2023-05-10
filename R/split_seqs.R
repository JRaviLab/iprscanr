.split_seqs <- function(fasta_path, outfolder) {

  tmp <- file.path(outfolder, "tmp")
  if (!(dir.exists(tmp))) {
    dir.create(tmp)
  }

  n_seqs <- nrow(fasta.index(fasta_path))
  i <- 0L
  batch_ct <- 0L
  # batch seqs into 30 to follow api reqs
  while (n_seqs > 0) {
    n_read <- ifelse(n_seqs < 30, n_seqs, 30L)
    batch <- readAAStringSet(fasta_path, nrec = n_read, skip = i)

    batch_name <- paste0("batch_", batch_ct, ".faa")
    out <- file.path(tmp, batch_name)
    writeXStringSet(batch, out)

    # substract seqs and add skipping index for next batch
    n_seqs <- n_seqs - n_read
    i <- i + n_read
    batch_ct <- batch_ct + 1L
  }
  # return the tmp outfolder path
  print('Exiting .split_seqs()')
  return(tmp)
}
