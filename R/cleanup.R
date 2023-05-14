# Colnames from R/colnames_molevol.R
# Cleanup functions from R/cleanup.R, upstream_scripts/01.4_ipr2lin.R
# IPR2DA from upstream_scripts/05a_ipr2da.R

##################
## Column Names ##
##################
# to be used post-lookup
ipr_cln_colnames <- c("DB.ID","TaxID","AccNum.noV","AccNum",
                      "SeqMD5Digest","SLength","Analysis","SignDesc",
                      "StartLoc","StopLoc","Score","Status","RunDate",
                      "IPRAcc","IPRDesc","FullAccNum","ProteinName",
                      "Length","SourceDB","Completeness","Lineage",
                      "Species","Name","ShortName","LookupTblDesc",
                      "ID","Label")


.rename_ipr <- function(iprpath="path/to/iprout.tsv"){
  #' @description An internal function to rename IPR output
  #' @param data raw IPRscan output data (tsv)

  ipr_colnames <- c("AccNum", "SeqMD5Digest", "SLength", "Analysis",
                    "DB.ID", "SignDesc", "StartLoc", "StopLoc", "Score",
                    "Status", "RunDate", "IPRAcc", "IPRDesc")
  col_classes <- c(
    "character", "character", "numeric", "character",
    "character", "character", "numeric", "numeric", "numeric",
    "character", "character", "character", "character")

  iprdata <- read_tsv(iprpath, col_names = T)
  colnames(iprdata) <- ipr_colnames
  #readr::cols(iprdata) <- col_classes
  return(iprdata)
}


# ipr2da function
ipr2da <- function(iprpath="path/to/iprout.tsv",
                   analysis=c("Pfam"))
                              #"SMART","Phobius",
                              #"Gene3D", "TMHMM", "SignalP_GRAM_POSITIVE",
                              #"SUPERFAMILY", "MobiDBLite", "TIGRFAM",
                              #"PANTHER", "Coils"))
{
  # For pre-cleaned, combined ipr_joined.tsv
  ipr_in <- .rename_ipr(iprpath)

  prefix <- basename(iprpath) |> str_replace(".tsv", "")
  # split dataframe into unique proteins
  x <- split(x = ipr_in, f = ipr_in$AccNum)

  #future::plan(strategy = "multicore", .skip = T) | not suitable for RStudio

  # within each data.table
  domarch <- map(x, function(y) {
    #domarch <- future_map(x, function(y) {
    acc_row <- data.frame(AccNum = y$AccNum[1],  stringsAsFactors = F)
    DAs <- data.frame(matrix(nrow = 1, ncol = length(analysis) ))
    DA <- y |> group_by(Analysis) |> arrange(StartLoc)
    i = 1
    for(a in analysis) {
      a_da <- DA |> filter(Analysis == a)
      # if (a == "SignalP_EUK" || a == "SignalP_GRAM_NEGATIVE" || a == "SignalP_GRAM_POSITIVE") {
      #   var_shortname = "DB.ID" }
      # else {
      #   var_shortname = "ShortName" } # Only available post-lookup/cleanup
      var_shortname = "DB.ID"
      var_shortname_sym = sym(var_shortname)
      a_da <- a_da |>
        ungroup() |>
        select({{var_shortname_sym}}) |>
        filter(!is.na({{var_shortname_sym}})) |>
        filter(!is.null({{var_shortname_sym}})) |>
        pull(var_shortname) |>
        paste(collapse = "+")
      DAs[1,i] = a_da
      i=(i+1)
    }

    colnames(DAs) = paste("DomArch", analysis, sep = ".")
    return(cbind(acc_row, DAs))
  })

  # select relevant rows from IPRscan input to add to domarch
  # ipr_select <- ipr_in |>
    # # applicable when cleaned up | post acc2info
    # select(Name, AccNum, Species, TaxID, Lineage,
    #        Lineage_long_na, Lineage_long, Lineage_med, Lineage_short,
    #        ProteinName, SourceDB, Completeness, AccNum.noV) |>
    # distinct()

  # combine domarchs to one data frame, merge w/ acc2info
  domarch2 <- do.call(rbind.data.frame, domarch)

  # Add lineages & metadata | use when "ipr_select" is active
  # domarch2 <- domarch2 |>
  #   merge(ipr_select, by = 'AccNum', all.x = T)

  # save domarch_lins file
  write_tsv(domarch2, file = paste0(dirname(iprpath), "/",
                                    prefix, '.ipr_domarch.tsv'),
            append = F, na = 'NA')

  # return domarch2 dataframe to append to blast results if given
  return(domarch2)
}
