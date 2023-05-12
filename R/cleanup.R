# Colnames from R/colnames_molevol.R
# Cleanup functions from R/cleanup.R, upstream_scripts/01.4_ipr2lin.R
# IPR2DA from upstream_scripts/05a_ipr2da.R

ipr_cln_colnames <- c("DB.ID","TaxID","AccNum.noV","AccNum",
                      "SeqMD5Digest","SLength","Analysis","SignDesc",
                      "StartLoc","StopLoc","Score","Status","RunDate",
                      "IPRAcc","IPRDesc","FullAccNum","ProteinName",
                      "Length","SourceDB","Completeness","Lineage",
                      "Species","Name","ShortName","LookupTblDesc",
                      "ID","Label")



.rename_ipr <- function(iprdata="path/to/iprout.tsv"){
#' @description
#' An internal function to rename IPR output
#' @param data raw IPRscan output data (tsv)

  ipr_colnames <- c("AccNum", "SeqMD5Digest", "SLength", "Analysis",
                  "DB.ID", "SignDesc", "StartLoc", "StopLoc", "Score",
                  "Status", "RunDate", "IPRAcc", "IPRDesc")
  colnames(iprdata) <- ipr_colnames
  return(iprdata)
}
