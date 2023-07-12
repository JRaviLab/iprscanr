#' @export
#' @title InterProScan applications
#' @details All supported application values for InterProScan5
APPL <- c("NCBIfam", "SFLD", "Phobius", "SignalP", "SignalP_EUK", "SignalP_GRAM_POSITIVE", "SignalP_GRAM_NEGATIVE", "SuperFamily", "Panther", "Gene3d", "HAMAP", "PrositeProfiles", "PrositePatterns", "Coils", "SMART", "CDD", "PRINTS", "PfamA", "MobiDBLite", "PIRSF", "TMHMM", "AntiFam", "FunFam", "PIRSR")
#' @export
#' @title InterProScan5 sequence type values
#' @details InterProScan5 can accept Amino Acid (protein) or Nucleic Acid (DNA) sequences. The API values are 'p' and 'n', respectively
SEQTYPE <- c("n", "p")