#' Plot a histogram of InterProScan API runtimes submitted with iprscanR
#'
#' @param ipr_log_path path to the iprscan api log
#'
#' @return A histogram plot object
#' @export
#'
#' @examples
job_time_hist <- function (ipr_log_path) {
  # possibly rename the 'api_log' in submit_ipr to 'ipr_log'
  # api is ambiguous in terms of the logs actual data

  # read log tsv as tibble
  tb_ipr_log <- tibble::as_tibble(readr::read_tsv(ipr_log_path))

  # calculate and add a total time column
  tb_ipr_log <- tb_ipr_log %>% dplyr::mutate(t_total = (t_latest_poll - t_submit))

  hist_plot <- ggplot2::ggplot(tb_ipr_log, ggplot2::aes(x = t_total, fill = ..x..)) +
    ggplot2::scale_fill_gradient(low = "green", high = "red") +
    ggplot2::geom_histogram() +
    ggplot2::labs(title = "IprScan API runtimes histogram",
                  x = "Job Runtimes (minutes)", y = "Count", fill = "minutes")
  return(hist_plot)
}
