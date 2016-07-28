#' Plot the basic PASS/FAIL information from a batch of FASTQC summaries
#'
#' @description Make a basic plot of all PASS/FAIL information rom a batch of FASTQC summaries
#'
#' @details This calls \code{\link{fastqcSummary}} and makes a simple plot for all samples
#'
#' \code{fastqcSummaryPlot} will plot the output from the function \code{fastqcSummary} using default parameters.
#' This function returns a ggplot2 object which can then be operated on using the standard ggplot2 syntax
#'
#'
#' @param fqNames the filenames to extract the totals for.
#' @param qcDir the directory to look in for the FASTQC reports
#' @param xLabLen the number of characters used for the labels on the x axis. Can be used to shorten long filenames.
#'
#'
#' @import ggplot2
#' @import dplyr
#' @import reshape2
#'
#' @seealso fastqc
#'
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#'
#' @rdname fastqcSummaryPlot
#' @export
fastqcSummaryPlot <- function(fqNames, qcDir, xLabLen = 20){

  counts <- fastqcSummary(fqNames, qcDir)
  cats <- rownames(counts)
  counts <- dplyr::mutate(counts,
                          category=factor(cats, levels=cats[seq(length(cats), 1, by=-1)]))
  n <- ncol(counts)

  ggCounts <- reshape2::melt(counts,
                             id.vars="category",
                             measure.vars=colnames(counts)[-n],
                             variable.name = "File",
                             value.name = "Status") %>%
    dplyr::mutate(File = substr(File, 1, xLabLen))

  ggplot(ggCounts, aes(x=File, y=category, fill=Status)) +
    geom_tile(colour="black") +
    scale_fill_manual(values=c(FAIL="red", PASS="green", WARN="yellow")) +
    labs(x="Source File", y="QC Category") +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0))

}
