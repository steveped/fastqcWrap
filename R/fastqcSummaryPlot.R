#' Plot the basic PASS/FAIL information from a batch of FASTQC summaries
#'
#' @description Make a basic plot of all PASS/FAIL information rom a batch of FASTQC summaries
#'
#' @details This calls \code{\link{fastqcSummary}} and makes a simple plot for all samples
#'
#' \code{fastqcSummaryPlot} will plot the output from the function \code{fastqcSummary} using default parameters.
#' This function relies on ggplot2 using pre-defined parameters.
#' Only the parameters for the x-axis labels (i.e. the file names) can be tweaked via the ellipsis.
#'
#'
#' @param fqNames the filenames to extract the totals for.
#' @param qcDir the directory to look in for the FASTQC reports
#' @param main the title of the plot. Will be automatically generated if missing. Use \code{main = c()} to remove the plot title
#' @param xLabLen the number of characters used for the labels on the x axis. Can be used to shorten long filenames.
#' @param showGuide show the guide to PASS, FAIL, WARN. Defaults to FALSE
#' @param ... passed to \code{element_text()} for \code{axis.text.x} only
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
fastqcSummaryPlot <- function(fqNames, qcDir, main, xLabLen = 20, showGuide = FALSE, ...){

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

  if (missing(main)) main <- paste("FASTQC Summary\n", qcDir, sep="")

  ggplot(ggCounts, aes(x=File, y=category, fill=Status)) +
    geom_tile(colour="black") +
    scale_fill_manual(values=c(FAIL="red", PASS="green", WARN="yellow")) +
    theme(axis.text.x = element_text(...)) +
    labs(x="Source File", y="QC Category", title=main) +
    guides(fill = showGuide) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0))

}
