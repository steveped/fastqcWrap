#' @title Extract the PASS/FAIL information from a batch of FASTQC summaries
#'
#' @description Look through a batch of .zip files produced by FASTQC to obtain summary information across a range of files
#'
#' @details This will look in a number of .zip files, as produced by FASTQC,
#' and will extract the PASS/FAIL/WARN information across the key report sections.
#'
#' \code{fastqcSummary} will return the PASS/FAIL/WARN information for any number of FASTQC reports.
#'
#' \code{fastqcSummaryPlot} will plot the output from the function \code{fastqcSummary} using default parameters.
#' This function relies on ggplot2 using pre-defined parameters.
#' Only the parameters for the x-axis labels (i.e. the file names) can be tweaked via the ellipsis.
#'
#'
#' @param fqName the filename to extract the totals for.
#' This should be the name of a single fastq file.
#' @param fqNames the filenames to extract the totals for.
#' @param qcDir the directory to look in for the FASTQC reports
#' @param main the title of the plot. Will be automatically generated if missing. Use \code{main = c()} to remove the plot title
#' @param xLabLen the number of characters used for the labels on the x axis. Can be used to shorten long filenames.
#' @param showGuide show the guide to PASS, FAIL, WARN. Defaults to FALSE
#' @param ... passed to \code{element_text()} for \code{axis.text.x} only
#'
#' @return the PASS/FAIL summary information from the FASTQC report(s) either as a
#' data.frame (\code{fastqcSummary}) or as a plot (\code{fastqcSummaryPlot})
#'
#' @seealso fastqc
#'
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#' @rdname fastqcSummary
#' @export
fastqcSummary <- function(fqNames, qcDir){

  if (missing(qcDir)) stop("The directory must be specified!")
  if (!file.exists(qcDir)) stop("Cannot find the requested directory: ", qcDir)
  if(missing(fqNames)) {
    fqNames <- unique(gsub(".(zip|html)$", "", list.files(qcDir)))
  }
  else {
    allFiles <- list.files(qcDir, pattern=".zip")
    allMatches <- sapply(fqNames, FUN=grep, x=allFiles, value=TRUE)
    if (length(unlist(allMatches)) < length(fqNames)) stop("Some requested files couldn't be found")
    if (length(unlist(allMatches)) > length(fqNames)) stop ("Some requested names matched to multiple files")
  }

  counts <-  data.frame(lapply(fqNames, FUN=extractFastqcSummary, qcDir=qcDir))
  colnames(counts) <- gsub("_fastqc", "", fqNames)
  return(counts)

}

#'
extractFastqcSummary <- function(fqName, qcDir){

  fullName <- list.files(qcDir, fqName, full.names=TRUE)
  zipFile <- grep(".zip$", fullName, value=TRUE)
  if (length(zipFile)>1) stop("The given fqName matches more than one file: ", fqName)

  intDir <- gsub(".zip", "", basename(zipFile))
  datFile <- file.path(intDir, "summary.txt")
  data <- readLines(uz <- unz(zipFile, datFile), 12L)
  close(uz)

  out.vec <- unlist(strsplit(data, split="\t"))
  out <- out.vec[seq(1, length(out.vec), by=3)]
  names(out) <- out.vec[seq(2, length(out.vec), by=3)]

  return(out)

}
