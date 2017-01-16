#' Make a boxplot of Read Qualities
#'
#' @description Look through the .zip file produced by FASTQC to obtain summary of read qualities
#'
#' @details Look in a .zip file, as produced by FASTQC, and extract any read quality summaries
#'
#' @param qcDir the directory to look in for the FASTQC reports
#' @param fqNames Optional vector of filenames to extract the totals for.
#' @param suffix The suffix to remove from the file names
#' @param merge Merge across all files. Defaults to TRUE.
#' @param mergeFun The averaging function to use at each position when merging files.
#' Defaults to \code{median}, but can be anything similar such as \code{min} or \code{mean}.
#'
#'
#' @return A boxplot using the ggplot2 syntax
#'
#'
#' @import ggplot2
#' @import dplyr
#' @import magrittr
#'
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#' @rdname plotReadQualities
#' @export
plotReadQualities <- function(qcDir, fqNames, suffix = "_fastqc.html",  merge = TRUE, mergeFun = median){

  stopifnot(file.exists(qcDir))
  allFiles <- list.files(qcDir, pattern = suffix)
  if(length(allFiles) == 0) stop("No fastqc files could be found")

  if(missing(fqNames)) {
    fqNames <- allFiles
    fqNames <- gsub(suffix, "", fqNames)
  }
  else {
    # If a file, or subset are selected, check they exist
    fqNames <- gsub(suffix, "", fqNames)
    checkFiles <- vapply(fqNames, function(fq){
      sum(grepl(fq, allFiles)) == 1
    },
    logical(1))
    if (!all(checkFiles)) stop(names(checkFiles), " was either not found, or had an ambiguous name")
  }

  # Look for zip files with the required information
  zipFiles <- list.files(qcDir,pattern = ".zip")
  zip <- dplyr::if_else(length(zipFiles) > 0, TRUE, FALSE)
  if (!zip) stop("This function is currently only implemented for zipped files")
  zipData <- dplyr::data_frame(fqNames = fqNames,
                        zipName = vapply(fqNames,
                                         grep,
                                         x = zipFiles,
                                         value = TRUE,
                                         character(1)))

  # Get the quantiles
  allQuant <- split(zipData, f = zipData$fqNames) %>%
    lapply(function(x){
      intDir <- gsub(".zip", "", x$zipName)
      datFile <- file.path(intDir, "fastqc_data.txt")
      allData <- readLines(uz <- unz(file.path(qcDir,x$zipName), datFile), -1)
      close(uz)
      startLine <- match(">>Per base sequence quality\tpass", allData) + 1
      endLine <- grep("END_MODULE", allData)
      endLine <- min(endLine[endLine > startLine])
      allData <- allData[(startLine + 1):(endLine - 1)]
      allData <- stringr::str_split(allData, pattern = "\t") %>%
        lapply(function(x){
          as.data.frame(as.list(x)) %>%
            set_names(c("Base", "Mean", "Median", "LQ", "UQ", "P10", "P90"))
        })
      dplyr::bind_rows(allData) %>%
        dplyr::mutate(Base = factor(Base, levels = unique(Base)))
    })

  # Extract the plot data
  plotData <- names(allQuant) %>%
    lapply(function(x){
      dplyr::mutate(allQuant[[x]], Sample = x)
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(Mean = as.numeric(Mean),
           Median = as.numeric(Median),
           LQ = as.numeric(LQ),
           UQ = as.numeric(UQ),
           P10 = as.numeric(P10),
           P90= as.numeric(P90))

  # Find how many bins for the x-axis & the maximum score
  xBins <- length(levels(plotData$Base))
  maxQ <- max(plotData$P90)

  # Make the plot
  if (merge){
    # message("Not yet implemented for merged libraries")
    # Find the read number then calculate the values weighted by library size
    # For now though, just take averages
    plotData %>%
      dplyr::group_by(Base) %>%
      dplyr::summarise(Median = mergeFun(Median),
                       P10 = mergeFun(P10),
                       P90 = mergeFun(P90),
                       UQ = mergeFun(UQ),
                       LQ = mergeFun(LQ)) %>%
      dplyr::mutate(Sample = "Merged Files") %>%
      ggplot(aes(x = as.integer(Base), y = Median)) +
      annotate("rect", xmin = 0, xmax = xBins+1, ymin = 30, ymax = maxQ, fill = rgb(0, 0.9, 0.6), alpha = 0.3) +
      annotate("rect", xmin = 0, xmax = xBins+1, ymin = 20, ymax = 30, fill = rgb(0.9, 0.9, 0.7), alpha = 0.5) +
      annotate("rect", xmin = 0, xmax = xBins+1, ymin = 0, ymax = 20, fill = rgb(0.8, 0.4, 0.5), alpha = 0.5) +
      geom_crossbar(aes(ymin = LQ, ymax = UQ), fill = "yellow") +
      geom_linerange(aes(ymin = P10, ymax = LQ)) +
      geom_linerange(aes(ymin = UQ, ymax = P90)) +
      scale_x_continuous(breaks = seq_along(levels(plotData$Base)),
                         labels = levels(plotData$Base),
                         expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, maxQ), expand = c(0,0)) +
      facet_wrap(~Sample) +
      xlab("Position in read (bp)") +
      ylab("Quality Sccores") +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 90,
                                       hjust = 1))
  }
  else{
    ggplot(plotData, aes(x = as.integer(Base), y = Median)) +
      annotate("rect", xmin = 0, xmax = xBins+1, ymin = 30, ymax = maxQ, fill = rgb(0, 0.9, 0.6), alpha = 0.3) +
      annotate("rect", xmin = 0, xmax = xBins+1, ymin = 20, ymax = 30, fill = rgb(0.9, 0.9, 0.7), alpha = 0.5) +
      annotate("rect", xmin = 0, xmax = xBins+1, ymin = 0, ymax = 20, fill = rgb(0.8, 0.4, 0.5), alpha = 0.5) +
      geom_crossbar(aes(ymin = LQ, ymax = UQ), fill = "yellow") +
      geom_linerange(aes(ymin = P10, ymax = LQ)) +
      geom_linerange(aes(ymin = UQ, ymax = P90)) +
      scale_x_continuous(breaks = seq_along(levels(plotData$Base)),
                         labels = levels(plotData$Base),
                         expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, maxQ), expand = c(0,0)) +
      xlab("Position in read (bp)") +
      ylab("Quality Sccores") +
      facet_wrap(~Sample) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 90,
                                       hjust = 1))
  }


}
# f <- list.files("/data/Mouni/E-MTAB-4634/trimmed/fastQC", pattern = "R2.+html")
# plotReadQualities("/data/Mouni/E-MTAB-4634/trimmed/fastQC", fqNames = f, merge = TRUE)
