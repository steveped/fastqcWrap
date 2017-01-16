#' Reproduce the sequence duplication plots from FastQC
#'
#' @description Look through the .zip file produced by FASTQC to obtain a summary of sequence duplication levels
#'
#' @details Look in a .zip file, as produced by FASTQC to obtain a summary of sequence duplication levels
#'
#' @param qcDir the directory to look in for the FASTQC reports
#' @param fqNames Optional vector of filenames to extract the totals for.
#' @param suffix The suffix to remove from the file names
#' @param merge Merge across all files. Defaults to TRUE.
#' @param mergeFun The averaging function to use at each position when merging files.
#' Defaults to \code{mean}, but can be anything similar such as \code{min} or \code{median}.
#'
#'
#' @return Plots which can be extended using the ggplot2 syntax
#'
#' @import ggplot2
#' @import dplyr
#' @import magrittr
#' @import reshape2
#'
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#' @rdname plotReadQualities
#' @export
plotDuplicationLevels <- function(qcDir, fqNames, suffix = "_fastqc.html",  merge = TRUE, mergeFun = mean){

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

  # Get the values
  allDup <- split(zipData, f = zipData$fqNames) %>%
    lapply(function(x){
      intDir <- gsub(".zip", "", x$zipName)
      datFile <- file.path(intDir, "fastqc_data.txt")
      allData <- readLines(uz <- unz(file.path(qcDir,x$zipName), datFile), -1)
      close(uz)
      startLine <- grep("Percentage of deduplicated", allData) + 1
      endLine <- grep("END_MODULE", allData)
      endLine <- min(endLine[endLine > startLine])
      allData <- allData[(startLine):(endLine - 1)]
      allData <- stringr::str_split(allData, pattern = "\t") %>%
        lapply(function(x){
          as.data.frame(as.list(x)) %>%
            set_names(c("Duplication Level", "Percent Deduplicated sequences", "Percent Total sequences"))
        })
      dplyr::bind_rows(allData) %>%
        dplyr::mutate(`Duplication Level` = factor(`Duplication Level`, levels = unique(`Duplication Level`)))
    })

  plotData <- names(allDup) %>%
    lapply(function(x){
      dplyr::mutate(allDup[[x]], Sample = x)
    }) %>%
    dplyr::bind_rows() %>%
    mutate(`Percent Deduplicated sequences` = as.numeric(`Percent Deduplicated sequences`),
           `Percent Total sequences` = as.numeric(`Percent Total sequences`))

  # Find how many bins for the x-axis
  xBins <- length(levels(plotData$`Duplication Level`))

  # Make the plot
  if (merge){
    # Find the read number then calculate the values weighted by library size
    # For now though, just take averages
    plotData %<>%
      dplyr::group_by(`Duplication Level`) %>%
      dplyr::summarise(`Percent Deduplicated sequences` = mergeFun(`Percent Deduplicated sequences`),
                       `Percent Total sequences` = mergeFun(`Percent Total sequences`)) %>%
      dplyr::mutate(Sample = "Merged Files")
  }

  plotData %>%
    reshape2::melt(id.vars = c("Duplication Level", "Sample"),
                   variable.name = "Type",
                   value.name = "Percent") %>%
    mutate(Type = gsub("Percent", "\\%", Type)) %>%
    ggplot(aes(x = as.integer(`Duplication Level`), y = Percent, colour = Type)) +
    geom_line() +
    scale_x_continuous(breaks = seq_along(levels(plotData$`Duplication Level`)),
                       labels = levels(plotData$`Duplication Level`),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
    scale_colour_manual(values = c("red", "blue")) +
    facet_wrap(~Sample) +
    xlab("Sequence Duplication Level") +
    ylab("") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.position = c(1, 1),
          legend.justification = c(1, 1),
          legend.box.background = element_blank())



}
# f <- list.files("/data/Mouni/E-MTAB-4634/trimmed/fastQC", pattern = "R2.+html")
# plotDuplicationLevels("/data/Mouni/E-MTAB-4634/trimmed/fastQC", fqNames = f, merge = FALSE)
