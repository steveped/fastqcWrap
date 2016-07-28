#' Extract the Kmer information from fastqc reports
#'
#' @description Look through the .zip file produced by FASTQC to obtain the k-mer information
#'
#' @details This will look in a .zip file, as produced by FASTQC, and will extract information about over-represented k-mers
#' This is primarily designed to explore k-mers at the 5' end of reads.
#' Designed to return the k-mer information from a single report only.
#'
#' @param fqName the filename to extract the totals for.
#' This should be the name of a single fastq file.
#' @param qcDir the directory to look in for the FASTQC reports
#' @param maxWidth the maximum combined length of k-mers to retain, starting from the 5' end.
#' Any k-mers which extend beyond this value will be ignored
#'
#' @return A \code{list} with components \code{$summary} and \code{$status}
#'
#' The component \code{summary} contains a \code{data.frame} with columns corresponding to the original report.
#' It should be noted that beyond position 10 of the read, fastqc reports only give a range of bases.
#' These are denoted in the additional column maxShift,
#' indicating the k-mer may need to be shifted when generating a consensus sequence.
#' This should be assessed manually.
#'
#' The additional column Proportion indicates the proportion of reads containing the k-mer
#'
#' The component \code{$status} returns the PASS/FAIL/WARN status
#'
#' @import magrittr
#' @import dplyr
#' @import stringr
#'
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#' @rdname readKmer
#' @export
readKmer <- function(fqName, qcDir, maxWidth = 25){

  fullName <- list.files(qcDir, fqName, full.names = TRUE)
  stopifnot(file.exists(fullName))

  # get the zip file info
  zipFile <- grep(".zip$", fullName, value = TRUE)
  if (length(zipFile) > 1)
    stop("The given fqName matches more than one file: ",
         fqName)
  intDir <- gsub(".zip", "", basename(zipFile))
  datFile <- file.path(intDir, "fastqc_data.txt")
  # Read in the full summary
  data <- readLines(uz <- unz(zipFile, datFile), -1)
  close(uz)

  totReads <- grep("Total Sequences", data, value = TRUE) %>%
    stringr::str_replace("Total Sequences\t", "") %>%
    as.integer()

  # The lines in the file where the k-mer information is contained
  start <- grep("Kmer", data)
  rng <- seq(start + 2, length(data) - 1)

  nm <- as.vector(stringr::str_split_fixed(data[start + 1], pattern = "\t", 5))
  nm[1] <- "Sequence"
  nm[5] <- "Position"
  df <- stringr::str_split_fixed(data[rng], pattern = "\t", 5) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    set_names(nm)

  k <- unique(nchar(df$Sequence))
  minPos <- stringr::str_replace(df$Position, "-.+", "") %>% as.integer()
  maxPos <- stringr::str_replace(df$Position, ".+-", "") %>% as.integer()

  df$Count %<>% as.integer()
  df$PValue %<>% as.numeric()
  df$Position <- minPos
  df$maxShift <- maxPos - minPos
  df$Proportion <- df$Count / totReads
  df$Sample <- fqName
  df %<>% dplyr::arrange(Position, Count) %>% dplyr::filter(Position < maxWidth - k)

  list(summary = df,
       status = stringr::str_replace(data[start], ".+\t", "") %>% stringr::str_to_upper())

}
