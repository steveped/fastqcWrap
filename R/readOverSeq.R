#' Extract the Over-represented Sequence Information
#'
#' @description Look through the .zip file produced by FASTQC to obtain the over-represented sequence information
#'
#' @details Look in a .zip file, as produced by FASTQC, and extract any information about over-represented sequences
#'
#' @param fqName the filename to extract the totals for.
#' This should be the name of a single fastq file.
#' @param qcDir the directory to look in for the FASTQC reports
#'
#' @return A \code{list} with components \code{$summary} and \code{$status}
#'
#' The component \code{summary} contains a \code{data.frame} with columns corresponding to the original report.
#'
#' The component \code{$status} returns the PASS/FAIL/WARN status
#'
#' @import magrittr
#' @import stringr
#'
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#' @rdname readOverSeq
#' @export
readOverSeq <- function(fqName, qcDir){

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

  # Find where the overrepresented sequences module starts & ends
  startMod <- grep("Overrepresented sequences", data)
  endMod <- grep("Adapter Content", data) - 1

  # Get the status.
  # If this is pass, there will be no data in the module, so a blank df will be formed
  status <- stringr::str_replace(data[startMod], ">>Overrepresented sequences\t", "")
  df <- data.frame(Sequence = c(), Count = c(), Percentage = c(), `Possible Source`= c())

  if (status != "pass") {

    nm <- data[startMod + 1] %>%
      stringr::str_split_fixed(pattern = "\t", 4) %>%
      stringr::str_replace("#", "")

    df <- data[seq(startMod + 2, endMod - 1)] %>%
      stringr::str_split_fixed(pattern = "\t", 4) %>%
      as.data.frame(stringsAsFactors = FALSE) %>%
      magrittr::set_names(nm)

  }

  list(summary = df,
       status = str_to_upper(status))

}
