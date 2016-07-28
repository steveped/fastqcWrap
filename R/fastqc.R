#' A wrapper for the bash shell command fastqc.
#'
#' @description Takes fastq files then runs FASTQC on them
#'
#' @details This is a simple wrapper function for controlling & running \code{fastqc} from within R.
#' This can be very useful for controlling & documenting an entire pipeline from within knitr to produce a simple report
#'
#' The output components \code{$fqcFiles} & \code{$outDir} can be passed directly to the
#' functions \code{\link{fqcSummary}} and \code{\link{fqcSummaryPlot}} as the appropriate arguments.
#'
#' @param files The \code{fastq}, \code{sam} or \code{bam} files to perform \code{fastqc} on. Must be specified with the complete path
#' @param outDir the output directory for writing the trimmed fastq file.
#'  The names wil be automatically the same as the original file
#' @param threads the number of threads to use for the process.
#'  This will be auto-detected on a Linux system & will default to half of the found number of cores.
#' @param args any remaining arguments to the function \code{fastqc}.
#'  See the list available by typing \code{fastqc -h} into a terminal
#' @param exec The path to the executable. If not specified this will be automatically detected.
#'
#' @return A \code{list} with the components \code{$complete}, \code{fqcFiles}, \code{$outDir}, \code{version}, \code{command}.
#'
#' The component \code{$complete} is the complete message as output by \code{fastqc}, in order of completion
#' The component \code{$fqcFiles} is the edited message from above, sorted in lexicographical order.

#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#' @seealso \code{\link{cat}}, \code{\link{system}}
#'
#' @import parallel

#' @export
fastqc <- function(files, outDir, threads=floor(detectCores()/2), args=c(), exec){

  # Find the executable
  if (missing(exec)){
    exec <- gsub("fastqc: ", "", system2("whereis", "fastqc", stdout = TRUE))
    exec <- strsplit(exec, split = " ")[[1]][1] # Just choose the first location
    if (exec == "fastqc:") stop("Executable for 'fastqc' not found")
    message("Executable found: ", exec)
  }
  stopifnot(file.exists(exec)) # Protects against manual mis-specification
  v <- system(paste(exec, "-v"), intern=TRUE) # Get the version

  # Check the given set of files
  fileCheck <- vapply(files, FUN=file.exists, logical(1))
  if (sum(!fileCheck) >= 1) stop("\nCouldn't find files: \n", paste(files[!fileCheck], collapse="\n"))
  # Check for a common filetype
  fileCheck <- grep("(.fastq|fastq.gz|.fq|.fq.gz|.bam|.sam)$", basename(files), invert=TRUE)
  if (length(fileCheck) >=1 ) stop("\nFiles can only be fastq, sam or bam files: \n", paste(files[fileCheck], collapse="\n"))

  fType <- unique(gsub(".+(\\.[A-Za-z]+$)", "\\1", basename(files)))
  if (length(fType)>1) stop("\nThis command is currently only configured to support a single file type",
                            "\nYou have supplied files of type", fType)
  if (fType == ".gz") {
    fType <- unique(gsub(".+(\\.[A-Za-z]+)(\\.gz$)", "\\1\\2", basename(files)))
  }

  # Make any required output directories
  if (missing(outDir)) stop("\nThe variable outDir must be provided to tell the process where to send the output")
  if (!file.exists(outDir)) dir.create(outDir, recursive=TRUE)
  message("Files will be written to ", outDir)

  # Set the arguments
  o <- paste("-o", outDir)
  threads <- max(1, threads) # at least one core...
  message(exec, " will run using ", threads, " threads.")
  t <- paste("-t", threads)
  f <- paste(files, collapse=" ")
  cm <- paste(exec, o, t, args, f)
  message("Running the command:\n\t", cm, "\n")
  log <- system(cm, intern=TRUE)
  reports <- paste(gsub(fType, "", gsub("Analysis complete for ", "", log)), "_fastqc", sep="")
  return(list(complete=log, fqcFiles=reports[order(reports)], outDir=outDir, version=v, command=cm))

}

#f <- list.files("/home/steveped/Documents/Barry/iClip/iCLIP_fastq/01b_trimmed_fastq", full.names=TRUE, pattern="fastq")[1:2]
#outDir<-"/home/steveped/Documents/Barry/iClip/QC/01b_trimmed_fastq"
#test <- fastqc(f, outDir)
