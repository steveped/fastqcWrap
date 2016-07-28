#' Convert a vector of k-mers into a DNAStringSet
#'
#' @description Convert a vector of k-mers into a DNAStringSet
#'
#' @details Take a vector of k-mers and their repective positions within a set of reads, then shift accordingly.
#' The DNA_ALPHBET value N is added by default to enable easy integration with other functions in the package Biostrings.
#'
#' To find a possible consensus sequence, simply pass the output to Biobase::consensusString()
#'
#' @param seq a vector of k-mers
#' @param pos a vector of positions within a set of reads
#' @param pad the DNA base to be added when shifting k-mers
#' @param ... not used
#'
#' @return a DNAStringSet
#'
#' @import magrittr
#' @import Biostrings
#' @import stringr
#'
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#' @rdname kmer2StringSet
#' @export
kmer2StringSet <- function(seq, pos, ..., pad = "N"){

  k <- unique(nchar(seq))
  if (length(x)!= 1) stop("There appear to be k-mers of varying lengths. Lengths must be equal for sequences.\n")
  stopifnot(length(pos) == length(seq))

  stringSet <- str_pad(seq, width = pos - 1 + k, side = "left", pad = pad) %>%
    str_pad(width = max(pos) + k - 1, side = "right", pad = pad) %>%
    DNAStringSet()

  stringSet

}
