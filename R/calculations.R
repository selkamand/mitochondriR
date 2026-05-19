#' Parse idxstats file
#'
#' parse idxstats file.
#'
#' @param path to idxstats file produced by `samtools idxstats`
#'
#' @returns data.frame with columns "reference","length", "n_mapped" and "n_unmapped".
#' @export
#'
#' @examples
#' path = system.file(package="mitochondrir", "idxstats")
#' parse_idxstats(path)
parse_idxstats <- function(path){
  assertions::assert_file_exists(path)

  cols = c("reference","length", "n_mapped", "n_unmapped")

  # parse idxstats
  read.csv(
    file = path,
    sep = "\t",
    header = false,
    col.names = cols,
    colclasses =
      c("character","numeric", "numeric", "numeric")
  )
}

#' Calculate mtdna burden from idxstats file
#'
#' @param path Path to idxstats file produced by `samtools idxstats <bam>`
#' @param read_length Length of each read. Used to estimate depth but if not known, feel free to leave at default.
#' It will change the absolute values of depths of coverage but not affect the calculation of mitochondrial burden.
#' @param mtname name of mitochondrial genome in reference genome
#' @param autosomes name of autosomal sequences in reference genome
#' @param verbose verbose mode
#' @inheritParams calculate_mtdna_tumour
#'
#' @returns number of mitochondrial genomes per cell. if \code{per_cell=false} returns number of mitochondrial genomes per single copy of the autosome.
#' @export
#'
#' @examples
#' path = system.file(package="mitochondrir", "idxstats")
#' parse_idxstats(path)
mtdna_from_idxstats <- function(path, background_ploidy = 2, tumour_ploidy = 2, tumour_purity = 1, read_length = 150, mtname =  "chrm", autosomes = paste0("chr", 1:22), per_cell = TRUE, verbose = true){

  # assertions
  assertions::assert_string(mtname)
  assertions::assert_character(autosomes)

  # parse idxstats file
  df <- parse_idxstats(path)

  # check bam includes all expected mitochondrial and autosome names
  assertions::assert_includes(df$reference, required = mtname)
  assertions::assert_includes(df$reference, required = autosomes)

  # compute per chromosome coverage
  df$depth <- df$n_mapped * read_length/df$length

  # split idxstats df into autosomal vs mitochondrial df
  df_autosomes <- df[df$reference %in% autosomes, ,drop=false]
  df_mt <- df[df$reference %in% mtname, ,drop=false]

  # compute total autosome median coverage
  autosome_total_length <- sum(df_autosomes$length)
  autosome_total_reads <- sum(df_autosomes$n_mapped)

  # compute median coverage per autosomal chromosome
  autosome_depth_median_of_averages <- stats::median(df_autosomes$depth)

  # compute average autosome depth of coverage
  autosome_depth_average <- autosome_total_reads*read_length/autosome_total_length

  # compute average mitochondrial depth
  mitochondrial_depth_average <- df_mt$depth

  if(verbose){
    message("total autosome reads: ", round(autosome_total_reads, digits = 1))
    message("total autosome length: ", round(autosome_total_length, digits = 1))
    message("average autosome depth: ", round(autosome_depth_average, digits = 1), "x")
    # message("median autosome depth (median of average depth per autosome): ", round(autosome_depth_median, digits = 1), "x")
    message("median autosome depth (median of average depth per autosome): ", round(autosome_depth_median_of_averages, digits = 1), "x")
    message("average mtdna depth: ", round(mitochondrial_depth_average, digits = 1), "x")
  }

  # calculate mtdna
  calculate_mtdna_tumour(
    mitochondrial_depth = mitochondrial_depth_average,
    autosome_depth = autosome_depth_median_of_averages,
    background_ploidy = background_ploidy,
    tumour_ploidy = tumour_ploidy,
    tumour_purity = tumour_purity
  )
}


#' Calculate typical number of mitochondria per cell
#'
#' from a bulk whole-genome tumour sample, calculate number of mitochondrial genomes per cell.
#' this function includes options to correct for tumour ploidy and purity,
#' but by default will assume all cells in sample are diploid
#'
#' @param mitochondrial_depth depth of coverage of the mitochondrial genome
#' @param autosome_depth depth of coverage of the autosomal genome (do not normalise for ploidy yet! this function performs the normalisation for you)
#' @param background_ploidy expected ploidy of healthy cells in bulk wgs sample (default: 2 for diploid background)
#' @param tumour_ploidy median ploidy of tumour cells in bulk sample (default: 2 for diploid tumours)
#' @param tumour_purity proportion of total cells in sample that are tumour cells (must be between 0 and 1). used to correct for tumour ploidy.
#' @param per_cell return mitochondrial genomes per cell. If FALSE, return mitochondrial genomes per autosome (single copy)
#' @returns number of mitochondrial genomes per cell. if \code{per_cell=false} returns number of mitochondrial genomes per single copy of the autosome.
#'
#' @export
#'
#' @examples
#' calculate_mtdna_tumour(
#'   mitochondrial_depth = 100,
#'   autosome_depth = 10,
#'   background_ploidy = 2,
#'   tumour_ploidy = 2,
#'   tumour_purity = 1
#' )
calculate_mtdna_tumour <- function(mitochondrial_depth, autosome_depth, background_ploidy = 2, tumour_ploidy = 2, tumour_purity = 1, per_cell=TRUE){
  typical_number_of_autosome_copies_per_cell <- background_ploidy * (1-tumour_purity) + tumour_ploidy * tumour_purity
  mtgenomes_per_autosome <- typical_number_of_autosome_copies_per_cell * mitochondrial_depth / autosome_depth

  if(per_cell)
    return(mtgenomes_per_autosome * typical_number_of_autosome_copies_per_cell)
  else
    return(mtgenomes_per_autosome)
}

#' calculate mitochondrial dna burden
#'
#' @param mitochondrial_depth depth of coverage of the mitochondrial genome
#' @param autosome_depth depth of coverage of the autosomal genome (do not normalise for ploidy yet! this function performs the normalisation for you)
#' @param ploidy autosomal ploidy. By default assumes cells are diploid (ploidy = 2)
#' @param per_cell return mitochondrial genomes per cell. If FALSE, return mitochondrial genomes per autosome (single copy)
#'
#' @returns Number of mitochondrial genomes per cell. If \code{per_cell=false} returns number of mitochondrial genomes per single copy of the autosome.
#' @export
#'
#' @examples
#' calculate_mtdna(
#'   mitochondrial_depth = 100,
#'   autosome_depth = 10,
#'   ploidy = 2
#' )
calculate_mtdna <- function(mitochondrial_depth, autosome_depth, ploidy = 2, per_cell=TRUE){
  if(per_cell)
    ploidy * ploidy * mitochondrial_depth / autosome_depth
  else
    ploidy * mitochondrial_depth / autosome_depth
}

