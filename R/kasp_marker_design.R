#' Design KASP markers based on causal variants.
#' @param vcf_file Path to the vcf file containing identified variants from
#' re-sequence data.
#' @param marker_ID Designated name of variant for marker design. Name must be
#' contained in `vcf_file`.
#' @param chr A character value representation the chromosome description in the
#' `genome_file`. Providing this helps to save memory in R.
#' @param genome_file Path to reference genome file in fasta format, either
#' compressed (.gz) or decompressed.
#' @param plot_draw A logical value indicating whether to plot 100 bp upstream
#' and downstream KASP sequence alignment to reference genome.
#' @param plot_file Path to save the sequence alignment if `plot_draw = TRUE`.
#' @param region_name A n optional character value for assigned region name.
#' @param vcf_geno_code A character vector for genotype coding for samples
#' in imported `vcf_file`. The order of the genotype codes in this vector matters:
#' `c('homo_alt_allele', 'heterozygous', 'homo_ref_allele', 'missing_genotype_call')`.
#' @param maf A numeric value between 0 and 1 representing the minor allele
#' frequency for variant subset.
#'
#' @examples
#' \donttest{
#'# Example to design a KASP marker on a substitution variant
#' path <- tempdir()
#' path1 <- "https://raw.githubusercontent.com/awkena/panGB/main/Chr02.fa.gz"
#' path2 <-  system.file("extdata", "Sobic.002G302700_SNP_snpeff.vcf",
#'                      package = "panGenomeBreedr",
#'                      mustWork = TRUE)
#' ma1 <- kasp_marker_design(vcf_file = path2,
#'                           genome_file = path1,
#'                           marker_ID = "SNP_Chr02_69200443",
#'                           chr = "Chr02",
#'                           plot_draw = TRUE,
#'                           plot_file = path,
#'                           region_name = "ma1",
#'                           maf = 0.05)
#'
#' # View marker alignment output from temp folder
#' path3 <- file.path(path, list.files(path = path, "alignment_"))
#' # system(paste0('open "', path3, '"'))
#'
#' on.exit(unlink(path))
#' }

#'
#' @details
#' This function provides the intertek sequence to be used for marker development
#' for the selected casual variants. It provides all the information of allele and
#' location to fill Intertek form.It needs a vcf file of variants calls, and a genome
#' sequence of the target crop in fasta format.
#'
#'
#' @returns A data frame containing all information required for KASP marker
#' design and a DNA sequence alignment to the reference genome.
#'
#' @importFrom Biostrings fasta.index readDNAStringSet DNAStringSet
#' @importFrom GenomicRanges GRanges
#' @import VariantAnnotation
#' @import msa
#' @importFrom IRanges IRanges
#' @importFrom stats na.omit
#'
#' @export

kasp_marker_design <- function(vcf_file, #path and file name of tbi index
                               marker_ID, # what variant want to design the marker for
                               chr = NULL, # Chromosome with the variant
                               genome_file, #genome file path
                               plot_draw = TRUE, #Do you want a plot drawn? Yes or No
                               plot_file = tempdir(),
                               region_name = 'loc_1', #assigned region name
                               vcf_geno_code = c('1|1', '0|1', '0|0', '.|.'),
                               maf = 0.05) # Minor allele frequency
{

  # Extracting genomic region from vcf file
  vcf_gene <- VariantAnnotation::readVcf(file = vcf_file)

  if (vcf_gene@fixed@nrows >= 1) {

    #converting to dataframe the variant names
    variant_table <- as.data.frame(names(vcf_gene))

    #add chromsome
    variant_table$chrom <- as.character(vcf_gene@rowRanges@seqnames)

    #adding position
    variant_table$pos <- vcf_gene@rowRanges@ranges@start

    # Extracting reference allele to column
    variant_table$reference <- as.character(VariantAnnotation::ref(vcf_gene))

    # Extracting alternate allele to column
    variant_table$alternate <- as.character(unlist(VariantAnnotation::alt(vcf_gene)))

    # trying type of variation
    # variant_table$type <- NA

    tryCatch({

      variant_table$type <-  ifelse(VariantAnnotation::isInsertion(vcf_gene), "Insertion",
                                    ifelse(VariantAnnotation::isSubstitution(vcf_gene), "Substitution",
                                           ifelse(VariantAnnotation::isDeletion(vcf_gene), "Deletion", "Others")))

    })

    # getting allele frequency
    geno <- as.data.frame(VariantAnnotation::geno(vcf_gene))
    geno <- geno[geno$group_name == 'GT', -c(1:2)] # getting only the alleles

    #converting to numeric to calculate MAF
    geno <- ifelse(geno == vcf_geno_code[1], 2,
                   ifelse(geno == vcf_geno_code[2], 1,
                          ifelse(geno == vcf_geno_code[3], 0, NA)))

    # converting all to numeric
    tryCatch({
      suppressWarnings(geno <- apply(geno, MARGIN = 2, as.numeric))
    })

    # calculating frequencies
    variant_table$MAF <- round(rowSums(geno, na.rm = TRUE)/(ncol(geno)*2), digits = 3)

    # Renaming column
    names(variant_table)[1] <- c("ID")
  }

  # Reading genome by chromosome subset or whole genome
  # Indexing for genome file
  dd <-  Biostrings::fasta.index(genome_file, seqtype = "DNA")

  # Read only the sequence for the chromosome with the variant if provided,
  # else read the whole genome sequence
  if (!is.null(chr) && chr %in% dd$desc) {

    indx <-  which(dd$desc == chr)

    genome <-  Biostrings::readDNAStringSet(dd[indx, ])

  } else {

    genome <- Biostrings::readDNAStringSet(genome_file)

  }

  # genome <- Biostrings::readDNAStringSet(genome_file)

  #### extracting information from variant table
  variant <- variant_table[variant_table$ID == marker_ID,] #variant information
  pos <- as.numeric(variant$pos) #position of variant
  type <- variant$type #type of variant
  w_var <- nchar(variant$reference) # length of variant
  chromosome <- variant$chrom

  # snp variants only for adjacent highly polymorphic regions
  variants_snp_only <- variant_table[variant_table$type == 'Substitution',] # getting SNPs

  ###### Setting boundaries for markers
  ####### Upstream
  if (variant$type == 'Deletion') {
    up_bound_end <- pos # end of upstream boundary
    up_bound_start <- up_bound_end-99 #start of upstream boundary
    down_bound_start <- pos+w_var #start of downstream boundary
    down_bound_end <- down_bound_start+99 # end of downstream boundary

  } else if (variant$type == 'Substitution') {
    up_bound_end <- pos-1 # end of upstream boundary
    up_bound_start <- up_bound_end-99 #start of upstream boundary
    down_bound_start <- pos+1 #start of downstream boundary
    down_bound_end <- down_bound_start+99 # end of downstream boundary

  } else if (variant$type == 'Insertion') {
    up_bound_end <- pos # end of upstream boundary
    up_bound_start <- up_bound_end-99 #start of upstream boundary
    down_bound_start <- pos+1 #start of downstream boundary
    down_bound_end <- down_bound_start+99 # end of downstream boundary
  }

  # Process upstream and downstream variants
  up_down_variant <- function(variant_snp_dat,
                              MAF = maf,
                              start,
                              end,
                              type = c('up', 'down')) {

    type <- match.arg(type)

    # variants dataframe
    variants_100bp <- data.frame(
      ID = NA, # marker id
      reference = NA, #reference allele
      alternate = NA, #alternate allele
      pos = NA, # position
      type = NA, #type of variant
      MAF = NA, # minor allele frequency
      ontology = NA, #ontology term
      effect = NA, # effect
      gene_name = NA, #gene name
      one_letter = NA # one letter code
    )

    # variants
    if (type == 'up') {

      variants_100bp <- variant_snp_dat[variant_snp_dat$pos >= start &
                                          variant_snp_dat$pos < end,]

    } else {

      variants_100bp <- variant_snp_dat[variant_snp_dat$pos > start &
                                          variant_snp_dat$pos <= end,]
    }

    # Filtering by MAF
    variants_100bp <- variants_100bp[variants_100bp$MAF > MAF,] #MAF filtering

    #ordering by MAF
    variants_100bp <- variants_100bp[order(variants_100bp$MAF, decreasing = TRUE),]

    #selecting the two variants with higher MAF
    variants_100bp <- stats::na.omit(variants_100bp[c(1,2),])

    #ordering by pos
    variants_100bp <- variants_100bp[order(variants_100bp$pos, decreasing = FALSE),]

    return(variants_100bp)

  }

  variants_100bp_up <- up_down_variant(variant_snp_dat = variants_snp_only,
                                       MAF = maf,
                                       start = up_bound_start,
                                       end = up_bound_end,
                                       type = 'up')
  # What happens if it returns no data?
  variants_100bp_down <- up_down_variant(variant_snp_dat = variants_snp_only,
                                         MAF = maf,
                                         start = down_bound_start,
                                         end = down_bound_end,
                                         type = 'down')

  # Creating iupacode table for substitutions
  iupacode <- data.frame(one_letter = c('R', 'Y','S','W','K','M', 'R', 'Y','S','W','K','M'),
                         allele_comb1 = c('AG','CT','GC','AT','GT','AC','GA','TC','CG','TA','TG','CA'))

  # Creating collapse column of alleles for substitution
  variants_100bp_down$collapse <- paste(variants_100bp_down$reference,
                                        variants_100bp_down$alternate, sep = "")
  variants_100bp_up$collapse <- paste(variants_100bp_up$reference,
                                      variants_100bp_up$alternate, sep = "")

  # Substitution collapse allele for the the one letter iupac code
  iupac_conv <- function(iupacode = iupacode, variants_100bp) {

    if (nrow(variants_100bp) > 0) {
      i <- 1
      while (i <= nrow(variants_100bp)) {
        variants_100bp$one_letter[i] <-
          iupacode[rowSums(iupacode == variants_100bp$collapse[i]) > 0, ]$one_letter
        i <- i+1
      }
    }

    return(variants_100bp)
  }

  #upstream
  variants_100bp_up <- iupac_conv(iupacode = iupacode,
                                  variants_100bp = variants_100bp_up)
  #downstream
  variants_100bp_down <- iupac_conv(iupacode = iupacode,
                                    variants_100bp = variants_100bp_down)

  ### Getting the actual upstream and downstream sequences

  get_actual_seq <- function(variants_100bp,
                             genome = genome,
                             bound_start,
                             bound_end) {

    ####### start stop list
    start_stop_list <- c(bound_start, variants_100bp$pos - 1, bound_end)

    # creating empty sequence
    stream_sequence <- Biostrings::DNAStringSet()

    j <- 1
    for ( i in 1:(length(start_stop_list)-1)) {

      #adjust start
      if (i == 1) {

        tmp_start <- start_stop_list[i] # first iteration

      } else {

        tmp_start <- start_stop_list[i] + 1 # adding a base

      }

      #adjust stop
      if (i == length(start_stop_list)-1) {

        tmp_stop <- start_stop_list[i+1] # if is the last iteration

      } else {

        tmp_stop <- start_stop_list[i+1]-1 # removing a base

      }

      # Extracting genomic region from reference genome
      region <- GenomicRanges::GRanges(chromosome, IRanges::IRanges(tmp_start, tmp_stop))

      # Creating empty base to add
      base_to_add <- Biostrings::DNAStringSet()

      if (j <= nrow(variants_100bp)) {

        #getting one letter code for the polymorphism upstream
        base_to_add <- Biostrings::DNAStringSet(variants_100bp$one_letter[j])

        j <- j+1 #counter

      }

      # getting sequence for the region
      sequence <- Biostrings::DNAStringSet(Biostrings::getSeq(genome, region))

      # merging sequence and base to add
      tmp_seq <- c(sequence, base_to_add)

      # unlisting and making new sequence
      new_seq <- Biostrings::DNAStringSet(unlist(tmp_seq))

      #saving to upstream sequence, merging with upstream sequence
      stream_sequence <- Biostrings::DNAStringSet(unlist(c(stream_sequence, new_seq)))
    }

    return(stream_sequence)
  }

  # Upstream
  upstream_sequence <- get_actual_seq(variants_100bp = variants_100bp_up,
                                      genome = genome,
                                      bound_start = up_bound_start,
                                      bound_end = up_bound_end)

  # Downdtream
  downstream_sequence <- get_actual_seq(variants_100bp = variants_100bp_down,
                                        genome = genome,
                                        bound_start = down_bound_start,
                                        bound_end = down_bound_end)



  ## Create output sequence annotated
  if (variant$type == 'Substitution') {

    reference_allele <- variant$reference # reference allele
    alternate_allele <- variant$alternate # alternate allele

    # intertek formatted sequence
    intertek_sequence <- paste(upstream_sequence, "[",
                               variant$reference,"/", variant$alternate, "]",
                               downstream_sequence, sep = "")

  } else if (variant$type == 'Deletion') {

    #variant without anchoring base
    varc <- substr(variant$reference, start = 2, stop = nchar(variant$reference))

    #reference allele
    reference_allele <- varc

    #alternate allele
    alternate_allele <- '-'

    # intertek formatted sequence
    intertek_sequence <- paste(upstream_sequence,"[", varc, "/-", "]",
                               downstream_sequence, sep = "")

  } else if (variant$type == 'Insertion') {

    #variant without anchoring base
    varc <- substr(variant$alternate, start = 2, stop = nchar(variant$alternate))

    #reference allele
    reference_allele <- '-'

    #alternate allele
    alternate_allele <- varc

    # intertek formatted sequence
    intertek_sequence <- paste(upstream_sequence, "[-/", varc, "]",
                               downstream_sequence, sep = "")

  }

  # Creating dataframe for export
  #Empty dataframe
  marker_data <- data.frame(SNP_Name = variant$ID,
                            SNP = variant$type,
                            Marker_Name = region_name,
                            Chromosome = chromosome,
                            Chromosome_Position = variant$pos,
                            Sequence = intertek_sequence,
                            ReferenceAllele = reference_allele,
                            AlternativeAllele = alternate_allele)

  # assign("marker_data", marker_data, envir = .GlobalEnv) # creating dataframe

  # Plot alignment of upstream and downstream sequences to the reference seq
  if (plot_draw == TRUE) {

    # getting genomic range to extract
    region_control <- GenomicRanges::GRanges(chromosome, IRanges::IRanges(up_bound_start, down_bound_end))

    # getting reference sequence
    sequence_control <- VariantAnnotation::getSeq(genome, region_control)

    # merging upstream and downstream sequence
    tmp <- c(upstream_sequence, downstream_sequence)

    # merging all sequences for plotting
    sequences <- c(sequence_control, tmp)

    # converting to DNAstring set
    sequences <- Biostrings::DNAStringSet(sequences)

    # adding sequence names
    names(sequences) <- c('reference', 'upstream', 'downstream')

    # alignment
    alg <- msa::msaClustalOmega(sequences, order = 'input')

    #converting to dnastringset
    alignment <- Biostrings::DNAStringSet(alg)

    # Extract the aligned sequences
    aligned_strings <- as.character(alignment)

    # Re-implementing alignment plot in ggplot
    y <- NULL
    group <- NULL
    x <- NULL

    x <- data.frame(
      x = 1,
      y = c(0.5, 0.4, 0.3),
      seq = aligned_strings,
      group = c('reference', 'upstream', 'downstream'))

    # Define the colors you want for each group
    group_colors <- c("reference" = "black",
                      "upstream" = "red",
                      "downstream" = "blue")

    pp <- ggplot2::ggplot(data = x, ggplot2::aes(x = x, y = y, label = seq,
                                                 color = group)) +
      ggplot2::geom_text(size = 4.5, hjust = 'outward', family = 'mono',
                         key_glyph = 'rect') +

      ggplot2::labs(color = paste("Alignment for marker", marker_ID, sep = " ")) +
      ggplot2::scale_x_continuous(breaks = c(-2, 0, 2), expand = c(0, 2)) +
      ggplot2::scale_y_continuous(breaks = c(-2, 0, 2), expand = c(0, 2)) +
      ggplot2::scale_color_manual(values = group_colors,
                                  breaks = c('reference', 'upstream', 'downstream'),
                                  label = c('Reference', 'Upstream', 'Downstream')) +
      ggplot2::guides(color = ggplot2::guide_legend(label.position = "right",
                                                    title.position = 'top',
                                                    title.hjust = 0.5)) +
      ggplot2::theme(legend.title = ggplot2::element_text(size  = 12,
                                                          color = 'black'),
                     panel.background = ggplot2::element_rect('white'),
                     axis.line = ggplot2::element_blank(),
                     axis.text = ggplot2::element_blank(),
                     axis.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     legend.text = ggplot2::element_text(size  = 12,
                                                         color = 'black'),
                     legend.position = 'inside',
                     legend.position.inside = c(0.5, 0.6),
                     legend.direction = 'horizontal')

    # Save ggplot object as a pdf device
    ggplot2::ggsave(filename = paste0('alignment_', marker_ID, '.pdf'),
                    plot = pp,
                    device = "pdf",
                    path = file.path(plot_file),
                    units = "in",
                    width = 24,
                    height = 9)
  }

  return(marker_data)

}
