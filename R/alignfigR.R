#' alignfigR.
#'
#' Creating figures multiple sequence alignments with ggplot2.
#' @name alignfigR
#' @docType package
#' @import ggplot2
#' @importFrom grDevices colors
NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to alignfigR!")
}

null_color           <- "grey85" ## Default missing and ambiguous
y1 <- y2 <- x1 <- x2 <- NULL
default_plot_step    <- 1


#' Read a multiple sequence alignment file.
#'
#' This function parses a FASTA file containing molecular sequence data and returns a named-array of the sequences.
#' @param file     File name. NOTE: Only FASTA file are supported!
#' @return seq_array, a named-array of the parsed sequence data
#' @examples
#' fasta_file <- system.file("extdata", "example.fasta", package = "alignfigR")
#' plot_frame <- read_alignment(file = fasta_file)
#' @export
read_alignment <- function(file){

    raw_data <- readLines( file, warn = FALSE )
    seq_vector <- c()
    seq_name <- ""
    for (line in raw_data){

        # New sequence record? Reset numbering
        if ( grepl("^>", line) ){
            seq_name <- sub("^>", "", line)
            seq_vector[seq_name] <- ""
        }
        else {
            temp_seq <- gsub(" ","",line)
            temp_seq <- gsub("\n","",temp_seq)
            seq_vector[seq_name] <- paste( seq_vector[seq_name], temp_seq, sep="" )
        }

    }
    # Is this an alignment?
    seq_list <- strsplit(seq_vector, split = "")
    lengths <- sapply(seq_list, length)
    if ( sum(lengths != lengths[1]) != 0 )
        stop("Your provided file is not an alignment. Please provide an alignment file in FASTA format to use alignfigR.")

    # Return sequence data parsed into named list
    seq_list
}







#' Define color palette.
#'
#' This function sets up, either using default or user-specified options, the color-coding
#' scheme used to plot sequences.
#' @param inpalette  Either a user-specified named-array of colors or flag for default options (dna, rna, protein, random, ...)
#' @param uniques    Unique characters found in alignment. Used to create the random color scheme.
#' @return Color palette named-array.
#' @examples
#' palette <- define_palette("DNA", c("A", "G", "T"))
#' palette <- define_palette("protein")
#' @export
define_palette <- function( inpalette, uniques )
{
    palette <- c()
    ambigc <- null_color
    ambig <- c("?" = ambigc , "-" = ambigc, "*" = ambigc)

    # Random colors (also called if nothing specified)
    if (  tolower(inpalette) == "random" || is.na(inpalette) )
    {
        for (m in names(ambig)){
            uniques <- uniques[!uniques == m]
        }
        subcolors <- colors()[colors() != ambigc] ## Ensure null_color is not in the random scheme
        palette <- sample( subcolors, length(uniques) )

        names(palette) <- uniques

        palette <- c( palette, ambig )
    }

    # Default alphabet colors
    else if (tolower(inpalette) == "rna" || tolower(inpalette) == "dna" || tolower(inpalette) == "protein")
    {
        if (tolower(inpalette) == "rna" || tolower(inpalette) == "dna")
        {
            missing_names <- c("M", "R", "W", "S", "Y", "K", "V", "H", "D", "B", "X", "N")
            palette <- c("A" = "mediumblue", "C" = "orangered1", "G" = "limegreen", "T" = "khaki1", "U" = "khaki1")
        }
        else {
            missing_names <- c("B", "X", "Z")
            palette <- c("A" = "limegreen", "G" = "limegreen",
                         "C" = "lightgreen",
                         "D" = "darkgreen", "E" = "darkgreen", "N" = "darkgreen", "Q" = "darkgreen",
                         "I" = "lightblue1", "L" = "lightblue1", "M" = "lightblue1", "V" = "lightblue1",
                         "F" = "lavender", "W" = "lavender", "Y" = "lavender",
                         "H" = "navy",
                         "K" = "orange", "R" = "orange",
                         "P" = "salmon",
                         "S" = "red", "T" = "red")
        }

        # Missing color default
        missing_color <- null_color
        missing_palette <- rep(missing_color, length(missing_names))
        names(missing_palette) <- missing_names

        palette <- c( palette, missing_palette, ambig )
    }


    ##
    ##
    ## Add more color schemes maybe?
    ##
    ##

    # Assign user-provided colors
    else {
        palette <- inpalette
    }
    palette

}




#' Extract subset of sequence alignment
#'
#' This function builds a data frame to plot an alignment from a specified subset of the full alignment.
#' @param seqs       Sequence list, as parsed by the function `read_alignment`
#' @param plot_step       Step size for alignment block. Default, 1.
#' @param tlist      Array of taxa intended to restrict figure to or to exclude from figure.
#' @param clist      Array of columns intended to restrict figure to or to exclude from figure.
#' @param texcl       Boolean indicating if taxa in tlist should be excluded. Default, False
#' @param cexcl       Boolean indicating if columns in clist should be excluded. Default, False
#' @param cincl       Boolean indicating if columns in clist should be included. Default, False
#' @return plot_frame, a data frame to be plotted
#' @examples
#' fasta_file <- system.file("extdata", "example.fasta", package = "alignfigR")
#' plot_frame <- read_alignment(file = fasta_file)
#' subset_seq_list <- extract_subalign(plot_frame, tlist = c("Cow", "Human", "Whale"), texcl = TRUE)
#' subset_seq_list <- extract_subalign(plot_frame, clist = 1:25)
#' @export
extract_subalign <- function(seqs, plot_step = 1, tlist = c(), clist = c(), texcl = FALSE, cexcl = FALSE, cincl = FALSE)
{
    # Create subset of seqs containing only the desired taxa to plot
    if (length(tlist) == 0){
        sub_seqs_raw <- seqs

    }else if (texcl){
        sub_seqs_raw <- seqs[!(names(seqs) %in% tlist)] # Exclude sequences in the provided list
    }else
    {
        sub_seqs_raw <- seqs[tlist]
    }

    # Futher subset the sequences to contain only the desired columns
    if (length(clist) == 0){
        clist <- 1:length(sub_seqs_raw[[1]])
    }
    if (cexcl){
        sub_seqs <- lapply(sub_seqs_raw, `[`, (-1*clist))
    }
    if (cincl){
    sub_seqs <- lapply(sub_seqs_raw, `[`, (1*clist))
   } else
    {
        sub_seqs <- lapply(sub_seqs_raw, `[`)
    }

    # Create the data frame to plot
    sub_seqs <- rev(sub_seqs) # For proper plotting direction
    each_length <- length(sub_seqs[[1]])
    seqnames <- c(t(replicate(each_length, names(sub_seqs))))
    seqletters <- unlist(sub_seqs)
    y1 <- c(t(replicate(each_length, 1:length(sub_seqs))))
    y2 <- y1 + plot_step
    x1 <- rep(1:each_length, length(sub_seqs))
    x2 <- x1 + plot_step

    plot_frame <- data.frame( "x1"  = x1,
                              "y1"  = y1,
                              "x2"  = x2,
                              "y2"  = y2,
                              "name" = seqnames,
                              "seq"  = seqletters)
    rownames(plot_frame) <- NULL
    plot_frame
}














#' Plot a multiple sequence alignment
#'
#' This function uses ggplot (in particular, w/ geom_rect) to plot a sequence alignment
#' @param seq_list         Sequence list parsed using the function `read_alignment`
#' @param palette          Named-array mapping sequence to color or a pre-defined color scheme (random, rainbow, etc.)
#' @param taxa             Array of taxa (the actual labels, not order) intended to restrict figure to or to exclude from figure
#' @param taxon_labels Boolean argument indicating that the Y-axis should be labeled with taxon names. Default: False
#' @param exclude_taxa     Boolean argument indicating that taxa should be excluded from plot. Default: False
#' @param columns          Array of columns (indexed from 1) intended to restrict figure to or to exclude from figure.
#' @param exclude_columns  Boolean argument indicating that columns should be excluded from plot. Default: False
#' @param include_columns  Boolean argument indicating that columns should be included from plot. Default: False
#' @param legend_title     String determining title of legend. Default: "Character"
#' @return ggplot object which may be saved or edited as desired
#' @examples
#' fasta_file <- system.file("extdata", "example.fasta", package = "alignfigR")
#' plot_frame <- read_alignment(file = fasta_file)
#' align_plot <- plot_alignment(plot_frame, "DNA")
#' align_plot <- plot_alignment(plot_frame, "protein")
#' align_plot <- plot_alignment(plot_frame, taxa = c("Cow", "Whale"), columns = c(1:25))
#' align_plot <- plot_alignment(plot_frame, taxa = c("Whale"), exclude_taxa = TRUE)
#' align_plot <- plot_alignment(plot_frame, legend_title = "") ## Remove the title
#' @export
plot_alignment <- function(seq_list, palette = NA, taxa = c(), taxon_labels = FALSE, columns = c(), exclude_taxa = FALSE, exclude_columns = FALSE, include_columns = FALSE, legend_title = "Character")
{
    # Extract desired alignment subset
    plot_frame <- extract_subalign(seq_list, default_plot_step, taxa, columns, exclude_taxa, exclude_columns, include_columns)

    # Determine alignment characters for palette construction
    unique_chars <- unique(plot_frame$seq)
    pal <- define_palette(palette, unique_chars)

    # Sort sequence columns so legend is alphabetical
    plot_frame$seq <- factor(plot_frame$seq, levels = sort(levels(plot_frame$seq)))

    # Plot
    theme_set(theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank()))
    if ( taxon_labels == FALSE){
    p <- ggplot() +
        geom_rect(plot_frame, mapping=aes(xmin=x1-1, xmax=x2-1, ymin=y1-1, ymax=y2-1, fill = seq), linetype=0) + scale_fill_manual(values=pal, name = legend_title)
    }
    else {
      p <- ggplot() + geom_rect(plot_frame, mapping=aes(xmin=x1-1, xmax=x2-1, ymin =
                                y1-1, ymax=y2-1, fill = seq), linetype=0) +                                         scale_fill_manual(values=pal, name = legend_title) +                                 scale_y_discrete(limits = names(seq_list))
    }
    p
}














