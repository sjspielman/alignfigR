#' alignfigR.
#'
#' Visualizing multiple sequence alignments with ggplot.
#' CAUTION: no sanity checking is performed to check if your provided data is an alignment or just sequences (but they will plot either way :) )
#' @name alignfigR
#' @docType package
#' @import ggplot2
NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to alignfigR!")
}
plot_step <- 1




#' Read a multiple sequence alignment file.
#'
#' This function parses a FASTA file containing molecular sequence data and returns a named-array of the sequences.
#' @param file     File name. NOTE: Only FASTA file are supported!
#' @return seq_array, a named-array of the parsed sequence data
#' @examples
#' plot_frame <- read_alignment(file = "/path/to/sequences/data.fasta")
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

 
    # Return sequence data parsed into named-array
    seq_vector 
}









#' Extract subset of sequence alignment
#'
#' This function builds a data frame to plot an alignment from a specified subset of the full alignment.
#' @param seqs        Sequence array
#' @param s       Step size for alignment block. Default, 1.
#' @param tlist   List of taxa (the actual labels, not order) to either restrict to or exclude from plot.
#' @param texl    Exclude those taxa?
#' @param clist   List of columns (indexed from 1) to either restrict to or exclude from plot.
#' @param cexl    Exclude those columns?
#' @return plot_frame, a data frame to be plotted
#' @examples
#' GIVE EXAMPLES
#' @export
extract_subalign <- function(seqs, s, tlist, clist, texcl, cexcl)
{
 
    # Initialize data frame with sequence information, coordinates for plotting  
    plot_frame <- data.frame( "x1"  = numeric(), 
                              "y1"  = numeric(),
                              "x2"  = numeric(),
                              "y2"  = numeric(), 
                              "name" = factor(), 
                              "seq"  = factor() )   

    # Remove or keep the columns via indexing
    if (!(cexcl){
        clist <- clist * -1
    }
    seq_index <- 1
    
    # Decide which sequences to keep, and which columns to keep within each sequence
    for (seq_name in names(seqs)){
        
        # Split the sequence for grabbing columns easily
        current_seq <- strsplit(seqs[seq_name], split="")[[1]]

    
        # Add sequence to plot_frame, if we want to keep it, selecting only desired columns
        if ( (texcl && (seq_name %in% tlist)) ||  (texcl == F && !(seq_name %in% tlist)) ){
            if (length(clist) == 0)
            {  
                retain_seq <- current_seq
            }
            else
            {
                retain_seq <- current_seq[clist]
            }
            pos_index <- 1
            for ( pos in retain_seq ){
                temp_frame <- data.frame( "x1" = pos_index, "y1" = seq_index, "x2" = pos_index + s, "y2" = seq_index + s, "name" = seq_name, "seq" = pos )
                plot_frame <- rbind(plot_frame, temp_frame)
                pos_index <- pos_index + 1
            }
            seq_index <- seq_index +  1 
        }
    }
              
    # Return data frame
    plot_frame    
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
#' @export
define_palette <- function( inpalette, uniques )
{
    palette <- c()
    ambigc <- "grey85"
    ambig <- c("?" = ambigc , "-" = ambigc, "*" = ambigc)

    # Random colors (also called if nothing specified)
    if (  tolower(inpalette) == "random" || is.na(inpalette) )
    {
        for (m in names(ambig)){
            uniques <- uniques[!uniques == m]
        }
        palette <- sample( colors(), length(uniques) )
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
        missing_color = "grey85"
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












#' Plot a multiple sequence alignment
#'
#' This function uses ggplot (in particular, w/ geom_rect) to plot a sequence alignment
#' @param seq_vector        Sequence vector parsed using read_alignment
#' @param palette           Named-array mapping sequence to color or a pre-defined color scheme (random, rainbow, etc.)
#' @param taxon_list        Array of taxa (the actual labels, not order) to either restrict to or exclude from plot.
#' @param taxon_exclude     Boolean argument indicating that taxa listed in taxon_list should be excluded from plot. Default: False
#' @param column_list       Array of columns (indexed from 1) to either restrict to or exclude from plot.
#' @param column_action     Boolean argument indicating that columns listed in column_list should be excluded from plot. Default: False
#' @return ggplot object which may be saved or edited as desired
#' @examples
#' align_plot <- plot_alignment(seq_vector, palette)
#' align_plot <- plot_alignment(seq_vector, palette, taxon_list = c("i_hate_this_organism"), column_list = c(1:25) )
#' align_plot <- plot_alignment(seq_vector, palette)
#' @export
plot_alignment <- function(seq_vector, palette = NA, step = 1, taxon_list = c(), column_list = c(), taxon_exclude = F, column_action = F)
{
    # Extract desired alignment subset
    plot_frame <- extract_subalign(seq_vector, step, taxon_list, column_list, taxon_exclude, column_exclude)


    # Determine alignment characters for palette construction
    unique_chars <- unique(plot_frame$seq)
    pal <- define_palette(palette, unique_chars)


    # Sort sequence columns so legend is alphabetical
    plot_frame$seq <- factor(plot_frame$seq, levels = sort(levels(plot_frame$seq)))
    
    # Set axis breaks. CURRENTLY TURNED OFF.
    #xbreaks <- calc_axis_breaks( max(plot_frame$x1) )
    #ybreaks <- calc_axis_breaks( length(unique(plot_frame$name)) )
    xbreaks <- c()
    ybreaks <- c()
    
    # Plot
    p <- ggplot() + geom_rect(plot_frame, mapping=aes(xmin=x1-1, xmax=x2-1, ymin=y1-1, ymax=y2-1, fill = seq), linetype=0) +
         scale_fill_manual(values=pal)  + 
         scale_x_continuous(breaks = xbreaks) + scale_y_continuous(breaks = ybreaks) + 
         theme(
            legend.title     = element_blank(),
            axis.title.x     = element_blank(), 
            axis.title.y     = element_blank(),
            panel.background = element_blank(), 
            panel.border     = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background  = element_blank())
             
    p
}


















     
            


################ CURRENTLY NOT USED. ##################
# Construct axis breaks for alignment plot
#
# Determine appropriate tick placement given the dataset size
# @param size     Number of sequences OR length of alignment (or longest sequence) to plot
# @return array axis break points
# @examples
# ybreaks <- calc_axis_breaks(number_sequences)
# xbreaks <- calc_axis_breaks(alignment_length)
# @export
calc_axis_breaks <- function(size)
{
    if (size <= 10){
        breaks <- c(seq(1,size), size)
    }
    else if (size > 10 & size <= 100){
        breaks <- c(1, seq(10, size, 10))
        if (size - breaks[length(breaks)] < 5){
            breaks[length(breaks)] <- size
        } 
        else{
            breaks <- c(breaks, size)
        }
    }
    else {
        breaks <- c(1, seq(50, size, 50))
        if (size - breaks[length(breaks)] < 25){
            breaks[length(breaks)] <- size
        } 
        else{
            breaks <- c(breaks, size)
        }
    }
    breaks <- sort(unique(breaks))
    breaks
    
}


                    
  






