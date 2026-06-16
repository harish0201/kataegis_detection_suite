
################ GGPLOT2 ############################
library(data.table)
library(ggplot2)
library(scales)

plot_rainfall_genome_ggplot2 <- function(maf_file, sample_id,
                                     log_scale = TRUE,
                                     highlight_kataegis = NULL,
                                     show_legend = TRUE,
                                     chrom_sizes = NULL,
                                     cex_axis = 0.65,
                                     include_other_mutations = TRUE,
                                     facet_by_chromosome = FALSE) {

  # --- Helper functions ---
  order_chromosomes <- function(chr_vec) {
    raw <- as.character(chr_vec)
    clean <- sub("^chr", "", raw, ignore.case = TRUE)
    numeric_val <- suppressWarnings(as.numeric(clean))
    is_num <- !is.na(numeric_val)
    special_rank <- rep(9999, length(clean))
    special_rank[toupper(clean) %in% c("X")] <- 1000
    special_rank[toupper(clean) %in% c("Y")] <- 1001
    special_rank[toupper(clean) %in% c("MT", "M", "MITO")] <- 1002
    sort_key1 <- ifelse(is_num, numeric_val, special_rank)
    sort_key2 <- ifelse(is_num, 0, as.numeric(factor(clean)))
    ord <- order(sort_key1, sort_key2, clean)
    unique(raw[ord])
  }

  get_mutation_type <- function(ref, alt) {
    mut <- paste0(ref, ">", alt)
    map <- c("C>A"="C>A","C>G"="C>G","C>T"="C>T",
             "T>A"="T>A","T>C"="T>C","T>G"="T>G")
    type <- map[mut]; type[is.na(type)] <- "Other"
    type
  }
  
  geom_mean <- function(x) exp(mean(log(x))) # Helper for geometric mean
  
  # --- Load and prep MAF ---
  maf <- data.table::fread(maf_file)
  maf <- maf[Tumor_Sample_Barcode == sample_id & Variant_Type == "SNP"]
  if (nrow(maf) == 0) stop("No SNPs found for that sample.")

  if (!all(c("Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2") %in% names(maf)))
    stop("MAF must have Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2 columns.")

  # --- Chromosome ordering and cumulative positions ---
  # Define desired chromosome order
  if (!is.null(chrom_sizes)) {
    cs <- data.table::as.data.table(chrom_sizes)
    stopifnot(all(c("Chromosome","Size") %in% names(cs)))
    chr_levels <- order_chromosomes(cs$Chromosome)
  } else {
    chr_levels <- order_chromosomes(unique(maf$Chromosome))
  }
  maf[, Chromosome := factor(as.character(Chromosome), levels = chr_levels)]

  # Calculate cumulative positions (needed for genome-wide view and kataegis mapping)
  if (!is.null(chrom_sizes)) {
    cs <- data.table::as.data.table(chrom_sizes)[, .(Chromosome, Size)]
    cs <- cs[order(match(Chromosome, chr_levels))]
    missing_chr <- setdiff(chr_levels, cs$Chromosome)
    if (length(missing_chr) > 0)
      cs <- rbind(cs, data.table::data.table(Chromosome = missing_chr, Size = 0))
    cs <- cs[match(chr_levels, Chromosome)]
    cs[, cum_start := data.table::shift(cumsum(as.numeric(Size)), fill = 0)]
    chr_lengths <- cs[, .(Chromosome, Size, cum_start)]
  } else {
    tmp <- maf[, .(Size = max(Start_Position, na.rm = TRUE)), by = Chromosome]
    tmp <- tmp[match(chr_levels, tmp$Chromosome)]
    tmp$Size[is.na(tmp$Size)] <- 0
    tmp[, cum_start := data.table::shift(cumsum(as.numeric(Size)), fill = 0)]
    chr_lengths <- tmp[, .(Chromosome, Size, cum_start)]
  }

  maf <- merge(maf, chr_lengths[, .(Chromosome, cum_start)], by = "Chromosome", all.x = TRUE, sort = FALSE)
  maf[, genome_pos := Start_Position + cum_start]
  data.table::setorder(maf, genome_pos)

  # --- Intermutation distance & Mutation Type ---
  maf[, dist_to_prev := c(NA, diff(genome_pos))]
  maf <- maf[!is.na(dist_to_prev)]
  maf[, dist_adj := dist_to_prev + 1] # dist_adj is IMD+1
  maf[, Mutation_Type := get_mutation_type(Reference_Allele, Tumor_Seq_Allele2)]
  
  # --- Setup for ggplot ---
  
  # Define mutation colors
  color_map <- c("C>A"="#4DBBD5","C>G"="#E64B35","C>T"="#00A087",
                 "T>A"="#3C5488","T>C"="#F39B7F","T>G"="#8491B4", "Other"="#999999")
  
  # Filter 'Other' mutations if requested
  if (!include_other_mutations) {
    maf <- maf[Mutation_Type != "Other"]
  }
  
  # Determine X-Axis and Faceting based on toggle
  if (facet_by_chromosome) {
      x_axis_col <- "Start_Position"
      x_label <- "Position (bp)"
      plot_title <- paste("Chromosome Rainfall Plot:", sample_id)
      plot_facet <- facet_wrap(~ Chromosome, scales = "free_x", ncol = 4)
  } else {
      x_axis_col <- "genome_pos"
      x_label <- "" # X-axis labels are handled by custom axis ticks/labels
      plot_title <- paste("Genome-wide Rainfall Plot:", sample_id)
      plot_facet <- NULL
  }
  
  # --- Base Plot ---
  p <- ggplot(maf, aes_string(x = x_axis_col, y = "dist_adj", fill = "Mutation_Type")) +
    
    # Points Layer
    geom_point(pch = 21, color = "black", size = 1.2, alpha = 0.8) +
    
    # Y-axis (Log Scale with clean labels)
    scale_y_continuous(
        name = "Intermutation Distance (bp)",
        trans = if (log_scale) "log10" else "identity",
        # Use trans_breaks to get 1, 10, 100, ... labels
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        # Use math_format to correctly label powers of 10
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    
    # Color Scale
    scale_fill_manual(
        values = color_map,
        name = "Mutation Type",
        breaks = names(color_map)[names(color_map) %in% maf$Mutation_Type] 
    ) +

    # Facet Layer (conditional)
    plot_facet +

    # Theme and Labels
    labs(x = x_label, title = plot_title) +
    
    theme_bw() +
    theme(
        legend.position = if (show_legend) "top" else "none",
        panel.grid.minor = element_blank()
    )

  # --- Custom Elements (only for genome-wide view) ---
  if (!facet_by_chromosome) {
    # Genome-wide X-axis labels (Chromosome midpoints)
    chr_mid <- chr_lengths[, .(Chromosome, midpoint = cum_start + Size / 2)]
    p <- p + scale_x_continuous(
        breaks = chr_mid$midpoint,
        labels = chr_mid$Chromosome,
        expand = expansion(mult = c(0.01, 0.01)) # Reduce margin at plot edges
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))

    # Vertical chromosome dividers
    vlines <- unique(chr_lengths$cum_start)
    p <- p + geom_vline(xintercept = vlines[-1], color = "gray85", linetype = "dotted")
  }

  # Horizontal line at 1 kbp (IMD=1000, plotted at log10(1001))
  p <- p + geom_hline(yintercept = 1001, color = "gray60", linetype = "dashed", size = 0.5)
  
  # --- Kataegis Arrows (using annotation_custom for precise placement) ---
  if (!is.null(highlight_kataegis)) {
    
    # Calculate geometric mean IMD for each highlight region
    kataegis_arrows <- lapply(highlight_kataegis, function(region) {
      if (!is.list(region) || is.null(region$sample) || region$sample != sample_id) return(NULL)
      
      chrname <- as.character(region$chromosome)
      if (!(chrname %in% chr_lengths$Chromosome)) return(NULL)
      
      # Determine genomic start/end for filtering
      chr_cum <- chr_lengths[Chromosome %in% chrname]$cum_start
      start_g <- region$start + chr_cum
      end_g <- region$end + chr_cum
      center_g <- (start_g + end_g) / 2
      
      # Identify mutations within the region (using genome_pos)
      in_region_maf <- maf[genome_pos >= start_g & genome_pos <= end_g]
      
      if (nrow(in_region_maf) > 0) {
        # Calculate geometric mean of IMD+1 values
        mid_imd <- geom_mean(in_region_maf$dist_adj)
        
        # Determine X-position for the arrow (local or genome-wide)
        x_pos <- if (facet_by_chromosome) (region$start + region$end) / 2 else center_g
        
        # Mark points with black border (by modifying the 'color' aesthetic for the range)
        p <<- p + geom_point(data = in_region_maf, 
                             aes_string(x = x_axis_col, y = "dist_adj"), 
                             pch = 21, color = "black", size = 1.2, alpha = 1)
        
        # Prepare data for the arrow
        return(data.table::data.table(
            x = x_pos,
            y_start = 10^(log10(min(maf$dist_adj)) * 1.5), # Start point near x-axis
            y_end = mid_imd,
            Chromosome = chrname
        ))
      }
      return(NULL)
    })
    
    kataegis_arrows <- rbindlist(kataegis_arrows)
    
    if (nrow(kataegis_arrows) > 0) {
      p <- p + geom_segment(
          data = kataegis_arrows,
          aes(x = x, xend = x, y = y_start, yend = y_end),
          color = "black", 
          size = 0.5,
          arrow = arrow(length = unit(0.2, "cm"), type = "open"),
          inherit.aes = FALSE
      )
    }
  }

  return(p)
}
#########################################################################################


# Note: This function requires the 'data.table' package to be loaded.
library(data.table)
plot_rainfall_genome_base <- function(maf_file, sample_id,
                                     chromosomes = NULL, # Default is NULL (Genome-wide)
                                     log_scale = TRUE,
                                     highlight_kataegis = NULL,
                                     ncol = 3,           # Only used if 'chromosomes' is supplied
                                     show_legend = TRUE,
                                     chrom_sizes = NULL,
                                     cex_axis = 0.65, 
                                     include_other_mutations = TRUE) {

  # --- Helper functions ---
  order_chromosomes <- function(chr_vec) {
    raw <- as.character(chr_vec)
    clean <- sub("^chr", "", raw, ignore.case = TRUE)
    numeric_val <- suppressWarnings(as.numeric(clean))
    is_num <- !is.na(numeric_val)
    special_rank <- rep(9999, length(clean))
    special_rank[toupper(clean) %in% c("X")] <- 1000
    special_rank[toupper(clean) %in% c("Y")] <- 1001
    special_rank[toupper(clean) %in% c("MT", "M", "MITO")] <- 1002
    sort_key1 <- ifelse(is_num, numeric_val, special_rank)
    sort_key2 <- ifelse(is_num, 0, as.numeric(factor(clean)))
    ord <- order(sort_key1, sort_key2, clean)
    unique(raw[ord])
  }

  get_mutation_colors_vec <- function(ref, alt) {
    mut <- paste0(ref, ">", alt)
    map <- c("C>A"="#4DBBD5","C>G"="#E64B35","C>T"="#00A087",
             "T>A"="#3C5488","T>C"="#F39B7F","T>G"="#8491B4")
    clr <- map[mut]; clr[is.na(clr)] <- "#999999"
    clr
  }
  
  geom_mean <- function(x) exp(mean(log(x)))
  
  # --- 1. Load and Initial Prep MAF ---
  maf_full <- data.table::fread(maf_file)
  maf_full <- maf_full[Tumor_Sample_Barcode == sample_id & Variant_Type == "SNP"]
  if (nrow(maf_full) == 0) stop("No SNPs found for that sample.")

  if (!all(c("Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2") %in% names(maf_full)))
    stop("MAF must have Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2 columns.")

  # --- DETERMINE PLOT MODE (Single vs. Multi-panel) ---
  is_facetted <- !is.null(chromosomes)
  
  # --- 2. Multi-Panel (Facetted) Mode ---
  if (is_facetted) {
    
    # Setup plotting environment for multiple panels
    n_chr <- length(chromosomes)
    maf_full <- maf_full[Chromosome %in% chromosomes]
    if (nrow(maf_full) == 0) stop("No SNPs found in the specified chromosomes.")

    nrow <- ceiling(n_chr / ncol)
    old_par <- par(no.readonly = TRUE) 
    on.exit(par(old_par)) 

    par(mfrow = c(nrow, ncol), 
        mar = c(3.5, 3.5, 1.5, 0.5), # Margins for each panel
        oma = c(0, 0, 2, 0))        # Outer margins

    # Loop over each selected chromosome
    for (i in 1:n_chr) {
      
      chr <- chromosomes[i]
      maf <- maf_full[Chromosome %in% chr]
      
      # --- PLOTTING LOGIC FOR ONE CHROMOSOME ---
      
      # Check for plotting sufficiency
      if (nrow(maf) <= 1) {
        x_max <- if (!is.null(chrom_sizes)) {
          size_dt <- data.table::as.data.table(chrom_sizes)
          size_dt[Chromosome %in% chr, Size][1]
        } else {
          1
        }
        plot(1, type="n", main=paste("Chr", chr), xaxt="n", yaxt="n", xlab="", ylab="", 
             xlim=c(0, x_max * 1.05))
        text(x_max * 0.5, 1, "Not enough SNPs")
        next
      }
      
      # Data prep specific to the current chromosome
      data.table::setorder(maf, Start_Position)
      maf[, dist_to_prev := c(NA, diff(Start_Position))]
      maf <- maf[!is.na(dist_to_prev)]
      maf[, dist_adj := dist_to_prev + 1] 

      # Colors and filtering
      fill_colors <- get_mutation_colors_vec(maf$Reference_Allele, maf$Tumor_Seq_Allele2)
      pt_border <- rep(NA_character_, length(fill_colors))

      if (!include_other_mutations) {
        other_color <- "#999999"
        not_other_idx <- which(fill_colors != other_color)
        maf <- maf[not_other_idx]
        fill_colors <- fill_colors[not_other_idx]
        pt_border <- pt_border[not_other_idx]
      }
      
      if (nrow(maf) == 0) {
        x_max <- if (!is.null(chrom_sizes)) {
          size_dt <- data.table::as.data.table(chrom_sizes)
          size_dt[Chromosome %in% chr, Size][1]
        } else {
          1
        }
        plot(1, type="n", main=paste("Chr", chr), xaxt="n", yaxt="n", xlab="", ylab="",
             xlim=c(0, x_max * 1.05))
        text(x_max * 0.5, 1, "No data after filtering")
        next
      }

      # Axis prep
      ylim_max <- max(maf$dist_adj, na.rm = TRUE)
      ylim_min <- min(maf$dist_adj[maf$dist_adj > 0], na.rm = TRUE)
      if (log_scale && ylim_min <= 1) ylim_min <- 1 
      
      # X-axis limits using chrom_sizes
      x_max <- if (!is.null(chrom_sizes)) {
        size_dt <- data.table::as.data.table(chrom_sizes)
        size_dt[Chromosome %in% chr, Size][1]
      } else {
        max(maf$Start_Position, na.rm = TRUE) 
      }
      x_max_plot <- x_max * 1.05

      # Empty plot
      plot(maf$Start_Position, maf$dist_adj,
           log = if (log_scale) "y" else "",
           type = "n",
           xlab = "Position (bp)",
           ylab = "Intermutation Distance (bp)",
           main = paste("Chr", chr),
           xaxt = "n",
           yaxt = "n",
           ylim = c(ylim_min, ylim_max),
           xlim = c(0, x_max_plot))

      # Axes and Lines
      axis(1, labels = TRUE, las = 1, cex.axis = cex_axis) 
      if (log_scale) {
        log_pos_valid <- 10^(0:ceiling(log10(ylim_max)))
        log_pos_valid <- log_pos_valid[log_pos_valid >= ylim_min & log_pos_valid <= ylim_max]
        labels <- formatC(log_pos_valid, format = "G")
        labels[log_pos_valid == 1] <- "1" 
        axis(2, at = log_pos_valid, labels = labels, las = 1, cex.axis = cex_axis)
        abline(h = 1001, col = "gray60", lty = 2, lwd = 1)
      } else {
        axis(2, las = 1, cex.axis = cex_axis)
        abline(h = 1000, col = "gray60", lty = 2, lwd = 1)
      }

      # Kataegis Arrows
      if (!is.null(highlight_kataegis)) {
        for (region in highlight_kataegis) {
          if (!is.list(region) || region$sample != sample_id || as.character(region$chromosome) != chr) next
          
          start_p <- region$start
          end_p <- region$end
          center_p <- (start_p + end_p) / 2
          in_idx <- which(maf$Start_Position >= start_p & maf$Start_Position <= end_p)
          
          if (length(in_idx) > 0) {
            pt_border[in_idx] <- "black"
            dist_in_region <- maf$dist_adj[in_idx]
            mid_imd <- geom_mean(dist_in_region)
            y_start <- ylim_min * (if (log_scale) 1.2 else 1.05) 
            y_end <- mid_imd 
            arrows(x0 = center_p, y0 = y_start, x1 = center_p, y1 = y_end, 
                   length = 0.08, angle = 30, code = 2, col = "red", lwd = 2) 
          }
        }
      }

      # Plot points
      points(maf$Start_Position, maf$dist_adj, pch = 21,
             bg = fill_colors, col = pt_border, cex = 0.7)

      # Legend (only in the first panel)
      show_leg <- show_legend && i == 1
      if (show_leg) {
        legend_labels <- c("C>A","C>G","C>T","T>A","T>C","T>G")
        legend_fills <- c("#4DBBD5","#E64B35","#00A087","#3C5488","#F39B7F","#8491B4")
        if (include_other_mutations) {
            legend_labels <- c(legend_labels, "Other")
            legend_fills <- c(legend_fills, "#999999")
        }
        legend("topright", legend = legend_labels, pt.bg = legend_fills, 
               col = "black", pch = 21, cex = 0.7, ncol = 1, bg = "white")
      }
    }
    
  # --- 3. Genome-Wide (Default) Mode ---
  } else {
    
    # Restore default par settings for single plot
    old_par <- par(no.readonly = TRUE) 
    on.exit(par(old_par)) 
    
    maf <- maf_full # Use full MAF
    
    # Chromosome ordering
    if (!is.null(chrom_sizes)) {
      cs <- data.table::as.data.table(chrom_sizes)
      desired_order <- order_chromosomes(cs$Chromosome)
      chr_levels <- desired_order
    } else {
      chr_levels <- order_chromosomes(unique(maf$Chromosome))
    }
    maf[, Chromosome := factor(as.character(Chromosome), levels = chr_levels)]

    # Compute cumulative positions
    if (!is.null(chrom_sizes)) {
      cs <- data.table::as.data.table(chrom_sizes)[, .(Chromosome, Size)]
      cs <- cs[order(match(Chromosome, chr_levels))]
      missing_chr <- setdiff(chr_levels, cs$Chromosome)
      if (length(missing_chr) > 0)
        cs <- rbind(cs, data.table::data.table(Chromosome = missing_chr, Size = 0))
      cs <- cs[match(chr_levels, Chromosome)]
      cs[, cum_start := data.table::shift(cumsum(as.numeric(Size)), fill = 0)]
      chr_lengths <- cs[, .(Chromosome, Size, cum_start)]
    } else {
      tmp <- maf[, .(Size = max(Start_Position, na.rm = TRUE)), by = Chromosome]
      tmp <- tmp[match(chr_levels, tmp$Chromosome)]
      tmp$Size[is.na(tmp$Size)] <- 0
      tmp[, cum_start := data.table::shift(cumsum(as.numeric(Size)), fill = 0)]
      chr_lengths <- tmp[, .(Chromosome, Size, cum_start)]
    }

    maf <- merge(maf, chr_lengths[, .(Chromosome, cum_start)], by = "Chromosome", all.x = TRUE, sort = FALSE)
    maf[, genome_pos := Start_Position + cum_start]
    data.table::setorder(maf, genome_pos)

    # Intermutation distance
    maf[, dist_to_prev := c(NA, diff(genome_pos))]
    maf <- maf[!is.na(dist_to_prev)]
    maf[, dist_adj := dist_to_prev + 1] 

    # Colors and borders
    fill_colors <- get_mutation_colors_vec(maf$Reference_Allele, maf$Tumor_Seq_Allele2)
    pt_border <- rep(NA_character_, length(fill_colors))
      
    # Axis prep
    chr_mid <- chr_lengths[, .(Chromosome, midpoint = cum_start + Size / 2)]
    ylim_max <- max(maf$dist_adj, na.rm = TRUE)
    ylim_min <- min(maf$dist_adj[maf$dist_adj > 0], na.rm = TRUE)
    if (log_scale && ylim_min <= 1) ylim_min <- 1 

    # Empty plot
    plot(maf$genome_pos, maf$dist_adj,
          log = if (log_scale) "y" else "",
          type = "n",
          xlab = "",
          ylab = if (log_scale) "Intermutation Distance (bp)" else "Intermutation Distance (bp)",
          main = paste("Genome-wide Rainfall Plot:", sample_id),
          xaxt = "n",
          yaxt = "n", 
          ylim = c(ylim_min, ylim_max))

    # 1. CUSTOM Y-AXIS (Log10 Labels: 1, 10, 100, ...)
    if (log_scale) {
      log_pos_valid <- 10^(0:ceiling(log10(ylim_max)))
      log_pos_valid <- log_pos_valid[log_pos_valid >= ylim_min & log_pos_valid <= ylim_max]
      labels <- formatC(log_pos_valid, format = "G")
      labels[log_pos_valid == 1] <- "1" 
      axis(2, at = log_pos_valid, labels = labels, las = 1, cex.axis = cex_axis)
      abline(h = 1001, col = "gray60", lty = 2, lwd = 1)
    } else {
      axis(2, las = 1, cex.axis = cex_axis)
      abline(h = 1000, col = "gray60", lty = 2, lwd = 1)
    }
    
    # X-Axis and Dividers
    axis(1, at = chr_mid$midpoint, labels = chr_mid$Chromosome, las = 2, cex.axis = cex_axis)
    vlines <- unique(chr_lengths$cum_start)
    if (length(vlines) > 1) abline(v = vlines[-1], col = "gray85", lty = 3)

    # --- Add kataegis up-arrows ---
    if (!is.null(highlight_kataegis)) {
      for (region in highlight_kataegis) {
        if (!is.list(region) || is.null(region$sample) || region$sample != sample_id) next
        chrname <- as.character(region$chromosome)
        if (!(chrname %in% chr_lengths$Chromosome)) next

        chr_cum <- chr_lengths[Chromosome %in% chrname]$cum_start
        start_g <- region$start + chr_cum
        end_g <- region$end + chr_cum
        center_g <- (start_g + end_g) / 2

        in_idx <- which(maf$genome_pos >= start_g & maf$genome_pos <= end_g)
        
        if (length(in_idx) > 0) {
          pt_border[in_idx] <- "black"
          dist_in_region <- maf$dist_adj[in_idx]
          mid_imd <- geom_mean(dist_in_region)

          y_start <- ylim_min * (if (log_scale) 1.01 else 1.001)
          y_end <- mid_imd
          
          arrows(x0 = center_g, y0 = y_start, x1 = center_g, y1 = y_end,
                 length = 0.1, angle = 30, code = 2, col = "black", lwd = 0.75)
        }
      }
    }

    # --- Plot points (Genome-wide) ---
    if (!include_other_mutations) {
      other_color <- "#999999"
      not_other_idx <- which(fill_colors != other_color)
      plot_maf <- maf[not_other_idx]
      plot_fill_colors <- fill_colors[not_other_idx]
      plot_pt_border <- pt_border[not_other_idx]
    } else {
      plot_maf <- maf
      plot_fill_colors <- fill_colors
      plot_pt_border <- pt_border
    }
    
    points(plot_maf$genome_pos, plot_maf$dist_adj, pch = 21,
            bg = plot_fill_colors, col = plot_pt_border, cex = 0.7)

    # --- Legend (Genome-wide) ---
    if (show_legend) {
      legend_labels <- c("C>A","C>G","C>T","T>A","T>C","T>G")
      legend_fills <- c("#4DBBD5","#E64B35","#00A087","#3C5488","#F39B7F","#8491B4")
      if (include_other_mutations) {
        legend_labels <- c(legend_labels, "Other")
        legend_fills <- c(legend_fills, "#999999")
      }
      legend("topright", legend = legend_labels, pt.bg = legend_fills, col = "black",
             pch = 21, cex = 0.7, ncol = 1, bg = "white")
    }
  }
}

##########
plot_rainfall <- function(maf_file, sample_id,
                                     chromosomes = NULL, 
                                     plot_mode = c("genome", "combined", "facetted"), # NEW CONTROL
                                     log_scale = TRUE,
                                     highlight_kataegis = NULL,
                                     ncol = 3,
                                     show_legend = TRUE,
                                     chrom_sizes = NULL,
                                     cex_axis = 0.65, 
                                     include_other_mutations = TRUE) {

  plot_mode <- match.arg(plot_mode) # Validate and match the plot_mode argument
  
  # --- Helper functions ---
  order_chromosomes <- function(chr_vec) {
    raw <- as.character(chr_vec)
    clean <- sub("^chr", "", raw, ignore.case = TRUE)
    numeric_val <- suppressWarnings(as.numeric(clean))
    is_num <- !is.na(numeric_val)
    special_rank <- rep(9999, length(clean))
    special_rank[toupper(clean) %in% c("X")] <- 1000
    special_rank[toupper(clean) %in% c("Y")] <- 1001
    special_rank[toupper(clean) %in% c("MT", "M", "MITO")] <- 1002
    sort_key1 <- ifelse(is_num, numeric_val, special_rank)
    sort_key2 <- ifelse(is_num, 0, as.numeric(factor(clean)))
    ord <- order(sort_key1, sort_key2, clean)
    unique(raw[ord])
  }

  get_mutation_colors_vec <- function(ref, alt) {
    mut <- paste0(ref, ">", alt)
    map <- c("C>A"="#4DBBD5","C>G"="#E64B35","C>T"="#00A087",
             "T>A"="#3C5488","T>C"="#F39B7F","T>G"="#8491B4")
    clr <- map[mut]; clr[is.na(clr)] <- "#999999"
    clr
  }
  
  geom_mean <- function(x) exp(mean(log(x)))
  
  # --- 1. Load and Initial Prep MAF ---
  maf_full <- data.table::fread(maf_file)
  # Handle NULL or "All" sample_id
  if (is.null(sample_id) || length(sample_id) == 0) {
    # If multiple samples present, pick the first one or stop
    available_samples <- unique(maf_full$Tumor_Sample_Barcode)
    if (length(available_samples) == 1) {
      sample_id <- available_samples[1]
    } else {
      stop(paste("Multiple samples in MAF. Please select one. Available:",
                paste(available_samples, collapse = ", ")))
    }
  }

  maf_full <- maf_full[Tumor_Sample_Barcode %in% sample_id & Variant_Type == "SNP"]

  if (nrow(maf_full) == 0) stop("No SNPs found for that sample.")

  if (!all(c("Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2") %in% names(maf_full)))
    stop("MAF must have Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2 columns.")

  # Save/Restore default par settings
  old_par <- par(no.readonly = TRUE) 
  on.exit(par(old_par)) 
  
  # --- 2. Facetted (Multi-Panel) Mode ---
  if (plot_mode == "facetted") {
      
      if (is.null(chromosomes) || length(chromosomes) < 1) {
          stop("When plot_mode='facetted', a list of chromosomes must be supplied.")
      }

      maf_full_subset <- maf_full[Chromosome %in% chromosomes]
      n_chr <- length(chromosomes)
      if (nrow(maf_full_subset) == 0) stop("No SNPs found in the specified chromosomes.")

      nrow <- ceiling(n_chr / ncol)
      par(mfrow = c(nrow, ncol), 
          mar = c(3.5, 3.5, 1.5, 0.5), # Margins for each panel
          oma = c(0, 0, 2, 0))        # Outer margins

      # Loop over each selected chromosome
      for (i in 1:n_chr) {
        
        chr <- chromosomes[i]
        maf <- maf_full_subset[Chromosome %in% chr]
        
        # --- PLOTTING LOGIC FOR ONE CHROMOSOME ---
        
        # Check for plotting sufficiency
        if (nrow(maf) <= 1) {
          x_max <- if (!is.null(chrom_sizes)) {
            size_dt <- data.table::as.data.table(chrom_sizes)
            size_dt[Chromosome %in% chr, Size][1]
          } else {
            1
          }
          plot(1, type="n", main=paste("Chr", chr), xaxt="n", yaxt="n", xlab="", ylab="", 
               xlim=c(0, x_max * 1.05))
          text(x_max * 0.5, 1, "Not enough SNPs")
          next
        }
        
        # Data prep specific to the current chromosome
        data.table::setorder(maf, Start_Position)
        maf[, dist_to_prev := c(NA, diff(Start_Position))]
        maf <- maf[!is.na(dist_to_prev)]
        maf[, dist_adj := dist_to_prev + 1] 

        # Colors and filtering
        fill_colors <- get_mutation_colors_vec(maf$Reference_Allele, maf$Tumor_Seq_Allele2)
        pt_border <- rep(NA_character_, length(fill_colors))

        if (!include_other_mutations) {
          other_color <- "#999999"
          not_other_idx <- which(fill_colors != other_color)
          maf <- maf[not_other_idx]
          fill_colors <- fill_colors[not_other_idx]
          pt_border <- pt_border[not_other_idx]
        }
        
        if (nrow(maf) == 0) {
          x_max <- if (!is.null(chrom_sizes)) {
            size_dt <- data.table::as.data.table(chrom_sizes)
            size_dt[Chromosome %in% chr, Size][1]
          } else {
            1
          }
          plot(1, type="n", main=paste("Chr", chr), xaxt="n", yaxt="n", xlab="", ylab="",
               xlim=c(0, x_max * 1.05))
          text(x_max * 0.5, 1, "No data after filtering")
          next
        }

        # Axis prep
        ylim_max <- max(maf$dist_adj, na.rm = TRUE)
        ylim_min <- min(maf$dist_adj[maf$dist_adj > 0], na.rm = TRUE)
        if (log_scale && ylim_min <= 1) ylim_min <- 1 
        
        # X-axis limits using chrom_sizes
        x_max <- if (!is.null(chrom_sizes)) {
          size_dt <- data.table::as.data.table(chrom_sizes)
          size_dt[Chromosome %in% chr, Size][1]
        } else {
          max(maf$Start_Position, na.rm = TRUE) 
        }
        x_max_plot <- x_max * 1.05

        # Empty plot
        plot(maf$Start_Position, maf$dist_adj,
             log = if (log_scale) "y" else "",
             type = "n",
             xlab = "Position (bp)",
             ylab = "Intermutation Distance (bp)",
             main = paste("Chr", chr),
             xaxt = "n",
             yaxt = "n",
             ylim = c(ylim_min, ylim_max),
             xlim = c(0, x_max_plot),
             xaxs = "i")

        # Grid Lines
        x_ticks <- axTicks(1)
        abline(v = x_ticks, col = "gray88", lty = 3, lwd = 0.5)
        
        if (log_scale) {
            y_minor <- 2^(0:ceiling(log2(ylim_max)))
            abline(h = y_minor, col = "gray88", lty = 3, lwd = 0.5)
        }

        # X-axis
        axis(1, labels = TRUE, las = 1, cex.axis = cex_axis) 

        # Custom Log Base 2 Y-Axis
        if (log_scale) {
          log_pos_valid <- 2^(0:ceiling(log2(ylim_max)))
          log_pos_valid <- log_pos_valid[log_pos_valid >= ylim_min & log_pos_valid <= ylim_max]
          labels <- formatC(log_pos_valid, format = "G")
          labels[log_pos_valid == 1] <- "1" 
          axis(2, at = log_pos_valid, labels = labels, las = 1, cex.axis = cex_axis)
          abline(h = 1001, col = "gray60", lty = 2, lwd = 1)
        } else {
          axis(2, las = 1, cex.axis = cex_axis)
          abline(h = 1000, col = "gray60", lty = 2, lwd = 1)
        }

        # Kataegis Arrows
        if (!is.null(highlight_kataegis)) {
          for (region in highlight_kataegis) {
            if (!is.list(region) || !isTRUE(region$sample == sample_id) || !isTRUE(as.character(region$chromosome) == chr)) next
            
            start_p <- region$start
            end_p <- region$end
            center_p <- (start_p + end_p) / 2
            in_idx <- which(maf$Start_Position >= start_p & maf$Start_Position <= end_p)
            
            if (length(in_idx) > 0) {
              pt_border[in_idx] <- "black"
              dist_in_region <- maf$dist_adj[in_idx]
              mid_imd <- geom_mean(dist_in_region)
              y_start <- ylim_min * (if (log_scale) 1.2 else 1.05) 
              y_end <- mid_imd 
              arrows(x0 = center_p, y0 = y_start, x1 = center_p, y1 = y_end, 
                     length = 0.08, angle = 30, code = 2, col = "red", lwd = 2) 
            }
          }
        }

        # Plot points
        points(maf$Start_Position, maf$dist_adj, pch = 21,
               bg = fill_colors, col = pt_border, cex = 0.7)

        # Legend (only in the first panel)
        show_leg <- show_legend && i == 1
        if (show_leg) {
          legend_labels <- c("C>A","C>G","C>T","T>A","T>C","T>G")
          legend_fills <- c("#4DBBD5","#E64B35","#00A087","#3C5488","#F39B7F","#8491B4")
          if (include_other_mutations) {
              legend_labels <- c(legend_labels, "Other")
              legend_fills <- c(legend_fills, "#999999")
          }
          legend("topright", legend = legend_labels, pt.bg = legend_fills, 
                 col = "black", pch = 21, cex = 0.7, ncol = 1, bg = "white")
        }
      }
    
  # --- 3. Genome-Wide or Combined (Single Plot Area) ---
  } else {
    
    # 1. Select the MAF subset
    if (plot_mode == "combined") {
        if (is.null(chromosomes) || length(chromosomes) < 1) {
            stop("When plot_mode='combined', a list of chromosomes must be supplied.")
        }
        maf <- maf_full[Chromosome %in% chromosomes]
        if (nrow(maf) == 0) stop("No SNPs found in the specified chromosomes for combined plot.")
    } else { # plot_mode == "genome"
        maf <- maf_full
    }
    
    # Chromosome ordering
    if (!is.null(chrom_sizes)) {
      cs <- data.table::as.data.table(chrom_sizes)
      # Order only the chromosomes present in the data/requested list
      chr_levels_data <- unique(maf$Chromosome)
      desired_order <- order_chromosomes(chr_levels_data)
      chr_levels <- desired_order[desired_order %in% chr_levels_data]
    } else {
      chr_levels <- order_chromosomes(unique(maf$Chromosome))
    }
    maf[, Chromosome := factor(as.character(Chromosome), levels = chr_levels)]

    # Compute cumulative positions
    if (!is.null(chrom_sizes)) {
      cs <- data.table::as.data.table(chrom_sizes)[, .(Chromosome, Size)]
      cs <- cs[Chromosome %in% chr_levels] # Subset to relevant chromosomes
      cs <- cs[order(match(Chromosome, chr_levels))]
      
      # Handle chromosomes missing from the MAF but present in chrom_sizes (less relevant here)
      missing_chr <- setdiff(chr_levels, cs$Chromosome)
      if (length(missing_chr) > 0)
        cs <- rbind(cs, data.table::data.table(Chromosome = missing_chr, Size = 0))
        
      cs <- cs[match(chr_levels, Chromosome)]
      cs[, cum_start := data.table::shift(cumsum(as.numeric(Size)), fill = 0)]
      chr_lengths <- cs[, .(Chromosome, Size, cum_start)]
    } else {
      tmp <- maf[, .(Size = max(Start_Position, na.rm = TRUE)), by = Chromosome]
      tmp <- tmp[match(chr_levels, tmp$Chromosome)]
      tmp$Size[is.na(tmp$Size)] <- 0
      tmp[, cum_start := data.table::shift(cumsum(as.numeric(Size)), fill = 0)]
      chr_lengths <- tmp[, .(Chromosome, Size, cum_start)]
    }

    maf <- merge(maf, chr_lengths[, .(Chromosome, cum_start)], by = "Chromosome", all.x = TRUE, sort = FALSE)
    maf[, genome_pos := Start_Position + cum_start]
    data.table::setorder(maf, genome_pos)

    # Intermutation distance
    maf[, dist_to_prev := c(NA, diff(genome_pos))]
    maf <- maf[!is.na(dist_to_prev)]
    maf[, dist_adj := dist_to_prev + 1] 

    # Colors and borders
    fill_colors <- get_mutation_colors_vec(maf$Reference_Allele, maf$Tumor_Seq_Allele2)
    pt_border <- rep(NA_character_, length(fill_colors))
      
    # Axis prep
    chr_mid <- chr_lengths[, .(Chromosome, midpoint = cum_start + Size / 2)]
    ylim_max <- max(maf$dist_adj, na.rm = TRUE)
    ylim_min <- min(maf$dist_adj[maf$dist_adj > 0], na.rm = TRUE)
    if (log_scale && ylim_min <= 1) ylim_min <- 1 

    # Empty plot
    plot(maf$genome_pos, maf$dist_adj,
          log = if (log_scale) "y" else "",
          type = "n",
          xlab = "",
          ylab = if (log_scale) "Intermutation Distance (bp)" else "Intermutation Distance (bp)",
          main = paste(ifelse(plot_mode == "combined", "Combined Chromosome", "Genome-wide"), 
                       "Rainfall Plot:", sample_id),
          xaxt = "n",
          yaxt = "n", 
          ylim = c(ylim_min, ylim_max),
          xaxs = "i") 

    # --- Grid Lines (Genome-wide/Combined) ---
    vlines <- unique(chr_lengths$cum_start)
    if (length(vlines) > 1) abline(v = vlines, col = "gray85", lty = 3) # Chromosome dividers
    
    if (log_scale) {
      y_minor <- 2^(0:ceiling(log2(ylim_max)))
      abline(h = y_minor, col = "gray88", lty = 3, lwd = 0.5)
    }

    # --- Custom Log Base 2 Y-Axis (Genome-wide/Combined) ---
    if (log_scale) {
      log_pos_valid <- 2^(0:ceiling(log2(ylim_max)))
      log_pos_valid <- log_pos_valid[log_pos_valid >= ylim_min & log_pos_valid <= ylim_max]
      
      labels <- formatC(log_pos_valid, format = "G")
      labels[log_pos_valid == 1] <- "1" 
      
      axis(2, at = log_pos_valid, labels = labels, las = 1, cex.axis = cex_axis)
      abline(h = 1001, col = "gray60", lty = 2, lwd = 1) 
    } else {
      axis(2, las = 1, cex.axis = cex_axis)
      abline(h = 1000, col = "gray60", lty = 2, lwd = 1)
    }
    
    # X-Axis Labels (Chromosome midpoints)
    axis(1, at = chr_mid$midpoint, labels = chr_mid$Chromosome, las = 2, cex.axis = cex_axis)

    # --- Add kataegis up-arrows ---
    if (!is.null(highlight_kataegis)) {
      for (region in highlight_kataegis) {
        if (!is.list(region) || is.null(region$sample) || region$sample != sample_id) next
        chrname <- as.character(region$chromosome)
        if (!(chrname %in% chr_lengths$Chromosome)) next

        chr_cum <- chr_lengths[Chromosome %in% chrname]$cum_start
        start_g <- region$start + chr_cum
        end_g <- region$end + chr_cum
        center_g <- (start_g + end_g) / 2

        in_idx <- which(maf$genome_pos >= start_g & maf$genome_pos <= end_g)
        
        if (length(in_idx) > 0) {
          pt_border[in_idx] <- "black"
          dist_in_region <- maf$dist_adj[in_idx]
          mid_imd <- geom_mean(dist_in_region)

          y_start <- ylim_min * (if (log_scale) 1.2 else 1.05)
          y_end <- mid_imd
          
          arrows(x0 = center_g, y0 = y_start, x1 = center_g, y1 = y_end,
                 length = 0.08, angle = 30, code = 2, col = "black", lwd = 0.75)
        }
      }
    }

    # --- Plot points ---
    if (!include_other_mutations) {
      other_color <- "#999999"
      not_other_idx <- which(fill_colors != other_color)
      plot_maf <- maf[not_other_idx]
      plot_fill_colors <- fill_colors[not_other_idx]
      plot_pt_border <- pt_border[not_other_idx]
    } else {
      plot_maf <- maf
      plot_fill_colors <- fill_colors
      plot_pt_border <- pt_border
    }
    
    points(plot_maf$genome_pos, plot_maf$dist_adj, pch = 21,
            bg = plot_fill_colors, col = plot_pt_border, cex = 0.7)

    # --- Legend ---
    if (show_legend) {
      legend_labels <- c("C>A","C>G","C>T","T>A","T>C","T>G")
      legend_fills <- c("#4DBBD5","#E64B35","#00A087","#3C5488","#F39B7F","#8491B4")
      if (include_other_mutations) {
        legend_labels <- c(legend_labels, "Other")
        legend_fills <- c(legend_fills, "#999999")
      }
      legend("topright", legend = legend_labels, pt.bg = legend_fills, col = "black",
             pch = 21, cex = 0.7, ncol = 1, bg = "white")
    }
  }
}
