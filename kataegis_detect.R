library(data.table)
library(dplyr)

# ==============================================================
# Script: identify_kataegis_sites.R
# Purpose: identify kataegis sites in a species agnostic manner.
# Inspiration: maftools detect_kataegis/changepoints approach, and katdetectr's IMD function
# Author: Harish Kothandaraman
# Date: 2026-06-11
# ==============================================================

# Function to calculate average intermutation distance
calc_avg_intermut_dist <- function(positions) {
  if (length(positions) < 2) return(Inf)
  distances <- diff(positions)
  return(mean(distances))
}

# Function to model sample mutation rate
modelSampleRate <- function(IMDs) {
  lambda <- log(2) / median(IMDs)
  return(lambda)
}

# Function for calculating the nth root of x
nthroot <- function(x, n) {
  y <- x^(1 / n)
  return(y)
}

# Function that defines the IMD cutoff specific for each chromosome
IMDcutoffFun <- function(chr_maf, chr_true_length = NULL) {
  chr_maf_sorted <- chr_maf[order(Start_Position)]
  IMDs <- diff(chr_maf_sorted$Start_Position)
  if (length(IMDs) == 0) return(1000)

  sampleRate <- modelSampleRate(IMDs)

  # Use true chromosome length if provided, else fall back to observed span
  width <- if (!is.null(chr_true_length) && chr_true_length > 0) {
    chr_true_length
  } else {
    max(chr_maf_sorted$Start_Position) - min(chr_maf_sorted$Start_Position)
  }                                          # FIX: was missing closing }

  totalVariants <- nrow(chr_maf_sorted)
  IMDcutoff <- -log(1 - nthroot(0.01 / width,
    ifelse(totalVariants != 0, totalVariants - 1, 1))) / sampleRate
  return(min(IMDcutoff, 1000))
}

# Function to detect kataegis from MAF file using deque algorithm
detect_kataegis <- function(maf_file,
                            sample_id        = NULL,
                            min_mutations    = 6,
                            max_avg_distance = 1000,
                            use_dynamic_cutoff = FALSE,
                            chrom_sizes      = NULL) {

  maf <- fread(maf_file)

  if (!is.null(sample_id)) {
    maf <- maf[Tumor_Sample_Barcode == sample_id]
  }

  maf <- maf[Variant_Type == "SNP",
             .(Chromosome, Start_Position, End_Position,
               Reference_Allele, Tumor_Seq_Allele2,
               Tumor_Sample_Barcode, Variant_Classification)]

  setorder(maf, Chromosome, Start_Position)

  kataegis_regions <- list()
  region_id <- 1

  for (samp in unique(maf$Tumor_Sample_Barcode)) {
    sample_maf <- maf[Tumor_Sample_Barcode == samp]

    for (chr in unique(sample_maf$Chromosome)) {
      chr_maf <- sample_maf[Chromosome %in% chr]

      if (nrow(chr_maf) < min_mutations) next

      if (use_dynamic_cutoff) {
        true_len <- NULL
        if (!is.null(chrom_sizes)) {
          cs  <- data.table::as.data.table(chrom_sizes)
          idx <- which(cs$Chromosome %in% chr)
          if (length(idx) > 0) true_len <- cs$Size[idx[1]]
        }
        chr_cutoff <- IMDcutoffFun(chr_maf, chr_true_length = true_len)
        cat("Sample:", samp, "Chromosome:", chr,
            "Dynamic cutoff:", round(chr_cutoff, 2), "bp\n")
      } else {
        chr_cutoff <- max_avg_distance
      }

      # Goretzky deque algorithm
      n_muts    <- nrow(chr_maf)
      start_idx <- 1
      end_idx   <- min_mutations

      while (end_idx <= n_muts) {
        queue_indices   <- start_idx:end_idx
        queue_positions <- chr_maf$Start_Position[queue_indices]
        avg_dist        <- calc_avg_intermut_dist(queue_positions)

        if (avg_dist > chr_cutoff) {
          # Slide window forward
          start_idx <- start_idx + 1
          end_idx   <- end_idx   + 1
        } else {
          # Extend cluster while avg IMD stays within cutoff
          while (avg_dist <= chr_cutoff && end_idx <= n_muts) {
            end_idx <- end_idx + 1
            if (end_idx <= n_muts) {
              queue_indices   <- start_idx:end_idx
              queue_positions <- chr_maf$Start_Position[queue_indices]
              avg_dist        <- calc_avg_intermut_dist(queue_positions)
            }
          }

          # Save region (exclude last mutation that broke threshold)
          final_indices <- start_idx:(end_idx - 1)
          cluster_muts  <- chr_maf[final_indices]

          cluster_muts[, dist_to_next := c(diff(Start_Position), NA)]

          start_pos      <- min(cluster_muts$Start_Position)
          end_pos        <- max(cluster_muts$Start_Position)
          span           <- end_pos - start_pos
          final_avg_dist <- calc_avg_intermut_dist(cluster_muts$Start_Position)

          cluster_muts[, mut_type := paste0(Reference_Allele, ">", Tumor_Seq_Allele2)]
          mut_counts <- table(cluster_muts$mut_type)

          kataegis_regions[[region_id]] <- list(
            region_id         = region_id,
            sample            = samp,
            chromosome        = chr,
            start             = start_pos,
            end               = end_pos,
            n_mutations       = nrow(cluster_muts),
            span_bp           = span,
            avg_intermut_dist = final_avg_dist,
            imd_cutoff_used   = chr_cutoff,
            mutations         = cluster_muts,
            mutation_spectrum = mut_counts
          )

          region_id <- region_id + 1

          # Reinitialize queue (allows potential overlap)
          start_idx <- end_idx
          end_idx   <- start_idx + min_mutations
        }
      }
    }
  }

  return(kataegis_regions)
}

# Function to create detailed dataframe with mutation spectrum
kataegis_to_dataframe <- function(kataegis_list, include_spectrum = TRUE, group_other = TRUE) {
  if (length(kataegis_list) == 0) {
    cat("No kataegis regions detected.\n")
    return(NULL)
  }

  df <- data.frame(                              # FIX: closing ) was missing
    region_id         = sapply(kataegis_list, function(x) x$region_id),
    sample            = sapply(kataegis_list, function(x) x$sample),
    chromosome        = sapply(kataegis_list, function(x) x$chromosome),
    start             = sapply(kataegis_list, function(x) x$start),
    end               = sapply(kataegis_list, function(x) x$end),
    n_mutations       = sapply(kataegis_list, function(x) x$n_mutations),
    span_bp           = sapply(kataegis_list, function(x) x$span_bp),
    avg_intermut_dist = sapply(kataegis_list, function(x) x$avg_intermut_dist)
  )

  if (include_spectrum) {
    main_mut_types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    all_mut_types  <- unique(unlist(lapply(kataegis_list, function(x) names(x$mutation_spectrum))))

    if (group_other) {
      for (mut_type in main_mut_types) {
        df[[mut_type]] <- sapply(kataegis_list, function(x) {
          if (mut_type %in% names(x$mutation_spectrum))
            return(as.numeric(x$mutation_spectrum[mut_type]))
          else
            return(0)
        })
      }

      df$Other <- sapply(kataegis_list, function(x) {
        sum_main <- sum(sapply(main_mut_types, function(m) {
          if (m %in% names(x$mutation_spectrum))
            return(unname(as.numeric(x$mutation_spectrum[m])))
          else
            return(0)
        }))
        return(x$n_mutations - sum_main)
      })

    } else {
      for (mut_type in all_mut_types) {
        df[[mut_type]] <- sapply(kataegis_list, function(x) {
          if (mut_type %in% names(x$mutation_spectrum))
            return(as.numeric(x$mutation_spectrum[mut_type]))
          else
            return(0)
        })
      }
    }

    if ("C>T" %in% names(df))
      df$pct_C_to_T <- round((df[["C>T"]] / df$n_mutations) * 100, 2)

    if ("C>G" %in% names(df))
      df$pct_C_to_G <- round((df[["C>G"]] / df$n_mutations) * 100, 2)
  }

  return(df)
}

# Function to get detailed mutations for a specific kataegis region
get_kataegis_mutations <- function(kataegis_list, region_id) {
  # FIX: use match() so region_id lookup is robust after subsetting
  pos <- match(region_id, sapply(kataegis_list, function(x) x$region_id))
  if (is.na(pos)) {
    cat("Region", region_id, "not found.\n")
    return(NULL)
  }
  region <- kataegis_list[[pos]]

  cat("Kataegis Region", region_id, "\n")
  cat("Sample:", region$sample, "\n")
  cat("Location:", region$chromosome, ":", region$start, "-", region$end, "\n")
  cat("Number of mutations:", region$n_mutations, "\n")
  cat("Average intermutation distance:", round(region$avg_intermut_dist, 2), "bp\n")
  cat("\nMutation Spectrum:\n")
  print(region$mutation_spectrum)
  cat("\n")

  return(as.data.frame(region$mutations))
}

# Function to export kataegis regions as BED format (0-based)
export_kataegis_bed <- function(kataegis_list, output_file) {
  if (length(kataegis_list) == 0) {
    cat("No kataegis regions to export.\n")
    return(NULL)
  }

  bed_df <- data.frame(
    chr   = sapply(kataegis_list, function(x) x$chromosome),
    start = sapply(kataegis_list, function(x) x$start - 1),  # BED is 0-based
    end   = sapply(kataegis_list, function(x) x$end),
    name  = sapply(kataegis_list, function(x) paste0(x$sample, "_region", x$region_id)),
    score = sapply(kataegis_list, function(x) x$n_mutations)
  )

  write.table(bed_df, output_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)

  cat("Exported", nrow(bed_df), "regions to", output_file, "\n")
}

# Function to subset detect_kataegis results by a BED file
subset_kataegis_by_bed <- function(kataegis_list, bed_file, overlap_type = "any") {
  if (length(kataegis_list) == 0) {
    cat("No kataegis regions to subset.\n")
    return(list())
  }

  bed <- fread(bed_file, header = FALSE)
  colnames(bed)[1:3] <- c("Chromosome", "Start_Position", "End_Position")
  bed[, Start_Position := Start_Position + 1]  # convert 0-based to 1-based

  kat_dt <- data.table(
    region_id      = sapply(kataegis_list, function(x) x$region_id),
    sample         = sapply(kataegis_list, function(x) x$sample),
    Chromosome     = sapply(kataegis_list, function(x) x$chromosome),
    Start_Position = sapply(kataegis_list, function(x) x$start),
    End_Position   = sapply(kataegis_list, function(x) x$end)
  )

  setkey(kat_dt, Chromosome, Start_Position, End_Position)
  setkey(bed,    Chromosome, Start_Position, End_Position)

  if (overlap_type == "any") {
    overlap <- foverlaps(kat_dt, bed, type = "any",    nomatch = NULL)
  } else if (overlap_type == "complete") {
    overlap <- foverlaps(kat_dt, bed, type = "within", nomatch = NULL)
  } else if (overlap_type == "start") {
    overlap <- foverlaps(kat_dt, bed, type = "any",    nomatch = NULL)
    overlap <- overlap[Start_Position >= i.Start_Position &
                       Start_Position <= i.End_Position]
  } else {
    stop("overlap_type must be \'any\', \'complete\', or \'start\'")
  }

  matched_ids <- unique(overlap$region_id)

  positions   <- match(matched_ids,
                       sapply(kataegis_list, function(x) x$region_id))
  subset_list <- kataegis_list[positions[!is.na(positions)]]

  cat("Found", length(subset_list), "kataegis regions overlapping with BED file\n")
  return(subset_list)
}
