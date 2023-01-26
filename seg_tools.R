# LIBRARIES ====================================================================
library('R.utils', quietly = TRUE)
library('magrittr', quietly = TRUE)
library('tidyverse', quietly = TRUE)
if (!require('BiocManager', quietly = TRUE)) {
    install.packages('BiocManager', quietly = TRUE)
}
if (!require('DNAcopy', quietly = TRUE)) {
    BiocManager::install('DNAcopy', quietly = TRUE)
}
if (!require('argparser', quietly = TRUE)) {
    install.packages('argparser', quietly = TRUE)
}
if (!require('IRanges', quietly = TRUE)) {
    BiocManager::install('IRanges', quietly = TRUE)
}
library('IRanges', quietly = TRUE)
library('argparser', quietly = TRUE)
library('DNAcopy', quietly = TRUE)
library('parallel', quietly = TRUE)

# FUNCTIONS ====================================================================

# SEGMETATION ------------------------------------------------------------------

# Get metadata
#   - A simple function to pull metadata (data set name, indexcov vs. mosdepth)
get_meta <- function(data_path) {
    set_name <- str_extract(data_path, "NWD[0-9]*")
    if (grepl('indexcov', data_path)) {
        type <- 'indexcov'
    } else if (grepl('mosdepth', data_path)) {
        type <- 'mosdepth'
    }
    return(list(
        set_name = set_name,
        type = type
    ))
}

# Load data
#   - Loads data from .bed.gz file. Interprets input file name to determine if
#     loading indexcov vs. mosdepth data. Untars, then gunzips file to a temp
#     directory, then fixes column names and returns data as a tibble.
#   - Inputs
#       - data_path: path to data file in .tar.gz format
#   - Outputs
#       - set_data: tibble of raw BED coverage data
load <- function(data_path) {
    set_name <- str_extract(data_path, "NWD[0-9]*")
    temp_dir <- paste0(tempdir(), '/', set_name, '/')
    if (grepl('indexcov', data_path)) {
        bed_file <- 'indexcov/indexcov-indexcov.bed.gz'
        header <- TRUE
    } else if (grepl('mosdepth', data_path)) {
        bed_file <- 'mosdepth/cnv.regions.bed.gz'
        header <- FALSE
    }
    set_file <- paste0(temp_dir, '/', str_sub(bed_file, 1, -4))
    untar(
        data_path,
        files = bed_file,
        exdir = temp_dir
        )
    gunzip(paste0(temp_dir, '/', bed_file), overwrite = TRUE)
    set_data <- read_tsv(
        set_file,
        col_names = header,
        col_types = 'ciid'
    )
    colnames(set_data) <- c('chromosome', 'start', 'end', 'depth')
    set_data$dataset <- set_name
    unlink(temp_dir)
    return(set_data)
}

# Filter by a specific gene
#   - Reads subset_expanded.bed file to pull chromosome or gene coordinates of
#       gene of interest, then filters set_data for the gene or coordinates.
#   - Inputs
#       - set_data: tibble, raw set data loaded from mosdepth or indexcov output data
#       - gene: string, gene name e.g. DNMT3A, can be lowercase
#       - meta_data: list, output from get_meta()
#   - Outputs
#       - set_data: tibble, set_data filtered by chromosome (indexcov) or gene location (mosdepth)
# TODO: filter by multiple genes
SUBSET_REGIONS <- '/storage1/fs1/timley/Active/aml_ppg/tmp/jonathanztang/breakpoint_reader/terra/mutect2/subset_expanded.bed'
filter_data <- function(set_data, gene, meta_data) {
    gene_window <- 20000
    subset_genes <- read.table(SUBSET_REGIONS, sep = '\t', header = FALSE, col.names = c('chromosome', 'start', 'end', 'name', 'zero', 'strand'))
    gene <- toupper(gene)
    subset_genes$name <- toupper(subset_genes$name)
    if (!gene %in% subset_genes$name) {
        stop(paste('Gene not found:', gene))
    }
    gene_info <- subset_genes[subset_genes$name == gene, ]
    if (meta_data$type == 'mosdepth') {
        filter_index <- (set_data$chromosome == gene_info$chromosome) &
            (set_data$start >= (gene_info$start - gene_window)) &
            (set_data$end <= (gene_info$end + gene_window))
    } else if (meta_data$type == 'indexcov') {
        filter_index <- set_data$chromosome == gene_info$chromosome
    }
    return(set_data = set_data[filter_index, ])
}

# Scale data by median of data set
scale_data <- function(set_data) {
    set_data$depth <- (set_data$depth / median(set_data$depth))
    return(set_data)
}

# Scale by pon median per bin
#   - Reads previousy generated PON from pharmU set indexcov or mosdepth data, then
#       scales data by the respective medians per bin to create scaled depth data.
#       Also filters any NA scaled depths.
#   - Inputs
#       - set_data: tibble, loaded depth data, can be filtered.
#       - type: string, type as determined by get_meta$type. Either indexcov or mosdepth.
#   - Outputs
#       - set_data: tibble, depth data scaled by pharmU median PON.
IND_PON <- '/storage1/fs1/timley/Active/aml_ppg/tmp/jonathanztang/breakpoint_reader/terra/outputs/pharmu/cnv_sample/indexcov_pon.csv'
MOS_PON <- '/storage1/fs1/timley/Active/aml_ppg/tmp/jonathanztang/breakpoint_reader/terra/outputs/pharmu/cnv_sample/mosdepth_pon.csv'
scale_pon <- function(set_data, type) {
    if (type == 'indexcov') {
        pon_path <- IND_PON
    } else if (type == 'mosdepth') {
        pon_path <- MOS_PON
    } else {
        stop('No type detected between indexcov or mosdepth.')
    }
    pon <- read_csv(pon_path, show_col_types = FALSE)
    set_pon <- left_join(set_data, pon, by = c('chromosome'='chromosome', 'start'='start', 'end'='end'), suffix = c('.set', '.pon'))
    set_data$depth <- set_pon$depth.set / set_pon$depth.pon
    set_data <- set_data %>% filter(!is.na(depth))
    return(set_data)
}

# Get log2
get_log <- function(set_data) {
    set_data$log2 <- log2(set_data$depth)
    return(set_data)
}

# Get segments
#   - Runs CBS segmentation. Creates an ordered variable list of chromosomes for sorting,
#       creates CNA from data, and segments via DNAcopy's segment() function. Uses values
#       from copyCat; undo_sd changed to 1.5 to increase sensitivity to low VAF deletions.
#   - Inputs
#       - set_data: tibble, loaded set data that has been scaled by median via scale_pon().
#       - meta_data: list, metadata obtained by get_meta().
#       - min_width: num, segment() input for minimum number of bins to call a deletion
#       - alpha: num, segment() input, refer to segment() documentation
#       - undo_sd: num, segment() input, probably minimum stdev difference to call deletion
#   - Outputs
#       - segs: DNAcopy, segment() output with $data (input), $out (segments), others.
get_segs <- function(set_data, meta_data, min_width = 3, alpha = 0.01, undo_sd = 1.5) {
    chromosomes <- set_data$chromosome %>%
        gsub("[a-z]", "", .) %>%
        ordered(levels = c(1:22, "X", "Y"))
    set_cna <- CNA(set_data$log2, chrom = chromosomes, maploc = set_data$start, data.type = 'logratio', sampleid = meta_data$set_name)
    set_cna <- smooth.CNA(set_cna)
    segs <- segment(set_cna, verbose = 0, alpha = alpha, min.width = min_width, undo.splits = "sdundo", undo.SD = undo_sd)
    return(segs)
}

# Filter segments
#   - Remove segments that contain centromeres, telomeres, etc. based on a reference if 25%+
#       of the segment overlaps with a reference interval. Standardizes chromosome notation,
#       selects reference gap regions based on present chromosomes, loads regions into
#       IRanges(), finds overlaps, intersects the ranges, and calculates percent overlap.
#       Then prepares a table to drop based on a minimum filter of 25%, filtering said rows
#       from the final segment output.
#   - Inputs
#       - set_segs: DNAcopy, raw segmentation output from get_segs()
#   - Outputs
#       - set_segs: DNAcopy, filtered segmentation data with segments overlapping with gap
#           regions more than 25% removed.
REMOVE_GAPS <- TRUE
GAPS_BED <- '/storage1/fs1/timley/Active/aml_ppg/tmp/jonathanztang/breakpoint_reader/terra/for_git/data/201_gaps.bed'
GAPS <- read.table(GAPS_BED) %>% as_tibble
names(GAPS) <- c('chr', 'start', 'end')
GAPS$start <- as.numeric(GAPS$start)
GAPS$end <- as.numeric(GAPS$end)

remove_gaps <- function(set_segs) {
    segments <- set_segs$output
    GAPS$chr <- GAPS$chr %>% str_replace('chr', '')
    ref_range <- GAPS %>% filter(chr == segments$chrom %>% unique)
    ref_range <- IRanges(start = ref_range$start, end = ref_range$end)
    seg_range <- IRanges(
        start = segments$loc.start,
        end = segments$loc.end,
        ID = segments$ID,
        chrom = segments$chrom,
        seg.mean = segments$seg.mean
        )
    hits <- findOverlaps(seg_range, ref_range)
    overlaps <- pintersect(seg_range[queryHits(hits)], ref_range[subjectHits(hits)])
    percent_overlap <- width(overlaps) / width(seg_range[queryHits(hits)])
    drop_seg <- tibble(
        seg_row = queryHits(hits),
        ref_row = subjectHits(hits),
        pct_ovr = percent_overlap,
        seg_siz = seg_range[queryHits(hits)] %>% width
    ) %>%
        group_by(seg_row) %>%
        summarize(pct = sum(pct_ovr), seg_siz = mean(seg_siz)) %>%
        filter(pct >= 0.25)
    if (nrow(drop_seg) == 0) {
        return(set_segs)
    } else {
        segments <- segments[-drop_seg$seg_row, ]
        set_segs$output <- segments
        return(set_segs)
    }
}

# DATA CHECKING ----------------------------------------------------------------
# Simple function to index raw_data rows that fall within gap_row intervals.
#   - Inputs
#       - gap_row: tibble, row of GAPS data.
#       - raw_data: tibble, unsegmented set_data data.
#   - Outputs
#       - idx: logical vector, index for matching data.
index_gap <- function(gap_row, raw_data) {
    idx <- (raw_data$chromosome == gap_row[1]) &
        (raw_data$start >= (gap_row[2] %>% as.numeric)) &
        (raw_data$start <= (gap_row[3] %>% as.numeric))
    return(idx)
}

# Multiprocessing function to vectorize interval checking
#   - Starts by selecting GAPS rows of relevant columns, then converting each row to a 
#       list entry. Then uses multiprocessing to perform index_gap to mark rows that 
#       fall within gaps, then use any() on all rows to finalize any row that has a gap.
#       Indexes raw_data against the opposite of gaps_idx to return rows that do not
#       fall within gaps.
#   - Inputs
#       - set_data: tibble, set data.
#   - Outputs
#       - out_data: tibble, set_data filtered with points within gaps removed.
filter_points <- function(set_data) {
    #set_data$chromosome <- set_data$chromosome %>% str_replace('chr', '')
    gaps <- GAPS[GAPS$chr %in% unique(set_data$chromosome),] %>% t %>% as.data.frame %>% as.list
    gaps_idx <- lapply(gaps, index_gap, raw_data = set_data) %>% bind_cols %>% apply(1, any)
    out_data <- set_data[!gaps_idx, ]
    return(out_data)
}

# Returns the median of a chunk
chunk_med <- function(chunk, set_data) {
    start_i <- min(chunk)
    end_i <- max(chunk)
    return(median(set_data$log2[start_i:end_i]))
}

# Interval median standard deviation
#   - Creates intervals of window n probes, then finds the median. Takes the standard deviation
#       across medians of all intervals to determine spread. Calculates a "waviness" statistic
#       across chromosomes 1:22 to determine how wavy a sample is.
#   - Inputs
#       - set_data: tibble, set data.
#       - win_size: number of bins to include in a window.
#   - Outputs
#       - tibble:
#           - dataset: string, name of dataset.
#           - medsd: num, standard deviation of interval medians
interval_medsd <- function(set_data, win_size = 100) {
    set_data <- filter_points(set_data) %>% filter(chromosome %in% paste0('chr', 1:22))
    set_name <- set_data$dataset %>% unique
    chunk_is <- split(set_data %>% row.names, ceiling(seq_along(set_data %>% row.names)/win_size))
    meds <- lapply(chunk_is, chunk_med, set_data = set_data) %>% unlist
    if (meds %>% mean >= 2) {
        return(tibble(
            dataset = set_name,
            medsd = 10))
    }
    out_med <- meds %>% sd
    return(tibble(
        dataset = set_name,
        medsd = out_med
        ))
}

# Batch interval median sd
# TODO: need data set names
#ims_batch <- function(set_list, win_size = 100) {
#    all_data <- mclapply(paste0(set_list, '/indexcov.tar.gz'), load)
#    all_data_2 <- mclapply(all_data, scale_pon, type = 'indexcov')
#    all_data <- mclapply(all_data_2, get_log)
#    ims <- mclapply(all_data, interval_medsd, win_size = win_size) %>% bind_rows
#    return(ims)
#}

# Calculates standard deviation of medians of intervals of window n size across multiple datasets.
# TODO: not actually needed in current workflow - useful if you have aggregate data from multiple sets.
#calc_wavy <- function(set_data, window = 100) {
#    filt_data <- filter_points(set_data)
#    medsd <- filt_data %>%
#        nest(msd = -"chromosome") %>%
#        mutate_at("msd", map, function(x) interval_medsd(x, window)) %>%
#        unnest(msd) %>%
#        arrange(desc(msd))
#    return(medsd)
#}

# WORKFLOW TOOLS ---------------------------------------------------------------
# Returns segment output
#   - Runs workflow for a single sample. Collects metadata, loads the data,
#       scales the data, obtains log2 of coverage depth, calculates waviness
#       statistic, filters data based on gene of interest (chromosome for
#       indexcov, gene window for mosdepth), performs segmentation, removes
#       segments that overlap gap regions, and returns filtered segments.
#   - Inputs
#       - data_path: string, path to coverage depth file.
#       - gene: string, name of gene of interest. Can be lower/mixed case.
#   - Outputs
#       - set_segs: DNAcopy, set segmentation data. Refer to get_segs().
seg_gene <- function(data_path, gene = 'DNMT3A') {
    meta_data <- get_meta(data_path)
    if (meta_data$set_name %>% is.na) {
        meta_data$set_name <- 'Sample'
    }
    set_data <- load(data_path)
    if (meta_data$type == 'mosdepth') {
        set_data <- scale_data(set_data)
    }
    set_data <- scale_pon(set_data, meta_data$type)
    set_data <- get_log(set_data)
    waviness <- interval_medsd(set_data)
    set_data <- set_data %>% filter_data(gene, meta_data)
    set_segs <- get_segs(set_data, meta_data)
    if (REMOVE_GAPS) {
        set_segs <- remove_gaps(set_segs)
    }
    set_segs$waviness <- waviness$medsd
    return(set_segs)
}

# Batch workflow for many samples, collects raw segment data
# TODO: remove sampling for actual run
# Batch workflow for many samples
#   - Using inputs, loads and performs segmentation on many samples. Uses
#       parallel processing to speed up the process. Recommend submitting
#       job with 32 GB of RAM to run 200 samples in chunks. Will error if
#       not enough memory is supplied, as samples are all read to memory.
#       Still much faster and more efficient than running each sample as
#       an individual job.
#   - Inputs
#       - inputs: list, refer to run_segs.R
#           - data_path: string, path to coverage depth file
#           - sample_file: string, path to sample file with data sets
#           - gene: string, name of gene of interest, can be mixed case
#           - type: string, indexcov or mosdepth
#   - Outputs
#       - seg_data: list, segmentation data where each element is a DNAcopy
#           output.
seg_gene_batch <- function(inputs) {
    folder_list <- read.table(inputs$sample_file)[[1]]
    folder_paths <- paste0(folder_list, '/', inputs$type, '.tar.gz')
    cluster <- makeForkCluster(detectCores() - 1)
    seg_data <- clusterMap(cluster,
        seg_gene,
        folder_paths,
        MoreArgs = list(
            gene = inputs$gene
            )
        )
    stopCluster(cluster)
    names(seg_data) <- names(seg_data) %>% str_split('/') %>% unlist %>% str_subset('NWD')
    return(seg_data)
}

# Filter segment output based on waviness value, for combined batch data
#   - Inputs
#       - seg_data: list, output data from seg_gene_batch()
#       - wav_thresh: num, waviness threshold
#   - Outputs
#       - seg_data: list, output data from seg_gene_batch() with filtered results removed
filter_wavy <- function(seg_data, wav_thresh = 0.042) {
    waviness <- lapply(seg_data, `[[`, 'waviness') %>% unlist
    wav_idx <- waviness <= wav_thresh
    return(seg_data = seg_data[wav_idx])
}

# Calls deletions based on thresholds, single dataset
#   - Inputs
#       - seg_data: tibble, seg_data$outputs?
#       - mean_thresh: num, mean threshold for log change
#       - bin_thresh: num, minimum numbre of bins to include
#   - Outputs
#       - deletions: tibble, table of deletions
call_dels <- function(seg_data, mean_thresh = -0.1, bin_thresh = 13) {
    deletions <- seg_data[seg_data$num.mark >= bin_thresh & seg_data$seg.mean <= mean_thresh,]# & seg_data$seg.mean >= -0.5,]
    return(deletions)
}

# Calls deletions based on thresholds in a batch dataset
#   - Inputs
#       - seg_data: list, seg_gene_batch() output
#       - mean_thresh/bin_thresh: mean and bin thresholds (refer to other fxs)
#   - Outputs
#       - list:
#           - data: list, raw seg_data from segmentation (to pass for analysis)
#           - deletions: tibble, table with deletions called from data
call_dels_batch <- function(seg_data, mean_thresh = -0.1, bin_thresh = 13) {
    segments <- seg_data %>% lapply(`[[`, 'output') %>% bind_rows
    deletions <- segments[segments$num.mark >= bin_thresh & segments$seg.mean <= mean_thresh,]# & segments$seg.mean >= -0.5,]
    waviness_raw <- seg_data %>% lapply(`[[`, 'waviness')
    waviness <- tibble(
        ID = names(waviness_raw),
        waviness = unlist(waviness_raw) %>% as.numeric
    )
    deletions <- left_join(deletions, waviness, by = 'ID')
    return(deletions)
}
 
# Workflow for a batch of data
#   - A one-stop shop that wraps the entire segmentation and deletion calling
#       pipeline. Runs batch segmentation, then calls deletions on segmentation data.
#       Returns both.
#   - Inputs
#       - inputs: input list, refer to run_segs.R
#           - data_path: string, path to coverage depth file
#           - sample_file: string, path to file with names of sample datasets
#           - gene: string, name of gene, can be mixed case
#           - type: string, indexcov or mosdepth
#           - out_dir: string, path to output directory
#           - mean_threshold: num, mean threshold to call deletions
#           - bin_threshold: num, minimum deletion size in bins
#           - wav_threshold: num, maximum waviness threshold for sample
#           - plot: logical, whether to plot outputs
#           - save_dels: logical, whether to save deletion outputs
#   - Outputs
#       - list:
#           - data: list, all segmentation data
#           - deletions: tibble, table of called deletions
del_gene_batch <- function(inputs) {
    gc()
    seg_data <- seg_gene_batch(inputs)# %>%
        #filter_wavy(inputs$wav_threshold)
    deletions <- call_dels_batch(seg_data, mean_thresh = inputs$mean_threshold, bin_thresh = inputs$bin_threshold)
    return(list(
        data = seg_data,
        deletions = deletions
    ))
}

# Plot raw segment data as alternative to plotSample from DNAcopy library
# Includes gene markers as vertical lines
GENE_REF <- '/storage1/fs1/timley/Active/aml_ppg/tmp/jonathanztang/breakpoint_reader/terra/for_git/data/gene_reference.bed'
plot_raw <- function(seg_data, gene, out_dir, cov_type, mean_thresh, bin_thresh) {
    points <- seg_data$data %>% as_tibble
    names(points) <- c('chrom', 'maploc', 'log2')
    segments <- seg_data$output %>% as_tibble
    gene_ref <- read.table(GENE_REF, sep = '\t', header = FALSE, col.names = c('chromosome', 'start', 'end', 'name', 'zero', 'strand'))
    names(gene_ref) <- c('chromosome', 'start', 'end', 'label', 'name', 'strand')
    gene <- toupper(gene)
    gene_info <- gene_ref[gene_ref$name == gene, ]
    set_name <- segments$ID %>% unique

    out_plot <- ggplot(data = points, aes(x = maploc, y = log2)) +
        geom_point(size = 0) +
        geom_segment(data = segments,
            aes(
                x = loc.start,
                xend = loc.end,
                y = seg.mean,
                yend = seg.mean
                ),
            size = 1,
            color = 'red'
            ) +
        ylim(-2, 2) +
        geom_segment(data = (GAPS %>% filter(chr == 'chr2')), aes(x = start, xend = end, y = -0.5, yend = -0.5), color = 'green') +
        geom_text(data = segments, aes(x = loc.start, y = seg.mean, label = ifelse((seg.mean <= mean_thresh & num.mark >= bin_thresh), '*', '')), color = 'pink', nudge_y = -0.075) +
        geom_vline(xintercept = gene_info$start, color = 'blue', size = 0.1) +
        geom_vline(xintercept = gene_info$end, color = 'blue', size = 0.1) +
        geom_hline(yintercept = 0, color = 'black', size = 0.1) +
        ggtitle(paste(set_name, gene, cov_type)) + xlab('Base Pair') + ylab('log2')
    ggsave(paste0(out_dir, '/', set_name, '_', cov_type, '_segments.png'), out_plot, device = 'png')
}

# Plot batch of segment data
# Inputs: raw output of del_gene_batch
plot_raw_batch <- function(deletion_data, inputs) {
    plot_out_dir <- paste0(inputs$out_dir, '/plots/')
    if (!file.exists(plot_out_dir)) {
        dir.create(plot_out_dir, recursive = TRUE)
    }
    cluster <- makeForkCluster(detectCores() - 1)
    clusterMap(cluster,
        plot_raw,
        deletion_data$data,
        MoreArgs = list(
            gene = inputs$gene,
            out_dir = plot_out_dir,
            cov_type = inputs$type,
            mean_thresh = inputs$mean_threshold,
            bin_thresh = inputs$bin_threshold
        ))
    stopCluster(cluster)
}

# SENSITIVITY/SPECIFICITY ANALYSIS ---------------------------------------------

# Quick and dirty thresholding
seg_thresh <- function(seg_data, mean_threshold = -0.1, size_threshold = 18) {
    out_data <- seg_data[(seg_data$seg.mean <= mean_threshold) & (seg_data$num.mark >= size_threshold), ]
    if (nrow(out_data) == 0) {
        seg_data$mean_thresh <- mean_threshold
        seg_data$size_thresh <- size_threshold
        return(seg_data)
    }
    out_data$mean_thresh <- mean_threshold
    out_data$size_thresh <- size_threshold
    return(out_data)
}

# False positive checking
fp <- function(seg_data) {
    sample_x <- seg_data$ID %>% unique %>% length
    return(sample_x/200)
}

# Threshold checking
#   - Batch with clusterMap based on different mean thresholds and size thresholds, compare against deletion testing
mean_ts_def <- seq(-0.5, 0, 0.05)
size_ts_def <- seq(1, 100, 1)
batch_thresh <- function(seg_list, mean_ts = mean_ts_def, size_ts = size_ts_def) {
    seg_data <- lapply(X = seg_list, FUN = `[[`, "output") %>% bind_rows
    all_options <- crossing(mean_ts, size_ts)
    mean_thresh <- all_options$mean_ts
    size_thresh <- all_options$size_ts
    cluster <- makeForkCluster(detectCores() - 1)
    thresh_data <- clusterMap(cluster,
        seg_thresh,
        mean_threshold = mean_thresh,
        size_threshold = size_thresh,
        MoreArgs = list(seg_data = seg_data)
    )
    stopCluster(cluster)
    fp_out <- vector('numeric', length(thresh_data))
    x <- 1
    for (i in thresh_data) {
        fp_out[x] <- fp(i)
        x <- x + 1
    }
    all_options$fp <- fp_out
    return(all_options)
}

min_thresh <- function(all_options, fp_thresh = 0.01) {
    return(all_options %>% filter(fp <= fp_thresh) %>% group_by(mean_ts) %>% summarize(min_size = min(size_ts)))
}

test_waviness <- function(results) {
    test_wav <- seq(0.01, 0.05, 0.0025)
    for (i in test_wav) {
        test_out <- filter_wavy(results$data, i)
        print(i)
        print(1-length(test_out)/length(results$data))
        print(batch_thresh(test_out) %>% min_thresh)
    }
}

# TODO: Deletion simulation, calling other script

# OTHER SCRIPTING/SANDBOX ------------------------------------------------------

# Gap interval generation
#orig_gaps <- read.table('super_gaps.bed') %>% as_tibble
#names(orig_gaps) <- c('chr', 'start', 'end')
#orig_gaps$start <- orig_gaps$start %>% as.double()
#orig_gaps$end <- orig_gaps$end %>% as.double()
#big_gaps <- orig_gaps
#big_gaps$width <- big_gaps$end - big_gaps$start
#big_gaps$exp <- big_gaps$width * 0.5
#chr_max_start <- big_gaps %>% group_by(chr) %>% summarize(max = max(start))
#
#for (row_i in 1:nrow(big_gaps)) {
#    check_row <- big_gaps[row_i, ]
#    if (check_row$start == 0) {
#        big_gaps[row_i, 'end'] <- big_gaps[row_i, 'end'] + big_gaps[row_i, 'width']
#    } else if (check_row$start %in% chr_max_start$max) {
#        big_gaps[row_i, 'start'] <- big_gaps[row_i, 'start'] - big_gaps[row_i, 'width']
#    } else {
#        big_gaps[row_i, 'start'] <- big_gaps[row_i, 'start'] - big_gaps[row_i, 'exp']
#        big_gaps[row_i, 'end'] <- big_gaps[row_i, 'end'] + big_gaps[row_i, 'exp']
#    }
#}
#
#big_gaps$start <- big_gaps$start %>% round(0)
#big_gaps$end <- big_gaps$end %>% round(0)
#big_gaps <- big_gaps[, c('chr', 'start', 'end')]
#
#write.table(big_gaps, 'super_gaps_200.bed', quote = F, row.names = F, col.names = F, sep = '\t')
# Make sure to combine overlapping intervals with bedtools
# sort -k1,1 -k2,2n super_gaps_200.bed > super_gaps_200_sorted.bed
# bedtools merge -d 1 -i super_gaps_200_sorted.bed > 200_gaps.bed

# Generate new sample from 200 samples reviewed with exclusions
#sample_list <- read.table('../analysis/plots/sample_200/samples.txt')[[1]]
#exclusion_list <- read.table('../analysis/plots/sample_200/exclusions.txt')[[1]] %>% str_replace('\\*', '')
#good_samples <- sample_list[!sample_list %in% exclusion_list][1:100]
#write.table(good_samples, '../analysis/sample.txt', quote = F, row.names = F, col.names = F)