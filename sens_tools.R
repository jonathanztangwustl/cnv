seg_tool_path <- '/storage1/fs1/timley/Active/aml_ppg/tmp/jonathanztang/breakpoint_reader/terra/scripts/cnv/seg_tools.R'
source(seg_tool_path)

# FUNCTIONS ====================================================================
DEL_START_IND <- 25116672
DEL_SIZES_IND <- c(250000, 500000, 750000, 1000000, 2500000, 5000000, 10000000)
DEL_PCTS <- c(1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7)

# Simulate deletion
sim_del <- function(set_data, chromosome, del_start, del_end, del_pct) {
    del_idx <- (set_data$start >= del_start) & (set_data$end <= del_end) & (set_data$chromosome == chromosome)
    set_data[del_idx, 'depth'] <- set_data[del_idx, 'depth'] * del_pct
    return(set_data)
}

# Run segmentation with simulated deletion
seg_del <- function(set_data, meta_data, chromosome, gene, del_start, del_end, del_pct, mean_threshold, bin_threshold, wav_threshold) {
    # Identical to seg_gene except with a simulated deletion and pre-loaded data
    set_data <- sim_del(set_data, chromosome, del_start, del_end, del_pct)
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
    set_segs$del_start <- del_start
    set_segs$del_end <- del_end
    set_segs$del_pct <- del_pct
    return(set_segs)
}

# Return what proportion of bins in the deletion are within a called deletion
check_del <- function(seg_row) {
    deletions <- seg_row$deletions
    del_start <- seg_row$del_start
    del_end <- seg_row$del_end
    if (nrow(deletions) == 0) {
        return(0)
    }
    bin_size <- (deletions$loc.end[1] - deletions$loc.start[1]) / (deletions$num.mark[1] - 1)
    probes <- seq(del_start, del_end, bin_size)
    in_seg <- vector('list', nrow(deletions))
    for (del_i in 1:nrow(deletions)) {
        in_seg[[del_i]] <- between(probes, deletions$loc.start[del_i], deletions$loc.end[del_i])
        names(in_seg[[del_i]]) <- probes
    }
    return((bind_rows(in_seg) %>% colSums %>% sum) / length(probes))
}

# Run a batch of segmentations for deletions for one sample
seg_del_batch <- function(data_path, chromosome = 'chr2', gene = 'DNMT3A', mean_threshold = -0.1, bin_threshold = 13) {
    meta_data <- get_meta(data_path)
    set_data <- load(data_path)
    del_pcts <- DEL_PCTS
    del_start <- DEL_START_IND
    del_sizes <- DEL_SIZES_IND
    all_options <- crossing(del_sizes, del_pcts)
    del_sizes <- all_options$del_sizes
    del_pcts <- all_options$del_pcts
    cluster <- makeForkCluster(detectCores())
    seg_list <- clusterMap(cluster,
        seg_del,
        del_end = del_start + del_sizes,
        del_pct = del_pcts,
        MoreArgs = list(
            set_data = set_data,
            meta_data = meta_data,
            chromosome = chromosome,
            gene = gene,
            del_start = del_start,
            mean_thres = mean_threshold,
            bin_thres = bin_threshold
            )
        )
    stopCluster(cluster)
    seg_table <- tibble(
        data = lapply(seg_list, `[[`, 'data') %>% lapply(as_tibble),
        segments = lapply(seg_list, `[[`, 'output') %>% lapply(as_tibble),
        deletions = lapply(seg_list, `[[`, 'output') %>% lapply(call_dels, mean_thresh = mean_threshold, bin_thresh = bin_threshold),
        waviness = lapply(seg_list, `[[`, 'waviness') %>% unlist(recursive = F),
        del_start = lapply(seg_list, `[[`, 'del_start') %>% unlist(recursive = F),
        del_end = lapply(seg_list, `[[`, 'del_end') %>% unlist(recursive = F),
        del_pct = lapply(seg_list, `[[`, 'del_pct') %>% unlist(recursive = F)
    )
    seg_table$del_size <- seg_table$del_end - seg_table$del_start
    seg_table$check_pct <- apply(seg_table, 1, check_del)
    return(seg_table)
}

# Return table with limited data for writing output sensitivity data
write_sens <- function(seg_table) {
    sens_data <- tibble(
        del_size = seg_table$del_size,
        del_pct = seg_table$del_pct,
        check_pct = seg_table$check_pct,
        waviness = seg_table$waviness
    )
    return(sens_data)
}

# Plot segment data for one row of a simulated deletion segmentation result table
plot_del_seg <- function(seg_table_row, out_dir, gene = 'DNMT3A', cov_type = 'indexcov', mean_thresh = -0.1, bin_thresh = 13) {
    options(scipen = 999)
    if (!file.exists(out_dir)) {
        dir.create(out_dir)
    }
    points <- seg_table_row$data %>% as_tibble
    names(points) <- c('chrom', 'maploc', 'log2')
    segments <- seg_table_row$segments %>% as_tibble
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
        geom_text(data = segments, aes(x = loc.start, y = seg.mean, label = ifelse((seg.mean <= mean_thresh & num.mark >= bin_thresh), '*', '')), color = 'blue', nudge_y = -0.075) +
        geom_vline(xintercept = gene_info$start, color = 'blue', size = 0.1) +
        geom_vline(xintercept = gene_info$end, color = 'blue', size = 0.1) +
        geom_hline(yintercept = 0, color = 'black', size = 0.1) +
        ggtitle(paste(set_name, gene, cov_type)) + xlab('Base Pair') + ylab('log2')
    ggsave(paste0(out_dir, '/', set_name, '_', cov_type, '_', seg_table_row$del_size, '_', seg_table_row$del_pct, '.png'), out_plot, device = 'png')
}

# Plot entire simulated deletion segmentation result table
plot_batch_del_seg <- function(seg_table, out_dir, gene = 'DNMT3A', cov_type = 'indexcov', mean_thresh = -0.1, bin_thresh = 16) {
    options(scipen = 999)
    cluster <- makeForkCluster(detectCores())
    parApply(cluster,
        seg_table, 1, plot_del_seg,
        out_dir = out_dir,
        gene = gene,
        cov_type = cov_type,
        mean_thresh = mean_thresh,
        bin_thresh = bin_thresh
        )
    stopCluster(cluster)
}
