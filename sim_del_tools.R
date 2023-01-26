seg_tool_path <- '/storage1/fs1/timley/Active/aml_ppg/tmp/jonathanztang/breakpoint_reader/terra/scripts/cnv/seg_tools.R'
source(seg_tool_path)

# FUNCTIONS ====================================================================
# Simulate deletion
#   - Multiplies segment between del_start and del_end by del_pct within set_data
#   - chr2, 25113821 to 25455977 for DNMT3A
sim_del <- function(set_data, chromosome, del_start, del_end, del_pct) {
    del_idx <- (set_data$start >= del_start) & (set_data$end <= del_end) & (set_data$chromosome == chromosome)
    set_data[del_idx, 'depth'] <- set_data[del_idx, 'depth'] * del_pct
    return(set_data)
}

# Check segments for simulated deletions
#   - Collects each bin of the simulated deletion and determines if it passes a threshold.
#       If enough bins pass a second threshold, TRUE is returned; else FALSE. Marks any values
#       not in a segment as 20 to fail.
# TODO: tweak to work with called segments
#MEAN_THRESHOLD <- -0.1
#BIN_THRESHOLD <- 0.75
check_dels <- function(set_data, out_segs, del_start, del_end, mean_thres, bin_thres) {
    bin_size <- set_data[1, ]$end - set_data[1, ]$start
    out_segs$loc.end <- out_segs$loc.end + bin_size
    del_idx <- ((set_data$start >= del_start) & (set_data$end <= del_end))
    del_bins <- tibble(
        start = set_data$start[del_idx],
        end = set_data$end[del_idx],
        seg_value = vector('numeric', sum(del_idx))
    )
    for (bin_i in 1:nrow(del_bins)) {
        bin_mean <- out_segs$seg.mean[(out_segs$loc.start <= del_bins$start[bin_i]) & (out_segs$loc.end >= del_bins$end[bin_i])]
        if (length(bin_mean) > 0) {
            del_bins[bin_i, 'seg_value'] <- bin_mean
        } else if (length(bin_mean) == 0) {
            del_bins[bin_i, 'seg_value'] <- 20
        }
    }
    bin_pass <- del_bins$seg_value <= mean_thres
    if ((sum(bin_pass)/length(bin_pass)) >= bin_thres) {
        return(match = TRUE)
    } else {
        return(match = FALSE)
    }
}

# Segment deletion
#   - Run segmentation on simulated deletion data
seg_del <- function(set_data, meta_data, chromosome, del_start, del_end, del_pct, mean_thres, bin_thres) {
    if (is.na(meta_data$set_name)) {
        meta_data$set_name <- 'Sample'
    }
    if (meta_data$type == 'mosdepth') {
        set_data <- scale_data(set_data)
    }
    set_data <- scale_pon(set_data, meta_data$type) %>%
        sim_del(chromosome, del_start, del_end, del_pct) %>%
        get_log
    set_segs <- get_segs(set_data, meta_data)
    out_segs <- set_segs$output %>% as_tibble
    is_del <- check_dels(set_data, out_segs, del_start, del_end, mean_thres, bin_thres)
    return(list(
        out_segs = out_segs,
        del_start = del_start,
        del_size = del_end - del_start,
        del_pct = del_pct,
        is_del = is_del
    ))
}

# Run segments for one data set with parallelization
DEL_START_MOS <- 25213821
DEL_START_IND <- 25116672
DEL_SIZES_MOS <- c(100, 500, 1000, 5000, 10000, 50000, 100000)
DEL_SIZES_IND <- c(50000, 100000, 500000, 1000000, 5000000)
DEL_PCTS <- c(1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.4, 0.3, 0.2, 0.1, 0)
batch_dels <- function(data_path, gene = 'DNMT3A', mean_thres = -0.1, bin_thres = 0.9) {
    meta_data <- get_meta(data_path)
    set_data <- load(data_path) %>% filter_data(gene, meta_data)
    chromosome <- set_data$chromosome %>% unique
    del_pcts <- DEL_PCTS
    if (meta_data$type == 'indexcov') {
        del_start <- DEL_START_IND
        del_sizes <- DEL_SIZES_IND
    } else if (meta_data$type == 'mosdepth') {
        del_start <- DEL_START_MOS
        del_sizes <- DEL_SIZES_MOS
    } else {
        stop('Error detecting data type: indexcov or mosdepth.')
    }
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
            del_start = del_start,
            mean_thres = mean_thres,
            bin_thres = bin_thres
            )
        )
    stopCluster(cluster)
    seg_table <- do.call(rbind, seg_list) %>% as_tibble
    seg_table <- unnest(seg_table, colnames(seg_table)[-1]) %>% add_column(dataset = meta_data$set_name, mean_thres = mean_thres, bin_thres = bin_thres, .before = 2)
    return(seg_table)
}

plot_agg <- function(data_dir) {
    file_list <- list.files(data_dir) %>% str_subset('NWD')
    #meta_list <- file_list %>% str_split('_', simplify = TRUE) %>% as_tibble(.name_repair = 'unique')
    #names(meta_list) <- c('dataset', 'type', 'mean_thres', 'bin_thres')
    #meta_list$bin_thres <- str_replace(meta_list$bin_thres, '.tsv', '')
    #meta_list$mean_thres <- meta_list$mean_thres %>% as.numeric
    #meta_list$bin_thres <- meta_list$bin_thres %>% as.numeric

    # TODO: Actual false positive data

    all_results <- vector('list', length(file_list))
    for (i in 1:length(all_results)) {
        result <- read.table(paste0(data_dir, '/', file_list[i]), header = TRUE) %>% as_tibble
        result$dataset <- file_list[i] %>% str_extract('NWD[0-9]*')
        all_results[[i]] <- result
    }
    all_results <- all_results %>% bind_rows
    #all_results$is_fp <- all_results$is_del & (all_results$del_pct == 0)

    mean_data_filt <- all_results %>%
        filter(check_pct > 0) %>%# & bin_thres == 0.75) %>%
        group_by(dataset, del_size) %>%#, mean_thres) %>%
        summarize(max_sens = max(del_pct)) %>%
        ungroup %>%
        group_by(del_size) %>%#, mean_thres) %>%
        summarize(
            mean = mean(max_sens),
            min = min(max_sens),
            max = max(max_sens),
            std = sd(max_sens)
        )
    mean_fp <- all_results %>% filter(is_del & bin_thres == 0.75) %>%
        group_by(del_size, mean_thres) %>%
        summarize(fp = sum(is_fp) / n())
    mean_data <- mean_data %>% full_join(mean_fp)
    mean_melt <- mean_data %>%
        select(del_size, mean_thres, mean, std, fp) %>%
        melt(id = c('del_size', 'mean_thres', 'std', 'fp'))
    mean_plot <- ggplot(
        data = mean_melt,
        aes(x = del_size,
            y = value,
            group = mean_thres,
            color = mean_thres
            )
        ) +
        geom_line() +
        geom_point() +
        geom_errorbar(aes(ymin = value - std, ymax = value + std)) +
        geom_line(data = mean_melt, aes(x = del_size, y = fp * 100, group = mean_thres, color = mean_thres)) +
        scale_y_continuous(
            name = '1 - VAF (%)',
            sec.axis = sec_axis(~., name = 'False Positive (%)'),
            limits = c(0, 1)
        ) +
        ggtitle(paste0('Sensitivity: mean threshold')) + xlab('Deletion size (bp)') + ylab('1 - VAF (%)')
    ggsave(paste0('../sensitivity_mean_threshold.png'), plot = mean_plot)

    bin_data <- all_results %>%
        filter(is_del & mean_thres == -0.1) %>%
        group_by(dataset, del_size, bin_thres) %>%
        summarize(max_sens = max(del_pct)) %>%
        ungroup %>%
        group_by(del_size, bin_thres) %>%
        summarize(
            mean = mean(max_sens),
            min = min(max_sens),
            max = max(max_sens),
            std = sd(max_sens)
        )
    bin_fp <- all_results %>% filter(is_del & mean_thres == -0.1) %>%
        group_by(del_size, bin_thres) %>%
        summarize(fp = sum(is_fp) / n())
    bin_data <- bin_data %>% full_join(bin_fp)
    bin_melt <- bin_data %>%
        select(del_size, bin_thres, mean, std, fp) %>%
        melt(id = c('del_size', 'bin_thres', 'std', 'fp'))
    bin_plot <- ggplot(
        data = bin_melt,
        aes(x = del_size,
            y = value,
            group = bin_thres,
            color = bin_thres
            )
        ) +
        geom_line() +
        geom_point() +
        geom_errorbar(aes(ymin = value - std, ymax = value + std)) +
        geom_line(data = bin_melt, aes(x = del_size, y = fp * 100, group = bin_thres, color = bin_thres)) +
        scale_y_continuous(
            name = '1 - VAF (%)',
            sec.axis = sec_axis(~., name = 'False Positive (%)'),
            limits = c(0, 1)
        ) +
        ggtitle(paste0('Sensitivity: bin threshold')) + xlab('Deletion size (bp)') + ylab('1 - VAF (%)')
    ggsave(paste0('../sensitivity_bin_threshold.png'), plot = bin_plot)

}

# TODO: TEMPORARY to pull raw seg data
#set.seed(1)
#samples <- list.files(data_dir) %>% str_subset('NWD') %>% sample(100)
#test <- vector('list', length(samples))
#names(test) <- samples
#for (sample in list.files(data_dir) %>% str_subset('NWD')) {
#    test[[sample]] <- batch_dels(paste0(sample, '/mosdepth.tar.gz'))
#}