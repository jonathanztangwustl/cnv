# Short script to concatenate results for one set of data files.
# Will expand to combine results for multiple sets of runs.
# Just run in directory of interest. Will capture all results in the folder,
# unless 'old' is in the directory.

library('tidyverse')

# I'm tired of writing this every time.
write_wrap <- function(input, output) {
    write.table(
        input,
        output,
        sep = '\t',
        quote = F,
        row.names = FALSE
        )
}

# The cluster sometimes doesn't refresh from storage, so sometimes
# there are ghost results. Sending ls -l usually cleans that up, and
# if it gets logged it'll help troubleshoot if it's just a cache
# issue.
system('ls -l')

# Exclude files in "old" folders. Only pull .txt files to speed up.
text_files <- list.files(recursive = T, pattern = '*.txt') %>%
    grep(pattern = 'old', invert = T, value = T)

del_files <- text_files %>% str_subset('_all_deletions')
gene_files <- text_files %>% str_subset('_gene_deletions')
stat_files <- text_files %>% str_subset('_stats')

all_dels <- vector('list', length(del_files))
gene_dels <- vector('list', length(gene_files))
stats <- vector('list', length(stat_files))

for (i in 1:length(del_files)) {
    del <- read.table(del_files[i], header = T)
    if (nrow(del) > 0) {
        all_dels[[i]] <- del
    }
}
all_dels <- all_dels %>% bind_rows
write_wrap(all_dels, 'all_deletions.txt')

for (i in 1:length(gene_files)) {
    del <- read.table(gene_files[i], header = T)
    if (nrow(del) > 0) {
        gene_dels[[i]] <- del
    }
}
gene_dels <- gene_dels %>% bind_rows
write_wrap(gene_dels, 'gene_deletions.txt')

for (i in 1:length(stat_files)) {
    stats[[i]] <- read.table(stat_files[i], header = T)
}
stats <- stats %>% bind_rows %>% as_tibble
stats <- add_row(stats, colSums(stats) %>% t %>% as_tibble)
stats <- add_column(stats, run = c(str_replace(stat_files, '_stats.txt', ''), 'Total'), .before = 1)
stats$pct <- stats$filt_deletion_samples / stats$filt_samples
write_wrap(stats, 'stats.txt')
