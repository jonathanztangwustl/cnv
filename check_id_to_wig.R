
## RAW WIG FOR COMPARISON ------------------------------------------------------
#raw_wig <- read_lines('normal.wig')
#raw_wig_id <- raw_wig %>% startsWith('fixedStep') %>% which
#raw_wig_id_end <- (raw_wig_id - 1) %>% .[-1] %>% c(., raw_wig %>% length)
#raw_wig_names <- raw_wig %>%
#    str_subset('fixedStep') %>%
#    str_extract('chr[0-9, X, Y].?\\s') %>%
#    str_replace(' ', '')
#
## Identify raw_wig chromosomes
#raw_wig_chrs <- vector('list', length(raw_wig_id))
#for (chr in 1:length(raw_wig_chrs)) {
#    raw_wig_chrs[[chr]] <- raw_wig[(raw_wig_id[[chr]] + 1):raw_wig_id_end[[chr]]]
#}
#names(raw_wig_chrs) <- raw_wig_names
#
## Identify offsets, as readcounter is a few bins larger than indexcov
#raw_wig_size <- raw_wig_chrs %>% lapply(length) %>% unlist
#offset_size <- (raw_wig_size - ind_size) %>% enframe
#write.table(offset_size, 'offset.txt', quote = F, row.names = F)
#

# ALIGNMENT CHECKING -----------------------------------------------------------
#rctr_wig <- read_lines('normal.wig')
#indx_wig <- read_lines('test.wig')
#
#wig_name_extract <- function(wig) {
#    wig_names <- wig %>%
#        str_subset('fixedStep') %>%
#        str_extract('chr[0-9, X, Y].?\\s') %>%
#        str_replace(' ', '')
#    return(wig_names)
#}
#
#wig_to_list <- function(wig) {
#    wig_id <- wig %>% startsWith('fixedStep') %>% which
#    wig_id_end <- (wig_id - 1) %>% .[-1] %>% c(., wig %>% length)
#    wig_names <- wig_name_extract(wig)
#    wig_chrs <- vector('list', length(wig_id))
#    for (chr in 1:length(wig_chrs)) {
#        wig_chrs[[chr]] <- wig[(wig_id[[chr]] + 1):wig_id_end[[chr]]]
#    }
#    names(wig_chrs) <- wig_names
#    return(wig_chrs)
#}
#
#rctr <- wig_to_list(rctr_wig)
#indx <- wig_to_list(indx_wig)
#rchr6 <- rctr[[6]] %>% as.numeric
#indx6 <- indx[[6]] %>% as.numeric
#
#rat <- rchr6/indx6
#test <- tibble(
#    x = 1:length(rchr6),
#    readcounter = scale(rchr6),
#    indexcov = scale(indx6)
#) %>% pivot_longer(
#    !x,
#    names_to = 'source'
#)
#
#plot <- ggplot(data = test, aes(x = x, y = value, group = source, color = source)) + geom_point()
#ggsave('test.png', plot)