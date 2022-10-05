# Script to de-isotope and de-adduct XCMS-picked peaks
# Called by Control.Rmd

# Include isIsoAdduct code here so we can safely source() the functions script
# rather than providing each of the functions because bplapply doesn't copy over
# the working environment
# given_file <- filterFile(xdata_filled, "180821_Smp_MS12C115m_A.mzML")
isIsoAdduct <- function(file_peaks, xdata_filled, addiso_masses){
  source("scripts/functions.R")
  library(data.table)
  setDTthreads(1)
  library(dplyr, warn.conflicts = FALSE)
  suppressPackageStartupMessages(library(xcms))
  
  file_addisos <- addiso_masses %>%
    {rbind(cbind(., type="self", dir=-1), cbind(., type="seek", dir=1))} %>%
    mutate(mz_diff=mz_diff*dir) %>%
    select(-dir) %>%
    group_by(addiso, type) %>%
    group_split()
  addiso_range <- range(sapply(file_addisos, function(x)x$mz_diff))
  
  given_file <- filterFile(xdata_filled, unique(file_peaks$filename))
  file_data_table <- as.data.table(given_file)[,c("rt", "mz", "i")]
  
  # Split up the large data.table into smaller ones to help speed up binary search
  split_vec <- round(1:nrow(file_peaks)/(nrow(file_peaks)/100))
  split_file_peaks <- split(file_peaks[order(file_peaks$rt),], split_vec)
  split_data_table <- lapply(split_file_peaks, function(file_peaks_i){
    file_data_table[rt%between%c(min(file_peaks_i$rtmin), max(file_peaks_i$rtmax))]
  })
  
  # test_peak <- file_peaks %>% slice_sample(n=1) # Debug
  # test_peak <- file_peaks %>% filter(mz%between%pmppm(118.0865, 5))
  # test_peak <- file_peaks %>% filter(feature=="FT0487")
  # feature_i <- test_peak$feature
  # mz_i <- test_peak$mz
  # rtmin_i <- test_peak$rtmin
  # rtmax_i <- test_peak$rtmax
  # split_table_i <- file_data_table
  
  addiso_all <- mapply(function(split_peaks_i, split_table_i){
    mapply(function(feature_i, mz_i, rtmin_i, rtmax_i){
      M_eic <- split_table_i[mz%between%pmppm(mz_i, 5)][rt%between%c(rtmin_i, rtmax_i)]
      setnames(M_eic, c("rt", "M_mz", "M_int"))
      if(nrow(M_eic)<5){
        return(data.frame(addiso="M+/-H", type="both", cor=1, area=0, feature=feature_i))
      }
      e <- lapply(file_addisos, function(row_data){
        v <- split_table_i[
          mz%between%pmppm(row_data$mz_diff+mz_i, 5)][
            rt%between%c(rtmin_i, rtmax_i)]
        v$addiso <- row_data$addiso
        v$type <- row_data$type
        v
      })
      addiso_data <- merge(rbindlist(e), M_eic, by="rt")[!(is.na(M_int)|is.na(i))]
      if(nrow(addiso_data)<5){
        return(data.frame(addiso="M+/-H", type="both", cor=1, feature=feature_i, 
                          area=trapz(M_eic$rt, M_eic$M_int)))
      }
      addiso_data[, .(cor=cor(i, M_int, use="pairwise"), 
                      area=trapz(rt, i)), by=c("addiso", "type")] %>%
        add_row(addiso="M+/-H", type="both", cor=1, 
                area=trapz(M_eic$rt, M_eic$M_int)) %>%
        mutate(feature=feature_i)
    }, SIMPLIFY = FALSE, feature_i=split_peaks_i$feature, mz_i=split_peaks_i$mz, 
    rtmin_i=split_peaks_i$rtmin, rtmax_i=split_peaks_i$rtmax) %>% 
      bind_rows()
  }, SIMPLIFY = FALSE, split_file_peaks, split_data_table) %>% 
    bind_rows() %>%
    mutate(filename = unique(file_peaks$filename))
  addiso_all
}

# Split filled_peaks by file so that we only load one file at a time (in parallel)
# file_peaks <- split(filled_peaks, filled_peaks$filename)[[1]] # Debug
start_time <- Sys.time()
isoadduct_df <- split(filled_peaks, filled_peaks$filename) %>% 
  bplapply(isIsoAdduct, xdata_filled, addiso_masses) %>%
  bind_rows() %>%
  arrange(filename, feature, addiso, type)
print(Sys.time()-start_time)
unique(warnings())
beepr::beep()

M_areas <- isoadduct_df %>%
  filter(addiso=="M+/-H") %>%
  select(filename, feature, M_area=area)
isoadduct_by_feature <- isoadduct_df %>%
  filter(area>0) %>%
  filter(addiso!="M+/-H") %>%
  left_join(M_areas, by = c("feature", "filename")) %>%
  group_by(feature, addiso, type) %>%
  filter(n()>5) %>%
  summarise(med_shape=median(cor, na.rm = TRUE), area_cor=cor(area, M_area, use = "pairwise"))

ggplot(isoadduct_by_feature) + geom_histogram(aes(x=area_cor), binwidth = 0.02) + xlim(-1, 1)
ggplot(isoadduct_by_feature) + geom_histogram(aes(x=med_shape), binwidth = 0.02) + xlim(-1, 1)

is_addiso <- isoadduct_by_feature %>%
  filter(type=="self") %>%
  filter(med_shape>shape_remove_threshold & 
           area_cor>area_remove_threshold) %>%
  anti_join(not_addisos) %>%
  ungroup()
has_addiso <- isoadduct_by_feature %>%
  filter(type=="seek") %>%
  filter(med_shape>shape_remove_threshold & area_cor>area_remove_threshold)

feature_envelopes <- has_addiso %>%
  select(feature, addiso, type) %>%
  left_join(M_areas) %>%
  left_join(isoadduct_df) %>%
  select(-cor, -type) %>%
  pivot_wider(names_from = addiso, values_from = area)
