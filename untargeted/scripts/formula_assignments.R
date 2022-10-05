
# Load addisos - this should be done earlier in the script but that's not guaranteed
if(!exists(addiso_masses)){
  warning("Using manually written addisos, may need updating")
  addiso_masses_all <- tribble(
    ~adduct_or_iso, ~polarity, ~addiso, ~nat_abund, ~mz_diff,
    "iso",          "both",    "C13",   0.0110,     1.003355, 
    "iso",          "both",    "X2C13", 0.0110,     2*1.003355, 
    "iso",          "both",    "S34",   0.0421,     1.995796, 
    "iso",          "both",    "S33",   0.0075,     0.999387,
    "iso",          "both",    "N15",   0.0037,     0.997035, 
    "iso",          "both",    "O18",   0.0020,     2.004244,
    "iso",          "both",    "Cl37",  0.2423,     1.99705,
    "iso",          "both",    "Br81",  0.4931,     1.997954,
    "adduct",       "both",    "M-H2O", NA,         -18.0106, 
    "adduct",       "pos",     "M+Na",  NA,          22.98922-1.007276 # Goes from [M+H] to [M+Na]
    # "adduct",       "pos",     "M+NH4", NA,          18.0338-1.007276,
    # "adduct",       "pos",     "M+K",   NA,          38.963708-1.007276,
    # "adduct",       "neg",     "M+Cl",  NA,          34.969402+1.007276,
    # "adduct",       "neg",     "M+Ac",  NA,          59.013851+1.007276,
    # "adduct",       "neg",     "M+Br",  NA,          78.918885+1.007276
  )
  addiso_masses <- addiso_masses_all %>% filter(polarity%in%c(!!polarity, "both"))
}


# Run SIRIUS ----
message("Running SIRIUS formulas...")
mzs <- feature_data$mzmed
sirius_cmd <- paste0('sirius --noCite decomp',
                     ' --mass=', paste(mzs, collapse = ","),
                     ' --ppm=5',
                     ' --elements=CHNOPS')
if(nchar(sirius_cmd)>8096){
  # Windows doesn't support commands longer than 8,096 characters
  n_chunks <- ceiling(nchar(sirius_cmd)/8096)
  mz_chunks <- split(mzs, cut(seq_along(mzs), n_chunks, labels = FALSE))
  sirius_formulas_i <- lapply(mz_chunks, function(mzs_i){
    sirius_cmd_i <- paste0('sirius --noCite decomp',
                         ' --mass=', paste(mzs_i, collapse = ","),
                         ' --ppm=5',
                         ' --elements=CHNOPS')
    sirius_output <- system(sirius_cmd_i, intern = TRUE)
    sirius_output <- sirius_output[grep(pattern = "m/z", sirius_output):length(sirius_output)]
    
    output_lengths <- sapply(sirius_output, nchar, USE.NAMES = FALSE)
    long_groups <- cumsum(output_lengths!=8095)-as.numeric(output_lengths!=8095)
    full_lines <- split(sirius_output, long_groups) %>%
      sapply(paste0, USE.NAMES = FALSE, collapse="")
    sirius_output <- grep(pattern = "^\\d", full_lines, value = TRUE)
    sirius_output <- gsub(sirius_output, pattern = "\\t$", replacement="\\\tNA")
    
    sirius_formulas <- read.table(text=sirius_output)
  })
  sirius_formulas <- do.call(sirius_formulas_i, what = "rbind")
} else {
  sirius_output <- system(sirius_cmd, intern = TRUE)
  sirius_output <- sirius_output[grep(pattern = "m/z", sirius_output):length(sirius_output)]
  
  output_lengths <- sapply(sirius_output, nchar, USE.NAMES = FALSE)
  long_groups <- cumsum(output_lengths!=8095)-as.numeric(output_lengths!=8095)
  full_lines <- split(sirius_output, long_groups) %>%
    sapply(paste0, USE.NAMES = FALSE, collapse="")
  sirius_output <- grep(pattern = "^\\d", full_lines, value = TRUE)
  sirius_output <- gsub(sirius_output, pattern = "\\t$", replacement="\\\tNA")
  
  sirius_formulas <- read.table(text=sirius_output)
}

sirius_formulas <- sirius_formulas[,2] %>%
  strsplit(split = ",")
names(sirius_formulas) <- feature_data$feature
sirius_formulas <- imap_dfr(sirius_formulas, data.frame) %>%
  select(feature=`.y..i..`, formula=`.x..i..`)


# Run Rdisop ----
message("Running Rdisop...")
iso_masses <- addiso_masses %>%
  filter(adduct_or_iso=="iso") %>%
  add_row(addiso="M", mz_diff=0, adduct_or_iso="iso", .before=1)

rdisop_formulas <- feature_envelopes %>%
  group_by(feature) %>%
  summarise(across(-c(filename, M_area), function(x){
    if(sum(!is.na(x))<5){
      return(NA)
    }
    coef(lm(x~M_area))[2]
  })) %>%
  mutate(M=1) %>%
  pivot_longer(-feature, names_to = "addiso", values_to = "rel_abund") %>%
  filter(!is.na(rel_abund)) %>%
  right_join(iso_masses) %>%
  left_join(feature_data) %>%
  mutate(addiso_mz=mzmed+mz_diff) %>%
  select(feature, addiso_mz, rel_abund) %>%
  arrange(feature, addiso_mz) %>%
  nest(envelope=-feature) %>%
  mutate(Rdisout=map(envelope, ~Rdisop::decomposeIsotopes(.x$addiso_mz, .x$rel_abund, ppm = 5)$formula)) %>%
  unnest_longer(Rdisout, values_to = "formula") %>%
  rowwise() %>%
  mutate(form_diff_H=ifelse(polarity=="pos", minH(formula), addH(formula))) %>%
  filter(formula%chin%database_formulae || form_diff_H%chin%database_formulae) %>%
  ungroup()

# rdisout <- Rdisop::decomposeIsotopes(c(76.07629, 77.07965, 77.07333, 78.08054),
#                           c(1, 0.03488086, 0.003816676, 0.002296900))


# Check isotope matches ----
message("Running isotope checks...")
manual_check <- feature_envelopes %>%
  pivot_longer(-c(feature, filename, M_area), names_to = "addiso", values_to = "value") %>%
  filter(!is.na(value)) %>%
  nest(envelope=-c(feature, addiso)) %>%
  mutate(lmout=map(envelope, ~broom::tidy(lm(.x$value~.x$M_area)))) %>%
  unnest(lmout) %>%
  filter(term!="(Intercept)") %>%
  select(feature, addiso, estimate, std.error, p.value) %>%
  filter(!addiso=="X2C13") %>%
  inner_join(addiso_masses %>% filter(adduct_or_iso=="iso")) %>%
  mutate(est_n=((estimate)*(1-nat_abund))/nat_abund) %>%
  mutate(elem=str_extract(addiso, "^[A-z]+")) %>%
  select(feature, elem, est_n)

# feature_envelopes %>%
#   filter(feature=="FT0203") %>%
#   left_join(feature_data) %>%
#   # with(summary(lm(N15~M_area+0)))
#   with(plot(O18~M_area))

feature_formulas <- rdisop_formulas %>%
  select(feature, formula) %>%
  inner_join(sirius_formulas) %>%
  mutate(formula_list=as.list(formula2elements(formula))) %>%
  unnest_longer(formula_list) %>%
  select(feature, formula, elem=formula_list_id, count=formula_list) %>%
  left_join(manual_check) %>%
  group_by(feature, formula) %>%
  summarize(euc_dist=sqrt(sum((est_n-count)^2, na.rm = TRUE))) %>%
  group_by(feature) %>%
  arrange(euc_dist) %>%
  slice(1) %>%
  select(feature, formula)
