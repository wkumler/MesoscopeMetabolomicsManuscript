# Grab MSMS data ----
print("Reading in MS2 data...")
msdata <- grabMSdata(MSMS_files, grab_what = c("MS1", "MS2"))
msdata$MS1$rt <- msdata$MS1$rt*60
msdata$MS2$rt <- msdata$MS2$rt*60

print("Determining which have MS2 frags...")
has_msms <- feature_data %>%
  select(feature, mzmed, rtmed) %>%
  split(.$feature) %>%
  pbsapply(function(x){
    msdata$MS2[premz%between%pmppm(x$mzmed, ppm = 5)] %>%
      filter(rt%between%(x$rtmed+c(-30, 30))) %>%
      nrow() %>%
      as.logical()
  })
data.frame(feature=names(has_msms), has_msms) %>%
  filter(has_msms) %>%
  left_join(stan_annotations, by=c("feature")) %>%
  # filter(!is.na(compound_name)) %>%
  select(feature, compound_name)
feature_data <- feature_data %>%
  filter(feature%in%names(has_msms)[has_msms])



# Create .mgf files for SIRIUS to read ----
if(dir.exists(sirius_project_dir)){
  unlink(sirius_project_dir, recursive = TRUE)
}
dir.create(sirius_project_dir)
dir.create(paste0(sirius_project_dir, "//raw_files"))
message("Writing .mgf files...")
pbsapply(feature_data$feature, function(feature_num){
  output_file <- paste0(sirius_project_dir, "\\raw_files\\", feature_num, ".mgf")
  feature_msdata <- feature_data[feature_data$feature==feature_num, ]
  ms1 <- rbind(c(feature_msdata$mzmed, feature_msdata$avg_area),
               c(feature_msdata$mzmed+1.003355, feature_msdata$C13),
               c(feature_msdata$mzmed+1.003355*2, feature_msdata$X2C13),
               c(feature_msdata$mzmed+1.995796, feature_msdata$S34),
               c(feature_msdata$mzmed+0.997035, feature_msdata$N15),
               c(feature_msdata$mzmed+2.004244, feature_msdata$O18))
  ms1 <- ms1[ms1[,2]!=0, , drop=FALSE]
  ms2 <- msdata$MS2[premz%between%pmppm(feature_msdata$mzmed)&
                      rt%between%(feature_msdata$rtmed+c(-30, 30))]
  if(nrow(ms2)==0){
    ms2 <- as.data.frame(ms1) %>%
      `names<-`(c("fragmz", "int")) %>%
      mutate(voltage=0)
  }
  mgf_maker(feature_msdata = feature_msdata, ms1 = ms1, 
            ms2 = ms2, output_file = output_file)
})

# Run SIRIUS ----
message("Running SIRIUS...")
if(polarity=="pos"){
  ion_type <- "[M+H]+"
} else {
  ion_type <- "[M-H]-"
}

# Are you fucking kidding me, you have to create a project before you can create
# a project? 
dir.create(paste0(sirius_project_dir, "//output_files/0_FT0000_FEATURE_1"),
           recursive = TRUE)
writeLines("index	0\nname	FT0000_FEATURE_1\nionMass	1.0\nionType	[M + ?]+", 
           con = paste0(sirius_project_dir, "//output_files/0_FT0000_FEATURE_1/compound.info"))
writeLines(paste0(">compound FEATURE_1\n>parentmass 1.0\n>ionization [M + ?]+",
                  ">instrumentation Unknown\n>source file:Unknown\n>quality ",
                  "UNKNOWN\n#index 1\n\n>ms1peaks\n1.0 1.0"),
           con = paste0(sirius_project_dir, "//output_files/0_FT0000_FEATURE_1/spectrum.ms"))


sirius_cmd <- paste0('"C://Program Files//sirius-gui//',
                     'sirius.exe"',
                     ' --recompute', 
                     ' --loglevel="SEVERE"',
                     ' -i "', normalizePath(sirius_project_dir), '//raw_files"',
                     ' -o "', normalizePath(sirius_project_dir), '\\output_files"',
                     ' formula',
                     ' --database PUBCHEM',
                     ' --profile orbitrap',
                     ' --ions-enforced ', ion_type,
                     ' -c 50',
                     # ' zodiac', # Removed because it hangs on badly connected cmpds
                     ' fingerid',
                     ' --database bio',
                     ' canopus')
# sirius_cmd <- paste0('"C://Program Files//sirius-gui//',
#                      'sirius.exe"',
#                      ' -i "', normalizePath(sirius_project_dir), '\\raw_files"',
#                      ' -o "', normalizePath(sirius_project_dir), '\\output_files"',
#                      ' formula')

message(sirius_cmd)
system(sirius_cmd)



# Get classifications ----
canopus_tree <- read_tsv(
  paste0(sirius_project_dir, "/output_files/canopus.tsv"), col_types = 
    cols(relativeIndex = col_double(), absoluteIndex = col_double(),
         id = col_character(), name = col_character(),
         parentId = col_character(), description = col_character())) %>%
  pull(name)
classy_classes <- read_tsv(paste0(sirius_project_dir, "/output_files/canopus_summary.tsv"),
                           col_types = cols(.default = "c")) %>%
  mutate(feature=str_extract(name, "FT\\d+")) %>%
  select(feature, formula=molecularFormula, classes=`all classifications`, name) %>%
  separate_rows(classes, sep = "; ") %>%
  arrange(feature)

classy_confs <- classy_classes %>%
  split(.$name) %>%
  map_dfr(function(class_df){
    file_path <- paste0(sirius_project_dir, "/output_files/", 
                        unique(class_df$name), "/canopus/")
    filenames <- list.files(file_path, pattern = "fpt", full.names = TRUE)[1] # FIX THIS LATER, ACCOUNT FOR MULTIPLE FILES
    
    data.frame(classes=canopus_tree, conf=scan(filenames, quiet = TRUE)) %>%
      right_join(class_df, confs, by="classes") %>%
      select(classes, conf, feature, formula)
  }) %>%
  arrange(feature, desc(conf))


# Get structures ----
structure_ids <- read_tsv(paste0(sirius_project_dir, "/output_files/compound_identifications.tsv"),
                          col_types = cols(.default = "c")) %>%
  mutate(feature=str_extract(id, "FT\\d+")) %>%
  select(feature, name, confidence="CSI:FingerIDScore")
