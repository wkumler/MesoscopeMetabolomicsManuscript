
# Functions ----
stan_guesser <- function(isotope_choice, mix_choice, match_choice, area_choice, rt_choice){
  if(all(is.na(c(isotope_choice, mix_choice, match_choice, area_choice, rt_choice)))){
    return("No peaks found")
  }
  if(!is.na(isotope_choice)){
    return(isotope_choice)
  }
  if(!is.na(match_choice)) {
    return(match_choice)
  }
  if(is.na(mix_choice)){
    mix_options <- NA
  } else if(!nchar(mix_choice)){
    mix_options <- NA
  } else {
    mix_options <- unlist(strsplit(mix_choice, split = "; "))
  }
  if(length(mix_options)==1&!all(is.na(mix_options))){
    return(mix_options)
  }
  
  if(length(mix_options)>1&!all(is.na(mix_options))){
    if(rt_choice==area_choice){
      return(rt_choice)
    }
    if(rt_choice%in%mix_options&!area_choice%in%mix_options){
      return(rt_choice)
    }
    if(area_choice%in%mix_options&!rt_choice%in%mix_options){
      return(area_choice)
    }
    if(area_choice%in%mix_options&rt_choice%in%mix_options){
      return(paste(c(area_choice, rt_choice), collapse = "; "))
    }
  }
  return(rt_choice)
}
reAnno <- function(anno_df, cmpd_name, which_new=NULL, new_feature_num=NULL){
  if(is.null(new_feature_num) & is.null(which_new)){
    stop(paste0("For ", cmpd_name, ", both `new_feature_name` and ",
                "`which_new` must not both be NULL"))
  }
  if(!is.null(which_new) && !which_new%in%c("largest", "first", "last", "second")){
    stop(paste0("For ", cmpd_name, ", ", which_new, " is not a known keyword"))
  }
  if(!is.null(new_feature_num)){
    new_df <- anno_df %>% rows_update(
      data.frame(compound_name=cmpd_name, feature=new_feature_num, isotope_validated="Manual"), 
      by = "compound_name"
    )
    return(new_df)
  }
  if(!is.null(which_new)){
    chosen_feature <- anno_df %>% filter(compound_name==cmpd_name) %>% pull(feature)
    if(length(chosen_feature)==0){
      stop(paste0(cmpd_name, " not found in supplied anno_df"))
    }
    ft_mz <- mean(filter(M_area_peaks, feature==chosen_feature)$mz)
    mz_peaks <- M_area_peaks %>% filter(mz%between%pmppm(ft_mz))
    possible_feats <- mz_peaks %>% group_by(feature) %>% 
      summarize(avgarea=mean(M_area), med_rtmin=median(rtmin), 
                med_rtmax=median(rtmax), .groups = "keep") %>%
      mutate(med_rt=sum(med_rtmin, med_rtmax)/2) %>%
      ungroup()
    if(which_new%in%c("first", "second", "last")){
      possible_feats <- possible_feats %>% arrange(med_rt)
      new_feat <- if(which_new=="first"){
        possible_feats %>% slice(1) %>% pull(feature)
      } else if(which_new=="second") {
        possible_feats %>% slice(2) %>% pull(feature)
      } else if(which_new=="last"){
        possible_feats %>% slice(n()) %>% pull(feature)
      } else {
        stop(paste0("which_new cannot be ", which_new))
      }
    } else if(which_new%in%c("largest")){
      new_feat <- possible_feats %>% arrange(avgarea) %>% 
        slice(n()) %>% pull(feature)
    } else {
      stop(paste0("which_new cannot be ", which_new))
    }
    new_df <- anno_df %>% rows_update(
      data.frame(compound_name=cmpd_name, feature=new_feat, 
                 isotope_validated="Manual"), by = "compound_name"
    )
    return(new_df)
  }
}
rmAnno <- function(anno_df, cmpd_name, new_feature_num){
  reAnno(anno_df, cmpd_name, new_feature_num = NA)
}
addAnno <- function(anno_df, new_cmpd_name, new_feature_num){
  anno_df %>% add_row(
    compound_name=new_cmpd_name, feature=new_feature_num, 
    isotope_validated="Manual"
  )
}

# Assignment loop ----

stan_annotations <- given_stans %>%
  split(.$compound_name) %>%
  pblapply(function(stan_data){
    stan_data <- slice(stan_data, 1)
    possible_stan_features <- feature_data %>%
      filter(mzmed%between%pmppm(as.numeric(stan_data["mz"]), 10)) %>%
      arrange(rtmed)
    if(!nrow(possible_stan_features)){
      return(rep(NA, 5))
    }
    
    # If the compound has an isotopologue, use RT matching
    # If the compound IS an isotopologue, same deal
    isotopologue <- stan_data$compound_name %>%
      paste0(", [15N\\d|13C\\d|2H\\d]") %>%
      grep(x = given_stans$compound_name, value = TRUE)
    if(length(isotopologue)==1){
      isotopo_data <- given_stans %>% 
        filter(compound_name==isotopologue)
      isotopo_features <- feature_data %>%
        filter(mzmed%between%pmppm(as.numeric(isotopo_data$mz), 6))
      if(nrow(isotopo_features)==1){
        isotope_choice <- possible_stan_features %>% 
          mutate(rtdiff=abs(rtmed-isotopo_features$rtmed)) %>%
          arrange(rtdiff) %>%
          slice(1) %>%
          pull(feature)
      } else {
        isotope_choice <- NA
      }
    } else if(grepl(pattern = ", 15N|13C|2H", x = stan_data$compound_name)) {
      isotopologue <- gsub(", 15N.*|, 13C.*|, 2H.*", "", stan_data$compound_name)
      isotopo_data <- given_stans %>% 
        filter(compound_name==isotopologue|
                 compound_name==gsub("^D", "", isotopologue)|
                 compound_name==gsub("^L", "", isotopologue))
      isotopo_features <- feature_data %>%
        filter(mzmed%between%pmppm(as.numeric(isotopo_data$mz), 5))
      isotope_choice <- possible_stan_features %>% 
        left_join(isotopo_features, by=character()) %>%
        select(-starts_with("mzmed")) %>%
        mutate(rtdiff=abs(rtmed.x-rtmed.y)) %>%
        arrange(rtdiff) %>%
        slice(1) %>%
        pull(feature.x)
    } else {
      isotope_choice <- NA
    }
    
    # If there's a peak that differs between the mixes
    if(is.na(stan_data$mix)){
      mix_choice <- NA
    } else {
      stan_peaks <- M_area_peaks %>%
        filter(feature%in%possible_stan_features$feature) %>%
        filter(str_detect(filename, "Mix")) %>%
        filter(grepl("Std", filename)) %>%
        filter(!grepl("H2OinMatrix", filename)) %>%
        mutate(date_run=strptime(str_extract(filename, "^\\d+"), format = "%y%m%d")) %>%
        filter(date_run>=stan_data$date_added) %>%
        select(feature, mz, rt, M_area, filename) %>%
        mutate(correct_mix=grepl(filename, pattern = stan_data["mix"])) %>%
        mutate(stan_type=str_extract(filename, "InH2O|InMatrix"))
      
      if(!nrow(stan_peaks)){
        mix_choice <- NA
      } else {
        mix_peaks <- stan_peaks %>% 
          group_by(feature, correct_mix, stan_type) %>%
          summarise(avgarea=mean(M_area), sdarea=sd(M_area)) %>%
          right_join(expand.grid(
            unique(.$feature),
            unique(.$correct_mix),
            unique(.$stan_type)
          ) %>% `names<-`(c("feature", "correct_mix", "stan_type")),
          by=c("feature", "correct_mix", "stan_type")) %>%
          mutate(avgarea=ifelse(is.na(avgarea), 0, avgarea)) %>%
          mutate(sdarea=ifelse(is.na(sdarea), 0, sdarea))
        
        if(length(unique(mix_peaks$feature))==1){
          mix_choice <- unique(mix_peaks$feature)
        } else {
          mix_choice <- mix_peaks %>%
            split(interaction(.$feature, .$stan_type)) %>%
            lapply(function(v){
              if(nrow(v)==1){
                if(v$correct_mix==TRUE){
                  data.frame(feature=unique(v$feature), 
                             stan_type=unique(v$stan_type),
                             diff_degree=Inf)
                } else {
                  data.frame(feature=unique(v$feature), 
                             stan_type=unique(v$stan_type),
                             diff_degree=-Inf)
                }
              } else {
                diff <- (v$avgarea[v$correct_mix] - v$avgarea[!v$correct_mix])/mean(v$sdarea)
                data.frame(feature=unique(v$feature), 
                           stan_type=unique(v$stan_type),
                           diff_degree=diff)
              }
            }) %>%
            do.call(what = rbind) %>% 
            `rownames<-`(NULL) %>%
            group_by(feature) %>%
            summarise(correct_mix_peak=mean(diff_degree)) %>%
            filter(correct_mix_peak>1) %>%
            arrange(desc(correct_mix_peak)) %>%
            pull(feature) %>%
            paste(collapse = "; ")
        }
      }
    }
    if(!is.na(mix_choice))if(!nchar(mix_choice))mix_choice<-NA

    # If there's one peak much closer in RT to expected than the others
    expected_rt <- stan_data$rt*60
    rt_choice <- possible_stan_features %>% 
      mutate(rtdiff=abs(rtmed-expected_rt)) %>%
      arrange(rtdiff) %>%
      slice(1) %>%
      pull(feature)
    
    # If there's same number of features as expected peaks, assume 1:1 and order by RT
    possible_other_stans <- given_stans %>%
      filter(mz%between%pmppm(stan_data$mz, 5)) %>%
      arrange(rt)
    if(nrow(possible_stan_features)==nrow(possible_other_stans)){
      match_choice <- possible_stan_features %>%
        arrange(rtmed) %>%
        cbind(possible_other_stans) %>%
        select(feature, compound_name) %>%
        filter(compound_name==stan_data$compound_name) %>%
        pull(feature)
    } else {
      match_choice <- NA
    }
    
    stan_peaks <- M_area_peaks %>%
      filter(feature%in%possible_stan_features$feature) %>%
      select(feature, mz, rt, M_area, filename) %>%
      filter(grepl("Std", filename))
    
    if(nrow(stan_peaks)==0){
      area_choice <- NA
    } else if(nrow(possible_stan_features)==nrow(possible_other_stans)&
              nrow(possible_stan_features)>1) {
      area_choice <- NA
    } else {
      area_choice <- stan_peaks %>% 
        group_by(feature) %>%
        summarise(areamed=median(M_area)) %>%
        arrange(desc(areamed)) %>%
        slice(1) %>%
        pull(feature)
    }
    
    feature <- stan_guesser(isotope_choice, mix_choice, match_choice, area_choice, rt_choice)
    
    return(data.frame(
      feature=feature,
      isotope_validated=isotope_choice,
      rt_matchup=match_choice,
      mix_matched=mix_choice,
      closer_rt=rt_choice,
      area_choice=area_choice
    ))
  }) %>%
  do.call(what=rbind) %>%
  mutate(compound_name=rownames(.)) %>%
  select(compound_name, everything())



# Initial dual-assignment check ----
# Check on features that are dual-assigned - these should be resolved by manual assignment
stan_annotations[(duplicated(stan_annotations$feature, fromLast=TRUE)|
                   duplicated(stan_annotations$feature))&
                   !is.na(stan_annotations$feature), ] %>%
  split(.$feature)


# Also check on compounds that are dual-assigned, also need manual curation
stan_annotations %>%
  filter(str_detect(feature, "; ")) %>%
  mutate(mz=str_extract(feature, "FT.*(?=; )")) %>%
  left_join(feature_data, by=c(mz="feature")) %>%
  select(compound_name, feature, mzmed)



# Manual falkor pos ----
if(cruise=="FK" & polarity=="pos"){
  stan_annotations <- stan_annotations %>%
    rmAnno("L-Alanine") %>%
    rmAnno("beta-Alanine") %>%
    addAnno("Alanine/beta-Alanine", new_feature_num = "FT0054") %>%
    reAnno("Sarcosine", new_feature_num = "FT0055") %>%
    reAnno("Homarine", new_feature_num = "FT0269") %>%
    rmAnno("O-Acetyl-L-serine") %>%
    rmAnno("beta-Glutamic acid") %>%
    rmAnno("Betonicine") %>%
    rmAnno("Turicine") %>%
    addAnno("Betonicine/Turicine", new_feature_num = "FT0398") %>%
    reAnno("Allopurinol", new_feature_num = "FT0258") %>%
    reAnno("Hypoxanthine", new_feature_num = "FT0256") %>%
    reAnno("N6-Acetyl-L-lysine", new_feature_num = "FT0541") %>%
    reAnno("Dexpanthenol", new_feature_num = "FT0603") %>%
    reAnno("Hydroxyproline", new_feature_num = "FT0218")
}
# Manual mesotransect pos ----
if(cruise=="MT" & polarity=="pos") {
  # anno_stan_manual <- stan_annotations %>%
  stan_annotations <- stan_annotations %>%
    rmAnno("L-Valine") %>%
    rmAnno("L-Isoleucine") %>%
    rmAnno("Muramic acid") %>%
    rmAnno("L-Threonine") %>%
    rmAnno("L-Homoserine") %>%
    addAnno("Threonine/Homoserine", new_feature_num = "FT0216") %>%
    rmAnno("4-Aminobutyric acid") %>%
    reAnno("Allopurinol", new_feature_num = "FT0359") %>%
    reAnno("N-Methyltaurine", new_feature_num = "FT0385") %>%
    reAnno(cmpd_name = "N6-Acetyl-L-lysine", new_feature_num = "FT0759") %>%
    reAnno("Nicotinic acid", new_feature_num = "FT0252") %>%
    reAnno("S-Adenosylmethionine", new_feature_num = "FT1537") %>%
    rmAnno("O-Acetyl-L-serine") %>%
    rmAnno("Betonicine") %>%
    rmAnno("Turicine") %>%
    addAnno("Betonicine/Turicine", new_feature_num = "FT0539") %>%
    # reAnno("Carnitine", new_feature_num = "FT0554") %>%
    reAnno("Choline", new_feature_num = "FT0139") %>%
    reAnno("L-Serine", new_feature_num = "FT0149") %>%
    addAnno("beta-Alaninebetaine", new_feature_num = "FT0314")
  warning("L-Serine has been poorly integrated and cannot be assigned to a distinct feature")
}

# Manual mesocenter pos ----
if(cruise=="MC" & polarity=="pos"){
  stan_annotations <- stan_annotations %>%
    rmAnno("L-Valine") %>%
    rmAnno("Dimethylglycine") %>%
    rmAnno("O-Acetyl-L-serine") %>%
    rmAnno("L-Isoleucine") %>%
    reAnno("L-Isoleucine, 15N", which_new = "first") %>%
    reAnno("5-Hydroxyectoine", which_new = "largest") %>%
    rmAnno("Muramic acid") %>%
    rmAnno("Pyro from glutamate?") %>%
    rmAnno("L-Threonine")
}

# Manual falkor neg ----
if(cruise=="FK" & polarity=="neg"){
  stan_annotations <- stan_annotations %>%
    rmAnno("Trehalose") %>%
    rmAnno("Maltose") %>%
    rmAnno("Lactose") %>%
    reAnno("Trehalose, 13C12", which_new = "last")
}
# Manual mesocenter neg ----
if(cruise=="MC" & polarity=="neg"){
  stan_annotations <- stan_annotations %>%
    rmAnno("Trehalose") %>%
    reAnno("Trehalose, 13C12", which_new = "last")
}
# Manual mesotransect neg ----
if(cruise=="MT" & polarity=="neg"){
  stan_annotations <- stan_annotations %>%
    rmAnno("Lactose") %>%
    rmAnno("Maltose") %>%
    reAnno("Trehalose, 13C12", new_feature_num = "FT418") %>%
    reAnno("5-Oxoproline", new_feature_num = "FT140") %>%
    rmAnno("Succinic semialdehyde")
}

# Post-manual assignment check ----
message("Features assigned to two compounds:")
stan_annotations[(duplicated(stan_annotations$feature, fromLast=TRUE)|
                    duplicated(stan_annotations$feature))&
                   !is.na(stan_annotations$feature), ] %>%
  left_join(given_stans %>% select(compound_name, mz)) %>%
  split(.$feature) %>%
  print()

message("Compounds assigned to two features:")
stan_annotations %>%
  filter(str_detect(feature, "; ")) %>%
  mutate(mz=str_extract(feature, "FT.*(?=; )")) %>%
  left_join(feature_data, by=c(mz="feature")) %>%
  select(compound_name, feature, mzmed) %>%
  print()
