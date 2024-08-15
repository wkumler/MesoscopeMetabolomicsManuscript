# MESOSCOPE Metabolomics Manuscript

This repository contains code and figures associated with the manuscript discussing the influence of mesoscale eddies on the marine microbial metabolome as well as various associated metadata. This manuscript was submitted to Frontiers in Marine Science on 8/15/2024.

## Abstract
Mesoscale eddies significantly alter open ocean environments such as those found in the subtropical gyres that cover a large fraction of the global ocean. Previous studies have explored eddy effects on biogeochemistry and microbial community composition but not on the molecular composition of particulate organic matter. This study reports the absolute concentration of 67 metabolites and relative abundances for 640 molecular features to understand how mesoscale eddies impact the metabolome of the North Pacific Subtropical Gyre during two cruises in 2017 and 2018. We find that many metabolites track biomass trends, but metabolites like isethionic acid, homarine, and trigonelline linked to eukaryotic phytoplankton were enriched at the deep chlorophyll maximum of the cyclonic features, while degradation products such as arsenobetaine were enriched in anticyclones. In every analysis, metabolites with the strongest responses were detected via untargeted mass spectrometry, indicating that the molecules most sensitive to environmental perturbation were not among the characterized metabolome. By analyzing depth variability (accounting for 20-40% of metabolomic variability across ~150 meters) and the vertical displacement of isopycnal surfaces (explaining 10-20% of variability across a sea level anomaly range of 40 centimeters and a spatial distance of 300 kilometers), this analysis constrains the importance of mesoscale eddies in shaping the chemical composition of particulate matter in the largest biomes on the planet.

## Repository structure

  - manuscript/
    - **manuscript.Rmd** R Markdown document containing code for analysis and written text
    - **figures/** Main text figures for the manuscript (see gallery below)
  - metadata/
    - **metadata_control.Rmd** R Markdown document that downloads and collates metadata into filled_file_metadata.csv for use in the analysis
    - **raw_data/** Raw data downloaded from the SCOPE website (https://scope.soest.hawaii.edu/data/)
    - **allmeso_metadata.pdf** Simple visualization of metadata columns across cruise track
    - **clean_stans.csv** Cleaned copy of the [Ingalls Lab Standards Sheet](https://github.com/IngallsLabUW/Ingalls_Standards)
  - targeted/
    - **targeted_control.Rmd** Script for running the targeted analyses
      - Requires metadata/filled_file_metadata.csv, metadata/clean_stans.csv
      - Requires Skyline files and SkylineRunner.exe (MSDIAL output must be manually exported)
      - Produces made_data/longdata.csv
    - **raw_data/** (ignored by Git)
      - Falkor/, Mesocenters/, and Mesotransect/ directories with Skyline and MSDIAL outputs
    - **made_data/**
      - **longdata.csv** Clean version of all targeted data and values after BMIS and Quant for use in the manuscript analysis
  - untargeted/
    - **untargeted_control.Rmd** Script for running the untargeted analyses
    - ***_output** Folders containing intermediate data products from XCMS for Falkor (FK_) and MESOSCOPE (MS_)
    - **all_feats.csv** Molecular feature information (mzmed, rtmed, compound/structure/class)
    - **all_peaks.csv** Peak area information (feature, filename, area)
  - mzMLs/ (ignored by Git)
    - neg/, pos/, and tricho_pos/ directories containing the corresponding [mzML files from Metabolomics Workbench](http://dx.doi.org/10.21228/M82719)

## Deliverables

  - [Exact (nM) concentrations for 68 particulate metabolites](https://github.com/wkumler/MesoscopeMetabolomicsManuscript/blob/master/targeted/Falkor_Targeted_nM.xlsx) in 24 Falkor samples
  - [Exact (nM) concentrations for 60 particulate metabolites](https://github.com/wkumler/MesoscopeMetabolomicsManuscript/blob/master/targeted/Falkor_Targeted_nM.xlsx) in 99 MESOSCOPE eddy transect samples and 30 MESOSCOPE eddy center samples
  - [Peak areas for 640 molecular features](https://github.com/wkumler/MesoscopeMetabolomicsManuscript/blob/master/untargeted/all_peaks.csv) (de-adducted and de-isotoped, normalized to BMIS) across Falkor and MESOSCOPE (analyzed separately)

## Main text figures

![](https://github.com/wkumler/MesoscopeMetabolomicsManuscript/blob/master/manuscript/figures/fig_1_MapForWill_v3.jpg)

*Figure 1: Sampling during MESO-SCOPE and the Hawaiian Eddy Experiment (HEE). Cruise bounds are shown in the large central map with Station ALOHA colored as a point in red near the Hawaiian Islands. Yellow stars denote sampling locations and station numbers relative to the sea level anomaly contours in the background during both the June 2017 MESO-SCOPE cruise (left) and the April 2018 HEE cruise (right).*

![](https://github.com/wkumler/MesoscopeMetabolomicsManuscript/blob/master/manuscript/figures/fig_2_nmds_and_med_metab.tif)

*Figure 2: Distribution of metabolites in multivariate space across adjacent eddies of opposite polarity during the MESO-SCOPE transect, broken down by depth. Top panels depict non-metric multidimensional scaling (NMDS) plots with individual samples colored based on their corrected sea level anomaly. NMDS stress values (s) have been reported in the bottom left corner, while PERMANOVA R$^2$ and p-values are reported in the top left. SLA trends are visible in the DCM and 175 meter samples, with dark blue circles consistently discriminating from the dark red circles along the first multidimensional axis. Bottom panels depict the direction and magnitude of this effect by plotting the grand mean value of three biological triplicates in color and the raw values in black behind.*

![](https://github.com/wkumler/MesoscopeMetabolomicsManuscript/blob/master/manuscript/figures/fig_3_targ_gp_w_sla_frac.tif)

*Figure 3: Differences in known metabolite concentration across the pair of adjacent eddies in the MESO-SCOPE transect separated by depth. Top panels depict concentrations of known compounds measured, where bar height corresponds to median triplicate concentration for each metabolite with the top 9 shown and the 44 other identified metabolites summed in grey. TMAO = trimethylamine N-oxide, HO-Ile = hydroxyisoleucine, GBT = glycine betaine, DCM = deep chlorophyll maximum. Lower panels depict the fraction of total particulate carbon in the known metabolites, with the mean value of three biological triplicates in color and the raw values in black behind. Colors correspond to corrected sea level anomaly (Corr. SLA), with dark red indicating anticyclonic (positive SLA) and dark blue indicating cyclonic (negative SLA) eddy state.*

![](https://github.com/wkumler/MesoscopeMetabolomicsManuscript/blob/master/manuscript/figures/fig_4_MC_nmds_gp.tif)

*Figure 4: NMDS plot of high-resolution depth sampling around the deep chlorophyll maximum (DCM, ~115 meters) at the two eddy centers during MESO-SCOPE. Red upward-pointed triangles are from the anticyclone and blue downward-pointed ones are from the cyclone. Shading intensity reflects the depth above or below the DCM. PERMANOVA estimates of the variance explained by depth and sea level anomaly (SLA) are noted in the upper left corner and the NMDS stress value is reported in the bottom left.*

![](https://github.com/wkumler/MesoscopeMetabolomicsManuscript/blob/master/manuscript/figures/fig_5_kclust_volcano_gp.tif)

*Figure 5: Distribution of metabolites in the high-resolution depth samples from the centers of each MESO-SCOPE eddy. The upper row of plots shows k-means clusters where points denote the average z-scored peak area for both known and unknown metabolites across the samples and are colored by the eddy from which they were taken. Clusters have been ordered by number of metabolites in each group and the total is denoted in the panel titles. Both depth trends (mostly a net decrease in metabolites with depth) and eddy effects (cyclonic enrichment in clusters 1 and 2, anticyclone enrichment in cluster 4) are observable. The lower plot shows the individual known and unknown metabolites where points correspond to the FDR-corrected p-value estimated by the nonparametric Mann-Whitney U test and the log$_2$ fold-change calculated with the average peak area in the cyclone divided by the average peak area in the anticyclone. Colors have been assigned using the k-means clusters and shapes have been assigned based on the status of the mass feature as either a known metabolite that was matched to an authentic standard or an unidentified metabolite. The dashed line across the figure represents the 0.05 level of significance as a visual cue for metabolites above which the differences between the eddies are unlikely to be due to chance. DMSP = dimethylsulfoniopropionate, DMS-Ac = dimethylsulfonioacetate*

![](https://github.com/wkumler/MesoscopeMetabolomicsManuscript/blob/master/manuscript/figures/fig_6_fk_nmdsplot.tif)

*Figure 6: Non-metric multidimensional scaling (NMDS) plot from the Hawaiian Eddy Experiment cruise data, in which points correspond to individual samples and have been colored and shaped by their source eddy status and shaded by the depth from which they were collected. The NMDS stress value has been reported in the bottom left corner, while PERMANOVA R$^2$ and p-values are reported in the top left. Samples from 25 meters deep are visibly distinct from the deep chlorophyll maximum (DCM) samples and an SLA signal is visible in the 25 meter samples only.*
