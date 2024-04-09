# MESOSCOPE Metabolomics Manuscript

This repository contains code and figures associated with the manuscript discussing the influence of mesoscale eddies on the marine microbial metabolome as well as various associated metadata.

## Abstract
Mesoscale eddies significantly alter open ocean environments such as the North Pacific Subtropical Gyre. Previous studies have explored eddy effects on biogeochemical measurements and community composition but not on the molecular composition of particulate organic matter. This study reports the exact concentration of 67 metabolites and relative abundances for 640 molecular features to understand how mesoscale eddies impact the NPSG’s metabolome during two cruises in 2017 and 2018. We find that many metabolites track biomass trends, but biosignals like isethionic acid, homarine, and trigonelline linked to picoeukaryotes were enriched at the deep chlorophyll maximum of the cyclonic features, while degradation products such as arsenobetaine were enriched in anticyclones. Additionally, we identify taurine betaine (N,N,N-trimethyltaurine) in open ocean samples for the first time. By analyzing depth variability (accounting for 20-40% of metabolomic variability across ~150 meters) and isopycnal disruption (explaining 10-20% of variability across a sea level anomaly range of 40 centimeters), this analysis constrains the importance of mesoscale eddies in modeling the chemical composition of particulate matter in the largest biomes on the planet.

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

## Figure gallery

![](https://github.com/wkumler/MesoscopeMetabolomicsManuscript/blob/master/manuscript/figures/nmds_and_med_metab.png)

*Figure 1: Plots of the metabolome across the eddy dipole, broken down by depth. The top row of plots are non-metric multidimensional scaling (NMDS) plots, in which points correspond to individual samples and have been colored by their corrected sea level height anomaly. NMDS stress values (s) have been reported in the bottom left corner, while PERMANOVA R^2 and p-values are reported in the top left. SLA trends are visible in the DCM and deep samples, with dark blue circles consistently discriminating from the dark red circles along the first multidimensional axis. The bottom row of plots show the direction and magnitude of this effect by plotting the grand median of the z-score normalized metabolite peak areas in the stations taken across the eddy dipole with the upper whisker extent, median bar, and lower whisker extent corresponding to the values measured in three biological triplicates.*

![](https://github.com/wkumler/MesoscopeMetabolomicsManuscript/blob/master/manuscript/figures/targ_gp_w_sla_frac.png)

*Figure 2: Upper row of stacked barplots shows known compounds measured across the eddy dipole and separated by depth. Bar height corresponds to median triplicate concentration for each metabolite, with the top 9 shown and the 44 other identified metabolites summed in grey. TMAO = trimethylamine N-oxide, HO-Ile = hydroxyisoleucine, GBT = glycine betaine, DCM = deep chlorophyll maximum. Lower row of boxplots shows the fraction of total particulate carbon (PC) in the known metabolites across the eddy dipole, with the upper whisker extent, median bar, and lower whisker extent corresponding to the values measured in three biological triplicates. Colors correspond to corrected sea level anomaly (Corr. SLA) in centimeters, with dark red indicating anticyclonic (positive SLA) and dark blue indicating cyclonic (negative SLA) eddy state.*

![](https://github.com/wkumler/MesoscopeMetabolomicsManuscript/blob/master/manuscript/figures/comb_pcnm_plot.png)

*Figure 2.5: Type I linear regressions of particulate carbon against total nM carbon in known metabolites both as a whole and divided by depth. Points have been colored according to the SLA and have been placed on independent axes. Outlier metabolite points have been removed from the surface at stations 6, 10, and 11 as well as a DCM sample from station 12, while two outlier PC points have been interpolated from beam attenuation for surface values at stations 10 and 11. For the full regression, see Supplemental Figure 3.*

![](https://github.com/wkumler/MesoscopeMetabolomicsManuscript/blob/master/manuscript/figures/MC_nmds_gp.png)

*Figure 3: NMDS plot of high-resolution depth sampling around the deep chlorophyll maximum (DCM, ~110 meters) at the eddy center of each polarity. Points are colored by the eddy from which the associated sample was collected and shaded by the depth. Visible separation between the cyclonic and anticyclonic eddy is evident along the first NMDS axis.*

![](https://github.com/wkumler/MesoscopeMetabolomicsManuscript/blob/master/manuscript/figures/kclust_volcano_gp.png)

*Figure 4: Distribution of metabolites in the high-resolution depth samples from the centers of each eddy. The upper row of plots shows k-means clusters where points denote the average z-scored peak area for both known and unknown metabolites across the samples and are colored by the eddy from which they were taken. Clusters have been ordered by number of metabolites in each group and the total is denoted in the panel titles. Both depth trends (mostly a net decrease in metabolites with depth) and eddy effects (cyclonic enrichment in clusters 1 and 2, anticyclone enrichment in cluster 4) are observable. The lower plot shows the individual known and unknown metabolites where points correspond to the FDR-corrected p-value estimated by the nonparametric Mann-Whitney U test and the log_2 fold-change calculated with the average peak area in the cyclone divided by the average peak area in the anticyclone. Colors have been assigned using the k-means clusters and shapes have been assigned based on the status of the mass feature as either a known metabolite that was matched to an authentic standard or an unidentified metabolite.*

![](https://github.com/wkumler/MesoscopeMetabolomicsManuscript/blob/master/manuscript/figures/fk_nmdsplot.png)

*Figure 5: Non-metric multidimensional scaling (NMDS) plot from the Falkor cruise data, in which points correspond to individual samples and have been colored by their corrected sea level height anomaly and shaped by the depth from which they were collected. NMDS stress values (s) have been reported in the bottom left corner, while PERMANOVA R^2 and p-values are reported in the top left. Surface samples are clearly distinct from the DCM samples and an SLA signal is visible at the surface.*
