# MESOSCOPE Metabolomics Manuscript

This repository contains code and figures associated with the manuscript discussing the influence of mesoscale eddies on the marine microbial metabolome as well as various associated metadata.

## Abstract
Mesoscale eddies significantly alter open ocean environments such as the North Pacific Subtropical Gyre. Previous studies have explored eddy effects on biogeochemical measurements and community composition but not on the molecular composition of particulate organic matter. This study reports the exact concentration of 67 metabolites and relative abundances for 640 molecular features to understand how mesoscale eddies impact the NPSG’s metabolome during two cruises in 2017 and 2018. We find that many metabolites track biomass trends, but biosignals like isethionic acid, homarine, and trigonelline linked to picoeukaryotes were enriched at the deep chlorophyll maximum of the cyclonic features, while degradation products such as arsenobetaine were enriched in anticyclones. Additionally, we identify taurine betaine (N,N,N-trimethyltaurine) in open ocean samples for the first time. By analyzing depth variability (accounting for 20-40% of metabolomic variability across ~150 meters) and isopycnal disruption (explaining 10-20% of variability across a sea level anomaly range of 40 centimeters), this analysis constrains the importance of mesoscale eddies in modeling the chemical composition of particulate matter in the largest biomes on the planet.

## Repository structure

## Deliverables

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
