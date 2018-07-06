# **_Readme file for Butterfly R script descriptions_**

1. **1_prep_UKBMS_data**: 
     - Take raw UKBMS data and clean the data to remove negative values & pilot years etc.
     - Creates 3 files: 
        - b_data: the cleaned UKBMS data with the following headings: SPECIES, BROOD, SITE, YEAR, & SINDEX
        - site_data_spp: estimate the number of years of data for each site
        - years_data_spp: estimate the number of survey data for each year

2. **2_calc_growth_rates**:
     - Takes b_data & calculates the growth rate for each species using a moving window (1980-2016)
     - Creates 1 file:
        - final_data_all_spp: growth rate data for all species with the following headings: site, year, gr, SINDEX, mid.year, start.year, end.year, name, & rec_id

3. **3_calc_pop_synchrony**: 
     - Calculates the population syncrhony for each pair-wise site comparison for each species using a moving window (1980-2016)
     - Also calculates the abundance of each species (i.e. the mean SINDEX) using a moving window (1980-2016)
     - Creates 3 files:
        - final_pair_data_all_species: synchrony (lag0) data for all species with the following headings: site1, site2, lag0, numYears, mid.year, start.year, end.year, & spp
        - final_summ_stats_all_spp: number of good pair-wise comparisons and number of unique sites with the following headings: total_comps, uni_sites, mid.year, start.year, end.year, & spp
        - abundance_data: abundance data for all species with the following headings: start.year, mid.year, end.year, abundance, & species

4. **4_temporal_trend_synchrony**:
     - Calculates the temporal trend in population synchrony
     - Create pair_attr file to begin with which includes synchrony data plus site attribute data 
     - Final model: lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp)
     - Creates 5 files:
        - pair_attr: synchrony & site attribute data with the following headings: spp, site1, site2, lag0, numYears, mid.year, start.year, end.year, site_a_EAST, site_a_NORTH, site_b_EAST, site_b_NORTH, mean_northing, distance, renk_hab_sim, COMMON_NAME, HABITAT, mobility.wil, & pair.id
        - results_final_all_spp: Synchrony model results for ALL species (i.e. one value for each year), with the following headings: Estimate, SD, t, parameter, rescaled_est, rescaled_sd, & rescaled_ci
        - results_final_sp: synchrony model results for EACH species (i.e. one value per species per year), with the following headings: sp, Estimate, SD, t, parameter, rescaled_est, rescaled_sd, rescaled_ci, COMMON_NAME, & HABITAT	
        - results_final_hab: synchrony model results for each HABITAT (i.e. one value per habitat per year), with the following headings: Estimate, SD, t, habitat, parameter, rescaled_est, rescaled_sd, & rescaled_ci
        - results_final_mobility: synchrony model results for each MOBILITY type (i.e. one very per mobility type per year), with the following headings: Estimate, SD, t, mobility.score, parameter, rescaled_est, rescaled_sd, & rescaled_ci
