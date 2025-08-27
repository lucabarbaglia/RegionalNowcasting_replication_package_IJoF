---
title: "Replication package for  the paper 'Nowcasting economic activity in European regions using a mixed-frequency dynamic factor model'"
output: html_document
---

This file details the workflow how to obtain the results of the paper "Nowcasting economic activity in European regions using a mixed-frequency dynamic factor model". 

*Creation date*: August 27th, 2025

*Contact information for the replication package*:

- Michael Pfarrhofer (Michael.Pfarrhofer\@wu.ac.at)

- Luca Barbaglia (luca.barbaglia\@ec.europa.eu)

*Software version*: R version 4.5.1 (2025-06-13) -- "Great Square Root"

The out-of-sample forecast exercise and in-sample results were obtained on a high-performance cluster which allows for parallel estimation across model specifications and holdouts. 
The corresponding files (see below) are used to produce the main outputs which are stored conveniently in "*.rda"-data archives ("collect_df", "preds_sub", "diag_df", "insample_df"). 
These can be used in Step 2 to produce all tables and figures contained in the paper.


**Step 0**, "00_data-prep" folder:

- "gva_eurostat.R", R-source file which downloads the most recent data from the Eurostat API

- "mf_data_OLD.rda", data archive which contains the data for a previous vintage

- "mf_data_new.rda", raw data used for all outputs contained in the paper. This archive results as output from running "gva_eurostat.R", and corresponds to "mf_data_2025-05.rda" (data as downloaded in May 2025, used for all estimations in the context of the paper) in the "01_estimation" folder


**Step 1**, "01_estimation" folder:

- "est_fcfunc.R", contains the source code for estimating the model; this file is called from both the out-of-sample estimation script and the in-sample estimation script (see below)

- "aux_func.R", contains several auxiliary functions required for estimation

- "mf_data_2025-05.rda" is the raw data file mentioned above

- "bestmodels.rda" is produced ex-post after the forecast evaluation, and selects the specifications required to be estimated in an in-sample context


The estimation files for the high-performance cluster are:

- "01a_main_fc_cluster_ffbs.R", which produces the vintages for out-of-sample design matrices, estimates the model and computes predictive losses ("est_fcfunc.R" file mentioned above)

- "01a_main_is_cluster_ffbs.R", which produces the in-sample results using all available information


The model-specific, raw outputs are stored in a folder "results", which is not contained in this replication package for storage reasons. However, the files:

- "01b_collect.R"

- "01b_collect_insample.R"

collect the out-of-sample and in-sample raw outputs in an efficient data archive after estimating everything in parallel. 

The outputs of running these two files after estimating all specifications on the cluster yields the main files required for producing all results shown in the paper:

- "collect_df.rda" contains the predictive losses across all considered specifications and different periods during the holdout

- "preds_sub.rda" is a subset of the collected file "preds_df.rda" (which is too big to be included in the package and available upon request) contains selected moments of the predictive distributions associated with the specifications in "collect_df.rda"

- "diag_df.rda" contains VAR state equation diagnostics and other functions of the estimated parameters

- "insample_df.rda" contains selected moments for the in-sample results shown in the paper when estimating the models with our full-sample.

These are the inputs required for running the output files in the main directory of this replication package, which are the source files for producing most outputs (figures and tables) contained in the paper. 


**Step 2**, main local output files:

- "02_main_fcout_local.R"; outputs:

-- Tables 1 (main paper), A.1--A.3 (appendix)

-- Figures 1, 3, 6, 7 (main paper)

- "02_main_insout_local.R"; outputs:

-- Figures 4 (main paper), A.2 (appendix)

- "02_main_maps_local.R"; outputs: 

-- Figures 2, 5 (maps, main paper)

The folder "auxiliary" contains additional outputs.


The mapping between tables/figures in the paper and files in this replication package is the following:

- Table 1  file  CRPS_h0_csrtighttable.pdf

- Figure 1 files CRPS_boxplot_factors.pdf,CRPSraw_boxplot_factors.pdf 

- Figure 2 file map_CRPS_MFREGvsHM.pdf

- Figure 3 files CRPS_boxplot_years.pdf,CRPSraw_boxplot_years.pdf 

- Figure 4 file insample_all.pdf 

- Figure 5 file map_GVA_20192020_wide.pdf

- Figure 6 file nuts_realtime.pdf 

- Figure 7 files paras_eigen.pdf,paras_logdet.pdf 

- Table A.1 file CRPS_h0_csrloosetable.pdf 

- Table A.2 file CRPS_h0_csrnotable.pdf 

- Table A.3 file CRPS_h1_csrtighttable.pdf 

- Figure A.1 file insample_all_lat.pdf 

- Figure A.2 file paras_uncondvar.pdf 


