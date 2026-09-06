Stochastic Character Mapping and Branch-Model, Branch-site Model Analysis with codeml (PAML)

This repository contains the scripts used to perform stochastic character mapping and branch-model, branch-site model analyses using codeml (PAML).
The workflow is general and can be applied to any gene family or phylogeny. These scripts were used for the analyses in the associated manuscript (https://www.biorxiv.org/content/10.64898/2026.02.03.700598v1.full).

## scripts:   
(1) *trim_full_scm_to_gene_subset.R* — trim trait tree for specific gene subset 

(2) *scm_and_trim.R*  — 2 modes (stochastic character mapping or trimming)  

(3) *run_branch_models_subset_parallel.sh*  — run codeml branch-model in parallel  

(4) *sample_prepare_codeml_branch_site_runs.py* — script to generate codeml branch-site .ctl files for every stochastic character mapping (SCM) realization. Expects SCM tree files for every SCM realization already exist.

(5) *sample_run_branch_site_in_batches.sh* — script to bulk run the .ctl files prepared with *sample_prepare_codeml_branch_site_runs.py* in batches (suitable for a computer cluster).


## Citation

If you use this workflow in your research, please cite this repository:

Diksha Bhalla. *SCM and codeml branch-, branch-site model analysis workflow*. GitHub.
https://github.com/dbhalla5/scm_branchmodels_codeml
