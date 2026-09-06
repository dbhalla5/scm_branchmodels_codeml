Stochastic Character Mapping and Branch-Model, Branch-site Model Analysis with codeml (PAML)

This repository contains the scripts used to perform stochastic character mapping and branch-model, branch-site model analyses using codeml (PAML).
The workflow is general and can be applied to any gene family or phylogeny. These scripts were used for the analyses in the associated manuscript (https://www.biorxiv.org/content/10.64898/2026.02.03.700598v1.full).

scripts:   
trim_full_scm_to_gene_subset.R  —  trim trait tree for specific gene subset  
scm_and_trim.R  — 2 modes (stochastic character mapping or trimming)  
run_branch_models_subset_parallel.sh  — run codeml branch-model in parallel  


## Citation

If you use this workflow in your research, please cite this repository:

Diksha Bhalla. *SCM and codeml branch-, branch-site model analysis workflow*. GitHub.
https://github.com/dbhalla5/scm_branchmodels_codeml
