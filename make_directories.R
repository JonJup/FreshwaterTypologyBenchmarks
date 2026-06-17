## Settting up all necesarry folders in the project location 
## 


# LEVEL 1 -----------------------------------------------------------------

dir.create("parent")
dir.create("R")
dir.create("diatoms_folder")
dir.create("fish_folder")
dir.create("invertebrates_folder")
dir.create("macrophytes_folder")

dirs <- list.dirs()
if ("." %in% dirs) dirs <- dirs[-which(dirs == ".")]
for (i in seq_along(dirs)){
        
        setwd(dirs[i])
        dir.create("data/")
        dir.create("code/")
        dir.create("logs/")
        dir.create("data/001_unfitted_hmsc_models")
        dir.create("data/002_initialized_hmsc_models")
        dir.create("data/003_fitted_hmsc_models")
        dir.create("data/004_model_fit")
        dir.create("data/005_variation_partitioning")
        dir.create("data/006_simulated_data")
        dir.create("data/007_evaluations_empirical")
        dir.create("data/007_evaluations")
        dir.create("data/008_qrf")
        dir.create("data/misc/")
        dir.create("data/misc/spatial_scale")
        dir.create("data/misc/taxonomic_resolution")
        dir.create("data/misc/unscaled_environments/")
        dir.create("data/misc/taxa_counts")
        dir.create("data/misc/simulation_diagnostics")
        
        if (grepl("parent",dirs[i])){
                dir.create("output/")
                dir.create("output/figures/")
                dir.create("output/figures/manuscript")
                dir.create("data/results/")
                dir.create("data/catchments/")
                dir.create("data/catchments_w_environment")
        }
        setwd("..")
}
