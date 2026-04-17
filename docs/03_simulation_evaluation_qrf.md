``` mermaid

graph TD
%% 1. Define custom colors for each type
classDef script fill:#E1F5FE,stroke:#0288D1,stroke-width:2px,color:#000
classDef helper fill:#A2G5FX,stroke:#0288D1,stroke-width:2px,color:#000
classDef file fill:#E8F5E9,stroke:#388E3C,stroke-width:2px,color:#000
classDef folder fill:#FFF9C4,stroke:#FBC02D,stroke-width:2px,color:#000

%% 2. Define Scripts (Stadium shape)
S1(["📜 09_simulate_data.R"]):::script
S2(["📜 10_evaluate_established_typologies.R"]):::script
S3(["📜 11_evaluate_simulated_typologies.R"]):::script
S4(["📜 12_fit_qrf.R"]):::script
	
%% 3. Define Files (Standard rectangle)
F1["📄 *taxon*_scheme.rds"]:::file
F2["📄 environmental_validation.rds"]:::file
F3["📄 eu_hydro_dem_w_enz.rds"]:::file


%% 4. Define Folders (Trapezoid shape)
D1[/"📁 001_unfitted_hmsc_models/"/]:::folder
D2[/"📁 003_fitted_hmsc_models/"/]:::folder
D3[/"📁 004_model_fit/"/]:::folder
D4[/"📁 005_variation_partitioning/"/]:::folder
D5[/"📁 unscaled_environments/"/]:::folder
D6[/"📁 catchments/"/]:::folder
D7[/"📁 simulation_diagnostics/"/]:::folder
D8[/"📁 006_simulated_data/"/]:::folder
D9[/"📁 spatial_scale/"/]:::folder
D10[/"📁 scheme_types/"/]:::folder
D11[/"📁 007_evaluationsOriginals/"/]:::folder
D12[/"📁 007_evaluations/"/]:::folder
D13[/"📁 008_qrf/"/]:::folder
D14[/"📁 taxonomic_resolution/"/]:::folder
D15[/"📁 taxa_counts/"/]:::folder


%% 5. Define Helpers 
H1(["⚙️ normalized_partitioning_entropy.R"]):::helper
H2(["⚙️ adjust_clusters_lda.R"]):::helper
H3(["⚙️ prop_sample.R"]):::helper
H4(["⚙️ calculate_auc.R"]):::helper
H5(["⚙️ group_sites_by_cluster.R"]):::helper
H6(["⚙️ balance_clusters.R"]):::helper
H7(["⚙️ cv_d_squared_V2.R"]):::helper
H8(["⚙️ render_table.R"]):::helper


	
%% 5. Draw the connections
%% 5.1 Helper to script 
H1 --> S1 
H2 --> S1 
H3 --> S2 
H4 --> S2 
H5 --> S2 
H6 --> S2 
H7 --> S2 
H8 --> S3
H5 --> S3
H3 --> S3
H4 --> S3

 
%% 5.2 Folders to script 
D1 --> S1
D1 --> S2
D2 --> S1
D3 --> S1
D4 --> S1
D5 --> S1
D6 --> S1
D8 --> S2
D8 --> S3
D9 --> S3
D10 --> S2
D12 --> S4
D4 --> S4
D9 --> S4
D14 --> S4
D15 --> S4

%% 5.3 Files to script 
F1 --> S1
F2 --> S1
F3 --> S1

%% 5.4 Script to folder 
S1 --> D7
S1 --> D8
S2 --> D11
S3 --> D12
```
