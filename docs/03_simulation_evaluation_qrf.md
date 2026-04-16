``` mermaid

graph TD
    %% 1. Define custom colors for each type
    classDef script fill:#E1F5FE,stroke:#0288D1,stroke-width:2px,color:#000
    classDef helper fill:#A2G5FX,stroke:#0288D1,stroke-width:2px,color:#000
    classDef file fill:#E8F5E9,stroke:#388E3C,stroke-width:2px,color:#000
    classDef folder fill:#FFF9C4,stroke:#FBC02D,stroke-width:2px,color:#000

    %% 2. Define Scripts (Stadium shape)

	
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


	
	%% 5. Define Helpers 
	H1(["⚙️ normalized_partitioning_entropy.R"]):::helper
	H2(["⚙️ adjust_clusters_lda.R"]):::helper
	
	
%% 5. Draw the connections
%% 5.1 Helper to script 
H1 --> S1 
H2 --> S1 

%% 5.2 Folders to script 
D1 --> S1
D2 --> S1
D3 --> S1
D4 --> S1
D5 --> S1
D6 --> S1

%% 5.3 Files to script 
F1 --> S1
F2 --> S1
F3 --> S1

%% 5.4 Script to folder 
S1 --> D7
S1 --> D8

```