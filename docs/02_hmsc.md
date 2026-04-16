``` mermaid

graph TD
    %% 1. Define custom colors for each type
    classDef script fill:#E1F5FE,stroke:#0288D1,stroke-width:2px,color:#000
    classDef helper fill:#A2G5FX,stroke:#0288D1,stroke-width:2px,color:#000
    classDef file fill:#E8F5E9,stroke:#388E3C,stroke-width:2px,color:#000
    classDef folder fill:#FFF9C4,stroke:#FBC02D,stroke-width:2px,color:#000

    %% 2. Define Scripts (Stadium shape)
    S1(["⚙️ determine_spatial_scale.R"]):::helper
    S2(["⚙️ remove_collinearity.R"]):::helper
    S3(["📜 05_build_hmsc_hpc_models.R"]):::script
    subgraph ContainerB ["🐳 hmsc.sif"]
	    S4(["📜 hmsc_array.sh"]):::script
	end
	subgraph ContainerA ["🐳 r1_4.sif"]
		S5(["📜 03_model_fit.sh"]):::script
		S6(["📜 04_varPart.sh"]):::script
	end
	S7(["📜 07_evaluate_model_fit.R"]):::script
	S8(["📜 08_variation_partitioning.R"]):::script
	S9(["⚙️ gelman_check.R"]):::helper 
    S10(["⚙️ c_score.R"]):::helper 
    S11(["⚙️ spec_env_test.r"]):::helper 
    S12(["⚙️ site_occupancy.R"]):::helper 
	
    %% 3. Define Files (Standard rectangle)
    F1["📄 *taxon*_scheme.rds"]:::file
    F2["📄 *taxon*_w_environment.rds"]:::file

    %% 4. Define Folders (Trapezoid shape)
    D1[/"📁 taxonomic_resolution/"/]:::folder
    D2[/"📁 scheme_types/"/]:::folder
    D3[/"📁 taxa_counts/"/]:::folder
    D4[/"📁 spatial_scale/"/]:::folder
    D5[/"📁 unscaled_environments/"/]:::folder
    D6[/"📁 001_unfitted_hmsc_models/"/]:::folder
    D7[/"📁 002_initialized_hmsc_models/"/]:::folder
	D8[/"📁 003_fitted_hmsc_models/"/]:::folder
	D9[/"📁 004_model_fit/"/]:::folder
	D10[/"📁 005_variation_partitioning/"/]:::folder
	
    %% 5. Draw the connections
    S1 --> S3
    S2 --> S3
    F1 --> S3
    F2 --> S3
	D7 --> S4	
	D8 --> S7
    D6 --> S7
    D8 --> S7	
    D6 --> S8
    D4 --> S7
    D8 --> S8
	D9 --> S8
    
    
	S7 --> S5    
    S8 --> S6
    S3 --> D1
    S3 --> D2
    S3 --> D3
    S3 --> D4
    S3 --> D5
    S3 --> D6
    S3 --> D7
	S4 --> D8
	S7 --> D9
	S8 --> D10
	S9 -->  S7
S10 --> S7 
S11 --> S7 
S12 --> S7 
	
	
```