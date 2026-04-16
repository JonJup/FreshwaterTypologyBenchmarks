``` mermaid

graph TD

    %% 1. Define custom colors for each type
	    classDef script fill:#E1F5FE,stroke:#0288D1,stroke-width:2px,color:#000
	    classDef helper fill:#A2G5FX,stroke:#0288D1,stroke-width:2px,color:#000
	    classDef file fill:#E8F5E9,stroke:#388E3C,stroke-width:2px,color:#000
	    classDef folder fill:#FFF9C4,stroke:#FBC02D,stroke-width:2px,color:#000
	    
	%% 2. Define Scripts (Stadium shape)
    S1(["📜 01_add_eu_hydro_to_biota.R"]):::script
    S2(["📜 02_combine_env.R"]):::script
    S3(["📜 03_add_env_biota.R"]):::script
    S4(["📜 04_define_schemes.R"]):::script
    
    %% 3. Define Helpers 
    H1(["⚙️ find_max_consequtive_sum.R"]):::helper
    
	%% 4. Define Files (Standard rectangle)
    F1["📄 *taxon*_w_catchment_id.rds"]:::file
    F2["📄 *EU_HYDRO_CATCHMENT*_w_variables.parquet"]:::file
    F3["📄 eu_hydro_dem_w_enz.rds"]:::file
    F4["📄 *taxon*_scheme.rds"]:::file
    F5["📄 *taxon*_w_environment.rds"]:::file
	
    %% 5. Define Folders (Trapezoid shape)
    D1[/"📁 MIDFIRE/"/]:::folder
	D2[/"📁 hydroDEM_parquet/"/]:::folder
	D3[/"📁 typologies/"/]:::folder
    
    D1 --> S1
     D2 --> S1
F2 --> S3
F1 --> S3
D3 --> S3 
F3 --> S4
F5 --> S4

H1 --> S4
     
     S1 --> F1
     S2 --> F2    
     S3 --> F5
S4 --> F4     

```