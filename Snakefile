#snakemake --cores 8 --sdm conda, start workflow + env
#conda activate snakemake-tutorial ; snakemake --dag | dot -Tsvg > dag.svg ; conda deactivate, svg du workflow

configfile:"config.yaml"
TIME = config["fastq"]["sname"]

rule all:
    input:
        "03_Output/Normalised/Seurat/Normalised_mouse_gastrulation.rds_Seurat.RDS",
        expand("03_Output/Data/{time}_data.rds", time=TIME)

rule data_atlas:
    output:
        "03_Output/Data/mouse_gastrulation_data.rds"
    params:
        mousegastrulation_samples = config["reference_sims"]["mousegastrulation_samples"].split(','),
    conda:
        "01_Container/preparation_sims.yaml"
    script:
        "02_Script/mouse_gastrulation_data.R"  

rule normalisation:
    input:
        input = "03_Output/Data/mouse_gastrulation_data.rds"
    output:
        "03_Output/Normalised/Seurat/Normalised_mouse_gastrulation.rds_Seurat.RDS"
    params:
        sample_id                 = config["reference_sims"]["sample_id"],
        mousegastrulation_samples = config["reference_sims"]["mousegastrulation_samples"].split(',')
    conda:
        "01_Container/preparation_sims.yaml"
    script:
        "02_Script/normalisation.R"  
 
rule seurat_to_sce:
    input:
        input  = expand("00_Seurat_object/{time}_filtered_seurat_object.rds",time=TIME)
    output:
        output = expand("03_Output/Data/{time}_data.rds", time=TIME)
    conda:
        "01_Container/seurat_to_sce.yaml"
        #"01_Container/preparation_sims.yaml"
    script:
        "02_Script/seurat_to_sce.R"

rule integration:
    input:
        input  = expand("00_Seurat_object/{time}_filtered_seurat_object.rds",time=["72h", "80h", "86h", "96h"])
    output:
        output = expand("03_Output/Data/{time}_data.rds", time=["72h", "80h", "86h", "96h"])
    conda:
        "01_Container/seurat_to_sce.yaml"
        #"01_Container/preparation_sims.yaml"
    script:
        "02_Script/seurat_to_sce.R"

