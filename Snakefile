#snakemake --cores 8 --sdm conda, start workflow + env
#conda activate snakemake-tutorial ; snakemake --dag | dot -Tsvg > dag.svg ; conda deactivate, svg du workflow

configfile:"config.yaml"
TIME = ["72h", "80h", "86h", "96h"]
DATA = ["mouse_gastrulation"] + TIME

rule all:
    input:
        expand("03_Output/Normalised/{data}_norm.rds", data=DATA)

rule data_atlas:
    output:
        "03_Output/Data/mouse_gastrulation_data.rds"
    params:
        mousegastrulation_samples = config["reference_sims"]["mousegastrulation_samples"].split(','),
    conda:
        "01_Container/preparation_sims.yaml"
    script:
        "02_Script/mouse_gastrulation_data.R"  
 
rule seurat_to_sce:
    input:
        input  = expand("00_Seurat_object/{time}_filtered_seurat_object.rds",time=TIME)
    output:
        output = expand("03_Output/Data/{time}_data.rds", time=TIME)
    conda:
        "01_Container/seurat_to_sce.yaml"
    script:
        "02_Script/seurat_to_sce.R"

rule normalisation:
    input:
        input  = expand("03_Output/Data/{data}_data.rds", data=DATA)
    output:
        output = expand("03_Output/Normalised/{data}_norm.rds", data=DATA)
    conda:
        "01_Container/preparation_sims.yaml"
    script:
        "02_Script/normalisation.R"  


rule integration:
    input:
        input  = expand("")
    output:
        output = expand("")
    conda:
        "01_Container/seurat_to_sce.yaml"
        #"01_Container/preparation_sims.yaml"
    script:
        "02_Script/seurat_to_sce.R"

