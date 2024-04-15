#snakemake --cores 8 --sdm conda, start workflow + env
#conda activate snakemake-tutorial ; snakemake --dag | dot -Tsvg > 05_Figure/dag.svg ; conda deactivate, svg du workflow

configfile:"config.yaml"
TIME = ["72h", "80h", "86h", "96h"]
DATA = ["mouse_gastrulation"] + TIME

rule all:
    input:
        expand("04_Output/Normalised/{data}_norm.rds", data=DATA)

rule data_atlas:
    output:
        "03_Input/Atlas/mouse_gastrulation_data.rds"
    params:
        mousegastrulation_samples = config["reference_sims"]["mousegastrulation_samples"].split(','),
    conda:
        "01_Container/preparation_sims.yaml"
    script:
        "02_Script/mouse_gastrulation_data.R"  

rule sce_to_seurat:
    input:
        input  = "03_Input/Atlas/mouse_gastrulation_data.rds"
    output:
        output = "03_Input/mouse_gastrulation_filtered_seurat_object.rds"
    conda:
        "01_Container/preparation_sims.yaml"
    script:
        "02_Script/sce_to_seurat.R"

rule normalisation:
    input:
        input = expand("03_Input/{data}_filtered_seurat_object.rds", data=DATA)       
    output:
        output = expand("04_Output/Normalised/{data}_norm.rds", data=DATA)
    conda:
        "01_Container/preparation_sims.yaml"
    script:
        "02_Script/normalisation.R"  


rule integration:
    input:
        input  = expand("04_Output/Normalised/{data}_norm.rds", data=DATA)
    output:
        output = "test.txt"
    conda:
        "01_Container/integration.yaml"
    script:
        "02_Script/integration.R"

#rule seurat_to_sce:
#    input:
#        input  = expand("00_Seurat_object/{time}_filtered_seurat_object.rds",time=TIME)
#    output:
#        #output = expand("03_Output/Data/{time}_data.rds", time=TIME)
#    conda:
#        "01_Container/seurat_to_sce.yaml"
#    script:
#        "02_Script/seurat_to_sce.R"
