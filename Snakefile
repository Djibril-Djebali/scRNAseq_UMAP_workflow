#snakemake --cores 8 --sdm conda, start workflow + env
#conda activate snakemake-tutorial ; snakemake --dag | dot -Tsvg > 05_Figure/dag.svg ; conda deactivate, svg du workflow

configfile:"config.yaml"
TIME = ["72h", "80h", "86h", "96h"]
DATA = ["mouse_gastrulation"] + TIME

rule all:
    input:
        "05_Figure/UMAP.pdf"
    

rule data_atlas:
    output:
        output = "03_Input/Atlas/mouse_gastrulation_data.rds"
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
        output = "03_Input/Atlas/Convert/Seurat/Convert_atlas.rds_Seurat.RDS"
    conda:
        "01_Container/preparation_sims.yaml"
    script:
        "02_Script/sce_to_seurat.R"

rule normalisation:
    input:
        input  = "03_Input/{data}_filtered_assigned_seurat_object.rds"    
    output:
        output = "04_Output/Normalised/{data}_norm.rds"
    conda:
        "01_Container/preparation_sims.yaml"
    script:
        "02_Script/normalisation.R"  


rule integration:
    input:
        input  = expand("04_Output/Normalised/{data}_norm.rds", data=DATA)
    output:
        output = "04_Output/Integrated/integrated.rds"
    conda:
        "01_Container/integration.yaml"
    script:
        "02_Script/integration.R"

rule new_file_path:
    input:
        "03_Input/Atlas/Convert/Seurat/Convert_atlas.rds_Seurat.RDS"
    output:
        "03_Input/mouse_gastrulation_filtered_assigned_seurat_object.rds"
    shell:
       """
       cp {input} duplicate.txt
       mv duplicate.txt {output}
       """

rule visualisation:
    input:
        input  = "04_Output/Integrated/integrated.rds"
    output:
        output = "05_Figure/UMAP.pdf"
    conda:
        "01_Container/integration.yaml"
    script:
        "02_Script/visualisation.R"
