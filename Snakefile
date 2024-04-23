#snakemake --cores 8 --sdm conda, start workflow + env
#conda activate snakemake-tutorial ; snakemake --dag | dot -Tsvg > 05_Figure/dag.svg ; conda deactivate, svg du workflow

configfile:"config.yaml"
TIME = config["fastq"]["sname"].split(",")
DATA = config["reference_sims"]["ref_name"].split(",") + TIME

rule all:
    input:
        "05_Figure/All_celltypes/UMAP_all_celltypes.pdf",
        "05_Figure/Single_celltypes/UMAP_gastruloid_single_celltypes.pdf",
        "05_Figure/Atlas_sample/UMAP_atlas_sample.pdf",
        "05_Figure/Time_simple/UMAP_time_sample.pdf",
        "05_Figure/Single_celltypes/UMAP_mouse_single_celltypes.pdf"
    
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
        input  = expand("03_Input/{data}_filtered_assigned_seurat_object.rds", data = DATA)    
    output:
        output = expand("04_Output/Normalised/{data}_norm.rds" , data = DATA)
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
        output1 = "05_Figure/All_celltypes/UMAP_all_celltypes.pdf",
        output2 = "05_Figure/Single_celltypes/UMAP_gastruloid_single_celltypes.pdf",
        output3 = "05_Figure/Atlas_sample/UMAP_atlas_sample.pdf",
        output4 = "05_Figure/Time_simple/UMAP_time_sample.pdf",
        output5 = "05_Figure/Single_celltypes/UMAP_mouse_single_celltypes.pdf"
    params:
        params = DATA
    conda:
        "01_Container/visualisation.yaml"
    script:
        "02_Script/visualisation_UMAP.R"

# rule allintegrate:
#     input:
#         input = "00_Bin/all_integrated.rds"
#     output:
#         output = "04_Output/Integrated/integrated.rds"
#     conda:
#         "00_Bin/test.yaml"
#     script:
#         "00_Bin/visualisation_UMAPtest.R"
