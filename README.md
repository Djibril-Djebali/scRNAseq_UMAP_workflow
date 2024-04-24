# Single-Cell RNA-seq Analysis Workflow for 10x GENOMICS data : data visualisation in an UMAP
## Author

Djibril Djebali (@Djibril-Djebali)

Thomas Vannier (@metavannier), https://centuri-livingsystems.org/t-vannier/

## About

This workflow performs a Snakemake pipeline to process 10x single-cell RNAseq data from RDS files to the visualisation of cell type in UMAPs compared to a reference atlas from the MouseGastrulationData R package. Correction for differences between datasets is done by a step of normalisation with Scran then the dataset is integrated with Seurat v4 following the tutorial "Mapping and annotating query datasets : Integration of 3 pancreatic islet cell datasets", next is performed a dimensionality reduction to create UMAPs with the package R ggplot2. 

Steps for the analysis:

Reference data extraction : 
Use of the MouseGastrulationData R package and the following documentation https://www.bioconductor.org/packages/release/data/experiment/vignettes/MouseGastrulationData/inst/doc/MouseGastrulationData.htmlto to extract 
the data of the following article [(doi:10.1016/j.stem.2023.04.018)](https://doi.org/10.1016/j.stem.2023.04.018).

Normalisation : 
Using the process from Pijuan-Sala et al. 2019 
[(doi:10.1038/s41586-019-0933-9)](https://doi.org/10.1038/s41586-019-0933-9) method with the Scran R package

Integration :
Using the method of integration of Seurat v4 https://satijalab.org/seurat/archive/v4.3/integration_mapping and then doing a dimensionality reduction following the same tutorial

Visualisation :
The visualisation of data was based from the figure of this article https://www.cell.com/cell-stem-cell/fulltext/S1934-5909(23)00170-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS1934590923001704%3Fshowall%3Dtrue

## Usage

You need to install [Singularity](https://github.com/hpcng/singularity/blob/master/INSTALL.md#install-golang) on your computer. This workflow also work in a slurm environment.

Each snakemake rules call a specific conda environment. In this way you can easily change/add tools for each step if necessary. 

### Step 1: Install workflow

You can use this workflow by downloading and extracting the latest release. If you intend to modify and further extend this workflow or want to work under version control, you can fork this repository.

We would be pleased if you use this workflow and participate in its improvement. If you use it in a paper, don't forget to give credits to the author by citing the URL of this repository and, if available, its [DOI](https://).

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files and repositories:
- 03_Input need the `"{data}_filtered_assigned_seurat_object.rds"` with {data} the name of your sample, that you will also specify in config.yaml.
- [config.yaml](/config.yaml) indicating the parameters to use.
- Comment the [Snakefile](/Snakefile) on the input line not expected for the pipeline.
- Build the singularity image of mambaforge:23.1.0-1

`singularity build mambaforge:23.1.0-1.sif docker://condaforge/mambaforge:23.1.0-1`
`mv mambaforge:23.1.0-1.sif 02_Container/`

### Step 3: Execute workflow

#### On your cumputer

- You need [Singularity v3.5.3](https://github.com/hpcng/singularity/blob/master/INSTALL.md#install-golang) installed on your computer or cluster.

- Load snakemake from a docker container and run the workflow from the working directory by using these commands:

`singularity run --bind ${PWD}:${PWD} docker://snakemake/snakemake:v6.3.0`

- Then execute the workflow locally via

`snakemake --use-conda --use-singularity --cores 12`

#### On a cluster

- Adapt the batch scripts run_slurm.sh to run your snakemake from the working directory

It will install snakemake with pip and run the workflow in the HPC:

`sbatch run_slurm.sh`
