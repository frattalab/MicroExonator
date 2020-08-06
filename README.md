# Introduction

MicroExonator is a fully-integrated computational pipeline that allows for systematic de novo discovery and quantification of microexons using raw RNA-seq data for any organism with a gene annotation. Compared to other available methods MicroExonator is more sensitive for discovering smaller microexons and it provides higher specificity for all lengths. Moreover, MicroExonator provides integrated downstream comparative analysis between cell types or tissues using [Whippet](https://github.com/timbitz/Whippet.jl) ([Sterne-Weiler et al. 2018](https://doi.org/10.1016/j.molcel.2018.08.018)). As a proof of principle MicroExonator  identified X novel microexons in Y RNA-seq samples from mouse early development to provide a systematic characterization based on time and tissue specificity.


# Changes in this fork

##### Additional control over which MEs to assess for differential splicing with Whippets
Given our interest in cryptic splicing events as a TDP43 LoF phenotype, we consider microexons with weaker splicing signals and sequence conservation to be potentially biologically relevant events. See additional config option in *Whippets* section to also consider microexons that fail the splice site signal filter (*Report/out_low_scored_ME.txt*) for differential inclusion assessment.

(exploratory)

##### Generate bedgraph files with standard reference genome coordinates from MicroExonator SAM/BAM files
This allows us to visualise microexon coverages in IGV with standard reference genomes, gene models etc. *A lot simpler* than reformatting BAM reference names and cigar strings (though would be great to be able to visualise reads/read alignments!)

(keep_bams branch - will merge when fully confident)

##### Job submit scripts for our SGE UCL cluster

`source submit_microexonator.sh <run_name>`

`source submit_whippets_microexonator.sh <run_name>`

##### On the list - Output data to a 'master results directory'
When I can be bothered to change the filepaths in all rules. Would save having to manually move outputs of one run before re-submitting (and running over and over in same directory saves needing duplicate conda installations!).


##### On the list - integrate automatic handling of paired end samples

See [open issue](https://github.com/hemberg-lab/MicroExonator/issues/4) at main repo - unlikely to come until paper has been published/ fully released. Would be handy to simplify processing our paired end data and improve discovery power.



Also other random updates to README and little typo fixes in scripts that I've noticed and reported. There is now a readthedocs hosted on main repo...


# Installation

Start by cloning MicroExonator

    git clone https://github.com/hemberg-lab/MicroExonator

Install [Miniconda 3](https://docs.conda.io/en/latest/miniconda.html)

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod +x Miniconda3-latest-Linux-x86_64.sh ./Miniconda3-latest-Linux-x86_64.sh


Finally, create an enviroment to run [snakemake](https://snakemake.readthedocs.io/en/stable/)

    conda create -n snakemake_env -c bioconda -c conda-forge snakemake

# Configuration

Before running MicroExonator you need to have certain files in the `MicroExonator/` directory. By default, all the downloaded files and the results will be stored in the `MicroExonator/` folder which means that several GBs of disk space may be required. First, you need to create a `config.yaml`, which should contain the path of the input files and certain paramethers. The format for the `config.yaml` file is shown below:

    #Make sure chromosome name format of header lines of Genome_fasta matches chromosome names in Gene_anontation_bed12
    Genome_fasta : /path/to/Genome.fa
    Gene_anontation_bed12 : /path/to/ensembl.bed12

    #Put "NA" if don't provide either of these two files
    GT_AG_U2_5 : /path/to/GT_AG_U2_5.good.matrix
    GT_AG_U2_3 : /path/to/GT_AG_U2_3.good.matrix

    conservation_bigwig : /path/to/conservation.bw  

    #make sure working_directory has trailing slash
    working_directory : /path/to/MicroExonator/

    ME_DB : /path/to/ME_DB.bed
    ME_len : 30
    Optimize_hard_drive : T

Here:

* `Genome_fasta` is a [multifasta](http://www.metagenomics.wiki/tools/fastq/multi-fasta-format) file containg the cromosomes.
* `Gene_anontation_bed12` is a BED file containing the transcript annotation. A collection of these files can be found at [UCSC genome browser](http://genome.ucsc.edu/cgi-bin/hgTables).
* `GT_AG_U2_5` and `GT_AG_U2_5` are splice site PWMs that can come from [SpliceRack](http://katahdin.cshl.edu//splice/splice_matrix_poster.cgi?database=spliceNew) (please remove the first lines so it looks like the files that we have in  `/PWM`). If you do not have the PWMs for the species you are interested in, you can leave it as `NA` and MicroExonator will generate the PWM internally based on annotated splice sites.
* `conservation_bigwig` is a bigwig file containing genome-wide conservation scores generated by Pylop or PhastCons, which can be downloaded from [UCSC genome browser](http://hgdownload.cse.ucsc.edu/downloads.html) for some species. If you do not have a `conservation_bigwig` file, you need to create an bigwig which has the value 0 for every position in your genome assembly.
* `working_directory` is the path to the MicroExonator folder that you are working with
* `ME_len` is microexon maximum lenght. The default value is 30.
* `ME_DB` is a path to a known Microexon database such as VAST DB (is this optional or no?)
* `Optimize_hard_drive` can be set as `T` or `F` (true of false). When it is set as `F`, fastq files will be downloaded or copied only once at `FASTQ/` directory. This can be inconvenient when you are analysing large amout of data (for instance when your input data is larger than 1TB), because the copy will not be deleted until MicroExonator finish completely. When `Optimize_hard_drive` is set as `T` instead, two independet local copies will be generated for the discovery and quantification moduled. As these are temporary copies, every copied fastq file will be deleted as soon as is mapped to the splice junction tags, which means the disk space usage will be set to the minimum while MicroExonator is runing.

To imput the RNA-seq data, you need to either create a `local_samples.tsv` (for fastq.gz files that are stored locally) and/or `NCBI_accession_list.txt`(for SRA accession names) which are automatically downloaded. You can find examples of these files inside the Examples folder.

Finnaly, if you are working on a high performace cluster, then it is very likely that you need to submit jobs to queueing systems such as lsf, qsub, SLURM, etc. To make MicroExonator work with these queueing systems, you need to create a `cluster.json` file. We currently provide in the Examples folder a `cluster.json` file to run MicroExonator with [lsf](https://www.ibm.com/support/knowledgecenter/en/SSETD4/product_welcome_platform_lsf.html). To adapt MicroExonator to other quequing systems please see the [SnakeMake documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html?highlight=cluster.json#cluster-configuration).



# Running

We highly recommend creating a screen before running MicroExonator

    screen -S session_name  #choose a meaning full name to you

To activate snakemake enviroment

    conda activate snakemake_env

In case you already have an older version of conda use instead:

    source activate snakemake_env


Then run

    snakemake -s MicroExonator.skm  --cluster-config cluster.json --cluster {cluster system params} --use-conda -k  -j {number of parallel jobs}

Before running it is recommended to check if SnakeMake can corretly generate all the steps given your input. To do this you can carry out a dry-run using the `-np` parameter:

    snakemake -s MicroExonator.skm  --cluster-config cluster.json --cluster {cluster system params} --use-conda -k  -j {number of parallel jobs} -np

The dry-run will display all the steps and commands that will be excecuted. If the dry-run cannot be initiated, make sure that you are running MicroExonator from inside the folder you cloned from this repository. Also make sure you have the right configuration inside `config.yaml`.

Notice that you should use `--cluster` only if you are working in a computer cluster that uses a queuing systems. We provide an example of `cluster.json` to work with lsf and in that case the cluster system params should be `"bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}"`. The number of parallel jobs can be a positive integer, the appropriate value depends on the capacity of your machine but for most users a value between 5 and 50 is appropriate.

If you want to process a large dataset, we recommend to set `Optimize_hard_drive` as `T`. If you tho this, you can also run the pipeline in two steps:

    snakemake -s MicroExonator.skm  --cluster-config cluster.json --cluster {cluster system params} --use-conda -k  -j {number of parallel jobs} discovery

Once the pipeline finish the discovery module, it can be resumed any time later as:

    snakemake -s MicroExonator.skm  --cluster-config cluster.json --cluster {cluster system params} --use-conda -k  -j {number of parallel jobs} quant

If you do not have space to have a let MicroExonator to have a complete copy of the data you want to process, you can limmit the number of jobs. For instance if you limit the number of jobs to 20, you can make sure no more than 20 fastq files will be stored at `FASTQ/` directory (if `Optimize_hard_drive` is set as `T`)

If you are working remotelly, the connection is likely to die before MicroExonator finish. However, as long as you are working within an screen, loggin off will not kill snakemake. To list your active screens you can do:

    screen -ls

To reattached and detach screens just use:

    screen -r session_name  # only detached screen can be reattached  
    screen -d session_name

## Alternative splicing analyses

Microexonator can couple the annotation and quantification of microexon with other downstream tools to asses alternative splicing changing. We have incorporated [whippet](https://github.com/timbitz/Whippet.jl) into a downstream snakemake workflow that can directly use the microexon annotation and quantification files to assess differential inclusion of microexons between multiple conditions. For this porpuse we recomend to install a custumised conda enviroment that has snkamemake and julia 0.6.1 installed. This can be don by using the following command from inside a MicroExonator folder:

`conda env create -f Whippet/julia_0.6.1.yaml`

Where `Whippet/julia_0.6.1.yam` is yaml file with the recipe to create a stable environment with julia 0.6.1 and snakemake.

Then, activate the newly created enviroment:

`source activate julia_0.6.1`

Enter julia's interactive mode:

`julia`

Install Whippet:

`Pkg.add("Whippet")`

Exit interactive julia session (`control + d`) and find Whippet's binary folder, that should be inside of your miniconda environment folder. Once you find the path to this folder, add it to `config.yaml` writing the following lines:

    whippet_bin_folder : path/to/miniconda/envs/julia_0.6.1/share/julia/site/v0.6/Whippet/bin

    #Would recommend generating Gene_anontation_bed12 from this GTF file so they match up
    Gene_anontation_GTF : path/to/gene_annotation.gtf
    whippet_delta: path/to/whippets_config.yaml

    #Optional key - if include value should be "all" or "high"
    whippet_which_MEs: <"all">/<"high">


Where `whippet_delta` correspond to a config file that has all the comparions betweeen samples that you would like to do. Each comparison can have an arbitrary name (here as `Control_vs_TreatmentA` and `Control_vs_TreatmentB`). For each comparison two groups of samples are compared; `A` and `B`. These correspond to comma-separated list of the replicates for the conditions you are interested to compare. The `whippet_delta` config file should be formatted like so:


    Control_vs_TreatmentA :
        A : SRR309144,SRR309143,SRR309142,SRR309141
        B : SRR309136,SRR309135,SRR309134,SRR309133
    Control_vs_TreatmentB :
        A : SRR309140,SRR309139
        B : SRR309138,SRR309137

`whippet_which_MEs` is an **optional key** that tells Microexonator which microexons to select for differential splicing analysis. If this option is **not included** or set to **"high"**, Whippets will only assess high confidence microexons that have passed all filters (`Report/out.high_quality.txt`). If set to **"low"**, Whippets will assess all microexons in `Report/out.high_quality.txt` and `Report/out_low_scored_ME.txt` for differential splicing.


Then we recomend to do a `dry-run` to check that all the inputs are in place:

    snakemake -s MicroExonator.skm  --cluster-config cluster.json --cluster {cluster system params} --use-conda -k  -j {number of parallel jobs} differential_inclusion -np


**Important**: If you arlready run MicroExonator quantification and discovery, and now you want to perform the these downstream alternative splicing analyses, you might need to skip microexonator steps to avoid them be re-run. To only perfom downstream alternative splicing analyses, you will need to include these extra key in config.yaml:

    downstream_only : T

If you do this, you will schedule processes that are related with `whippet` indexing, quantification and differential inclusion. For example if you have a total of 6 samples, you will see something like this:

    Job counts:
            count   jobs
            6       ME_psi_to_quant
            1       delta_ME_from_MicroExonator
            1       delta_ME_from_whippet
            1       differential_inclusion
            6       download_fastq
            1       get_GTF
            6       gzip_ME_psi_to_quant
            1       whippet_delta
            1       whippet_delta_ME
            1       whippet_index
            6       whippet_quant
            31

After you are sure everthing is in place, you can sumbit the runnig command:

    snakemake -s MicroExonator.skm  --cluster-config cluster.json --cluster {cluster system params} --use-conda -k  -j {number of parallel jobs} differential_inclusion

# Output

The main results of MicroExonator discovery and quantification modules can be found at the Results folder. All the detected microexons that passed though the quantitative filters can be found at `out.high_quality.txt`. This is a tabular separated file with 14 columns that contain the folowing information:


|Column| Description|
|:--|:-------|
|`ME` | Microexon Coordinates|
|`Transcript` | Transcript where the microexon was detected |
|`Total_coverage` | Total coverage across all microexon splice junctions|
|`Total_SJs` | Splice junctions where the micrexon was detected in|
|`ME_coverages` | Coma-separated coverage values for each microexon splice junction|
|`ME_length` | Microexon length|
|`ME_seq` | Microexon sequence|
|`ME_matches` | Microexon number of matches inside the intron|
|`U2_score` | U2 splicing score|
|`Mean_conservation` | Mean conservation values (if phylop score was provided)|
|`P_MEs` | Microexon confidence score |
|`Total_ME` | ME coordinates, U2 score and conservation for all microexon matches|
|`ME_P_value` | Value used for the final microexon filters |
|`ME_type` | Microexon type (IN, RESCUED or OUT)|

MicroExonator also reports microexons that do not meet the confidence filtering criteria. Detected microexons that are equal or shorter than 3 nt are reported at  `out_shorter_than_3_ME.txt`. Microexons that are longer than 3 nt, but did not have sufficiently low ME_P_value are reported at `out_low_scored_ME.txt`. Finally, microexons that had ME_P_values below the threshold, but they can also correspond to alternative splicing acceptors or donors (as the microexon sequence matches either at the begining or the end of an intron) are reported at `out.ambiguous.txt`.

On the other hand, microexon quantification is provided as `out_filtered_ME.PSI.txt` file. This file contains the following information:

|Column| Description|
|:--|:-------|
|`File` | Sample |
|`ME_coords` | Microexon Coordinates|
|`SJ_coords` | Splice junctions where the micrexon was detected ins|
|`ME_coverages` | Coma-separtated values corresponding to microexon converage on each splice junction|
|`SJ_coverages` | Coma-separtated values corresponding to converage on each splice junction (exon skiping)|
|`PSI` | Computed PSI values for a given microexon on a given sample|
|`CI_Lo` | Lower bound of the PSI confidence interval |
|`CI_Hi` | Upper bound of the PSI confidence interval |
|`Alt5` | Alternative donor coordinate|
|`Alt3` | Alternative aceptor coordinate|
|`Alt5_coverages` | Alternative donor coverage|
|`Alt3_coverages` | Alternative aceptor coverage |

Differentially included microexon analyses that can be obtained with Whippet, are reported at `Whippet/Delta` folder. MicroExonator performs these analyses using both PSI values calculated internally by the pipeline and PSI values directly calculated with Whippet. These results are reported under the same format than the `diff.gz` descrived at the [Whippet's GitHub page](https://github.com/timbitz/Whippet.jl#output-formats). However, to provide easier interpretation, we filter the Whippet splicing nodes that correspond to microexon inclusion events, these are reported as `.microexons` files, where `.diff.ME.microexons` files correspond to the output when MicroExonator PSI values are taken as input and `.diff.microexons` when Whippet PSI  values are taken as input.

# Troubleshooting

Check that every necesary key is writen correctly at `config.yaml`. Also, when you provide the gene annotation files, ensure that all the gene chromosomes are available at the genome fasta that you are providing. Any missing chromosome, will generate errors at different parts of the workflow.

After the release of conda 4.4 some additional configuration might be required to work with snakemake worflow that have fix conda enviroments, particulaly there might be some inteference of some conda enviroments when that pipeline call is done from an screen [interference with screen comand](https://stackoverflow.com/questions/50591901/screen-inside-the-conda-environment-doesnt-work). Thus, we recomend modify your `~/.bashrc` file and add the following line:

    . /path/to/miniconda/etc/profile.d/conda.sh

Where `/path/to/miniconda` is the directory where you installed Miniconda, which is at your home directory by default, but during Miniconda instalation can be set at any directory.

If you have any errors while you are running MicroExonator is useful to read the logs that are reported by the queuing system. Some errors may occur because when not enough memory has been allocated for a given step. Resources for each step can be this can be configured inside `cluster.json` (check example file at MicroExonator/Examples/Cluster_config/lsf/)

When some step fails, you can always resume MicroExonator run to avoid starting from scratch. To do this is important to check that there are not interrupted files. If snakemake finish or is manually interrupted (only once) it will flag interrupted files and delete them before you resume the pipeline. However, when downloading steps are interrupted, snakemake will not check whether the fastq files inside `FASTQ/` are complete or not. Thus, if you are not sure if fastq files inside  `FASTQ` are complete, please delete all fastq files inside `FASTQ/` before running MicroExonator again. Ensure avoiding re-running completed jobs you can include `rerun_incomplete_round1` or `rerun_incomplete_round2` as targets by adding either of these at the end of MicroExonator running command. Downloading errors are often due to network overload when too many fastq files are downloading at the same time, to limit the number of simultaneous downloads use get_data=x, where x is the desired maximum number of simultaneous downloading process executed by MicroExonator.

Once you do this, you run MicroExonator again as dry-mode (with `-np`), some jobs might be already completed, so you might have less remaining steps to finish the analysis.

# Contact

For questions, ideas, feature requests and potential bug reports please contact gp7@sanger.ac.uk.
