# Protamine Tools

Protamine Tools provides code used in portions of the analysis of spermatid
ATAC-seq and Cut&Tag data and its relationship with other transcription and
chromatin data sets, as reported in 
Rabbani et al. (in preparation), 
a collaboration between the 
[Sue Hammoud](https://hammoud.lab.medicine.umich.edu/)
and 
[Tom Wilson](https://wilsonte-umich.github.io/)
laboratories.

## Installation

Protamine Tools is implemented in the
[Michigan Data Interface](https://midataint.github.io/) (MDI),
a framework for developing, installing and running 
Stage 1 HPC **pipelines** and Stage 2 interactive web **apps**.
For this tool suite, we use and recommend a multi-suite 
installation (see the [MDI documentation](https://midataint.github.io/) 
for single-suite installations), which is accomplished by:
- cloning and installing the MDI framework
- adding this tool suite (and potentially others) to your MDI installation
- calling the `mdi` command line interface (CLI) to use tools from any installed suite

### Install the MDI framework

Change into a folder of your choosing and run the following commands to
obtain and install the MDI framework. To start, choose installation
option 1 to only install the Stage 1 HPC pipelines framework - you can return
to install the R Shiny apps interface later, which takes longer (see below).

```bash
git clone https://github.com/MiDataInt/mdi.git
cd mdi
./install.sh
```

### OPTIONAL: Add an mdi alias to .bashrc

The following commands will create a permanent named alias, called `protaminer`, 
to the `mdi` target script in your new installation.

```bash
./mdi alias --help
./mdi alias --alias protaminer # change the alias name if you'd like
```

You will need to open a new shell for the alias to become active.

If you prefer not to use an alias, add the installation directory to your PATH variable,
or `cd` into the directory prior to calling `./mdi`, instead of `protaminer` as in commands below.

### Test the command line interface (CLI)

For top-level help, call the CLI or one of its commands with no arguments
or with the `--help` option. 

```bash
protaminer # or ./mdi if you didn't create an alias
protaminer --help
protaminer <command>
protaminer <command> --help
```

### Add the Protamine Tools suite to your MDI installation

The following commands will install Protamine Tools into your MDI installation.

```bash
protaminer add --help
protaminer add -s HammoudWilson/protamine-tools 
```

Alternatively, perform the required steps manually:

```sh
nano config/suites.yml # to edit suites.yml as below
```

```yml
# mdi/config/suites.yml
suites:
    - HammoudWilson/protamine-tools  # add this tool suite
```

```sh
protaminer install --help
protaminer install # re-install to obtain the new tool suite
```

Finally, see the tools that are available in your installation:

```sh
protaminer list
```

### Build the required Conda environments

Protamine Tools pipelines use version-controlled 3rd-party software built 
into conda environments. Build those environments in your MDI installation as follows:

```sh
protaminer nascent conda --create # for preparing round spermatid Pro-seq data
protaminer atac conda --create    # for executing ATAC-seq data analysis
```

You must have `conda` available in your server environment. If you need to run
a command to make conda available, follow the instructions in
`.../mdi/config/stage1-pipelines.yml`, which is pre-configured to work on
the University of Michigan Great Lakes cluster.

HINT: on our shared server environment, we need to build the conda environments
on a worker node, as the conda build process exceeds memory limits on our login nodes.

### Obtain additional external resources and data sets

The following external resources and data sets are required to execute Protamine Tools pipelines.

Round spermatid Pro-seq nascent transcription data from
[Kaye et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38287033/)
must be [downloaded from NIH GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228452), specifically files:
- GSE228452_mm10_rs_proseq_n3spike_F.bw
- GSE228452_mm10_rs_proseq_n3spike_R.bw

Because the Pro-seq data are in mm10 genome coordinates, you must also
download the 
[mm10 to mm39 liftover chain from UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz).

Paths to the above files are passed to pipelines using options
`--bigwig-file-forward`, 
`--bigwig-file-reverse`, and 
`--liftOver-chain`.

### Create required directory paths

If you will use default values for genome options (recommended), 
you should manually create the following directory:

```sh
mkdir -p /path/to/mdi/resources/genomes # adjust the path as needed
```

## Execute pipelines using the CLI

### Pipeline and action workflow 

There are five named pipelines in the Protamine Tools suite, used in this order:

- `genome` = download, assemble, and bin a mouse-fly composite genome for alignment
- `nascent` = prepare round spermatid nascent RNA-seq data for analysis
- `atac` = analyze spermatid ATAC-seq data, prepare packages for use in interactive apps
- `cuttag` = analyze spermatid Cut&Tag data to extend ATAC-seq analysis
- `enrichment` = assemble BED files of genome regions for ATAC-seq enrichment analysis

Each pipeline has a series of named actions executed in the order,
and with the dependencies, in the following diagram:

![Protamine Tools Workflow](https://github.com/HammoudWilson/protamine-tools/blob/main/protaminer_workflow.png)

Details of the execution and options for each action are beyond the scope
of this README but can be listed using the CLI, e.g.

```bash
protaminer genome download --help # show the options help for one pipeline action
protaminer genome template --all-options # output a job file template for all pipeline actions
# etc.
```

Details of the work done by each pipeline and action are found in the scripts
organized into sub-folders of the top-level `pipelines` folder, 
of the same names as the pipelines and actions themselves. 
Each pipeline's `pipeline.yml` file describes its actions,
and each action's `Workflow.sh` script coordinates the work it accomplishes.

Note that Protamine Tools requires specifically configured genome resource files,
which you must obtain and assemble using `genome download`, etc. Some actions
of the `genome` pipeline are slow but only need to be performed once per genome,
independently of any sample data to be analyzed. Similarly, the `nascent` pipeline
only needs to be executed once on the Pro-seq data files listed above. 

The `cuttag` pipeline is unusual in that it depends on the output of the
`atac collate` pipeline action, but once `cuttag score` has been executed, you
must re-run `atac collate` and `atac tss` to update the data packages used for
data visualizaiton in the apps.

### Universally required options

Options `--output-dir/-O` and `--data-name/-N` are required by all pipeline actions.

Several of the pipelines create data package files to be loaded into the
R Shiny interactive visualization apps, in addition to various BAM, BED and
other formatted files, all placed into directory `<--output-dir>/<--data-name>`.

### Sample metadata file(s)

You must provide a path to an ATAC-seq sample metadata file as option 
`--metadata-file`. The sample file we used is included in the downloaded repository
as file `mdi/suites/definitive/protamine-tools/resources/ATAC_sample_file_v6.csv`.
It included various names and spermatid stage descriptions. If you create your 
own file, it must follow the same column format.

Addtionally, to run the `cuttag score` pipeline action, you must provide 
a path to an Cut&Tag sample metadata file as option `--metadata-file`. The sample
file we used is again included in the downloaded repository
as file `mdi/suites/definitive/protamine-tools/resources/CutAndTag_sample_file_v2.csv`.
Some columns are shared with ATAC-seq samples file while others, like `antibody target`,
are specific to Cut&Tag.

### Job files

Protamine Tools pipelines can be called entirely using the CLI. However, you 
are strongly encouraged to create YAML-format job configuration files that define the
parameters for your job and coordinate execution steps.

See 
the [protamine-tools/templates](https://github.com/HammoudWilson/protamine-tools/tree/main/templates) folder
for job file templates for all pipelines, which are the same as the job files we used except for file paths, 
and <https://midataint.github.io/mdi/docs/job_config_files.html>
for extended help on using job files. Job file templates can also be generated with 
command `protaminer <pipeline> template` as in the example above.

The MDI CLI and job files can run pipeline actions either 
inline in the calling shell or by submitting jobs to your server job scheduler,
which is recommended for most use cases. Thus, our most common usage pattern is:

```sh
protaminer inspect myJob.yml          # check the formatting of your job file
protaminer mkdir myJob.yml            # create any missing output directories
protaminer submit --dry-run myJob.yml # test the job file to see what will happen
protaminer submit myJob.yml           # submit the job to Slurm or your scheduler
protaminer myJob.yml status           # show the state of all submitted jobs
protaminer myJob.yml top              # monitor a running job
protaminer myJob.yml report           # show a job log report
protaminer myJob.yml ls               # show the contents of a job's output diretory
```

### Example status and log files

The [protamine-tools/logs](https://github.com/HammoudWilson/protamine-tools/tree/main/logs) folder
has files with the output of our execution of the different pipelines, as 
generated by `protaminer status` followed by `protaminer report`. 
These are the program logs as we ran them except that environment-specific values 
such as accounts and file paths are replaced with dummy values. 
Comparison to your logs may help identify the source of problems.

## Launch the interactive apps server

Protamine Tools uses interactive R Shiny apps for data visualization and plotting.
To install and launch the Protamine Tools interactive apps server, 
we recommend using the 
[MDI Desktop app](https://midataint.github.io/mdi-desktop-app),
which allows you to control both local and remote MDI web servers.

After following the instructions to run the Desktop on your local machine
or server, load a Protamine Tools data package file ending in `mdi.package.zip`,
into the app interface. The two most important data packages are:
- `xxx.atac.collate.mdi.package.zip` = sample and aggregate-level summary plots
- `xxx.atac.tss.mdi.package.zip` = plots based on transcription analysis, including a genome Track Browser

You can also install the apps framework into your MDI installation from the command line:

```bash
protaminer install --install-packages --n-cpu 4
```

An active R installation must be present on your server when installing
the apps framework. Note that you will still want to install and use the 
[MDI Desktop app](https://midataint.github.io/mdi-desktop-app)
to run the apps server even if you install the apps framework at the command line.
