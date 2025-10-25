# Protamine Tools

Protamine Tools provides code used in portions of the analysis of spermatid
ATAC-seq data and its relationship with other transcription and
chromatin data sets, as reported in Rabbani et al., a collaboration
between the 
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
to install the R Shiny apps interface later, which takes considerably longer (see below).

```bash
git clone https://github.com/MiDataInt/mdi.git
cd mdi
./install.sh
```

### OPTIONAL: Add an mdi alias to .bashrc

These commands will create a permanent named alias to the `mdi`
target script in your new installation.

```bash
./mdi alias --help
./mdi alias --alias mdi # change the alias name if you'd like 
```

You will need to open a new shell for the alias to become active.

If you prefer not to use an alias, 
you can add the installation directory to your PATH variable,
or `cd` into the directory prior to calling `./mdi`, instead of just
`mdi` in commands below.

### Test the command line interface (CLI)

For top-level help, call the CLI or one of its commands with no arguments
or with the `--help` option. 

```bash
mdi # or ./mdi if you didn't create an alias
mdi --help
mdi <command>
mdi <command> --help
```

### Add the Protamine Tools suite to your MDI installation

The following commands will install Protamine Tools into your MDI installation.

```bash
./mdi add --help
./mdi add -p -s HammoudWilson/protamine-tools 
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
mdi install --help
mdi install # re-install to obtain the new tool suite
```

### Build the required Conda environments

Protamine Tools pipelines use version-controlled 3rd-party software built 
into conda environments. Build those environment in your MDI installation as follows:

```sh
mdi nascent conda --create # for preparing round spermatid Pro-seq data
mdi atac conda --create    # for executing ATAC-seq data analysis
```

You must have `conda` available in your server environment. If you need to run
a command to make conda available, follow the instructions in
`.../mdi/config/stage1-pipelines.yml`, which is pre-configured to work on
the University of Michigan Great Lakes cluster.

## Execute pipelines using the CLI

### Pipeline and action workflow 

There are four named pipelines in the Protamine Tools suite, used in this order:

- `genome` = download, assemble, and bin a mouse-fly composite genome for alignment
- `nascent` = prepare round spermatid nascent RNA-seq data for analysis
- `atac` = analyze spermatid ATAC-seq data, prepare packages for use in interactive apps
- `enrichment` = assemble BED files of genome regions for ATAC-seq enrichment analysis

Each pipeline has a series of named actions executed in the order,
and with the dependencies, in the following diagram:

![Protamine Tools Workflow](https://github.com/HammoudWilson/protamine-tools/blob/main/protaminer_workflow.png)

Details of the execution and options for each action are beyond the scope
of this README but can be listed using the CLI, e.g.

```bash
mdi genome download --help # show the options help for one pipeline action
mdi genome template --all-options # output a job file template for all pipeline actions
# etc.
```

### Universally required options

Options `--output-dir/-O` and `--data-name/-N` are required by all pipeline actions.

Several of the pipelines create data package files to be loaded into the
R Shiny interactive visualization apps, in addition to various BAM, BED and
other formatted files, all placed into directory `<--output-dir>/<--data-name>`.

### Job files

Protamine Tools pipelines can be called entirely using the CLI. However, you 
are strongly encouraged to create YAML-format job configuration files that define the
parameters for your job and coordinate execution steps.

See 
[the templates folder](https://github.com/HammoudWilson/protamine-tools/tree/main/templates)
for job file templates for all pipelines and actions, and 
<https://midataint.github.io/mdi/docs/job_config_files.html>
for extended help on using job files. Job file templates can also be generated with 
command `mdi <pipeline> template` as in the example above.

The MDI CLI and job files can run pipeline actions either 
inline in the calling shell or by submitting jobs to your server job scheduler,
which is recommended for most use cases. Thus, our most common usage pattern is:

```sh
mdi inspect myJob.yml          # check the formatting of your job file
mdi mkdir myJob.yml            # create any missing output directories
mdi submit --dry-run myJob.yml # test the job file to see what will happen
mdi submit myJob.yml           # submit the job to Slurm or your scheduler
mdi myJob.yml status           # show the state of all submitted jobs
mdi myJob.yml top              # monitor a running job
mdi myJob.yml report           # show a job log report
mdi myJob.yml ls               # show the contents of a job's output diretory
```

## Launch the interactive apps server

Protamine Tools uses interactive R Shiny apps for data visualization and plotting.
To install and launch the Protamine Tools interactive apps server, 
we recommend using the 
[MDI Desktop app](https://midataint.github.io/mdi-desktop-app),
which allows you to control both local and remote MDI web servers.

After following the instructions to run the Desktop on your local machine
or server, load a Protamine Tools data package file ending in `mdi.package.zip`,
into the app interface. 

You can also install the apps framework into your MDI installation from the command line:

```bash
mdi install --install-packages --n-cpu 4
```

An active R installation must be present on your server when installing
the apps framework. Note that you will still want to install and use
the 
[MDI Desktop app](https://midataint.github.io/mdi-desktop-app)
to run the apps server even if you install the apps framework at the command line.
