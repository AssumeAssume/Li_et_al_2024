[![Snakemake](https://img.shields.io/badge/snakemake-≥3.6.0-brightgreen.svg?style=flat-square)](http://snakemake.bitbucket.org)
[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/mageck-vispr/badges/downloads.svg)](https://anaconda.org/bioconda/mageck-vispr)

Quality Control, Modeling and Visualization of CRISPR/Cas9 Screens with MAGeCK-VISPR
====================================================================================
MAGeCK-VISPR is a comprehensive quality control, analysis and visualization workflow
for CRISPR/Cas9 screens.
For recent changes, see the [change log](CHANGELOG.md).
The workflow combines

* the [MAGeCK](http://mageck.sourceforge.net) algorithm to identify essential
genes from CRISPR/Cas9 screens considering multiple conditions
* with [VISPR](https://bitbucket.org/liulab/vispr) to interactively explore results and quality control in a web-based frontend.

Please visit the [MAGeCK](https://bitbucket.org/liulab/mageck) and [VISPR](https://bitbucket.org/liulab/vispr) pages for details about these tools. In particular, visit [VISPR](https://bitbucket.org/liulab/vispr) if you only want to test the interactive visualization.


Join our [MAGeCK-VISPR google group](https://groups.google.com/d/forum/mageck) if you have any questions or comments. 


Table of Contents
-----------------

* Prerequisites
* Running MAGeCK-VISPR Docker Images
  - Running on AWS EC2 instance
* Installation
* Testing MAGeCK-VISPR with two examples
* Video Tutorials
* Authors
* License



Prerequisites
-------------

The input of the MAGeCK-VISPR workflow are raw FASTQ files.
If you already have [MAGeCK](http://liulab.dfci.harvard.edu/Mageck/)
results or only want to test the interactive visualization, you can directly use [VISPR](https://bitbucket.org/liulab/vispr) and don't need to install the whole workflow.
Further, it is recommended to execute this workflow on a linux server or cluster rather
than your laptop or workstation, since it is computationally expensive.


Running MAGeCK-VISPR Docker Images
------------

Don't want to install MAGeCK-VISPR by yourself and want to try the latest MAGeCK-VISPR? [MAGeCK-VISPR Docker Image](https://hub.docker.com/r/davidliwei/mageck-vispr) provides an easy and reproducible way to run MAGeCK-VISPR pipeline. 
The image is automatically built upon each commit in our bitbucket source code.

To run through Docker image, you need to first install [Docker](https://www.docker.com) on your own system. 
Then type (get docker image via internet)

    docker pull davidliwei/mageck-vispr
or (build docker image via Dockerfile)
    
    cd /your_path_to/mageck-vispr
    docker build -t davidliwei/mageck-vispr . 
and

    docker run -it -d -v /your_path_to/data_folder:/container_path_to/data_folder -p 5000:5000 --cpus="2" --privileged=true --name mageck-vispr davidliwei/mageck-vispr bash

Where

- -v: set up sharing directory between host and container
- -p: vispr uses port 5000 by default to run web server
- --cpus: optional, do not exceed your cpu cores of host machine

to run MAGeCK-VISPR Docker. For more options (e.g., data access), refer to our [MAGeCK Docker tutorial](https://sourceforge.net/p/mageck/wiki/demo/#tutorial-3-run-mageck-on-docker-image) for more details.

to enter the container on **macOS** or **Linux**

    docker exec -it mageck-vispr bash

to enter the container on **Windows**

    winpty docker exec -it mageck-vispr bash

- another way of playing with Docker on **Windows** is using the **Docker Hub**, which has a GUI for manipulating docker

to test vispr in Docker, type command

    vispr test --host 0.0.0.0

and open the page at http://127.0.0.1:5000 in your **host machine's** browser.

To stop the Docker service, use

    docker stop mageck-vispr


## Run VISPR Docker image on AWS and access it in the local browser

1. When run the container on AWS, please add the parameter --net="host"


    sudo docker run -it -d -v /data:/data -p 5000:5000 --cpus="16" --net="host" --privileged=true --name mageck-vispr davidliwei/mageck-vispr bash


2. When access AWS server by SSH command, please add the parameter -L5000:localhost:5000


    ssh -L5000:localhost:5000 -i your-key.pem ec2-user@address


3. Next, type ``http://127.0.0.1:5000`` on your local browser to access the server.

Installation
------------

To install MAGeCK-VISPR you have to use the **Python 3** variant of the
Miniconda Python distribution (http://conda.pydata.org/miniconda.html).
MAGeCK-VISPR cannot be installed on the Python 2 variant of Miniconda.
When installing Miniconda, make sure that you answer ``yes`` to this question:

    Do you wish the installer to prepend the Miniconda3 install location to PATH ...? [yes|no]

Also, make sure that you do not have set the ``PYTHONPATH`` environment variable, because it will interfere with the Miniconda setup.




You can create an isolated software environment for mageck-vispr by executing

    conda create -n mageck-vispr  python=3.7

(We found dependency errors on python 3.8 but 3.7 works well).

Then, activate the environment by running

    conda activate mageck-vispr

To install mageck-vispr into the current channel, use

    conda install -c bioconda -c conda-forge mageck-vispr

in a **terminal**.

**Important: conda has a known issue of taking very long time solving software dependencies in some cases.** 
If you had this issue on your computer, we recommend using [mamba](https://github.com/QuantStack/mamba) instead. Here are the steps:

First, install mamba using

    conda install -c conda-forge mamba
    
Then, use mamba instead of conda to install mageck-vispr:

    mamba install -c bioconda -c conda-forge mageck-vispr




To update mageck-vispr, run

    conda update mageck-vispr

from within the environment.

This environment can be deactivated via

    source deactivate

If you are using an old version of MacOS X and the `conda` command is not available after installation of Miniconda, you have to change your shell to `bash`. To do this permanently, issue

    chsh -s /bin/bash

Note: We recommend setting up **Bioconda** as described [here (step 2)](http://bioconda.github.io). **Make sure to include bioconda and conda-forge as the default channel.**

Usage
-----

All steps below have to be executed in a **terminal**.

### Testing MAGeCK-VISPR

We provide test data for exploring the interactive visualization interface of VISPR at the [VISPR homepage](https://bitbucket.org/liulab/vispr).

The MAGeCK-VISPR workflow is harder to test, because input data (raw FASTQ files) is usually big. However, we provide a subsampled test dataset based on [Yusa et al. (Nat. Biotechn. 2014)](http://www.ncbi.nlm.nih.gov/pubmed/24535568). While this test data won't yield reasonable results, it can be used to see how the workflow is configured and executed. There are two test datasets available:

* configuring the workflow from scratch (Step 2, see Usage): [download](https://bitbucket.org/liulab/mageck-vispr/downloads/esc.testdata.step2.tar.bz2)
* running a configured workflow (Step 4, see Usage): [download](https://bitbucket.org/liulab/mageck-vispr/downloads/esc.testdata.step4.tar.bz2)

Please choose from which step you want to test the MAGeCK-VISPR, download and extract the corresponding dataset using tar. For example,

     tar xvf esc.testdata.step2.tar.gz2

After that, start with the corresponding step from the Usage section below.

For usage of MAGeCK-VISPR on real data, please see all steps below.

### Step 0: Activate the mageck-vispr environment

If you have installed mageck-vispr as shown above, make sure to activate the corresponding conda environment before conducting any mageck-vispr related command via

    source activate mageck-vispr

Conda environments allow you to use mutliple versions of the same software without conflicts. E.g., mageck-vispr uses Python 3.5, which might be not the default Python version in your system.

### Step 1: Choosing a workflow directory

To analyze a particular set of raw FASTQ files representing a CRISPR/Cas9 screening experiment, MAGeCK-VISPR first requires you to select a directory where the workflow shall be executed and results will be stored.
Please choose a meaningful name and ensure that the directory is empty.

### Step 2: Initializing a new workflow

We assume that you have chosen a workflow directory (here `path/to/my/workflow`), and that you have a set of FASTQ files (`path/to/sample1.fastq path/to/sample2.fastq ...`) you want to analyze.
You can initialize the workflow with

    mageck-vispr init path/to/my/workflow --reads path/to/sample1.fastq path/to/sample2.fastq ...

which installs a `README`, a config file `config.yaml` and a
[Snakemake](https://bitbucket.org/johanneskoester/snakemake) workflow definition
(a `Snakefile`) to the given directory.
For further information about the mageck-vispr command issue

    mageck-vispr --help

for further information.

### Step 3: Configure the workflow

If you provided the FASTQ files to mageck-vispr in above command, the config file already contains a mapping between putative sample names and the file paths.
Now change into your workflow directory via

    cd path/to/my/workflow

and open the `config.yaml` file. Edit the config file to your needs. Especially, define experiments for use with MAGeCK.
Here, you can choose between providing treatment and control samples or a design matrix. See the [MAGeCK homepage](http://liulab.dfci.harvard.edu/Mageck/) for details.

### Step 4: Execute the workflow

Once configured, the workflow can be executed with [Snakemake](https://bitbucket.org/johanneskoester/snakemake).
First, it is advisable to invoke a dry-run of the workflow with

    snakemake -n

which will display all jobs that will be executed. If there are errors in your config file, you will be notified by Snakemake.
The actual execution of the workflow can be started with

    snakemake --cores 8

here using 8 CPU cores in parallel.
For more information on how to use Snakemake, refer to the
[Snakemake Documentation](https://bitbucket.org/johanneskoester/snakemake/wiki/Documentation#markdown-header-executing-snakefiles).

### Step 5: Visualize results with VISPR

Once workflow execution has finished, you can visualize the generated results
and quality controls by issueing

    vispr server results/*.vispr.yaml

in your workflow directory.
For more information and VISPR options, please visit the [VISPR](https://bitbucket.org/liulab/vispr) homepage.

Video Tutorials
-------------

We provide a set of YouTube video tutorials to help you go through the usage of MAGeCK-VISPR (thanks to Anya Zhang). The tutorial includes videos for installing MAGeCK-VISPR, exploring quality control and comparison results, and running the whole pipeline:

* [Installation](https://www.youtube.com/watch?v=IslJQyJb_jE)
* [Pipeline](https://www.youtube.com/watch?v=3maSxhy1JL0)
* [Quality control](https://www.youtube.com/watch?v=-Um-Mrf4_7s)
* [Results](https://www.youtube.com/watch?v=VI-8Z-fYxh4)



Authors
-------

* Johannes Köster <koester@jimmy.harvard.edu>
* Wei Li <wli2@childrensnational.org>

License
-------

Licensed under the MIT license (http://opensource.org/licenses/MIT). This project
may not be copied, modified, or distributed except according to those terms.
