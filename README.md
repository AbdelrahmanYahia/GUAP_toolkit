# **Genomics Unit Analysis Pipelines (GUAP) Toolkit**
# **!!!** This tools is still under development and **shouldn't be run!!**

The Genomics Unit Analysis Pipelines (GUAP) toolkit is a collection of workflows and pipelines designed for routine NGS (Next-Generation Sequencing) work in the field of genomics. GUAP is built to provide reproducible analysis while allowing fine-tuning of workflows based on sample conditions. It is fast, easy to use, and operates on the Linux operating system via a command-line interface.

## **Features**

- Reproducible analysis with customizable workflows
- Fast and easy-to-use command-line interface
- User-defined input directory, output directory, and analysis pipeline
- Command-line options to specify workflow parameters

## **Usage**

A basic GUAP command may look like this:

```bash
guap 16s -i indir -o outdir -m sample-metadat.csv -c classifier.qza
```

However, GUAP provides flexibility for more advanced specifications:

```bash
guap 16s -i samples -o out -m sample-metadata.tsv \
	-c silva/Silva_Classifier.qza --skip-QC \
	 -t 65 -mo 4 -n run1 -cm pooled  \
	--export-figs --condition-name condition \
	-ef 9 -er 10 --remove-primers
```

## **Supported Workflows**

GUAP currently supports the following workflows:

- 16s rRNA analysis
- RNAseq differential expression
- Whole Exome Sequencing (WES)
- Bacterial Whole Genome Sequencing (WGS)
- Metagenomics

Each workflow comes with default options and allows modification of individual steps, including the ability to change aligners, quantifiers, variant callers, and more. Each step (tool) in the workflow is customization and extra args can be added

## **Snakemake Integration**

The core engine behind GUAP's workflows is Snakemake, which provides advanced features such as parallelization, speed, efficiency, and the ability to resume workflows in case of failure. GUAP leverages Snakemake's capabilities to deliver powerful and customizable analysis pipelines.

Using output names and Python functions, workflows can be easily customized. GUAP also parses Snakemake's standard output to generate a progress bar, keeping users informed about the ongoing analysis in a clear way.

## **Sample Parsing and Configuration**

GUAP automatically parses sample naming patterns and generates sample tables. It also parses user arguments and generates configuration files for Snakemake workflows. This ensures smooth running of the workflows, maintaining compatibility with the generated config and sample tables.

## **Conda Environment and Reproducibility**

By default, GUAP utilizes Conda to manage dependencies and ensure reproducibility. Each workflow in GUAP has its own isolated Conda environment to overcome package conflicts and guarantee consistency across analyses. These Conda environments can be reused within GUAP or outside of it.

## **Additional Functionality**

GUAP incorporates various modules, classes, and functions written in Python to implement sample parsing logic and Snakemake workflow management. Moreover, each workflow includes custom scripts that automate certain aspects of the analysis. For example, in bacterial whole-genome sequencing, GUAP provides a script that utilizes Kraken2 to identify species in the sample and automatically downloads the reference from NCBI using Entrez.

## **Ease of Use and Problem Solving**

Great thought and effort have been invested in GUAP to ensure it is user-friendly, fast, reliable, and reproducible. GUAP aims to address common bioinformatics challenges faced by researchers, such as package installation, workflow creation, and testing different approaches.

During the development of GUAP, the creator aimed to automate tasks encountered bioinformatic researchs. GUAP incorporates solutions to various challenges, including testing different settings and tools to determine the best fit for the samples. With GUAP, users can modify the entire workflow with just a few arguments, providing incredible flexibility and adaptability.

With GUAP, there's no need to create sample sheets, modify configurations, or install packages manually. During installation, GUAP automatically installs most required packages in Conda. These environments can be reused within GUAP or independently.

## **Installation**

Installation instructions and detailed usage will be updated once the tools is finished

## **Contact**

For any inquiries or questions, please contact abdelrhmanahmedyahia@gmail.com
