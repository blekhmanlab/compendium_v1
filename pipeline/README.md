# Human Microbiome Compendium pipeline v1

The files in this directory illustrate the process used to download, analyze and integrate the samples used in version 1 of the Human Microbiome Compendium. **This will not be maintained or updated**—while future iterations will be developed as reusable applications, this is intended only as a supplement to our publication describing the data. To improve readability, some very general edits have been made: Primarily, many of the steps in our process were written in shell scripts intended for the Slurm HPC scheduler. Most were uninformative in this context and have been excluded. Two of them, `run_dada.sh` and `run_fasterq.sh`, were reformatted as conventional bash scripts, but the commands are exactly as executed.

## Downloading metadata

The `main.py` file has multiple functions that accommodate different steps in the process of pulling data out of the Sequence Read Archive and putting it into the samples database.

* **`load_xml()`:** This loads an XML file **exported from a BioSample search** and puts the sample metadata into the database.
  * The BioSample search for this project was `txid408170[Organism] AND biosample sra[filter] AND "public"[filter]`
* **`find_runs()`:** This step uses the sample metadata (*already in the database*) to find *runs* ("SRR" codes), which are the entities actually associated with sequencing data. This makes lots of API calls to NCBI.
* **`write_lists()`:** This function fetches the SRA projects and generates text files containing runs for each project. **These accession lists are the input for the processing pipeline.**

The `get_instrument.py` file contains a script to augment the initial data recorded by main.py. This calls the NCBI API to retrieve sample-level metadata about sequencing instruments. Three other files, `db.py` and `config.py`, contain helper functions used in the other Python scripts. `requirements.txt` lists the dependencies used by these scripts.

## Fetching, processing SRA files

1. The `run_fasterq.sh` script takes a project ID as a parameter and uses the SRA Toolkit command-line tool to download the SRA files and convert them to FASTQ files. A final step standardizes the file names of single- and paired-end projects.

2. `run_dada.sh` is triggered at the end of the previous step. This launches a Singularity container with DADA2 installed in it and runs one of two scripts:
   * `process_forwards.R`: Generates an ASV table and associated taxonomic inferences for single-end reads.
   * `process_paired_end.R`: Follows the same process as the single-end script, but adds an additional step to merge paired-end reads. This script includes an option to be re-run in "matchids" mode, which uses read-level metadata to attempt to find the mate pairs for FASTQ files in which the order of the forward reads is different than the order of the reverse reads.

3. The `consolidate_taxa.py` script loads the results from a single project at a time and consolidates ASV-specific read counts into read counts grouped by the taxonomic assignment of each ASV. The result is that each completed project has a new taxonomic table without any duplicate column names.

4. `combine_studies.py` takes the taxonomic tables from the previous step and consolidates them together into a single table.
