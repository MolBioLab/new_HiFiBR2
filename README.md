# new-HiFiBR2 

new-HiFiBR2 is an automatic analysis pipeline for double-strand-breaks (DSBs)-repaired junctions. This pipeline uses amplicon sequences which consist of repair junctions as input. Recently, our tool only performs NGS data. The final result is an excel table that describes the main features of each junction in a pool. 

new-HiFiBR2 is developed based on HiFiBR software (check out HiFiBR published [paper](https://doi.org/10.1016/bs.mie.2017.11.028) and [GitHub](https://github.com/Alexander-Brown13/Hi-FiBR)). Notably, the pipeline is modified in order to handle sequences generated from our study model. The below table outlines certain features in our study model compared to the recently common model in studying repair junction and corresponding additional analysis steps. 

Additional features | Our study model | Common study model | 
| ----------------- | --------------- | ------------------ | 
| Analyze junction preciesly with parameter for 2 cutting sites position | 2 distal DSBs (induced by ISce-I) | 1 DSB  | 
| Mapping reads with optimized parameter | Large size junction (with insertion / deletion >10 bp)  | Small size junction |
| Filter out sequences do not consist of full junction | Sequencing RE-cut products | Sequencing PCR products | 

This repo contains the code for analyzing experimental data in our published [paper]()


## Requirements
The new-HiFiBR2 currently runs on UNIX systems (tested on Ubuntu 16.04 LTS, x86_64 GNU / UNIX)
#### Required software for running pipeline:
- <span style="color:blue"> Ubuntu </span>
- <span style="color:blue"> python >= 3.6.15</span>. Most of the workflow scripts were built in the programming language Python.
- <span style="color:blue"> trimmomatic == 0.39</span>
- <span style="color:blue"> flash2 == 2.2.00</span>
- <span style="color:blue"> Geneious Prime 2022 or later</span>. Geneious is the software for sequence data analysis. This workflow uses Geneious Mapper for mapping reads with large size junction.
- <span style="color:blue"> samtools == 1.11.0</span>
- <span style="color:blue"> bamtools == 2.4.0</span> 

#### Required Python packages:
- pysam == 0.19.1
- pandas == 1.1.5
- numpy == 1.19.5
- scikit-learn == 0.24.2
- scipy == 1.5.3
- diff-match-patch
- tqdm == 4.46.0
- xlsxwriter == 3.0.3

## Installation
```
git clone https://github.com/MolBioLab/new_HiFiBR2

# If your machine already installed required tools and Python packages above 
(except Geneious Prime tool, because it's a purchase software), then pass 
below codes

# Setup conda environment
conda update conda 

# Create newHiFiBR2 environment with all required tools and Python packages
conda env create -f newHiFiBR2_env.yml
```
## Run new-HiFiBR2
### Input
- Config file
  - Config file for *test run* in defaults/config_test.ini. 
  - Config file *template* for your run in bin/config.ini.
- Sequence files. **Note:** Tool gets FASTQ format with .fastq, .fq, .fastq.gz, .fq.gz extension.
  - Get sequence file for *test run* in defaults/mini_simu_SE.fq.gz
- Reference sequence file with FASTA format <span style="color:blue"> (NOT multi-FASTA) </span>. This is the sequence of experimental subtrate (e.g. CD4+ subtrate in our study) from forward primer to reverse primer.
  - Ref seq file for *test run* in defaults/ref.fa.

### Excuting pipeline
- Create a new folder (e.g. test).
- Create a new configuration file and save in new folder. **Important:** You should fill in these variables:
  - <span style="color:blue"> trimmomatic_path </span>. Copy path to trimmomatic.jar file. If you're using trimmomatic in conda newHiFiBR2 environment, you could use the following path path/to/miniconda3/envs/newHiFiBR2/share/trimmomatic-0.39-2/trimmomatic.jar
  - <span style="color:blue"> ref_path</span>. Copy path to reference sequence FASTA file
- Copy sequence files into new folder. 
```
# User manual
$ bash new_HiFiBR2.sh 
Usage: new_HiFiBR2.sh [options]

-w|--working-dir		Working directory which contains .fq files and config file. 
						And also output files destination
-c|--config-file    	Name of config-file
-s|--software-path   	Path to new_HiFiBR2 software
-h|--help

# Test run
$ bash new_HiFiBR2.sh -w /path/to/test/folder -c config_test.ini -s /path/to/software/new_HiFiBR2

```

### Output
Features of each type of junctions in the pool of experimental sample. The below table lists the columns in the output excel file and the details (see an example in defaults/mini_simu_Padded_Final.xlsx). 

**Note:** If "Status" column is marked as "new", then these are advanced features in our pipeline. But if the status is blank, the HiFiBR pipeline has already analyzed these features and the explaination is refer to Table 1 in [published paper](https://doi.org/10.1016/bs.mie.2017.11.028).

| Column | Explaination | Status |
| ------ | ------------ | -------|
| Ref seq name | Name of reference sequence to which read was mapped | |
| Original CIGAR | Original CIGAR string (CIGAR from mapping result) | |
| Len of read | Length of read | |
| Recontructed CIGAR | Original CIGAR string; later may contain a reconstructed CIGAR string | |
| Len of matching left | Corrected length of sequence matching the left flank of the reference | |
| Len of matching right | Corrected length of sequence matching the right flank of the reference | |
| Distance: break to left matching | Distance between the break and the end of the left matching sequence | |
| Distance: break to right matching | Distance between the break and the end of the right matching sequence | |
| Nu del left side | Number of nucleotides deleted on the left side of the break | |
| Nu del right side | Number of nucleotides deleted on the left side of the break | |
| Total del nu | Total number of deleted nucleotides | |
| Start of ins | Start of inserted sequence | |
| End of ins | End of inserted sequence | |
| Len of ins | Length of inserted sequence | |
| Seq of ins | Sequence of inserted nucleotides |
| Repair event | Repair event class | |
| Reconstructed seq | Reconstructed read sequence | |
| Times of event (popu) | Number of times reconstructed read occurred in read population | |
| Microhomo seq | Potential microhomology that mediated deletion | |
| Len of microhomo | Length of microhomology | |
| Match_mismatch MH | Match (M) / Mismatch (u) of MH compare to reference | new |
| Times of event (file) | Number of times reconstructed read is observed in this file (in a row) | |
| % read contain mismatch | Percentage of reads that contain mismatches to reconstructed read | |
| Represented read | Header ID of a read that represents for the considering type of juction | new |
| No. Time of event (file) | No. group of reconstructed read | new |	
| Maybe same event (if any) | This type of junction maybe belongs to other junction | new |

