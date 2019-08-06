# SNP CRISPR
SNP-targeted CRISPR design pipeline for design of sgRNA's in Fly, Human, Mouse, Zebrafish, and Rat genomes.  
Online version available at: [https://www.flyrnai.org/tools/snp_crispr](https://www.flyrnai.org/tools/snp_crispr)

## Prerequisites
- Anaconda or Miniconda for environment management
- Species genome in fasta format

**Note:** To install without conda, see `environment.yml` for dependencies

## Install
- From the project directory run `conda env create -f environment.yml`
- Then `conda activate snp_crispr`
- FlyBase release 6.28 genome fasta file and blast database are included for testing. `dm.fasta.zip` must be unzipped before running. Remaining installation steps are only necessary for other species
- Download species genome fasta file and put in `fasta_files/<species>.fasta`
- Run `makeblastdb -in fasta_files/<species>.fasta -dbtype nucl -parse_seqids -title <species> -out blast_dbs/<species>`
- Species chromosome names -> fasta id mapping files included in `fasta_files/<species>_chr_ids.txt`

**Important Note:** `<species>` name/abbreviation in `fasta_files/<species>.fasta`, `fasta_files/<species>_chr_ids.txt`, and `blast_dbs/<species>` must all match

## Species Genome Fasta Download Links
Human - https://www.ncbi.nlm.nih.gov/genome/?term=txid9606  
Mouse - https://www.ncbi.nlm.nih.gov/genome/?term=txid10090  
Zebrafish - https://www.ncbi.nlm.nih.gov/genome/?term=txid7955  
Rat - https://www.ncbi.nlm.nih.gov/genome?term=txid10116  

## Usage
- Add SNP's or INDEL's of interest to input csv file, see `sample_input.csv` for format
- From the project directory run `./snp_crispr.sh <species> <input_file> <PAM> <all>`
- Both `-NGG` and `-NAG` PAM sequences supported
- Running with `-all` argument designs guides where all SNPs are targeted within the 23-mer
- Design results found in `results/designs.csv`
- To test the pipeline run `./snp_crispr.sh dm sample_input.csv -NGG` and compare the output in `results/designs.csv` with `results/expected_results.csv`

**Note:** The pipeline can be run using any species or genome assembly by adding a new chromosome name to fasta id mapping file in the same format as the examples in `fasta_files/<species>_chr_ids.txt` and then following the same installation instructions

## Example Commands
- Fly: `./snp_crispr.sh dm dm_snps.csv -NGG -all`
- Human: `./snp_crispr.sh hs hs_snps.csv -NAG -all`
- Mouse: `./snp_crispr.sh mm mm_snps.csv -NGG`

## Design Scores
**Housden Efficiency Score**  
- These scores were computed using a position matrix. Detailed information about the input dataset and the algorithm can be found in Housden et al. Sci Signal. 2015 https://www.ncbi.nlm.nih.gov/pubmed/26350902
- Scores range from 1.47-12.32 (higher is better, > 5 recommended)

**Off Target Score**  
- Calculated based on sgRNA sequence blast results.
- Scores range from 0-5441.73 (lower is better, < 1 recommended)

## Precomputed Human Clinically Associated SNP Targeted Designs
- Precomputed CRISPR design results provided in [clinical_hs_snp_designs.csv](results/clinical_hs_designs.csv)
- Clinically associated variant data from [Ensembl](https://ensembl.org/Homo_sapiens)

## Authors
**Jonathan Rodiger** - [jrodiger](https://github.com/jrodiger)  
**Verena Chung** - [vpchung](https://github.com/vpchung)

## License
This project is licensed under the MIT license - see [LICENSE](LICENSE) for details

## Acknowledgements
**Perrimon Lab** (https://perrimon.med.harvard.edu/)    
**DRSC** (https://fgr.hms.harvard.edu/)
