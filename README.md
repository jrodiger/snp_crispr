# SNP CRISPR
SNP-targeted CRISPR design pipeline for design of sgRNA's in Fly, Human, Mouse, Zebrafish, and Rat genomes. Online version available at: [https://www.flyrnai.org/tools/snp_crispr](https://www.flyrnai.org/tools/snp_crispr)

## Prerequisites
- Anaconda or Miniconda for environment management
**Note:** To install without conda, see `environment.yml` for dependencies
- Species genome in fasta format

## Install
- From project directory run `conda env create -f environment.yml`
- Then `conda activate snp_crispr`
- Put species genome fasta file in `fasta_files/<species>.fasta`
- Run `makeblastdb -in fasta_files/<species>.fasta -dbtype nucl -parse_seqids -title <species> -out blast_dbs/<species>`
- Species chromosome names -> fasta id mapping files included in `fasta_files/<species>_chr_ids.txt`
**Important Note:** `<species>` name/abbreviation in `fasta_files/<species>.fasta`, `fasta_files/<species>_chr_ids.txt`, and `blast_dbs/<species>` must all match

## Usage
- Add SNP's of interest to input csv file, see `sample_input.csv` for format
- From project directory run `./snp_crispr.sh <species> <input_file> <PAM> <all>`
- Both `-NGG` and `-NAG` PAM sequences supported
- Running with `-all` argument designs guides where all SNPs are targeted within the 23-mer
- Design results found in `results.csv`

## Example Commands
- Fly: `./snp_crispr.sh dm dm_snps.csv -NGG -all`
- Human: `./snp_crispr.sh hs hs_snps.csv -NAG -all`
- Mouse: `./snp_crispr.sh mm mm_snps.csv -NGG`

## Authors
**Jonathan Rodiger** - [jrodiger](https://github.com/jrodiger)
**Verena Chung** - [vpchung](https://github.com/vpchung)

## License
This project is licensed under the MIT license - see [LICENSE.txt](LICENSE.txt) for details

## Acknowledgements
**Perrimon Lab** (https://perrimon.med.harvard.edu/)    
**DRSC** (https://fgr.hms.harvard.edu/)
