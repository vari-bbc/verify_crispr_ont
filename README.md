# Verify CRISPR ONT

Analyze ONT long-read sequencing adaptive sampling data to verify whether 
CRISPR edits were successful.

## Usage

### Step 1: Configure the workflow

- Copy or symlink your ONT sequencing data (FASTQ file) into the `raw_data/` directory. 
  - For best results, these should be basecalled with the SUP model within dorado. 
  - Note that the name of the FASTQ file (up to the `.fastq.gz` suffix) must
  match the `project_name` specified in the `config.yaml` file - this will be 
  perpetuated throughout the workflow.

- Modify the `config.yaml` file to specify the paths to your reference genome, 
and any parameters for the analysis. 
  - Note that the reference sequence should already contain the expected CRISPR 
  edits for accurate verification - this can be done using 
  [reform](https://gencore.bio.nyu.edu/reform/).

### Step 2: Run the analysis

Before submitting the job, you can test your configuration by performing a dry-run via:

```bash
module load bbc2/snakemake/snakemake-8.25.2
snakemake -np
```

When everything looks good, you can submit the job to SLURM using the provided `run_snake.sh` script.

```bash
sbatch bin/run_snake.sh
```

## Output

The results will be saved in the `results` directory, including:

- `qc`: Quality control html reports for the sequencing data (length and quality of reads).
- `align`: BAM files containing the aligned reads to the reference genome.
- `mosdepth`: Coverage information for the target regions within an html report.
- `deepvariant`: VCF files containing the variant calls for the target regions.
- `sniffles`: VCF files containing structural variant calls *and* an IGV-style
html report summarizing the structural variants detected.



## License

This project is licensed under the GNU General Public License v3.0. 
See the LICENSE file for details.

Created by [Kaitlyn](https://github.com/kew24) 
within the [Bioinformatics and Biostatistics Core](https://github.com/vari-bbc) of
the [Van Andel Institute](https://www.vai.com/).
