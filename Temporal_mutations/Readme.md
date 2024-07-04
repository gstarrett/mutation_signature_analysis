# User Guide for temporal_mutations.py

The script determines the category of variants based on a combination of statistical analyses and predefined thresholds related to allele frequency, copy number, and purity. This categorization process is encapsulated within the `timing` function, which utilizes the `bootstrap` function to estimate confidence intervals around allele frequencies and then compares these estimates against specific criteria to assign a category to each variant. Here's a breakdown of how the script categorizes variants:

1. **Adjustment for Purity**: Before categorization, the script adjusts allele frequencies based on a given purity value. This adjustment accounts for the dilution effect of normal cells in tumor samples, aiming to provide a more accurate representation of the true allele frequencies in the tumor genome.

2. **Bootstrap Method for Confidence Interval Estimation**: The `bootstrap` function is used to estimate the 95% confidence interval for allele frequencies based on read coverage at that position. This method involves repeatedly sampling from the observed allele frequencies with replacement and calculating the mean and standard deviation of these samples. The resulting confidence interval provides a measure of the uncertainty associated with the estimated allele frequency.

3. **Categorization Criteria**: The `timing` function uses the adjusted allele frequency and the confidence interval estimated by the `bootstrap` function to categorize each variant. Variants are classified into one of four categories based on comparisons against predefined thresholds:
   - **EARLY_clonal**: If the lower bound of the confidence interval is greater than 0 and less than 0.9 divided by the copy number, and the upper bound is greater than 1.8 divided by the copy number.
   - **LATE_clonal**: If the lower bound of the confidence interval is greater than or equal to 0.9 divided by the copy number.
   - **LATE_subclonal**: If the lower bound of the confidence interval is greater than 0 and less than 0.45, regardless of the copy number.
   - **UNK**: If none of the above conditions are met.

## Prerequisites

- Ensure Python 3.x is installed on your system.
- Install required Python packages (`argparse`, `math`, `random`, `pybedtools`, `os`) via pip:

```bash
pip install pybedtools
```

Note: Some dependencies like `pybedtools` might require additional system-level tools to be installed.

## Running the Script

The script accepts command-line arguments for input files and optional parameters. Run the script from the terminal as follows:

```bash
python script_name.py --vcf path_to_vcf_file --seg path_to_seg_file [--prefix output_prefix] [--purity purity_value]
```

Replace `script_name.py` with the name of your Python script file. Replace `path_to_vcf_file` and `path_to_seg_file` with the actual paths to your VCF and SEG files. The `--prefix` option specifies the prefix for the output file names, and `--purity` allows you to specify the purity value for the analysis.

### Example Command

```bash
python script_name.py --vcf variants.vcf --seg segments.seg --prefix my_output --purity 0.98
```

This command will process `variants.vcf` and `segments.seg`, applying a purity value of 0.98, and save the output under the prefix `my_output`.

## Understanding the Output

The script generates a text file named `<output_prefix>.txt` in the current directory, where `<output_prefix>` is the value provided to the `--prefix` argument. This file contains the analyzed and categorized information about the genetic variants found at the intersections of VCF and SEG data.

Each line in the output file represents a variant and includes fields such as allele frequency, copy number, coverage, and a category indicating whether the variant is considered early clonal, late clonal, late subclonal, or unknown based on the analysis.
