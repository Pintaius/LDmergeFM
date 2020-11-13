
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LDmergeFM

This script contains an implementation in R of a weighted average
procedure to generate consensus locus-specific LD matrices from multiple
single-cohort [LDSTORE](http://www.christianbenner.com/) files. Together
with [FINEMAP](https://doi.org/10.1093/bioinformatics/btw018), it has
been used in the generation of the
[PGC3-SCZ](https://doi.org/10.1101/2020.09.12.20192922) fine-mapping
results.

## Input files

`LDmergeFM` requires the
[readr](https://cran.r-project.org/web/packages/readr/index.html),
[dplyr](https://cran.r-project.org/web/packages/dplyr/index.html),
[purrr](https://cran.r-project.org/web/packages/purrr/index.html),
[reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)
and [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html)
packages to be available for your local R installation. It is designed
to be run from the command-line as:

    Rscript --vanilla LDmergeFM.R $LOCUS

Where `$LOCUS` is the locus identifier for the LD matrix being
calculated, present in all of the input files:

| Filename                | Contents                                                                                                                                                    |
| ----------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `$LOCUS.ref`            | Two column whitespace-delimited file. **Column 1:** SNP name. **Column 2:** Effect allele. Equivalent to columns 1 and 4 of a `FINEMAP` .z file. No header. |
| `$COHORT_$LOCUS.fam`    | 1x `PLINK v1.07+` .fam file for each cohort being analysed (with case/control phenotype).                                                                   |
| `$COHORT_$LOCUS.cor.gz` | 1x `LDSTORE v1.1` .cor file (compressed output of *–table* flag) for each cohort being analysed.                                                            |

Single-cohort names (`$COHORT`) should be unique but can contain any
non-whitespace characters. Note the underscore ("\_") separation with
`$LOCUS`.

## Output files

| Filename             | Contents                                                                                                       |
| -------------------- | -------------------------------------------------------------------------------------------------------------- |
| `$LOCUS.ld`          | Square consensus LD matrix. SNPs are given on the same order as `$LOCUS.ref`.                                  |
| `$LOCUS.snps.log`    | SNPs used in the computation of the consensus LD matrix. Should match those on `$LOCUS.ref`.                   |
| `$LOCUS.samples.log` | Cohorts used in the computation of the consensus LD matrix. Should match all of those provided as input files. |

<!--## Testing

The `./test/` folder contains some simulated input/output files that can be used to conduct a reproducible run.-->

## Notes

`LDmergeFM` has not been tested with correlation table file from
`LDSTORE v2+`, please conduct a test run before using these in important
analyses.

`LDmergeFM` can work with an arbitrary number of input matrices but in
its current state is not optimised to take advantage of multicore
environments or matrix sparsity, and thus can be potentially
resource-hungry. If working on systems with resource quotas, please
check `.log` files to make sure the computation of the consensus LD
matrix has used all available data.

`LDmergeFM` is not ancestry-aware. If ancestry-specific consensus
matrices are needed (e.g. for trans-ancestry fine-mapping purposes) you
should run the script separately for each group of single-ancestry
inputs.

## Citation

If this script is helpful for your work, just reference the main
[PGC3-SCZ](https://doi.org/10.1101/2020.09.12.20192922) paper. If it
ends up being *very* helpful, please let me know so I can keep fighting
impostor syndrome one day at a time :relieved:.

## Contact

Please submit suggestions and bug-reports at
<https://github.com/pintaius/LDmergeFM/issues>.
