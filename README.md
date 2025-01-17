
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LDmergeFM

This script contains an implementation in R of a weighted average
procedure to generate consensus locus-specific LD matrices from multiple
single-cohort correlation files. Together with
[FINEMAP](https://doi.org/10.1093/bioinformatics/btw018), it has been
used in the generation of the
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

    Rscript --vanilla LDmergeFM.R $LOCUS $COR_FORMAT $ESS_FORMULA 

Where the argument `$LOCUS` is the locus identifier for the LD matrix
being calculated, present in all of the input files:

| Filename                | Contents                                                                                                                                                                                                                                                     |
|-------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `$LOCUS.ref`            | Two-column whitespace-delimited file. **Column 1:** SNP name. **Column 2:** Effect allele. Equivalent to columns 1 and 4 of a `FINEMAP` .z file. No header.                                                                                                  |
| `$COHORT_$LOCUS.fam`    | 1x `PLINK v1.07+` .fam file for each cohort being analysed (with case/control phenotype). Individuals with missing phenotypes not used to compute the pairwise correlations should be excluded from this file.                                               |
| `$COHORT_$LOCUS.cor.gz` | 1x `LDSTORE v1.1` .cor file (compressed output of *–table* flag) for each cohort being analysed. Output of the `PLINK v1.9+` *–r inter-chr gz* flag (`$COHORT_$LOCUS.ld.gz`) is also acceptable if the `$COR_FORMAT` argument is changed as described below. |

Single-cohort names (`$COHORT`) should be unique but can contain any
non-whitespace characters. The underscore ("\_") separation with
`$LOCUS` is mandatory.

## Changing cohort weights and correlation file format

The other two arguments of the script are optional:

`$COR_FORMAT` indicates whether the input correlations have been
computed with *“LDSTORE”* or *“PLINK”*, allowing the script to correctly
process these files. Defaults to *“LDSTORE”* if not explicit.

`$ESS_FORMULA` indicates how to compute the effective sample size used
as weight of each LD matrix. Options are *“METAL”* for the formula used
in [Willer et al. 2010](https://doi.org/10.1093/bioinformatics/btq340)
or *“NCP”* for the definition of [Matti
Pirinen](https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS3.html)
and [Vukcevic et al. 2011](https://doi.org/10.1002/gepi.20576). Defaults
to *“METAL”* if not explicit.

Note that if these last two arguments are used, they have to be used
**in the order above**. This implies that to change `$ESS_FORMULA` one
needs to be explicit and state the value of `$COR_FORMAT` as well (but
the converse is not true).

## Output files

| Filename             | Contents                                                                                                                                                                                                            |
|----------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `$LOCUS.ld`          | Square consensus LD matrix. SNPs are given on the same order as `$LOCUS.ref`.                                                                                                                                       |
| `$LOCUS.snps.log`    | SNPs used in the computation of the consensus LD matrix. Should match those on `$LOCUS.ref`.                                                                                                                        |
| `$LOCUS.samples.log` | Cohorts used in the computation of the consensus LD matrix. Should match all of those provided as input files.                                                                                                      |
| `$LOCUS.heatmap.png` | Basic illustration of the consensus LD structure at the locus. Intended for troubleshooting or to identify regions that could be problematic for fine-mapping. Only generated if R installation has PNG capability. |

## Testing

The [`./test/`](test/) folder contains some input/output files that can
be used to conduct a reproducible run. For illustration purposes, these
files include the region around exon 12 of the
[EDAR](https://www.genecards.org/cgi-bin/carddisp.pl?gene=EDAR) gene,
which contains some very strong linkage as previously discussed by
[Sabeti et al. 2007](https://dx.doi.org/10.1038%2Fnature06250).
Genotypes were derived from polymorphic SNPs from four subpopulations
(Europeans, Sub-Saharan Africans, East Asians and Native Americans) of
the public
[HGDP](ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/)
dataset. Please reference [Bergström et
al. 2020](https://dx.doi.org/10.1126/science.aay5012) if you find this
data useful for other purposes.

## Assumptions

`LDmergeFM` has been designed with fine-mapping in a meta-analytic
case-control GWAS setting in mind, so one of its implicit requirements
(in line with [FINEMAP](https://doi.org/10.1093/bioinformatics/btw018))
is that the reference allele for the correlations is the **same** in all
cohorts. Given that inconsistent criteria are currently used to decide
effect/reference alleles, it can help to set these explictly using the
`PLINK` *–a1-allele*/*–ref-allele* flag, which in fact can accept the
format of the `$LOCUS.ref` file.

For a similar reason, `LDmergeFM` expects that all SNPs in the
`$LOCUS.ref` file will be uniquely named and that each of them can be
found in the correlation file of **at least** one cohort. Duplicated or
missing SNPs might cause the script to fail silently, returning
erroneous output, so please ensure these are not present.

## Notes

`LDmergeFM` has not been tested with correlation table files from
`LDSTORE v2+`, please conduct a test run before using those in important
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

## Major version history

**2021-03-09** =&gt; *Added some internal checks for better error
reporting. Introduced arguments to accommodate other correlation file
formats and change the calculation for effective sample size weights if
desired. New basic heatmap output.*

**2020-11-13** =&gt; *Upload of initial version with essential
functionality.*

## Additional software

`FINEMAP`/`LDSTORE`: <http://www.christianbenner.com/>

`PLINK`: <https://www.cog-genomics.org/plink/1.9/>

## Citation

If this script is helpful for your work, just reference the main
[PGC3-SCZ](https://doi.org/10.1101/2020.09.12.20192922) paper. If it
ends up being *very* helpful, please let me know so I can keep fighting
impostor syndrome one day at a time :relieved:.

## Contact

Please submit suggestions and bug-reports at
<https://github.com/pintaius/LDmergeFM/issues>.
