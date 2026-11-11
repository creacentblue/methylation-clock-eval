# References and data attribution

Every item below was checked against the publisher/PMC page, GEO, or the metadata that `pyaging` ships with its clocks
(`results/clock_info.csv` stores the pyaging citation text and DOI per clock). Items I could not verify are not listed.

## Datasets (downloaded at run time; not redistributed in this repository)

| Accession | What | Source paper |
|---|---|---|
| GEO [GSE40279](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279) | 656 whole-blood 450K profiles, ages 19–101 | Hannum G, et al. Genome-wide methylation profiles reveal quantitative views of human aging rates. *Mol Cell* 49(2):359–367 (2013). doi:10.1016/j.molcel.2012.10.016 |
| GEO [GSE87571](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87571) | 732 GEO records (729 individuals), whole-blood 450K, ages 14–94, Northern Sweden Population Health Study (NSPHS) | Johansson Å, et al. Continuous aging of the human DNA methylome throughout the human lifespan. *PLoS ONE* 8(6):e67378 (2013). doi:10.1371/journal.pone.0067378 (421-sample batch). Ahsan M, et al. The relative contribution of DNA methylation and genetic variants on protein biomarkers for human diseases. *PLoS Genet* 13(9):e1007005 (2017). doi:10.1371/journal.pgen.1007005 (additional samples; genotype/cell-estimate supplement) |
| GEO [GSE35069](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35069) | Purified blood cell types (6 male donors), used as the cell-type reference | Reinius LE, et al. Differential DNA methylation in purified human blood cells: implications for cell lineage and studies on disease susceptibility. *PLoS ONE* 7(7):e41361 (2012). doi:10.1371/journal.pone.0041361 |

GEO data are public and de-identified; they and the Reinius reference keep their own terms of use. Please cite the source papers above if you use the data.

## Software and methods

- **pyaging** — de Lima Camillo LP. pyaging: a Python-based compendium of GPU-optimized aging clocks. *Bioinformatics* 40(4):btae200 (2024). doi:10.1093/bioinformatics/btae200 (version used here: 0.6.2).
- **Houseman constrained projection** — Houseman EA, et al. DNA methylation arrays as surrogate measures of cell mixture distribution. *BMC Bioinformatics* 13:86 (2012). doi:10.1186/1471-2105-13-86.
- **minfi** (the probe-selection recipe `pickCompProbes` that `scripts/02b_cell_reference.py` imitates; minfi itself is not used) — Aryee MJ, et al. Minfi: a flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays. *Bioinformatics* 30(10):1363–1369 (2014). doi:10.1093/bioinformatics/btu049.
- **Salas 2018 EPIC deconvolution** (pyaging's `DeconvoluteBloodEPIC*`, run only as a diagnostic) — Salas LA, et al. An optimized library for reference-based deconvolution of whole-blood biospecimens assayed using the Illumina HumanMethylationEPIC BeadArray. *Genome Biology* 19:64 (2018). doi:10.1186/s13059-018-1448-7 (citation as given in pyaging metadata).

## Clocks (citation text as given in pyaging metadata)

| pyaging name | Reference |
|---|---|
| horvath2013 | Horvath S. DNA methylation age of human tissues and cell types. *Genome Biology* 14, R115 (2013). doi:10.1186/gb-2013-14-10-r115 |
| hannum | Hannum G, et al. Genome-wide methylation profiles reveal quantitative views of human aging rates. *Molecular Cell* 49(2), 359–367 (2013). doi:10.1016/j.molcel.2012.10.016 |
| skinandblood | Horvath S, et al. Epigenetic clock for skin and blood cells applied to Hutchinson Gilford Progeria Syndrome and ex vivo studies. *Aging* 10(7), 1758–1775 (2018). doi:10.18632/aging.101508 |
| zhangen | Zhang Q, Vallerga CL, Walker RM, et al. Improved precision of epigenetic clock estimates across tissues and its implication for biological ageing. *Genome Medicine* 11, 54 (2019). doi:10.1186/s13073-019-0667-1 |
| lin | Lin Q, et al. DNA methylation levels at individual age-associated CpG sites can be indicative for life expectancy. *Aging* 8(2): 394–401 (2016). doi:10.18632/aging.100908 |
| intrinclock | Tomusiak A, et al. Development of an epigenetic clock resistant to changes in immune cell composition. *Communications Biology* 7, 934 (2024). doi:10.1038/s42003-024-06609-4 |
| dnamphenoage | Levine ME, et al. An epigenetic biomarker of aging for lifespan and healthspan. *Aging* 10(4): 573–591 (2018). doi:10.18632/aging.101414 |
| dunedinpace | Belsky DW, Caspi A, Corcoran DL, et al. DunedinPACE, a DNA methylation biomarker of the pace of aging. *eLife* 11, e73420 (2022). doi:10.7554/elife.73420 |
