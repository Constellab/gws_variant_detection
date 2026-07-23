# Dependencies

This document lists the dependencies of the `gws_variant_detection` brick: the Constellab bricks it requires, and the external bioinformatics tools/packages used by its tasks.

## Constellab bricks

Declared in [`settings.json`](./settings.json):

| Brick | Version |
|---|---|
| `gws_core` | >= 0.20.7 |

## Task environments

`settings.json` declares no pip or git packages: each task creates its own isolated conda/pip environment at runtime (via `ShellProxy`), defined by an env file in `base_env/env_files/`.

| Task | Env file | Key packages |
|---|---|---|
| NanoQualitycheck, NanoAdaptersElimination, VariantCallingStats, NanoFreebayesVariantsDetection, NanoMinipileupVariantsDetection (`quality_check/`, `trimming_porechop/`, `vcf_stats/`, `variants_freebayes/`, `variants_minipileup/`) | `base_env/env_files/Variants_env.yml` | `minimap2==2.26`, `samtools==1.19.2`, `bcftools==1.19`, `freebayes==1.3.6`, `porechop==0.2.4`, `nanoplot==1.42.0`, `fastp==0.23.4`, `whatshap==2.3`, `pandas==2.2.2` |
| NanoReadsMapping (`mapping_minimap2/`) | `base_env/env_files/Alignment_env.yml` | `python==3.10`, `minimap2==2.26`, `samtools==1.19.2`, `pysam==0.22.1`, `pandas==2.2.2`, `numpy==1.26.*`, `scikit-learn==1.4.2`, `plotly==5.24.1` |
| NanoReadsFiltering (`filtering_chopper/`) | `base_env/env_files/Filtering_env.yml` | `chopper==0.9.0`, `mappy==2.28`, `pysam==0.23.0`, `pandas==2.2.2`, `scikit-learn==1.4.2`, `plotly==5.24.1` |
| ReadsClassification (`reads_classification/`) | `base_env/env_files/Mappy_env.txt` (pip) | `mappy==2.28`, `pysam==0.22.0` |

`variants_minipileup/` additionally downloads the `minipileup` executable at runtime from a Gencovery-hosted archive (not a pip/conda dependency).
