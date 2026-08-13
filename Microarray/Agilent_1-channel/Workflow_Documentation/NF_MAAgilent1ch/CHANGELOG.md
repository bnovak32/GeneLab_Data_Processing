# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.5](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_MAAgilent1ch_1.0.5/Microarray/Agilent_1-channel/Workflow_Documentation/NF_MAAgilent1ch) - 2025-05-xx

### Added

- Add support for custom annotations, see [specification](examples/annotations/README.md)
- Add option to skip differential expression analysis (`--skipDE`) ([#104](https://github.com/nasa/GeneLab_Data_Processing/issues/104))
- Add a retry wrapper for functions that utilize internet resources, syncing with [NF_MAAffymetrix_1.0.3](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_MAAffymetrix_1.0.3/Microarray/Affymetrix/Workflow_Documentation/NF_MAAffymetrix)
- Workflow can now be run using an ISA archive by supplying parameter: 'isaArchivePath' (as either a local path or public web uri), syncing with [NF_MAAffymetrix_1.0.2](https://github.com/asaravia-butler/GeneLab_Data_Processing/tree/NF_MAAffymetrix_1.0.2/Microarray/Affymetrix/Workflow_Documentation/NF_MAAffymetrix)
- Add nextflow schema support for parameter validation and help text generation
- Add conda support for easier local development and debugging

### Changed

- Replace `RUNSHEET_FROM_GLDS` and `RUNSHEET_FROM_ISA` processes and their associated workflow logic with a new staging analysis subworkflow supporting both accession-based and input-file-based execution modes
- Rework publish directory behavior as part of the staging analysis subworkflow where `outdir` is now the base directory for `GLDS-NNN/` output directory if `--accession` is provided, or the base directory for `results/` output directory if `--runsheet` is provided
- Bump gl-microarray image from version 1.0.0 to 1.1.0 to match R package updates in the [GL-DPPD-7112-A pipeline document](../../Pipeline_GL-DPPD-7112_Versions/GL-DPPD-7112-A.md)
- Rename `annotation_config_path` as `array_annot_path` and `config.csv` as `design_info.csv` throughout the workflow and documentation to better reflect the purpose of the file and its contents
- Convert `generate_protocol.sh` to a Python script for automated handling of reference/annotation parameters
- Add `create_date` to design_info.csv and parse it in the protocol
- Move protocol creation from post-processing to main nextflow script to make passing needed values easier and more robust
- Rename module files from UPPERCASE.nf to lowercase.nf following Nextflow community convention
- Flatten directory-based modules (PROCESS_NAME/ with scripts under resources/usr/bin/) to single lowercase process_name.nf files directly under modules/
- Move process scripts from modules/PROCESS_NAME/resources/usr/bin/ to the top-level bin/ directory
- Fixes in `Agile1CMP.qmd`
  - Replace live biomaRt::getBM() queries in the QMD with direct downloads of Ensembl's FTP mart-dump tables; drop chunking/retry/Sys.sleep tied to those queries
  - Simplify group sample retrieval during differential expression group-wise statistics computation to use a more concise `filter/pull/sort` chain instead of `group_by/summarize/filter/pull`, addressing the deprecation warning in dplyr >= 1.1.0 where returning more than 1 row per `summarise()` group is deprecated
- Changes to post-processing workflow
  - Resolve output directory `GLDS-NNN/` or `results/` to match main workflow behavior
  - Replace dp_tools dependency in assay table update and md5sum table generation with standalone scripts
  - Rename `UPDATE_ISA_TABLES` and `update_curation_table.py` to `UPDATE_ASSAY_TABLE` and `update_assay_table.py` to better reflect their purpose
  - Add new PURGE_PROCESSING_INFO Nextflow module to strip full paths in nextflow_processing_info_GLmicroarray.txt before publishing
  - Add parameter validation and summary log from nf-schema

### Removed

- Packages `R.utils`, `purrr`, and `biomaRt` are no longer used in the processing code, and have been removed from software table generation

## [1.0.4](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_MAAgilent1ch_1.0.4/Microarray/Agilent_1-channel/Workflow_Documentation/NF_MAAgilent1ch) - 2024-10-02

### Added

- Add automatic generation of processed data protocol ([#85](https://github.com/nasa/GeneLab_Data_Processing/issues/85))

### Changed

- Small bug fixes in `Agile1CMP.qmd`
  - Check if `getBM()` returned results before concatenating it to dataframe to avoid error in `bind_rows()` ([#96](https://github.com/nasa/GeneLab_Data_Processing/issues/96))
  - When renaming column names, specify which columns to rename to avoid unintentional renaming ([#97](https://github.com/nasa/GeneLab_Data_Processing/issues/97))
  - When renaming factor names, prevent cases where a factor is partially renamed because it contains a substring that is another factor ([#100](https://github.com/nasa/GeneLab_Data_Processing/issues/100))
- Update software table generation to exclude `R.utils` from table if data files are not compressed ([#99](https://github.com/nasa/GeneLab_Data_Processing/issues/99))

## [1.0.3](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_MAAgilent1ch_1.0.3/Microarray/Agilent_1-channel/Workflow_Documentation/NF_MAAgilent1ch) - 2024-05-17

### Changed

- Fix cache location issues that arose in `quarto render` when using Nextflow v.23.10.1 ([#82](https://github.com/nasa/GeneLab_Data_Processing/issues/82))

## [1.0.2](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_MAAgilent1ch_1.0.2/Microarray/Agilent_1-channel/Workflow_Documentation/NF_MAAgilent1ch) - 2023-04-28

### Added

- Support for Arabidposis Thaliana datasets using the plants ensembl FTP server.

### Changed

- When encountering error about column reordering, the expected order is saved for debugging purposes.
- Post Processing Workflow: Assay Table Update now added '_array_' prefix to processed files instead of '_microarray_' prefix.

## [1.0.1](https://github.com/asaravia-butler/GeneLab_Data_Processing/tree/NF_MAAgilent1ch_1.0.1/Microarray/Agilent_1-channel/Workflow_Documentation/NF_MAAgilent1ch) - 2023-03-31

### Removed

- Deprecated column renaming code (abcd380)

### Fixed

- Bumped dp_tools from 1.3.0 to 1.3.1 to address 'ISO-8859-1' encoded ISA archive files (example: OSD-271-v2) (d518f40)
- Added handling for raw data that lacks the ProbeUID column (example: OSD-271-v2) (efbc237)

### Changed

- Reordering error message is now more informative (007e36c)

## [1.0.0](https://github.com/asaravia-butler/GeneLab_Data_Processing/tree/NF_MAAgilent1ch_1.0.0/Microarray/Agilent_1-channel/Workflow_Documentation/NF_MAAgilent1ch) - 2023-03-22

### Added

- First internal production ready release of the Agilent 1 Channel Microarray Processing Workflow