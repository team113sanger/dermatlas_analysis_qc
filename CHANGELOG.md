# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [0.5.2] 2025-02-20
### Added
- Liscence, mirrroring and initial release

## [0.5.1] 2025-02-14
### Changed
- Updated build to inherit from MAF 0.6.2


## [0.5.0] 2024-10-31
### Changed
- Changes to `somatic_variants_qc.sh` to accommodate changes to QC plotting that allow for germline variant plotting and plotting of genes without Hugo symbols

## [0.4.3] 2024-09-27
### Added
- Added a Docker build and continuous integration to the repo. Build is on top of MAF and inherits all of it's available perl/R dependencies. Current version is 0.5.3.

## [0.4.2] 2024-09-03
### Fixed
- Fixed calculation of VAF vs depth plot height in `somatic_variants_qc.sh` when number of samples is divible by 5

## [0.4.1] 2024-08-19
### Changed
- Add option for the user to provide a gene/transcript list when running `reformat_vcf2maf.pl`

## [0.4.0] Publishable Unit 3
### Added
- Add `mutations_per_MB.sh` and `source.me`

### Changed
- Large refactor to move contstants and modules to `source.me`
- Many scripts now have usage information and interfaces
- Other small changes to improve code quality

## [0.3.0] Publishable Unit 3
### Remove
- Remove `parse_metadata_biosamples.pl` and `somatic_variants_qc_release_v2.sh`

## [0.2.0] Publishable Unit 2
### Added
- Add `parse_metadata_biosamples.pl` and `somatic_variants_qc_release_v2.sh`

## [0.1.0] Publishable Unit 1
### Added
- Add initial Perl script `check_vcf_tn_pairs.pl`,
  `make_samplelists_from_manifest.pl` and
  `one_sample_per_patient_from_manifest.pl`.
- Add initial shell scripts `add_commonSNPs2vcf.sh`,
  `check_caveman_pindel_status.sh`, `select_vcf_pass_ontarget.sh` and
  `somatic_variants_qc.sh`.
