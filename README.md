# Breast cancer treatment geography in Korea, 2011–2024

> **Code only (updated 2026-11-05).** This repository contains the analysis code. The processed district-year data, result tables, adjacency graph, and boundary files were removed; they are available from the author on reasonable request. All source statistics are public.


Version **2.0.0** — corrected analysis dated **2026-11-04**. Code, aggregate data and selected results for `Seongdae-Kim-Hongik/breast-cancer-spatiotemporal-korea-v2`.

## Contents

- `data/`: analysis-ready aggregate district-year panel, expected counts, canonical graph, node order, configuration and variable metadata.
- `code/primary.R`: minimal primary-model entrypoint requiring INLA and jsonlite, without sf or a shapefile.
- `code/archived_execution/`: privacy-redacted copies of workflow, preparation and review scripts. These are provenance code, not a fully portable pipeline.
- `results/`: selected analysis tables and documented numerical comparisons. Adopted KH IV and PPC outputs are in `results/adopted/`; the approved S41 combination is under `results/ppc_s41_audit/S41_combine_snapshot/`.
- `provenance/`: explicit inference selection, final S41 validation, source-ID-to-repository-path mapping and file hashes. Use the mapping to resolve historical source-relative IDs in selection metadata; those IDs are not runnable local paths.

## Minimal primary workflow

From the repository root:

```sh
Rscript code/primary.R .
Rscript code/primary.R . --fit /path/to/new_empty_output
```

The first command validates inputs without fitting. The optional second command fits the primary model into a new directory and records diagnostics; it does not silently switch approximation methods. Dependencies are R, INLA and jsonlite. The original primary used R 4.6.0, INLA 25.10.19 and one INLA thread (`1:1`). Stored per-run session information provides further package provenance.

The supplied 3,206 rows, 24 covariates, expected counts (E), and node and time indices were checked against the stored clean primary model. The design matrix has rank 25; the maximum differences in the design matrix and E were both zero. The distributed data and entrypoint retain the validated file hashes. The primary entrypoint also passed input validation without fitting. These checks do not constitute an independent refit or a complete rerun of the pipeline from this repository. Platform, optimizer and package differences can affect estimates; bitwise reproducibility is not promised.

Archived scripts contain `${BRCA_WORK}` and `${HOME}` placeholders after privacy redaction. They may require omitted population/registry sources, geometry, multi-disease inputs and local helper configuration. Adapt copies in a separate directory and document changes. Python dependencies can include numpy, pandas, scipy, geopandas, matplotlib and pyarrow; R helpers can require sf, spdep, Matrix, digest, ggplot2, dplyr and patchwork. Do not run these scripts blindly or overwrite archived results.

## Corrected analysis and inference selection

Version 2 corrects Cheongju/Sejong female denominators and the officially verified Changwon 2023 population/age discrepancy, propagating changes into expected counts and population-derived variables. The canonical graph has 229 nodes and 593 undirected edges, with two connected components of 227 and 2 nodes. Missing polygon-contact pairs were restored while retaining documented island links. Historical territorial alignment and source CRS limitations remain; fully verified 2024 boundaries are not claimed.

The primary reference is the original clean final S3 fit, not the warned duplicate. KH IV uses validated adopted A; PPC uses clean retry 1. S41 uses the approved snapshot of 68 analyses and 69 status records, including 12 adopted A retries and conditional 2024 retry A. Initial warned or quarantined fits must not be substituted by filename matching. Consult `provenance/inference_routing.json` and `provenance/source_to_repository.json`.

PPC component summaries describe posterior-mean component variation with covariance explicitly retained, not variance explained or R-squared. Noise-perturbation pooling uses CrI-based log-normal standard errors and is not proper multiple imputation. Conditional 2024 checks retain outcome-calibrated E/full-panel preprocessing and are not prospective forecasts. No comparison-outcome NCO fits were run.

## Measurement and source coverage

The panel contains aggregated annual treated-person counts, not individual patient data or incident diagnoses. Counts may include the same person in multiple years and include men; female resident-population denominators are proxies, not measured person-time. Regional outcome coverage includes Health Insurance and Medical Aid. External national age counts (code 73/C50; eight age groups: 0–9, 10–19, 20–29, 30–39, 40–49, 50–59, 60–69, and 70 years or older) cover Health Insurance only, excluding Medical Aid and non-covered care. Age-pattern expected counts are calibrated each year to the combined-coverage regional total; this does not recover the unobserved Medical Aid age profile or quantify bias from coverage differences.

Licensed raw shapefiles, unrelated MASTER data, patient records, credentials, private review links and personal home/cloud-storage paths are not distributed. Original source/public-copy SHA256 mappings preserve provenance without disclosing original private absolute paths.

## Licensing and data rights

`LICENSE` preserves the existing repository's MIT license (Copyright 2026 Seongdae Kim) for project software and associated documentation. No new blanket grant is made for third-party data, source statistics or geometry. Data and derived tables retain applicable source rights and conditions; access or derivation does not imply that the author can relicense every source under MIT. See `DATA_RIGHTS.md` and the panel/source metadata. Raw source redistribution permissions are not asserted.

## Version history and DOI

Version 2.0.0 is the authoritative corrected analysis for the results distributed here. Earlier v1 commits and artifacts are superseded. Their history is retained separately in the original private repository as provenance and is not included in this clean v2 repository. They contain known incorrect or incomplete population-denominator and geographic definitions, as well as pre-correction results; they should not be used as the current scientific specification or combined with v2 outputs. Preserving that separate private history does not endorse the earlier definitions. Use the versioned data, explicit inference selection and manifests together.

The earlier restricted v1.0 record is Zenodo record 23527444 (DOI `10.5281/zenodo.23527444`), associated with concept record 23527443 (`10.5281/zenodo.23527443`). These identifiers refer to the earlier version and the concept, respectively; neither is presented as a new v2 version DOI. **The corrected v2.0.0 archive is available at [10.5281/zenodo.23938401](https://doi.org/10.5281/zenodo.23938401)** ([Zenodo record](https://zenodo.org/records/23938401)), published on 2026-11-04. The immutable `v2.0.0` tag identifies scientific commit `6174fd805ea1332f966aca1b96e189fc4ec26e85`. Subsequent DOI metadata updates on `main` do not alter that tag or the archived release files.
