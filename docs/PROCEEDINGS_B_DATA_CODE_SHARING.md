# Proceedings B data and code sharing notes

This file summarises the release steps needed to make this repository ready for a submission to *Proceedings of the Royal Society B*.

## What the journal expects

- The Methods section should contain enough detail to allow interpretation and replication.
- Submission requires a Data Availability statement.
- Deposited datasets should have stable persistent identifiers where possible.
- The exact version of the repository used for the manuscript should be archived, not left only as a moving GitHub branch.

The Royal Society currently collects end-section statements during submission in ScholarOne rather than requiring them to remain in the manuscript source file, so prepare the wording in advance when you submit.

## Recommended release workflow for this repository

1. Push the manuscript-ready repository to GitHub.
2. Create a tagged release for the exact manuscript version.
3. Archive that release in Zenodo or Dryad and obtain a DOI.
4. Use the DOI and repository URL in the submission system's Data Availability statement.
5. If requested by the journal, add the archived dataset DOI to the reference list.
6. Confirm that all coauthors approve the final sharing scope and chosen licences.

## Suggested Data Availability statement

All data and code required to reproduce the analyses reported in this study are available in the project GitHub repository at [GITHUB_URL] and in the archived release at [DOI]. The archived release includes extracted call WAV files, processed sequence and spectral tables, model-ready analysis tables, and manuscript summary outputs.

## Suggested Code Availability sentence

Custom analysis code for preprocessing, sequence-distance estimation, spectral-distance estimation, and Bayesian multilevel modelling is included in the archived repository release referenced above.

## Final checks before submission

- Replace `[GITHUB_URL]` and `[DOI]`.
- Confirm that the archived release is the exact version submitted to the journal.
- Confirm that the repository is public and downloadable without requesting access.
- Confirm that the released files are the ones actually needed to reproduce the current manuscript analyses.
- Choose explicit code and data licences before the final public release.

## Repo-specific note

Folder names such as `seqeunce data` are retained intentionally for backwards compatibility with the existing notebooks and scripts.
