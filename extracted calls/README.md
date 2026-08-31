# Partial audio collection

This Git history contains 728 WAV files, whereas the canonical processed table
contains 3,612 unique calls. The directory is therefore an incomplete,
illustrative subset and is not sufficient for rebuilding exact DTW or VAE
artifacts.

Cached reproduction does not read audio. Before an explicitly requested full
rebuild, `marmoset-repro preflight-audio` compares WAV basenames with all 3,612
canonical `call_id` values and stops before computation if any file is missing.
The complete audio archive should be distributed with the final data DOI rather
than represented as complete in this repository.
