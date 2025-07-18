# Tag-and-Probe Muon Efficiency Analysis
[![CMS Open Data](https://img.shields.io/badge/CMS-Open_Data-FF6D00?logo=cern)](https://opendata.cern.ch)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> Muon identification/isolation efficiency measurement using CMS NanoAOD

## Repository Structure
```bash
src/
├── MuonExtractor/      # Muon selection from NanoAOD
├── TagAndProbe/        # Efficiency calculations
data/                   # ROOT files (git-ignored)
plots/                  # Output plots (git-ignored)
```
# Compile and run
root -l -q src/TagAndProbe/tag_and_probe_analysis.C

## Citation 

If you use this work in your research, please cite:

```bibtex
@software{TagAndProbeAnalysis,
  author = {Nisha},
  title = {Tag-and-Probe Muon Efficiency Analysis},
  year = {2024},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/Nisha-Dalal/Tag-and-Probe-analysis}}
