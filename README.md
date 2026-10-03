# PhysioNeurophys

<!-- badges: start -->
[![r-universe](https://x-biosignal.r-universe.dev/badges/PhysioNeurophys)](https://x-biosignal.r-universe.dev/PhysioNeurophys)
<!-- badges: end -->

TMS and peripheral motor-neurophysiology analysis for the Physio ecosystem.

The package provides TMS pulse metadata, MEP amplitude and latency,
Boltzmann recruitment curves, H-reflex recruitment and Hmax/Mmax summaries,
F-wave persistence and latency statistics, paired-pulse SICI/ICF, cortical
silent-period timing, sample-exact TMS pulse removal, source-informed SSP-SIR
and SOUND cleaning, canonical TMS-evoked-potential summaries, and PREP2
prognosis.

TMS-EEG cleaning keeps the original assay, records every altered sample and
fixed cleaning matrix, and requires exact lead-field/channel identity. TEP
component labels are measurement windows, not proof of cortical source
identity; auditory and somatosensory confounds can remain.

## Installation

```r
# the containers build on Bioconductor, so its repositories are needed too
install.packages("BiocManager", repos = "https://cloud.r-project.org")
install.packages("PhysioNeurophys",
  repos = c("https://x-biosignal.r-universe.dev", BiocManager::repositories()))
```

## Part of the x-biosignal ecosystem

See the [x-biosignal](https://github.com/x-biosignal) organization and
[x-biosignal.r-universe.dev](https://x-biosignal.r-universe.dev) for the full package suite.
