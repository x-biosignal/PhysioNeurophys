# PhysioNeurophys

TMS and peripheral motor-neurophysiology analysis for the Physio
ecosystem.

The package provides TMS pulse metadata, MEP amplitude and latency,
Boltzmann recruitment curves, H-reflex recruitment and Hmax/Mmax
summaries, F-wave persistence and latency statistics, paired-pulse
SICI/ICF, cortical silent-period timing, sample-exact TMS pulse removal,
source-informed SSP-SIR and SOUND cleaning, canonical
TMS-evoked-potential summaries, and PREP2 prognosis.

TMS-EEG cleaning keeps the original assay, records every altered sample
and fixed cleaning matrix, and requires exact lead-field/channel
identity. TEP component labels are measurement windows, not proof of
cortical source identity; auditory and somatosensory confounds can
remain.

## Installation

``` r

install.packages("PhysioNeurophys",
  repos = c("https://x-biosignal.r-universe.dev", "https://cloud.r-project.org"))
```

## Part of the x-biosignal ecosystem

See the [x-biosignal](https://github.com/x-biosignal) organization and
[x-biosignal.r-universe.dev](https://x-biosignal.r-universe.dev) for the
full package suite.
