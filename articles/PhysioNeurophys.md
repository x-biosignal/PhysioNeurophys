# TMS motor-evoked potentials with PhysioNeurophys

PhysioNeurophys analyses transcranial magnetic stimulation (TMS) and
peripheral motor neurophysiology on `PhysioExperiment` objects:
motor-evoked-potential (MEP) amplitude and latency, Boltzmann
recruitment curves, H-reflex and F-wave metrics, paired-pulse SICI/ICF,
cortical silent period, TMS-EEG cleaning, and the PREP2 prognosis
algorithm.

This vignette uses the built-in simulator so it runs offline. TEP/MEP
labels are measurement windows, not proof of a cortical source.

``` r

library(PhysioNeurophys)
#> Loading required package: PhysioExperiment
```

## 1. Simulate epoched MEPs

[`make_mep()`](https://x-biosignal.github.io/PhysioNeurophys/reference/make_mep.md)
generates biphasic MEP waveforms with a known latency and a Boltzmann
amplitude-intensity relationship, and attaches the per-trial TMS pulse
metadata. Epochs start 20 ms before the pulse.

``` r

pe <- make_mep(n_trials = 12, seed = 1)
dim(SummarizedExperiment::assay(pe))   # samples x channels x trials
#> [1] 600   1  12
getTMSpulses(pe)[1:3, ]
#>   pulse intensity_pct_mso target coil rmt_pct_mso onset_sec
#> 1     1          30.00000   <NA> <NA>          NA        NA
#> 2     2          33.63636   <NA> <NA>          NA        NA
#> 3     3          37.27273   <NA> <NA>          NA        NA
```

## 2. MEP amplitude

Peak-to-peak amplitude is measured within a response window, per channel
and trial. Pass `epoch_start_ms = -20` to match the simulated epoch.

``` r

amp <- mepAmplitude(pe, epoch_start_ms = -20)
head(amp[, c("trial", "intensity_pct_mso", "p2p")])
#>   trial intensity_pct_mso       p2p
#> 1     1          30.00000 0.1306161
#> 2     2          33.63636 0.1487786
#> 3     3          37.27273 0.2217329
#> 4     4          40.90909 0.3351076
#> 5     5          44.54545 0.5334475
#> 6     6          48.18182 0.8296879
```

## 3. MEP onset latency

``` r

lat <- mepLatency(pe, epoch_start_ms = -20)
head(lat[, c("trial", "onset_sample", "latency_ms")])
#>   trial onset_sample latency_ms
#> 1     1           NA         NA
#> 2     2           NA         NA
#> 3     3           NA         NA
#> 4     4          206       21.0
#> 5     5          204       20.6
#> 6     6          203       20.4
```

## 4. Boltzmann recruitment curve

Amplitude rises sigmoidally with stimulator intensity.
[`recruitmentCurve()`](https://x-biosignal.github.io/PhysioNeurophys/reference/recruitmentCurve.md)
fits a three-parameter logistic and returns the fitted parameters and
fit quality.

``` r

rc <- recruitmentCurve(pe, epoch_start_ms = -20)
rc
#> MEP recruitment curve
#>   MEP_max: 2.0731
#>   S50: 50.097 %MSO
#>   k: 5.6192 %MSO
#>   Maximum slope: 0.092231
#>   R-squared: 0.99778
#>   Trials: 12; converged: TRUE
```

## Where to go next

[`?PhysioNeurophys`](https://x-biosignal.github.io/PhysioNeurophys/reference/PhysioNeurophys-package.md)
lists every entry point, including peripheral measures
([`hReflexRecruitment()`](https://x-biosignal.github.io/PhysioNeurophys/reference/hReflexRecruitment.md),
[`hMaxMMax()`](https://x-biosignal.github.io/PhysioNeurophys/reference/hMaxMMax.md),
[`fWaveDetect()`](https://x-biosignal.github.io/PhysioNeurophys/reference/fWaveDetect.md)),
paired-pulse inhibition/facilitation
([`pairedPulse()`](https://x-biosignal.github.io/PhysioNeurophys/reference/pairedPulse.md),
[`siciIcf()`](https://x-biosignal.github.io/PhysioNeurophys/reference/siciIcf.md)),
the cortical silent period
([`corticalSilentPeriod()`](https://x-biosignal.github.io/PhysioNeurophys/reference/corticalSilentPeriod.md)),
TMS-EEG pulse removal and cleaning
([`removeTMSpulse()`](https://x-biosignal.github.io/PhysioNeurophys/reference/removeTMSpulse.md),
[`tmsSSPSIR()`](https://x-biosignal.github.io/PhysioNeurophys/reference/tmsSSPSIR.md),
[`soundClean()`](https://x-biosignal.github.io/PhysioNeurophys/reference/soundClean.md),
[`tepAverage()`](https://x-biosignal.github.io/PhysioNeurophys/reference/tepAverage.md)),
and the PREP2 upper-limb prognosis
([`prep2()`](https://x-biosignal.github.io/PhysioNeurophys/reference/prep2.md)).
