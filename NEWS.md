# PhysioNeurophys 0.4.1

## Documentation

* A vignette carries one task end to end on synthetic or bundled data, offline,
  and is built and run by `R CMD check`.
* Runnable `@examples` added or corrected across 1 help pages. Each runs
  offline in seconds, writes nothing outside `tempdir()`, and is executed by
  `R CMD check`; anything needing a device, a download or an optional backend is
  fenced with the reason stated.

# PhysioNeurophys 0.4.0

- Added sample-exact TMS pulse and recharge removal with inspectable
  interpolation or blanking masks and diagnostics.
- Added channel-matched SSP-SIR and SOUND source-informed TMS-EEG cleaning.
- Added trial-wise baseline-corrected TEP averages and canonical signed peak
  summaries.

# PhysioNeurophys 0.3.0

- Added paired-pulse MEP ratios and conditioned-trial-weighted SICI/ICF
  summaries.
- Added sample-exact cortical silent-period onset, recovery, and duration
  analysis for tonic EMG.

# PhysioNeurophys 0.2.0

- Added H-reflex recruitment, Hmax/Mmax, latency, and threshold analysis.
- Added F-wave detection with persistence, latency distribution,
  chronodispersion, and mean F/M ratio.

# PhysioNeurophys 0.1.0

- Initial TMS motor-evoked-potential analysis release.
- Added MEP amplitude and latency analysis, Boltzmann recruitment curves,
  PREP2 prognosis, TMS metadata accessors, and `make_mep()`.
