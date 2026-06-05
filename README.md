[![R](https://img.shields.io/badge/-script-276DC3.svg?style=flat&logo=R)](https://cran.r-project.org)
[![Docker container build](https://github.com/KarlssonLaboratory/multiplexing_d.magna/actions/workflows/build-deploy-docker-container.yml/badge.svg)](https://github.com/KarlssonLaboratory/multiplexing_d.magna/actions/workflows/build-deploy-docker-container.yml)

# README

This project uses a small planktonic crustaceans (_Daphnia_), 0.2-6.0 mm. A _Daphnia_ is transparent are perfect for staining with various dyes able to penetrate into the tissues. Dyes can stain proteins, membranes, molecules etc and are detected with flourescece images. By treating the _Daphnia_ with different toxic compounds we can compare the differences of flourescence between the groups and understand which tissues, molecules are affected by the treatment.

## Raw data

During the submission process the dataset can be obtained by contacting Cedric Abele <Cedric.Abele@aces.su.se>. Upon publication the dataset will be included in the repo.

## Reproduce the results

> [!IMPORTANT]
> Make sure the dataset is located inside `data/`

```sh
# apptainer
apptainer exec \
  docker://ghcr.io/karlssonlaboratory/multiplexing_d.magna:a580599 \
  Rscript bin/LLM.R

# docker
docker run --rm \
  -v $(pwd):/data \
  ghcr.io/karlssonlaboratory/multiplexing_d.magna:a580599 \
  Rscript bin/LLM.R
```

## Reproducibility

A docker container holding the environment is produced within the repo and can be accessed [here](https://github.com/karlssonlaboratory/multiplexing_d.magna/pkgs/container/multiplexing_d.magna).

The environment setup can be inspected here: [Dockerfile](./Dockerfile).

## Experimental design

The experiment is setup in plates of 8 x 12 (96 wells), each well holds one _Daphnia_, with each row containing individuals with the same treatment, replicates. Each well translate into a fluorescence intensity values. These values are log2-transformed to correct for right-skewed data distribution and make the distribution normalised. In total, three replication of the plates were ran, giving us a batch effect.

We fitted a Linear Mixed-effect Model (LMM), as the model handles non-independent (dependent) observations, where each well are dependent on the plate which has confounding effects like imaging session (laser power, exposure etc), staining batch, environmental noise (person pipetted etc).

> So if plate 2 happened to have slightly brighter staining overall, every well on plate 2 inherits that brightness. The wells aren't giving us independent information about the underlying biology. This is called **clustering** or **pseudo-replication**.

Groups was considered as **fixed effects**, while experimental replicates were considered as **random effect** to account for inter-plate variability. Model validity and adherence to statistical assumptions (e.g., residual uniformity and homoscedasticity) were verified using simulated residual diagnostics (DHARMa R-package). Overall group effects were evaluated using Analysis of Deviance (Type II Wald chi-square tests) followed by Dunnett's post-hoc test to compare each treatment group against the control.

<!--

```
# data table sample:
                          Experiment_ID Well Row group   conc time survival
1 20260120-DM-MPLX-CCCPIII-24h_measures  E20   E   100 100.00   24      0.8
2    20260109-DM-MPLX-CCCP-24h_measures  C03   C  12.5   1.00   24      1.0
3 20260120-DM-MPLX-CCCPIII-24h_measures  F16   F   200 200.00   24      0.8
4    20260109-DM-MPLX-CCCP-24h_measures  D07   D    50   2.00   24      1.0
5    20260109-DM-MPLX-CCCP-24h_measures  C02   C  12.5   1.00   24      1.0
6    20260109-DM-MPLX-CCCP-24h_measures  E03   E   100   0.25   24      1.0
  channel intensity obs_id
1    DAPI  1665.287    352
2    DAPI  2466.174     28
3    DAPI  1898.649    367
4    DAPI  2600.560     70
5    DAPI  3040.146     25
6    DAPI  1986.162     82
```

```r
# formula in R
# fixed-effect  : response ~ predictor + other_terms
# random-effect : ( | )

fit_log_mixed <- lme4::lmer(
  formula = log(intensity) ~ group + (1 | Experiment_ID), 
  data = dye_data,

  # restricted maximum likelihood or maximum likelihood
  REML = TRUE
)
```

-->