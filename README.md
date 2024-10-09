# Simulation study: random sample generation from PKBD

We consider a simulation study to compare the performance of the `QuadratiK` package with the `Directional` package for generating random samples from the Poisson kernel-based distribution (PKBD) in `R`.

## Content description

-   `sim_shell.sh`, `sim_slurm.sh`: files for running simulation on cluster
-   `sim.R`: main functions
-   `parameters.R`: parameter setting
-   `merge.R`: merge simulation results and save final data set
-   `figures_and_tables.R`: create plots in the figures folder
-   `\figures`: folder of figures
-   `\results`: folder of results
-   `QuadratiK_1.1.2.tar.gz`, `Directional_6.8.tar.gz`, `Rnanoflann_0.0.3.tar.gz`: needed packages
-   `res_rpkbd.RData`: saved simulation results

## Parameter setting

We generated

-   $n = 500, 1000, 2000, 5000, 10000$ random observations ,
-   dimension $d = 2, 3, 5, 10$,
-   concentration parameter $\rho=0.6, 0.7, 0.8, 0.9, 0.95$,
-   as mean direction a vector with 1 as first element and 0 in the remaining, that is $\mathbf{\mu} = (1, 0, \ldots, 0)$,
-   for each configuration, we considered $N=100$ replications.

Samples are generated using the function `rpkbd` of the `Directional` package and with the function `rpkb` in `QuadratiK` using the three available methods `rejvmf`, `rejacg` and `rejpsaw`. To ensure fairness, we used identical seed values for both methods to facilitate a direct comparison.

Notice that `rpkbd` and `rejacg` implement the same algorithm.

## Results

In terms of sample quality, `QuadratiK` and `Directional` produced nearly indistinguishable distributions.\
The method `rejvmf` demonstrates the highest execution time and it is displayed only in the first figures (time_rho_d2.pdf, time_rho_d3.pdf, time_rho_d5.pdf, time_rho_d10.pdf). It is then excluded for a better visualization.

The `rpkbd` function requires less computational time when the dimension is low. However, the difference is in order of seconds, or less than a second. For higher dimension, when the value of $\rho$ increases the computational time needed by `rpkbd` increases as well, while the methods `rejacg` and `rejpsaw` are consistent with respect to $\rho$.

## Additional investigations

The `Directional` package provides a function for estimating, via maximum likelihood, the concentration parameter $\rho$ and the mean direction $\mu$ of the PKBD model. Note that, the clustering algorithm in `QuadratiK` can be used for the same purpose choosing the number of clusters as 1.

Hence, we also investigate the performance for estimating the parameters of PKBD using the function `pkbd.mle()` from the `Directional` package and the function `pkbc` from `QuadratiK`.

We considered two cases

-   1: data are generated using `rpkbd`
-   2: data are generated using `rpkb` with method `rejacg`

The computational time does not change if the data are generated with the `rpkbd` function or by `rpkb`. In both cases, the function `pkbc` takes more time in estimating the parameters. This is more evident for increasing sample size. We create an additional plot of Time vs n. The clustering algorithm used for estimating the parameters of a single PKBD is slower than the MLE implemented in the `Directional` package.

We now investigate the performance of the clustering algorithm and the MLE in estimating the model parameters. For evaluating the performance we compute

-   the difference between the estimated rho and the "true" value used for generating the data set.
-   the MSE (Mean squared error) between the estimated mu and the value used for generating the data set.

#### Results - estimation

Results in estimating $\rho$ depend if data points are generated using the `Directional` package or the `QuadratiK` package. However the performance in estimation is the same.\
In particular, if data are generated using `rpkbd` the difference between the estimated $\rho$ and the "true" value has greater variability, while if the function `rpkb` is used, the difference has less variation and it is closer to zero.

For estimating the location vector mu, the considered functions perform almost identically, independently of the method used for generate the data.\
For curiosity we compute also the angle between the estimated location and the original mean vector as additional performance measure for the estimation of mu. We obtain identical results.
