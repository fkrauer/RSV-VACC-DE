# Effectiveness and efficiency of immunisation strategies to prevent RSV among infants and older adults in Germany: a modelling study

Fabienne Krauer, Felix Guenther, Marina Treskova-Schwarzbach, Viktoria Schoenfeld, Mihaly Koltai,  Mark Jit, David Hodgson, Udo Schneider, Ole Wichmann, Thomas Harder, Frank G. Sandmann,  Stefan Flasche

This study is currently published as a preprint on [Medrxiv](https://www.medrxiv.org/content/10.1101/2024.06.20.24309248v1). You are free to use and adapt this model provided you cite the publication and this repository according to scientific citation standards. 

## Abstract
Background 

Recently, several novel RSV immunisation products that protect infants and older adults against RSV disease have been licensed in Europe. We estimated the effectiveness and efficiency of introducing these RSV immunisation strategies in Germany. 

Methods 

We used a Bayesian framework to fit an age-structured dynamic transmission model of RSV to sentinel surveillance and RSV-specific hospitalisation data in Germany from 2015-2019. The calibrated model was used to evaluate different RSV intervention strategies over five years: long-acting, single-dose monoclonal antibodies (mAbs) in high-risk infants aged 1-5 months; long-acting mAbs in all infants aged 1-5 months; seasonal vaccination of pregnant women and one-time seasonal vaccination of older adults (75+/65+/55+ years). We performed sensitivity analysis on vaccine uptake, seasonal vs. year-round maternal vaccination, and the effect of under-ascertainment for older adults. 

Results 

The model was able to match the various RSV datasets. Replacing the current short-acting mAB for high-risk infants with long-acting mAbs prevented 1.1% of RSV-specific hospitalisations in infants per year at the same uptake. Expanding the long-acting mAB programme to all infants prevented 39.3% of infant hospitalisations per year. Maternal vaccination required a larger number to be immunised to prevent one additional hospitalisation than a long-acting mAB for the same uptake. Vaccination of adults older than 75 years at an uptake of 40% in addition to Nirsevimab in all infants prevented an additional 4.5% of all RSV-hospitalisations over five years, with substantial uncertainty in the correction for under-ascertainment of the RSV burden. 

Conclusions 

Immunisation has the potential to reduce the RSV disease burden in Germany.


## Software
The main model code is written in Julia and fitted with [Turing](https://turing.ml/dev/). You can download Julia [here](https://julialang.org/downloads/). 

## Repository content and instructions

This repository contains all Julia and R code for this analysis as well as the input data and the final posterior samples (chains). The scripts, which you can run from the console are:

| script name | description |
| :--- | :--- |
| **scripts/fit.jl** | fits the model to the RSV data. Produces a date-time-stamped directory for each seed with the chain object, an arguments object, the initial thetas and the model object |
| **scripts/visualize.jl** | to quickly visualize the results of the selected fit(s). Produces a PDF file |
| **scripts/postprocessing.jl** | simulates the posterior predictive of the selected fit(s). Produces several csv files and figures |
| **scripts/vaccsim.jl** | simulates the theoretical disease dynamics under different vaccination scenarios, based on the posteriors from the selected fit(s). Produces a csv file |


The repo also contains a `src` directory with several functions, which are all sourced when the project is compiled:
| **src/models.jl** | Contains the turingmodel |
| **src/differential_equations/** | Contains the ODE models for the fitting and the vaccination simulation |


## Getting started

To get started with this project, clone the Github repository into your local directory. You then need to instantiate the project, which downloads all the necessary packages with the correct version numbers. This only has to be done before the first use. 

Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

## Reproducing the results

The structure and workflow of the code is as follows:

1. Initial data cleaning and reshaping. This is done with the R files 1-3. 
2. Fitting of the model to the different datasets. This is done with `fit.jl`. This scripts takes several arguments: *seed* is any seed no. you prefer. *initialize* samples from the prior for the initial draw. *verbose* returns more detailed information of what the script is doing at any time, we recommend to turn it on. The fitting script stores the model, the trace, the arguments and the state after adaptation in a date-time stamped folder within the `output` directory. To reproduce our results use the seeds 1-4 and the following commands: (replace "seed" with the individual seed you choose). If you are using VSCode, you can run several instances of Julia in separate terminals to run the scripts in parallel. Keep an eye on CPU and memoary usage though. Depending on your computer and the seed, a fit takes about 48 hours to complete. There is currently no option for multithreading. 

```
julia --project scripts/fit.jl seed --initialize --verbose
```

If you receive an error, you can print the error message to a log file using

```
julia --project scripts/fit.jl seed --initialize --verbose 2>&1 > log.log
```

<br/> 

3. Diagnostics/visualization. To inspect the trace plots and the marginal posterior distributions as well as the fit to the data sets, do
```sh
julia --project scripts/visualize.jl --verbose --prefix=prefix --ignore-commit
```
where prefix is the date or date-time prefix of the posterior directories you want to include in the diagnostics. To reproduce our results, use `-prefix=20-8-2023_8`. You can also type out the individual directory names serially instead of using the prefix:
```sh
julia --project scripts/visualize.jl --verbose output/20-8-2023_8-40-8-909 output/20-8-2023_8-40-3-073 --ignore-commit
```

<br/> 

4. Postprocessing. To generate the **posterior predictive** of the fit and resulting figures, do:
```sh
julia --project scripts/postprocessing.jl --verbose --prefix=prefix --ignore-commit
```
To reproduce our results, use the prefix as in point 3.

<br/> 

5. Vaccination strategies.

To simulate the **main vaccination strategies**, do:
```sh
julia --project scripts/vaccsim.jl --prefix=prefix --ignore-commit --verbose
```

To simulate the **sensitivity vaccination strategies**, do:
```sh
julia --project scripts/vaccsim.jl --prefix=prefix --ignore-commit --verbose --sensitivity --comparator ="0" --strategies="[1, 2, 3, 4, 5, 6, 7, 8]"
```
<br/> 


