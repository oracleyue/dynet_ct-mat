# Toolbox for Dynamic Network Reconstruction 

This is a MATLAB toolbox of *continuous-time methods* for dynamic network
reconstruction from data with low sampling frequencies.

(Note that this is not the original git repository for development. It is only used
for released versions on github. The development is maintained in a private gitLab
repository.)

The basic idea is to model causal interactions between measurable variables by a
dynamical network model, called *dynamical structure function* (DSF), and to infer
network structures (with feedback loops allowed) by identifying continuous-time DSF
models. Moreover, it takes into account particular properties of biological time
series:

-  limited lengths of time series, and
-  low sampling frequencies.

The current released version only deals with full-state measurements.

## File organisation

The folders in the project root are mainly two types:

- folders with name starting with capital letters are used for storing application/simualation data, outputs/results, which you may skip if you want to use this toolbox;
- folders with lower-case names (incl. starting with `@`) are part of toolbox.

The toolbox consists of the following essential directories:

- `functions`: all essential functions for CT DSF model simulation and network inference;
- `supports`: the third-party MATLAB toolbox/functions as dependency;

One has to first run `sSPID_Init.m` to setup the search paths for the toolbox.

And the data/build directories:

- `Data`: the default folder to store simulation data of network models for further inference;
- `Results`: the default folder to save inference results;
- `Backup` and `Test`: the folders for backup or temporary scripts, which shouldn't appear in released versions.

If you didn't see the folders `Data` and `Results`, you may create them by `mkdir Data Results`.

## Dependency

- `mftoolbox`: a MATLAB toolbox for computation of matrix functions, developed by
N. J. Higham and released on his [homepage](http://www.ma.man.ac.uk/~higham/mftoolbox/).

## Examples

- `sSysSim.m`: demo to simulate continuous-time DSFs models for time series data;
- `sSPIDNoInput.m`: demo to perform network inference for models driven by noises;
- `sSPIDPerfCurve.m`: demo to draw performance curves (P-R or ROC)

## References

1. [Systems Aliasing in Dynamic Network Reconstruction: Issues on Low Sampling
Frequencies ](https://arxiv.org/abs/1605.08590)
2. [Dynamic Network Reconstruction in Systems Biology: Methods and Algorithms](http://publications.uni.lu/handle/10993/35580)