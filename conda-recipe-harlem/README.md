# simulation-base-conda-src

This repository contains the source files necessary to build the
simulations-base conda package. Please note that this is a work in
progress and that some portions are incomplete.

Currently, the build only works on / produces packages for the linux64
target. Since arbalest is distributed in binary form, there is no point
in building for other targets. WSL may work, but is untested.

## Building

To build the conda package, you need a conda environment that has
conda-build available:

```
$ mkdir ~/conda
$ conda create conda-build -p ~/conda/condabuild
$ conda activate ~/conda/condabuild
```

Clone the git repo:

```
TODO SGREENSLADE: ~/conda/simulation-base
```

Make any needed modifications to the conda package source files. Be sure
to update meta.yaml with a new build number / package version if you are
intending to release the package to a conda channel. Reset the build
number to zero if you increment the package version.

Build the simulations-base package:

```
conda build ~/conda/simulation-base
```

Test the new package to verify that it has not broken anything. If it
tests OK, push it to the distribution channel.

```
TODO SGREENSLADE: Proper push instructions once we have a real web server.
```
