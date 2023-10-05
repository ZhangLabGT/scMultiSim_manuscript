# scMultiSim Manuscript

This repository contains code used in the scMultiSim manuscript, including datasets generation and benchmarks.

## Datasets

All **datasets** used in the manuscript are available under [this Dropbox share](https://www.dropbox.com/sh/sfkn5hweaejbrir/AAB9liDyL8QuXd7LAgUsihGfa?dl=0).
Please note that certain imformation that were not used during benchmarking,
such as the GIV and CIF values, were omitted from the datasets for the sake of space.
You may use the code under `datasets/` to regenerate them.

Please note that generating the datasets from scratch requires long computation time and
a large amount of memory.

- Each main dataset may take 30-120 minutes depending on your hardware, number of cells, and number of genes.
- Each main dataset (or other CCI dataset) requires up to 16 GB of memory. Therefore, you may need to adjust the number of parallel processes in `run_simulation.py` to avoid running out of memory.
- For other datasets without the CCI, you may expect several minutes of simulation time for each dataset.

## Installing scMultiSim

To use scMultiSim, please refer to the [scMultiSim](https://github.com/ZhangLabGT/scMultiSim) repository,
where you can find vignettes and documentation.

Use the following command to install the development version of scMultiSim:

```R
devtools::install_github("ZhangLabGT/scMultiSim@main")
```

For your reference, the datasets used in the manuscript were generated at commit `d25f1c8`.

## Installing other dependencies

To run the benchmarks, you will need to have R (we used 4.1.0) and python (we used 3.9.13).
A Linux environment is required (we used Manjaro Linux).

The following dependencies are required for the benchmarks:

- CCI
  - [Giotto R 1.1.2](https://giottosuite.readthedocs.io/en/master/gettingstarted.html)
  - [SpaTalk 1.0](https://github.com/ZJUFanLab/SpaTalk)
  - [SpaOTsc](https://github.com/zcang/SpaOTsc) (_Python_)
- GRN
  - [BEELINE commit 79775f0](https://murali-group.github.io/Beeline/BEELINE.html#getting-started)
    - You also need Docker to run BEELINE.
- Integration
  - [Seurat v4 (`feat/dictionary` branch)](https://satijalab.org/seurat/)
    - The `feat/dictionary` branch has been included in Seurat v5. The original branch may be removed in the future.
  - [rliger 1.1.0](https://cran.r-project.org/web/packages/rliger/index.html) can be installed from CRAN.
- Clustetring
  - [CIDR](https://github.com/VCCRI/CIDR)
  - [Monocle3](https://cole-trapnell-lab.github.io/monocle3/docs/installation/)
  - [scater](https://bioconductor.org/packages/release/bioc/html/scater.html) can be installed from BioConductor
  - [SC3](https://bioconductor.org/packages/release/bioc/html/SC3.html) can be installed from BioConductor
  - [TSCAN](https://www.bioconductor.org/packages/release/bioc/html/TSCAN.html) can be installed from BioConductor
- Trajectory Inference
  - [dyno ](https://dynverse.org/users/1-installation/)
  - [Slingshot](https://www.bioconductor.org/packages/release/bioc/html/slingshot.html) can be installed from BioConductor
- RNA Velocity
  - [scVelo 0.2.4](https://scvelo.readthedocs.io/en/stable/installation/) (_Python_)
  - [velocyto 0.17.17](https://velocyto.org/velocyto.py/install/index.html) (_Python_)
- Other
  - `uwot` and `mclust` can be installed from CRAN.
  - For the visualization code, you need to have `ggplot` in R and `seaborn` in Python.
