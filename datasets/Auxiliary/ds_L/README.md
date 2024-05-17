# About very large dataset (L)

scMultiSim's support for very large datasets is experimental.
**You should switch to the `large-mem` branch before running the code in this folder.**

The memory optimization stores most data on disk in HDF5 format.
Data is loaded into memory by chunks when needed.
This allows scMultiSim to handle datasets that are too large to fit into memory.
Therefore, the empriical runtime performance depends on the I/O speed of your disk and the number of threads used.

We simulated 1 million cells and 1000 genes on a MacBook Pro (M2 Max, 64GB RAM) in about three hours.
