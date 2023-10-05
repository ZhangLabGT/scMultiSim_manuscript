# Running BEELINE

After running the post-processing script in the `datasets/Main` directory,
all input files to BEELINE (`ExpressionData.csv`, `refNetwork.csv`, `PseudoTime.csv`) should be ready under the designated GRN directory (`GRN_DIR`).
We also can provide these files upon request.

BEELINE can be run using the following command:

```sh
python BLRunner.py --config <config_file>
```

Since the main datasets are very large, we ran BEELINE in multiple batches,
each only containing 12 datasets.

The `beeline_config.yaml` file contains the configuration used for one batch,
and other configuration files are the same except for the dataset names.

When running BEELINE, some methods may take too long to run.
We manually killed the container after 8 hours of running, and BEELINE would continue to run the next method.
