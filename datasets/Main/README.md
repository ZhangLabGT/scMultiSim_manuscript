# The Main Datasets

Run the following command to generate **one** main dataset:

```bash
Rscript --vanilla datasets/main.R <tree> <n>
```

Where

- `<tree>` is the population type. `1` means the linear tree (`ML`), `3` means the continuous tree (`MT`), and `5` means discrete (`MD`).
- `<n>` is the configuration number (1-48).

For example,

```bash
Rscript --vanilla datasets/main.R 1 1
```

generates dataset `ML1a`.

## Generating all main datasets

The `run_simulation.py` provides an example of running the tasks in parallel.

However, it is not recommended to generate all main datasets at once, as it may take a long time and a large amount of memory (see the main README of this repository).
