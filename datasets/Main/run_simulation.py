import os
import subprocess
from multiprocessing import Pool

PROC_COUNT = 7
TREES = [1, 3, 5]
CONFIGS = [x + 1 for x in range(48)]
TASKS = [(tree, conf) for tree in TREES for conf in CONFIGS]

def f(x):
    print(x)
    subprocess.run(["Rscript", "--vanilla", "sim.R", f"{x[0]}", f"{x[1]}"])
    return x

if __name__ == '__main__':
    with Pool(PROC_COUNT) as p:
        print(p.map(f, TASKS))
