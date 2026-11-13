# HPC tools

## Overview

It is a suite of codes for the evaluation of simulation performance in EC-Earth4.

## Documentation

Files needed:

- [ ]  `NODE.001_01`
- [ ]  `leginfo.yml`
- [ ]  `timing.log`
- [ ]  `paths.yml` (optional)

## Examples

Loaded libraries:

``` 
import numpy as np
import matplotlib.pyplot as plt
import lobster as lobs
```

Available commands:

- [x] single plot: `lobs.plot_sypd_vs_time(expname, leg)`
- [x] multi plots:  `lobs.multiplot_vs_time(expname, leg)`
- [x] multi experiments: `lobs.plot_sypd_vs_nptot([exp1,exp2,...], leg)`
- [x] save table: `lobs.save_table(expnames, leg)`

For more info, see jupiter notebook `test.ipynb`

![image](scaling.png)