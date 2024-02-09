# Create the proper conda environment

If it does not exist yet, create the conda environnement:

```
conda create --strict-channel-priority --solver libmamba --override-channels --channel conda-forge --channel bioconda --channel defaults --name r-styler r-styler r-argparse
```

# get the styler.R script in your environment and ensure it is executable

to make style.R executable, `chmod 755 styler.R`

# `conda activate r-styler`

# Run styler.R

:warning: styler.R will modify your target file using the following command; work in a
versioned (git) environment !

```
<path>/./styler.R --dry off <path>/<Rfile-to-lint>
```

# Reference

[https://styler.r-lib.org/index.html](https://styler.r-lib.org/index.html)
