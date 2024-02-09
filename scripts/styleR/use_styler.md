# Introduction

styleR will parse your R codes and fix formating issue at various levels, adapting correct
spacing, correct indentation, correct break lines and even variable naming.

The main advantage, compare to lintr is that styleR does the reformating job *in your
place*, which in some occasion may be a huge time saving.

The disavantage of this avantage is that styleR change your code *in your place*... so be
prepared.

I have tested it and was literally impressed by its ability to radically intervene on the
code. In my test case, styleR split a recursive if-else loop in tree independent
non-recursive loops. It was a beautiful simplification, and a real progress for the code
maintenance.

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

- Note that `--dry on` will **not** make any change in your code.
- For more info on the multiple parameters of styler (which can be integrated to Rstudio
and will by available in the plug-in menu), see the complete reference below.
# Reference

[https://styler.r-lib.org/index.html](https://styler.r-lib.org/index.html)
