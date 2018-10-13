# Contributing

This document is mostly inspired from the [IUC contributing policy](https://github.com/galaxyproject/tools-iuc/blob/master/CONTRIBUTING.md)
and describes how to contribute to this repository. Pull
requests containing bug fixes, updates, and extensions to the existing
tools and tool suites in this repository will be considered for
inclusion.

## How to Contribute

* Make sure you have a [GitHub account](https://github.com/signup/free)
* Make sure you have git [installed](https://help.github.com/articles/set-up-git)
* Fork the repository on [GitHub](https://github.com/ARTbio/tools-artbio/fork)
* Make the desired modifications - consider using a [feature branch](https://github.com/Kunena/Kunena-Forum/wiki/Create-a-new-branch-with-git-and-manage-branches).
* Make sure you have added the necessary tests for your changes and they pass.
* Make sure submitted tools meet IUC [Best Practices](https://galaxy-iuc-standards.readthedocs.io/en/latest/)
* Open a [pull request](https://help.github.com/articles/using-pull-requests)
  with these changes.

## What to contribute

* Wrappers for new tools
* Visualization Plugins
* Updates for tools
* Enhancements for tools (e.g. supporting new parameters)
* Bug fixes
* Documentation improvements
* New test cases in `<tests>` section of tools

### Abandoned Tools

* If this is of general interest, ARTbio may adopt tools that are not updated
or abandoned.
* If there are tools that you find useful but seem to be abandoned and not
  updated, you're welcome to create an issue recommending that the ARTbio adopts
  that tool.

## What not to contribute

* Tools already wrapped and currently maintained by other users, unless there are
clear indications for missing functionalities not supported in these tools.
* Wrappers without tests
* New datatypes
    * When possible, new datatypes should be added directly to the Galaxy
      codebase.

## Tests

Contributed tools must include test cases. They need not
necessarily cover all uses of the program, but should ensure that it is
generally working. The Galaxy Wiki has a
[page](https://wiki.galaxyproject.org/Admin/Tools/WritingTests) on writing
tests.

ARTbio strongly recommends testing with [planemo](https://github.com/galaxyproject/planemo/),
which provides a simple command line utility for testing functionality

```console
$ planemo test --install_galaxy my_tool.xml
```

## Requirements for Contributions

Before a PR will be accepted, the ARTbio has [some requirements](https://wiki.galaxyproject.org/Tools/BestPractices) on the
submitted code, which we will be happy to help you achieve if you need the assistance.

* Tools must contain tests and test-data
* The tools should pass linting by planemo (`planemo lint my_tool.xml`)
* If used, python codes must pass linting by [flake8](http://flake8.pycqa.org/en/latest/)
* The tools must pass tests by planemo (`planemo test --install_galaxy my_tool.xml`), which
are automatically performed by Travis CI in this repository
* If there is a relevant paper for the tool, it should be cited in a [citation](https://docs.galaxyproject.org/en/latest/dev/schema.html#tool-citations) block

