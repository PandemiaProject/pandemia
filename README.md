# Pandemia
![Integration](https://github.com/?/workflows/Integration/badge.svg?branch=master)
![Pytest](https://github.com/?/workflows/Pytest/badge.svg)
![Pylint](https://github.com/?/workflows/Pylint/badge.svg)
[![CodeFactor](https://www.codefactor.io/repository/github/?/badge?s=006dc8f386c6ea6d2a7a90377ff30fcf15328919)](https://www.codefactor.io/repository/github/?)

Pandemia is an individual based pandemic simulator.

![pandemia Logo](pandemia_logo.jpg)

## Overview
Scenarios are configured using YAML.

The code is mixed Python and C.

### Input Data
Input data are defined in the [Scenarios](Scenarios/) directory.

### Output Data
Output data are stored in a configured output directory.

## Requirements

 * python 3.10

## Usage

 * `pip install .`
 * `pandemia Scenarios/Global/global_config.yaml`

## Testing
To test:

    pip install .[test]
    pytest

## Docs
To generate documentation:

    pip install pdoc
    pdoc --html --overwrite --html-dir docs pandemia

## Citing This Work
Pandemia is based on the [ABMlux model](https://github.com/abm-covid-lux/abmlux), which was used for the article Thompson, J. and Wattam, S. "Estimating the impact of interventions against COVID-19: from lockdown to vaccination", 2021, PLOS ONE, https://doi.org/10.1371/journal.pone.0261330.

If you publish using technology from this repository, please cite [the above article](https://doi.org/10.1371/journal.pone.0261330), using this BibTeX:

    @article{10.1371/journal.pone.0261330,
        doi = {10.1371/journal.pone.0261330},
        author = {Thompson, James AND Wattam, Stephen},
        journal = {PLOS ONE},
        publisher = {Public Library of Science},
        title = {Estimating the impact of interventions against COVID-19: From lockdown to vaccination},
        year = {2021},
        month = {12},
        volume = {16},
        url = {https://doi.org/10.1371/journal.pone.0261330},
        pages = {1-51},
        number = {12},
    }

## License
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
