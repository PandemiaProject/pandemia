# Pandemia
<!-- ![Integration](https://github.com/?/workflows/Integration/badge.svg?branch=master)
![Pytest](https://github.com/?/workflows/Pytest/badge.svg)
![Pylint](https://github.com/?/workflows/Pylint/badge.svg)
[![CodeFactor](https://www.codefactor.io/repository/github/?/badge?s=006dc8f386c6ea6d2a7a90377ff30fcf15328919)](https://www.codefactor.io/repository/github/?) -->

[![End to End tests](https://github.com/PandemiaProject/pandemia/actions/workflows/end-to-end-tests.yml/badge.svg)](https://github.com/PandemiaProject/pandemia/actions/workflows/end-to-end-tests.yml)

Pandemia is an individual-based stochastic pandemic simulator, currently a work in progress. It is
able to simulate the spread of an infectious disease over multiple regions, for example the entire
world. The model is fast and scalable, able to simulate extremely large numbers of individuals, and
supports a wide range of features.

![pandemia Logo](pandemia_logo.jpg)

## Overview

Scenarios are configured using YAML.

The code is mixed Python and C.

### Input Data

Input data are defined in the [Scenarios](Scenarios/) directory.

### Output Data

Output data are stored in a configured output directory.

## Requirements

- Python 3.10
  - Python version other than 3.10 may produce errors.
- A GCC complier.

## Usage

To build C libraries

```
make
```

To install pandemia python package:

```
pip install .
```
To run the homogeneous mixing example:
```
pandemia Scenarios/Global/global_config.yaml
```
To run the heterogeneous mixing example:
```
pandemia Scenarios/Global_Grid/global_grid_config.yaml
```
To run the heterogeneous mixing example and save the world file:
```
pandemia Scenarios/Global_Grid/global_grid_config.yaml Scenarios/Global_Grid/global_grid_world.wld
```
To run the heterogeneous mixing example using the saved world file, thereby skipping the build phase:
```
pandemia Scenarios/Global_Grid/global_grid_config.yaml Scenarios/Global_Grid/global_grid_world.wld
```
The heterogeneous mixing example uses the following population grid data, available under a CC BY 4.0 license:

Center for International Earth Science Information Network - CIESIN - Columbia University. 2018.
Gridded Population of the World, Version 4 (GPWv4): Population Density, Revision 11. Palisades,
New York: NASA Socioeconomic Data and Applications Center (SEDAC). https://doi.org/10.7927/H49C6VHW.
Accessed 31 OCTOBER 2022.

To run the ABMlux scenario, an additional .abm file is required.

## Testing

To install additional dependencies required for testing:

```
pip install .[test]
```

### Unit tests (Future)

**Currently there are no unit tests**. When they've been written, they will be run using pytest:

```
pytest
```

### Integration tests

Integration tests (and other tests which take a long time to execute) should be marked with the `@pytest.mark.slow` decorator, eg:

```
@pytest.mark.slow
def test_long_processing_time():
    sleep(500)
```

These tests will **not** be run when `pytest` is called without arguments. (See [pytest.ini](pytest.ini) for details). To execute these tests, use the `-m slow` argument. eg:

```
pytest -m slow
```

### (Re)creating the "gold standard" outputs

The integration tests launch complete runs of pandemia and then compare the resulting output file with a set of "gold standard" files for each scenario. Occasionally (depending on the development of the relevant module) it may be necessary to recreate these. To recreate the gold standard outputs, use `pytest`'s `basetemp` dir option. **This can overwrite all the existing gold standard output files**

```
pytest -m slow --basetemp=./tests/recreate_gold_standard
```

This command can be combined with selecting individual tests if required.


### Test Coverage

To obtain a test


## Docs
To generate the documentation:
```
pip install pdoc
pdoc --html --overwrite --html-dir docs pandemia
```
## Citing This Work
Pandemia is based on [ABMlux](https://github.com/abm-covid-lux/abmlux), an epidemic model which was used in the article Thompson, J. and Wattam, S. "Estimating the impact of interventions against COVID-19: from lockdown to vaccination", 2021, PLOS ONE, https://doi.org/10.1371/journal.pone.0261330.

If you publish using technology from this repository, please cite the above article using this BibTeX:
``` 
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
```
## License
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
