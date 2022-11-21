# Pandemia
<!-- ![Integration](https://github.com/?/workflows/Integration/badge.svg?branch=master)
![Pytest](https://github.com/?/workflows/Pytest/badge.svg)
![Pylint](https://github.com/?/workflows/Pylint/badge.svg)
[![CodeFactor](https://www.codefactor.io/repository/github/?/badge?s=006dc8f386c6ea6d2a7a90377ff30fcf15328919)](https://www.codefactor.io/repository/github/?) -->

Pandemia is an individual-based stochastic pandemic simulator. It is currently a work in progress.

![pandemia Logo](pandemia_logo.jpg)

Pandemia is able to simulate the spread of an infectious disease across multiple regions. These
regions might, for example, represent all the countries of the world, or the administrative
divisions of a single country. The model is fast and scalable, able to simulate extremely large
numbers of individuals while supporting a wide range of features.

The code is mixed Python and C.

## Overview
The Pandemia simulator acts upon a **World**. A **World** consists of regions and a **Travel Matrix**. The
travel matrix represents the number of individuals travelling from each region to each other region
each day. Each region consists of individuals (referred to as agents), locations and activities. A
**World** is built using a **World Factory**.

The **World** is then converted into a **Vector World**. This is done by converting each **Region** into a 
**Vector Region**. A **Vector Region** is a vectorized version of a **Region**, in which data is formatted as
arrays of integers and floats, as opposed to Python lists and dictionaries. This facilitates
interface with libraries of functions written in C.

Building the **World** and **Clock** is the first step in building a **Simulator**. Once these objects are
built, Pandemia builds the simulation components. Each **Component** represents a pandemic submodel.
These include models of health, movement, hospitalization, travel, vaccination, diagnostics,
seasonality and input.

The input component allows the user to specify a **Policy**, consisting of interventions. Featured
interventions include vaccination, border closure, lockdown, testing and contact tracing,
quarantine and face masks.

Reporters collect output data for visualization and analysis.

A number of **World Factory** and **Component** examples are provided for the user. In particular, for each
**Component**, a default model is provided, as well as a void model in case the user does not wish for
this component to be active during a simulation. Among the **World Factory** examples are **Global** and
**Global Grid**. Both factories model all the countries in the world, using air travel data to
configure travel between countries. Whereas **Global** implements homogeneous mixing within each
country, **Global Grid** implements a simple model of hetergeneous mixing, based on average household
size and population density grids. **Global Grid** also allows the user to limit the simulation to a
subset of countries. For both of these world factories, the recommended scale factor is
0.0005.

Scenarios are configured using YAML. A scenario consists of a choice of world factory, and a choice
of submodel for each of the simulation components, together with configurations for each of these
selected objects and the reporters. Example scenarios can be found in the [Scenarios](Scenarios/)
directory. The homogeneous mixing scenerio **Global** uses the **Global** world factory, while the
heterogeneous mixing scenario **Global Grid** uses the **Global Grid** world factory.

### Input Data
Input data for each scenario are found in the [Scenarios/](Scenarios/) directory. For example, all
input data for the **Global** scenario is found in the [Scenarios/Global/data](Scenarios/Global/data).
All input data for the **Global Grid** scenario is found in the [Scenarios/Global_Grid/data](Scenarios/Global_Grid/data).

The **Global Grid** factory uses the following grid data, available under a CC BY 4.0 license:

Center for International Earth Science Information Network - CIESIN - Columbia University. 2018.
Gridded Population of the World, Version 4 (GPWv4): Population Density, Revision 11. Palisades,
New York: NASA Socioeconomic Data and Applications Center (SEDAC). https://doi.org/10.7927/H49C6VHW.
Accessed 31 OCTOBER 2022.

### Output Data
Output data are stored in a output directory, configured by the user in the reporters section of the
scenario configuration.

## Requirements

 * python 3.10

## Usage
To install:

    pip install -e .[test]

To run the homogeneous mixing scenario:

    pandemia Scenarios/Global/global_config.yaml

To run the heterogeneous mixing scenario:

    pandemia Scenarios/Global_Grid/global_grid_config.yaml

To run the heterogeneous mixing scenario and save after the world building phase:

    pandemia Scenarios/Global_Grid/global_grid_config.yaml Scenarios/Global_Grid/global_grid_world.wld

To run the heterogeneous mixing scenario using the save, thereby skipping the world building phase:

    pandemia Scenarios/Global_Grid/global_grid_config.yaml Scenarios/Global_Grid/global_grid_world.wld

To configure a new scenario, the user can simply edit one of the configs already provided.

## Testing
To test:

    pip install .[test]
    pytest

## Docs
To generate documentation:

    pip install pdoc
    pdoc --html --overwrite --html-dir docs pandemia

## Citing This Work
Pandemia is based on [ABMlux](https://github.com/abm-covid-lux/abmlux), an epidemic model which was used in the article Thompson, J. and Wattam, S. "Estimating the impact of interventions against COVID-19: from lockdown to vaccination", 2021, PLOS ONE, https://doi.org/10.1371/journal.pone.0261330.

If you publish using technology from this repository, please cite the above article using this BibTeX:

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
