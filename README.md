# Pandemia
<!-- ![Integration](https://github.com/?/workflows/Integration/badge.svg?branch=master)
![Pytest](https://github.com/?/workflows/Pytest/badge.svg)
![Pylint](https://github.com/?/workflows/Pylint/badge.svg)
[![CodeFactor](https://www.codefactor.io/repository/github/?/badge?s=006dc8f386c6ea6d2a7a90377ff30fcf15328919)](https://www.codefactor.io/repository/github/?) -->

![pandemia Logo](pandemia_logo.jpg)

Pandemia is an individual-based stochastic pandemic simulator. It is currently a work in progress.

Pandemia is able to simulate and visualize the spread of an infectious disease across multiple
geographical regions, for example the countries of the world or the administrative divisions of
a country. The model is fast and scalable, able to simulate extremely large numbers of individuals
while supporting a wide range of highly adaptable features.

This software can be used by researchers looking to assess the impact of policy in the context of a
public health emergency caused by an infectious disease of humans. The disease could, for example, be a
respiratory infectious disease spread by a coronavirus or an influenza virus. The emergence and
re-emergence of infectious diseases threatens the health and well-being of people all over the world,
and tools such as Pandemia can play a vital role in supporting pandemic preparedness and response.

This document contains an overview of the model, a quickstart guide and information for
contributors. The full user guide can be found [here](docs/user_guide.pdf).

The code is mixed Python and C.

## Overview
The Pandemia simulator acts upon a **World**. A **World** consists of a number of objects of type
**Region** and a **Travel Matrix**. The **Travel Matrix** represents the number of individuals travelling
from each **Region** to each other **Region** each day. Each **Region** consists of individuals (referred to
as agents), locations and activities. A **World** is built using a **World Factory**.

After the building, the **World** is then converted into a **Vector World**. This is done by
converting each **Region** into a **Vector Region**. A **Vector Region** is a vectorized version of
a **Region**, in which data is formatted as arrays of integers and floats, as opposed to Python
lists and dictionaries. This facilitates interface with libraries of functions written in C.

Building the **Clock** and **World** are the first steps in building a **Simulator**. Once these
objects are built, Pandemia then builds a number of objects of type **Component**. These components
include models of:

* Health
* Movement
* Hospitalization
* Travel
* Vaccination
* Diagnostics
* Seasonality
* Input

The models provided support, for example, multiple strains, multiple vaccines, a novel approach to
modelling the health of individuals and a complex representation of immunity. The input component
allows the user to specify a **Policy**, consisting of interventions. Featured interventions
include:

* Border Closures
* Lockdowns
* Testing and Contact Tracing
* Quarantine
* Face Masks
* Vaccination

Objects of type **Reporter** collect output data for visualization and analysis.

A number of **World Factory** and **Component** examples are provided for the user. In particular,
for each **Component**, a default model is provided, as well as a void model in case the user does
not wish for this component to be active during a simulation. Among the **World Factory** examples
are **Homogeneous** and **Heterogeneous**. Both these factories build all the countries in the world, using
air travel data to configure travel between countries. However, whereas **Homogeneous** implements
homogeneous mixing within each country, **Heterogeneous** implements a simple model of hetergeneous
mixing, based on average household size and population density grids. **Heterogeneous** also allows
the user to limit the simulation to a chosen subset of countries. In particular, the user can run
**Heterogeneous** on only a single country, if the user wishes. For both of these world factories,
the recommended scale factor is 0.0005.

Scenarios are configured using YAML. A scenario consists of a choice of world factory, and a choice
of submodel for each of the simulation components, together with configurations for each of these
objects and the reporters. Example scenarios can be found in the [Scenarios](Scenarios/)
directory.

The homogeneous mixing scenario uses the **Homogeneous** world factory. In this scenario,
individuals mix homogeneously within each region, with mixing between regions being determined using
air travel data. Colouring regions according to prevalance, the scenario can be visualized as
follows:

![pandemia Logo](pandemia_homogeneous.jpg)

The heterogeneous mixing scenario uses the **Heterogeneous** world factory. In this scenario,
individuals mix heterogeneously within each region. This uses data on average household size,
population distribution and a simple gravity model of mobility. Colouring grid squares according to
numbers infected, the scenario has the following visualization:

![pandemia Logo](pandemia_heterogeneous.jpg)

For each scenario, all parameter values are set in the corresponding YAML file. For example, the
**Homogeneous** config can be found [here](Scenarios/Homogeneous/homogeneous_config.yaml) while the **Heterogeneous**
config can be found [here](Scenarios/Heterogeneous/heterogeneous_config.yaml).

### Input Data
Input data for each scenario are found in the [Scenarios/](Scenarios/) directory. For example, all
input data for the **Homogeneous** scenario are found in [Scenarios/Homogeneous/data](Scenarios/Homogeneous/data).
All input data for the **Heterogeneous** scenario are found in [Scenarios/Heterogeneous/data](Scenarios/Heterogeneous/data).

The **Heterogeneous** world factory uses the following grid data, available under a CC BY 4.0 license:

Center for International Earth Science Information Network - CIESIN - Columbia University. 2018.
Gridded Population of the World, Version 4 (GPWv4): Population Density, Revision 11. Palisades,
New York: NASA Socioeconomic Data and Applications Center (SEDAC). https://doi.org/10.7927/H49C6VHW.
Accessed 31 OCTOBER 2022.

### Output Data
Output data are stored in a output directory, configured by the user in the reporters section of the
scenario configuration.

## Requirements

- Python 3.10
  - Python version other than 3.10 may produce errors. 
- A GCC complier 

## Quickstart

To build C libraries
```
make
```
To install pandemia python package:
```
pip install .
```
To run the homogeneous mixing scenario:
```
pandemia Scenarios/Homogeneous/homogeneous_config.yaml
```
To run the heterogeneous mixing scenario:
```
pandemia Scenarios/Heterogeneous/heterogeneous_config.yaml
```
To run the heterogeneous mixing scenario and save after the world building phase:
```
pandemia Scenarios/Heterogeneous/heterogeneous_config.yaml Scenarios/Heterogeneous/heterogeneous_world.wld
```
To run the heterogeneous mixing scenario using the save, thereby skipping the world building phase:
```
pandemia Scenarios/Heterogeneous/heterogeneous_config.yaml Scenarios/Heterogeneous/heterogeneous_world.wld
```
To configure a new scenario, the user should choose one of the configs already provided, as a
template, which they should then edit as necessary. These configs are annotated with comments
briefly explaining the meaning of each of the parameters. For more information, consult the user
manual.

The advanced user might even wish to write their own models, to be used instead of the default
models provided.

## Testing
To test:
```
pip install .[test]
pytest
```

## Documentation
Consult the [User Guide](docs/user_guide.pdf), or generate documentation using:
```
pip install pdoc
pdoc --html --overwrite --html-dir docs pandemia
```

## Contributors

## Acknowledgements

## Citing this work
Pandemia is based on [ABMlux](https://github.com/abm-covid-lux/abmlux), an epidemic model used in the article Thompson, J. and Wattam, S. "Estimating the impact of interventions against COVID-19: from lockdown to vaccination", 2021, PLOS ONE, https://doi.org/10.1371/journal.pone.0261330. In that article, the authors presented an agent-based model of the COVID-19 pandemic in Luxembourg, and used it simulate the impact of interventions over the first 6 months of the pandemic. Pandemia is a far-reaching generalization of the ABMlux model.

If you publish using technology from this repository, please cite the above article using this BibTeX:
```
@article{10.1371/journal.pone.0261330,
    doi = {10.1371/journal.pone.0261330},
    author = {Thompson, James AND Wattam, Stephen},
    journal = {PLoS One},
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
