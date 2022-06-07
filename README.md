# Pandemia
![Integration](https://github.com/?/workflows/Integration/badge.svg?branch=master)
![Pytest](https://github.com/?/workflows/Pytest/badge.svg)
![Pylint](https://github.com/?/workflows/Pylint/badge.svg)
[![CodeFactor](https://www.codefactor.io/repository/github/?/badge?s=006dc8f386c6ea6d2a7a90377ff30fcf15328919)](https://www.codefactor.io/repository/github/?)

Pandemia is an individual-based pandemic simulator.

![pandemia Logo](pandemia_logo.jpg)

## Overview
Scenarios are configured using YAML.

The code is mixed Python and C.

The ABMlux model, on which Pandemia is based, has been used for the preprint Thompson, J. and Wattam, S. "Estimating the impact of interventions against COVID-19: from lockdown to vaccination", 2021, https://doi.org/10.1101/2021.03.21.21254049.

### Input Data
Input data are defined in the [Scenarios](Scenarios/) directory.

### Output Data
Output data are stored in a configured output directory.

## Requirements

 * python 3.9

## Usage

 * `pip install .`
 * `pandemia Scenarios/Luxembourg/config.yaml`

## Testing
To test:

    pip install .[test]
    pytest

## Docs
To generate documentation:

    pip install pdoc
    pdoc --html --overwrite --html-dir docs pandemia

## Citing This Work
If you publish using technology from this repository, please [give us a citation](https://www.medrxiv.org/content/10.1101/2021.03.21.21254049v1), using this BibTeX:

    @article {Thompson2021.03.21.21254049,
        author = {Thompson, James and Wattam, Stephen},
        title = {Estimating the impact of interventions against COVID-19: from lockdown to vaccination},
        elocation-id = {2021.03.21.21254049},
        year = {2021},
        doi = {10.1101/2021.03.21.21254049},
        publisher = {Cold Spring Harbor Laboratory Press},
        URL = {https://www.medrxiv.org/content/early/2021/03/26/2021.03.21.21254049},
        eprint = {https://www.medrxiv.org/content/early/2021/03/26/2021.03.21.21254049.full.pdf},
        journal = {medRxiv}
    }

## License
<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

Read the full text for details, but basically this means:
 * No commercial exploitation (contact us for another license in this case);
 * You must re-publish the source if you modify the application.

We would like this work to be useful to non-profit and academic users without significant effort.  If the license is an impediment to you using the work, please get in touch with us to discuss other licensing options.
