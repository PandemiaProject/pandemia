---
title: 'Pandemia: A stochastic pandemic simulator'
tags:
  - Python
  - pandemic
  - epidemic
  - COVID-19
  - agent-based
authors:
  - name: James Thompson
    orcid: 0000-0000-0000-0000
    corresponding: true
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Stephen Wattam
    orcid: 0000-0000-0000-0000
    affiliation: 2
  - name: Andrew Smith
    orcid: 0000-0000-0000-0000
    affiliation: 3
  - name: Aoife Hughes
    orcid: 0000-0000-0000-0000
    affiliation: 4
affiliations:
 - name: Research Associate, Imperial College London, UK
   index: 1
 - name: xxx, WAP Academic Consulting, UK
   index: 2
 - name: xxx, The Alan Turing Institute, UK
   index: 3
 - name: xxx, The Alan Turing Institute, UK
   index: 4
date: 21 March 2023
bibliography: pandemia.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Emerging infectious diseases threaten the physical and mental health of people all over the world.

The ongoing COVID-19 pandemic, for example, is among the most disruptive global events in modern history. The SARS-CoV-2 virus has resulted in millions of deaths and hundreds of millions of confirmed cases @WHO. A future pandemic could be an order of magnitude worse @GPMB.

It is of vital importance that we continue to build our understanding of the epidemiology of pandemic pathogens, and predict the impact of public health interventions, to help policy makers respond to future pandemics with cost-effective strategies, that save lives while balancing the economic and social impact.

Such might be the devastating impact of a future pandemic, that the difference between an optimal strategy and a sub-optimal strategy could equate to millions of lives lost, and widespread economic and social disruption, that could otherwise have been avoided.

# Statement of need

While a large number of mathematical and computational epidemic models have been developed over the course of the COVID-19 pandemic, many are lacking in key features and limited in scope @LORIG2021.

A useful pandemic model must encompass a wide range of advanced features, including spatial dynamics, international travel, multiple variants, waning and cross immunity, and multiple vaccines of differing efficacy. As evidenced during the COVID-19 pandemic, these factors cannot be ignored @KANG202096 @FRANCHPARDO2020140033 @CHINAZIIETAL @CASCELLAETAL.

Moreover, there is a need for flexible model architecture that supports fast and scalable simulations and variable resolution at the level of the synthetic populations. While for many high-income countries, sufficient data are available to model those populations in some detail, this is simply not the case throughout much of the rest of the world @WORLDBANK.

The universal impact of a pandemic calls for scientific methods that are well-documented and transparent, based on open-source software accessible to researchers around the world, not requiring the expensive resources of a high-performance computer.

`Pandemia` is an individual-based stochastic pandemic simulator, able to simulate and visualize the spread of an infectious disease of humans across multiple geographical regions, that addresses these needs.

In an individual-based model, the simultaneous actions and interactions of multiple individuals are simulated in an attempt to re-create and predict the emergence of complex phenomena resulting from their collective behaviour. This bottom-up approach contrasts with the top-down equation-based approach of Kermack–McKendrick theory @KERMACKMCKENDRICK, an example of the latter being the SIR model. While individual-based models are more computationally intensive than equation-based models, they provide the closest approximations of social interaction and infectious disease dynamics.

The `Pandemia` software is written in Python and C. Python facilitates a flexible and user-friendly class-based API, while enabling the wrapping of low-level C functions for speed. These low-level functions provide efficient implementations of disease transmission, changes to individual health and immunity, mobility, travel dynamics, diagnostics and vaccination. The model architecture has been designed with parallelism in mind, making the model both fast and scalable, able to simulate large numbers of individuals while supporting a wide range of adaptable features.

`Pandemia` is based on the ABMlux model, developed in Luxembourg between July 2020 and February 2021 and funded by the National Research Fund of Luxembourg @ABMLUX. ABMlux was used in the article @THOMPSON2021 to estimate the impact of interventions against COVID-19 in Luxembourg during the first six months of the pandemic. `Pandemia` represents a far-reaching generalization of this earlier work.

In `Pandemia`, the policy maker is directly represented as a component in the simulation. The policy maker enacts interventions on a temporal and regional basis, including lockdowns, border closures, vaccination, testing, contact tracing and face masks. In the first instance, this allows the user to simulate predefined strategies and assess their impact. Using optimization algorithms, `Pandemia` can take a step beyond this, searching systematically for optimal strategies which minimize the total cost of an outbreak. These optimization tools are, at the time of writing, a work in progress and represent the next step in the development of the `Pandemia` software.

# Acknowledgements

We acknowledge contributions from ?

# References