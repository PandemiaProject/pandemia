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
    corresponding: true
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Stephen Wattam
    affiliation: 2
  - name: Andrew Smith
    affiliation: 3
  - name: Aoife Hughes
    affiliation: 4
affiliations:
 - name: Department of Infectious Disease Epidemiology, Imperial College London, UK
   index: 1
 - name: WAP Academic Consulting, UK
   index: 2
 - name: The Alan Turing Institute, UK
   index: 3
 - name: The Alan Turing Institute, UK
   index: 4
date: 21 March 2023
bibliography: paper.bib
---

# Summary

Emerging infectious diseases pose a significant threat to the physical and mental well-being of individuals across the globe. The ongoing COVID-19 pandemic, for example, is one of the most disruptive worldwide events in recent history. The SARS-CoV-2 virus has led to millions of deaths and hundreds of millions of confirmed cases [@WHO], with potential future pandemics posing even graver consequences [@GPMB].

It is of vital importance that we deepen our understanding of the epidemiology of pandemic pathogens, and predict the impact of public health interventions, to help policy makers respond to future pandemics with cost-effective strategies that not only save lives but also balance the economic and social implications.

Such might be the devastating impact of a future pandemic, that the difference between an optimal strategy and a sub-optimal strategy could equate to millions of lives lost, and widespread economic and social disruption, that could otherwise have been averted.

# Statement of need

While a large number of mathematical and computational epidemic models have been developed over the course of the COVID-19 pandemic, many are lacking in key features and limited in scope [@LORIG2021].

A useful pandemic model must encompass a wide range of advanced features, including spatial dynamics, international travel, multiple variants, waning and cross immunity, and multiple vaccines of differing efficacies. As evidenced during the COVID-19 pandemic, these factors cannot be ignored [@KANG202096] [@FRANCHPARDO2020140033] [@CHINAZIIETAL] [@CASCELLAETAL].

Moreover, there is a need for flexible model architecture capable of supporting fast and scalable simulations and variable resolution at the level of the synthetic populations. While ample data exists for modeling populations in many high-income countries, the same level of detail is often unattainable for a substantial portion of the global population [@WORLDBANK].

The global ramifications of a pandemic necessitate the adoption of robust, transparent, and well-documented scientific approaches, which rely on open-source software that is readily available to researchers worldwide, without mandating the use of costly high-performance computing resources.

`Pandemia` is an individual-based stochastic pandemic simulator, able to simulate and visualize the spread of an infectious disease of humans across multiple geographical regions, that addresses these needs.

In an individual-based model, the simultaneous actions and interactions of multiple individuals are simulated in an attempt to re-create and predict the emergence of complex phenomena resulting from their collective behaviour. This bottom-up approach contrasts with the top-down equation-based approach of Kermack–McKendrick theory [@KERMACKMCKENDRICK], an example of the latter being the SIR model. While individual-based models are more computationally intensive than equation-based models, they provide the closest approximations of social interaction and infectious disease dynamics.

The `Pandemia` software is written in Python and C. Python facilitates a flexible and user-friendly class-based API, while enabling the wrapping of low-level C functions for speed. These low-level functions provide efficient implementations of disease transmission, changes to individual health and immunity, mobility, travel dynamics, diagnostics and vaccination. The model architecture has been designed with parallelism in mind, making the model both fast and scalable, able to simulate large numbers of individuals while supporting a wide range of adaptable features.

`Pandemia` is based on the ABMlux model, developed in Luxembourg between July 2020 and February 2021 and funded by the National Research Fund of Luxembourg [@ABMLUX]. ABMlux was used in the article [@THOMPSON2021] to estimate the impact of interventions against COVID-19 in Luxembourg during the first six months of the pandemic. `Pandemia` represents a far-reaching generalization of this earlier work.

In `Pandemia`, the policy maker is directly represented as a component in the simulation. The policy maker enacts interventions on a temporal and regional basis, including lockdowns, border closures, vaccination, testing, contact tracing and face masks. In the first instance, this allows the user to simulate predefined strategies and assess their impact. Using optimization algorithms, `Pandemia` can take a step beyond this, searching systematically for optimal strategies which minimize the total cost of an outbreak. These optimization tools are, at the time of writing, a work in progress and represent the next step in the development of the `Pandemia` software.

# Acknowledgements

The Pandemia software was created by James Thompson in early 2022, based on the ABMlux software written by Stephen Wattam and James Thompson. Between October 2022 and February 2023, Andy Smith and Aoife Hughes contributed to the project as members of the Research Engineering Group at The Alan Turing Institute.

Since June 2022, James Thompson has been employed as a Research Associate at the Department of Infectious Disease Epidemiology at Imperial College London, having been previously employed by The Alan Turing Institute, between April 2021 and May 2022.

Stephen Wattam contributed to the ABMlux project via WAP Academic Consulting Ltd.

# References