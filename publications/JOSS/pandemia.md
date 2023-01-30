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
date: 30 January 2023
bibliography: pandemia.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Emerging infectious diseases threaten the health and well-being of people all over the world.

The ongoing COVID-19 pandemic, for example, is among the most disruptive global events in modern history. The SARS-CoV-2 virus has spread across world, resulting in approximately seven million confirmed deaths and a hundred times as many confirmed cases, with the true death toll being likely much higher. The next pandemic could be even worse.

It is of vital importance that we continue to build our understanding of the epidemiology of such viruses and predict the impact of interventions, to help policy makers formulate effective strategies that save lives while simultaneously balancing the economic and social impact.

Such is the sweeping impact of a pandemic, that the difference between an optimal strategy and a sub-optimal one could equate to millions of lives lost, and widespread economic disruption, that could otherwise have been avoided.

# Statement of need

`Pandemia` is an individual-based stochastic epidemic simulator, able to simulate and visualize the spread of an infectious disease across multiple geographical regions. It is written in Python and C. Python facilitates a flexibile and user-friendly class-based API, while enabling the wrapping of low-level C functions for speed. These low-level functions provide fast implementations for disease transmission, changes to individual health and immunity, mobility, travel dynamics, diagnostics and vaccination. The model architecture naturally supports parallelism, making the model both fast and scalable, able to simulate extremely large numbers of individuals while supporting a wide range of highly adaptable features.

In an individual-based model, the simultaneous actions and interactions of multiple individuals are simulated in an attempt to re-create and predict the emergence of complex phenomena resulting from their collective behaviour. This bottom-up approach contrasts with the top-down equation-based approach to epidemic modelling, an example of the latter being the SIR model. While individual-based models are more computationally intensive, they provide the most realistic descriptions of social interaction and infectious disease dynamics.

This software will be of great use to researchers looking to assess the impact of policy in the context of a global public health emergency caused by an infectious disease of humans.

`Pandemia` is based on the @ABMlux model, developed in Luxembourg between July 2020 and February 2021 and funded by the National Research Fund of Luxembourg. ABMlux was used in the article @Thompson:2021 to estimate the impact of interventions against COVID-19 in Luxembourg. `Pandemia` represents a far-reaching generalization of ABMlux, adding speed and scalability along with a host of other features, including support for multiple virus strains with cross and waning immunity, multiple vaccines, international travel, seasonality and more.

# Acknowledgements

We acknowledge contributions from ?

# References