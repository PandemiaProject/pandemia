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

It is of vital importance that we continue to build our understanding of the epidemiology of such viruses and predict the impact of interventions, to help policy makers formulate effective strategies that save lives while simultaneously balancing their economic and social impact.

Such is the sweeping impact of a pandemic, that the difference between an optimal strategy and a sub-optimal strategy could equate to millions of lives lost, and widespread economic disruption, that could otherwise have been avoided.

# Statement of need

`Pandemia` is an individual-based stochastic pandemic simulator, able to simulate and visualize the spread of an infectious disease of humans across multiple geographical regions. The software is aimed at researchers looking to assess the impact of policy in the context of such a public health emergency.

In an individual-based model, the simultaneous actions and interactions of multiple individuals are simulated in an attempt to re-create and predict the emergence of complex phenomena resulting from their collective behaviour. This bottom-up approach contrasts with the top-down equation-based approach to epidemic modelling, an example of the latter being the SIR model. While individual-based models are more computationally intensive, they provide the most realistic descriptions of social interaction and infectious disease dynamics.

`Pandemia` is written in Python and C. Python facilitates a flexibile and user-friendly class-based API, while enabling the wrapping of low-level C functions for speed. These low-level functions provide fast implementations of disease transmission, changes to individual health and immunity, mobility, travel dynamics, diagnostics and vaccination. The model architecture was designed with parallelism in mind, making the model both fast and scalable, able to simulate large numbers of individuals while supporting a wide range of highly adaptable features.

`Pandemia` is based on the @ABMlux model, developed in Luxembourg between July 2020 and February 2021 and funded by the National Research Fund of Luxembourg. ABMlux was used in the article @Thompson:2021 to estimate the impact of interventions against COVID-19 in Luxembourg during the first six months of the pandemic. `Pandemia` represents a far-reaching generalization of ABMlux, adding speed and scalability along with a host of other features, including multiple viral variants, cross and waning immunity, multiple vaccines, international travel, seasonality and more.

In `Pandemia`, the policy maker is directly represented as a component in the simulation. The policy maker enacts interventions on a temporal and regional basis, including lockdowns, border closures, vaccination, testing and contact tracing and face masks. In the first instance, this allows the user to simulate a predefined strategy and assess its impact. Using optimization algorithms, `Pandemia` can take a step beyond this, searching systematically for optimal strategies which minimize the total cost of an outbreak. These optimization tools are, at the time of writing, a work in progress and represent the next step in the development of `Pandemia`.

TODO:
- cite some reviews of agent-based covid or epi models?
- cite SIR model if mentioning it?

# Acknowledgements

We acknowledge contributions from ?

# References