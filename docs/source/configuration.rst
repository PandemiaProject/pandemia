.. Just include the readme from the root dir of the repo

Configuration Files
===================

The pandemia configuration file is a yaml file. An annotated example is shown below.

The example given is of the Heterogeneous Mixing scenerio with interventions. The user
can create a new example by making a copy of one of the yaml files provided, before
editing the values contained within. The configuration files are divided into sections,
reflecting the structure of the Pandemia model.

Basic Parameters
----------------

.. include:: ../../Scenarios/Heterogeneous/heterogeneous_interventions.yaml
   :start-after: # Basic Parameters #
   :end-before: # Clock #
   :code:

Clock
-----

.. include:: ../../Scenarios/Heterogeneous/heterogeneous_interventions.yaml
   :start-after: # Clock #
   :end-before: # World #
   :code:

World
-----

.. include:: ../../Scenarios/Heterogeneous/heterogeneous_interventions.yaml
   :start-after: # World #
   :end-before: # Seasonal Effects #
   :code:

Models Components
-----------------

.. include:: ../../Scenarios/Heterogeneous/heterogeneous_interventions.yaml
   :start-after: # Seasonal Effects #
   :end-before: # Output #
   :code:

Outputs
-------

.. include:: ../../Scenarios/Heterogeneous/heterogeneous_interventions.yaml
   :start-after: # Output #
   :end-before: # Logging #
   :code:

Logging
-------

.. include:: ../../Scenarios/Heterogeneous/heterogeneous_interventions.yaml
   :start-after: # Logging #
   :code:

