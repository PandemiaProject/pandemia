"""Reporters that output to CSV files"""

import os
import os.path
import csv
import numpy as np

from .reporters import Reporter

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class StrainCounts(Reporter):
    """Reporter that writes to a CSV file as it runs. This reporter records the number of people in
    each region infected with each strain, each day."""

    def __init__(self, telemetry_bus, config):
        super().__init__(telemetry_bus)

        self.filename = config['filename']

        self.subscribe("strain_counts.initialize", self.initialize)
        self.subscribe("strain_counts.update", self.update_counts)
        self.subscribe("simulation.end", self.stop_sim)

    def initialize(self, number_of_regions, number_of_strains, region_names, population_sizes):
        """Called when the simulation starts.  Writes headers and creates the file handle."""

        self.number_of_regions = number_of_regions
        self.number_of_strains = number_of_strains

        # Check dir exists and open handle
        dirname = os.path.dirname(self.filename)
        if dirname != '':
            os.makedirs(dirname, exist_ok=True)
        self.handle = open(self.filename, 'w', newline='')
        self.writer = csv.writer(self.handle)

        # Write header
        header = ["day", "iso8601"] + [str(region_names[r]) + ", Strain: " + str(s)
                                       for r in range(number_of_regions)
                                       for s in range(number_of_strains)]
        self.writer.writerow(header)

    def update_counts(self, clock, strain_counts):
        """Update the CSV, writing a single row for every clock tick"""

        row = [clock.day, clock.iso8601()]
        row += strain_counts.tolist()
        self.writer.writerow(row)

    def stop_sim(self):
        """Called when the simulation ends.  Closes the file handle."""

        if self.handle is not None:
            self.handle.close()
