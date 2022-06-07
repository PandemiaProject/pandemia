"""Reporters that output a plot"""

import matplotlib.pyplot as plt
from pandemia.reporters import Reporter
import numpy as np

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class PlotInfected(Reporter):
    """Reporter that returns a plot at the end of the simulation."""

    def __init__(self, telemetry_bus, config):
        super().__init__(telemetry_bus)

        self.show     = config['show']
        self.savefig  = config['savefig']
        self.filename = config['filename']

        self.num_x_ticks = 20

        self.days = []
        self.infected = []
        self.total_population = None
        self.scale_factor = None

        self.subscribe("strain_counts.initialize", self.initialize)
        self.subscribe("strain_counts.update", self.update_counts)
        self.subscribe("simulation.end", self.stop_sim)

    def initialize(self, number_of_regions, number_of_strains, region_names, population_sizes):
        """Called when the simulation starts."""

        self.total_population = sum(population_sizes)

    def update_counts(self, clock, strain_counts):
        """Update counts"""

        self.days.append(clock.day)
        self.infected.append(np.sum(strain_counts))

    def stop_sim(self):
        """Called when the simulation ends."""

        plt.figure(figsize=(12, 6))

        font = {'size' : 10}

        plt.rc('font', **font)
        plot_label = 'Infected'
        plt.plot(list(range(len(self.days))),
                 self.infected,'black', linewidth=1, alpha=1.0, label=plot_label)

        # self.plot_sir()

        plt.xlabel('Day')
        increment = int(max(len(self.days) // self.num_x_ticks, len(self.days)))
        plt.xticks(ticks=[i*increment for i in range((len(self.days) // increment))],
            labels=[self.days[i*increment] for i in range((len(self.days) // increment))])
        plt.grid(False)
        plt.xlim([0, len(self.days)])
        plt.gca().set_ylim(bottom=0)
        plt.legend(loc='upper right')
        fig = plt.gcf()
        if self.show:
            plt.show()
        if self.savefig:
            fig.savefig(self.filename, bbox_inches='tight')

    def plot_sir(self):
        """Plots solution to an SIR model"""

        step_size = 1
        gamma = 1/5
        beta = 0.48

        T = len(self.days)
        N = self.total_population
        S = np.zeros((T,), dtype=float)
        I = np.zeros((T,), dtype=float)
        R = np.zeros((T,), dtype=float)

        # Initial conditions
        R[0] = 0
        I[0] = self.infected[0]
        S[0] = N - I[0]

        # Simulate
        for t in range(T-1):
            S[t+1] = S[t] - (step_size * beta * (I[t] / N) * S[t])
            I[t+1] = I[t] + (step_size * beta * (I[t] / N) * S[t]) - (step_size * gamma * I[t])
            R[t+1] = R[t] + (step_size * gamma * I[t])

        # Plot output
        plot_label = 'Infected SIR'
        plt.plot(list(range(len(self.days))), I,'red', linewidth=1, alpha=1.0, label=plot_label)
