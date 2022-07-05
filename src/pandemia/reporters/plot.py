"""Reporters that output a plot"""

import matplotlib.pyplot as plt
from pandemia.reporters import Reporter
import numpy as np
import csv
from numpy import genfromtxt
import os

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
        # plt.plot(list(range(len(self.days))),
        #          self.infected, 'black', linewidth=1, alpha=1.0, label=plot_label)

        # vector_region = vector_world.vector_regions[0]
        # pop_by_age_group = np.zeros((17,), dtype=np.uint64)
        # initial_cases_by_age_group = np.zeros((17,), dtype=int)
        # for n in range(vector_region.number_of_agents):
        #     pop_by_age_group[vector_region.age_group[n]] += 1
        #     initial_cases_by_age_group[vector_region.age_group[n]] += vector_region.initial_cases[n]
        # self.plot_sir_age(pop_by_age_group, initial_cases_by_age_group)
        # num_initial_cases = np.sum(initial_cases_by_age_group)
        # self.infected.insert(0, num_initial_cases)

        plt.bar(list(range(len(self.days))), self.infected, label=plot_label)

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
            dirname = os.path.dirname(self.filename)
            if dirname != '':
                os.makedirs(dirname, exist_ok=True)
            fig.savefig(self.filename, bbox_inches='tight')

        # dirname = os.path.dirname('/tmp/infected.csv')
        # if dirname != '':
        #     os.makedirs(dirname, exist_ok=True)
        # handle = open('/tmp/infected.csv', 'w', newline='')
        # writer = csv.writer(handle)
        # for row in self.infected:
        #     writer.writerow([row])
        # handle.close()

    def plot_sir(self):
        """Plots solution to an SIR model"""

        step_size = 1
        gamma = 1/9
        beta = 0.27

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

        # dirname = os.path.dirname('/tmp/infected_sir.csv')
        # if dirname != '':
        #     os.makedirs(dirname, exist_ok=True)
        # handle = open('/tmp/infected_sir.csv', 'w', newline='')
        # writer = csv.writer(handle)
        # for row in I:
        #     writer.writerow([row])
        # handle.close()

    def plot_sir_age(self, pop_by_age_group, initial_cases_by_age_group):
        """Plots solution to an SIR model"""

        mat = genfromtxt("Scenarios/SIR/data/age_mixing_matrices/sir_region.csv",
                         delimiter=',', dtype=float)
        row_sums = mat.sum(axis=1)
        mat = mat / row_sums[:, np.newaxis]
        A = 17

        step_size = 1
        gamma = 1/ (8 + (23/24))
        beta = [1.5 * 0.27, 1.5 * 0.27, 1.5 * 0.27, 1.5 * 0.27, 1.5 * 0.27, 1.5 * 0.27, 1.5 * 0.27, 1.5 * 0.27, 1.5 * 0.27, 1.5 * 0.27, 1.5 * 0.27, 1.5 * 0.27, 1.5 * 0.27, 0.5 * 0.27, 0.5 * 0.27, 0.5 * 0.27, 0.5 * 0.27]

        T = len(self.days)
        N = pop_by_age_group
        S = np.zeros((T, A), dtype=float)
        I = np.zeros((T, A), dtype=float)
        R = np.zeros((T, A), dtype=float)

        # Initial conditions
        for b in range(A):
            R[0][b] = 0
            I[0][b] = initial_cases_by_age_group[b]
            S[0][b] = N[b] - I[0][b]

        # Simulate
        for t in range(T-1):
            for a in range(A):
                B = sum([mat[a][b] * ((I[t][b] * beta[b]) / (N[b])) for b in range(A)])
                S[t+1][a] = S[t][a] - (step_size * B * S[t][a])
                I[t+1][a] = I[t][a] + (step_size * B * S[t][a]) -\
                                      (step_size * gamma * I[t][a])
                R[t+1][a] = R[t][a] + (step_size * gamma * I[t][a])

        # Plot output
        plot_label = 'Infected SIR'
        plt.plot(list(range(T)), I.sum(axis=1),'red',
                 linewidth=1, alpha=1.0, label=plot_label)
