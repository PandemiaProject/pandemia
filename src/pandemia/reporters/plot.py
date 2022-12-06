"""Reporters that output a plot"""

from collections import defaultdict
import matplotlib.pyplot as plt
from . import Reporter
import numpy as np
import csv
import datetime

from numpy import genfromtxt
import os

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init

class PlotDeaths(Reporter):
    """Reporter that returns a plot at the end of the simulation. This reporter plots the number of
    deaths (both daily deaths and cumulative deaths) each day. If a filepath is given to historical
    data, this reporter will plot that too, for comparison."""

    def __init__(self, telemetry_bus, config):
        super().__init__(telemetry_bus)

        self.show                       = config['show']
        self.savefig                    = config['savefig']
        self.daily_deaths_filename      = config['daily_deaths_filename']
        self.cumulative_deaths_filename = config['cumulative_deaths_filename']
        self.historical_data_fp         = config['historical_data']
        self.deaths_by_country_filename = config['deaths_by_country_filename']

        self.num_x_ticks = 20

        self.days = []
        self.scale_factor = None
        self.day_to_date = {}
        self.daily_deaths = {}
        self.cumulative_deaths = {}

        self.subscribe("deaths.initialize", self.initialize)
        self.subscribe("deaths.update", self.update_counts)
        self.subscribe("simulation.end", self.stop_sim)

    def initialize(self, region_names_to_super_regions, regions_omitted):
        """Called when the simulation starts."""

        self.regions_omitted = regions_omitted
        self.region_names = list(region_names_to_super_regions.keys())
        self.super_regions_dict = region_names_to_super_regions
        self.super_regions = set()
        for region_name in self.super_regions_dict:
            if region_name not in self.regions_omitted:
                super_region = self.super_regions_dict[region_name]
                if super_region is not None:
                    self.super_regions.add(super_region)
        self.super_regions = list(self.super_regions)

    def update_counts(self, clock, daily_deaths, cumulative_deaths):
        """Update counts"""

        date = clock.iso8601()
        self.days.append(clock.day)
        self.day_to_date[clock.day] = date
        self.daily_deaths[date] = daily_deaths
        self.cumulative_deaths[date] = cumulative_deaths

    def stop_sim(self):
        """Called when the simulation ends."""

        historical_daily_deaths_dict      = {self.day_to_date[day]: 0 for day in self.days}
        historical_cumulative_deaths_dict = {self.day_to_date[day]: 0 for day in self.days}
        if self.historical_data_fp is not None:
            with open(self.historical_data_fp, newline='') as csvfile:
                next(csvfile)
                deaths_data = csv.reader(csvfile, delimiter=',')
                for row in deaths_data:
                    date =\
                        datetime.datetime.strptime(str(row[0]), '%d/%m/%Y').strftime('%d/%m/%Y%Z')
                    historical_daily_deaths_dict[date]      = float(row[4])
                    historical_cumulative_deaths_dict[date] = int(row[1])
        hist_daily_deaths = []
        hist_cumulative_deaths = []
        daily_deaths = []
        cumulative_deaths = []
        for day in range(len(self.days)):
            hist_daily_deaths.append(historical_daily_deaths_dict[self.day_to_date[day]])
            hist_cumulative_deaths.append(historical_cumulative_deaths_dict[self.day_to_date[day]])
            min_day = max(0, day - 10)
            max_day = min(len(self.days), day + 10)
            num_days = max_day - min_day
            average_deaths = 0
            for day_sample in range(min_day, max_day):
                for region_name in self.region_names:
                    if region_name not in self.regions_omitted:
                        average_deaths +=\
                            self.daily_deaths[self.day_to_date[day_sample]][region_name]
            average_deaths /= num_days
            daily_deaths.append(average_deaths)
            cumulative_deaths_count = 0
            for region_name in self.region_names:
                if region_name not in self.regions_omitted:
                    cumulative_deaths_count +=\
                        self.cumulative_deaths[self.day_to_date[day]][region_name]         
            cumulative_deaths.append(cumulative_deaths_count)

        if self.deaths_by_country_filename is not None:
            handle = open(self.deaths_by_country_filename, 'w', newline='')
            writer = csv.writer(handle)
            for region_name in self.region_names:
                if region_name not in self.regions_omitted:
                    total_deaths =\
                        self.cumulative_deaths[self.day_to_date[self.days[-1]]][region_name]
                    row = [region_name, total_deaths]
                    writer.writerow(row)
            handle.close()

        if len(self.super_regions) > 0:
            daily_deaths_by_super_region = defaultdict(list)
            average_deaths_by_super_region = {sr: 0 for sr in self.super_regions}
            for day in range(len(self.days)):
                min_day = max(0, day - 10)
                max_day = min(len(self.days), day + 10)
                num_days = max_day - min_day
                for day_sample in range(min_day, max_day):
                    for region_name in self.region_names:
                        if region_name not in self.regions_omitted:
                            super_region = self.super_regions_dict[region_name]
                            average_deaths_by_super_region[super_region] +=\
                                self.daily_deaths[self.day_to_date[day_sample]][region_name]
                for super_region in self.super_regions:
                    average_deaths_by_super_region[super_region] /= num_days
                    average = average_deaths_by_super_region[super_region]
                    daily_deaths_by_super_region[super_region].append(average)

        plt.figure(figsize=(12, 6))

        font = {'size' : 10}

        plt.rc('font', **font)

        plt.bar(list(range(len(self.days))), daily_deaths, width=1.0, label='Deaths - Simulation',
                alpha=0.3)
        plt.plot(list(range(len(self.days))),
                 hist_daily_deaths, 'red', linewidth=1, linestyle='dotted',
                 alpha=1.0, label='Deaths - Historical')
        if len(self.super_regions) > 0:
            for super_region in self.super_regions:
                plt.plot(list(range(len(self.days))), daily_deaths_by_super_region[super_region],
                         linewidth=1, alpha=1.0, label='Deaths - Simulation: ' + super_region)

        plt.xlabel('Date')
        increment = int(min(len(self.days) // self.num_x_ticks, len(self.days)))
        plt.xticks(ticks=[i*increment for i in range((len(self.days) // increment))],
            labels=[self.day_to_date[self.days[i*increment]] for i in
                    range((len(self.days) // increment))])
        plt.xticks(rotation=90)
        plt.grid(False)
        plt.xlim([0, len(self.days)])
        plt.gca().set_ylim(bottom=0)
        plt.legend(loc='upper right')
        fig = plt.gcf()
        if self.show:
            plt.show()
        if self.savefig:
            dirname = os.path.dirname(self.daily_deaths_filename)
            if dirname != '':
                os.makedirs(dirname, exist_ok=True)
            fig.savefig(self.daily_deaths_filename, bbox_inches='tight')

        plt.figure(figsize=(12, 6))

        font = {'size' : 10}

        plt.rc('font', **font)

        plt.bar(list(range(len(self.days))), cumulative_deaths, label='Deaths - Simulation')
        plt.plot(list(range(len(self.days))),
                 hist_cumulative_deaths, 'red', linewidth=1, alpha=1.0, label='Deaths - Historical')
        plt.xlabel('Date')
        increment = int(min(len(self.days) // self.num_x_ticks, len(self.days)))
        plt.xticks(ticks=[i*increment for i in range((len(self.days) // increment))],
            labels=[self.day_to_date[self.days[i*increment]] for i in
                    range((len(self.days) // increment))])
        plt.xticks(rotation=90)
        plt.grid(False)
        plt.xlim([0, len(self.days)])
        plt.gca().set_ylim(bottom=0)
        plt.legend(loc='upper right')
        fig = plt.gcf()
        if self.show:
            plt.show()
        if self.savefig:
            dirname = os.path.dirname(self.cumulative_deaths_filename)
            if dirname != '':
                os.makedirs(dirname, exist_ok=True)
            fig.savefig(self.cumulative_deaths_filename, bbox_inches='tight')

class PlotInfected(Reporter):
    """Reporter that returns a plot at the end of the simulation. This reporter plots the number of
    people infected each day (daily prevalence)."""

    def __init__(self, telemetry_bus, config):
        super().__init__(telemetry_bus)

        self.show               = config['show']
        self.savefig            = config['savefig']
        self.filename           = config['filename']
        self.historical_data_fp = config['historical_data']

        self.num_x_ticks = 20

        self.days = []
        self.infected = []
        self.total_population = None
        self.scale_factor = None
        self.day_to_date = {}

        self.subscribe("strain_counts.initialize", self.initialize)
        self.subscribe("strain_counts.update", self.update_counts)
        self.subscribe("simulation.end", self.stop_sim)

    def initialize(self, number_of_regions, number_of_strains, region_names, population_sizes):
        """Called when the simulation starts."""

        self.total_population = sum(population_sizes)

    def update_counts(self, clock, strain_counts):
        """Update counts"""

        self.days.append(clock.day)
        self.day_to_date[clock.day] = clock.iso8601()
        self.infected.append(np.sum(strain_counts))

    def stop_sim(self):
        """Called when the simulation ends."""

        plt.figure(figsize=(12, 6))

        font = {'size' : 10}

        plt.rc('font', **font)

        plt.bar(list(range(len(self.days))), self.infected, label='Infected - Simulation')

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
        S = np.zeros((T,), dtype=np.float64)
        I = np.zeros((T,), dtype=np.float64)
        R = np.zeros((T,), dtype=np.float64)

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
                         delimiter=',', dtype=np.float64)
        row_sums = mat.sum(axis=1)
        mat = mat / row_sums[:, np.newaxis]
        A = 17

        step_size = 1
        gamma = 1/ (8 + (23/24))
        beta = [1.5 * 0.27, 1.5 * 0.27, 1.5 * 0.27, 1.5 * 0.27, 1.5 * 0.27, 1.5 * 0.27, 1.5 * 0.27,
                1.5 * 0.27, 1.5 * 0.27, 1.5 * 0.27, 1.5 * 0.27, 1.5 * 0.27, 1.5 * 0.27, 0.5 * 0.27,
                0.5 * 0.27, 0.5 * 0.27, 0.5 * 0.27]

        T = len(self.days)
        N = pop_by_age_group
        S = np.zeros((T, A), dtype=np.float64)
        I = np.zeros((T, A), dtype=np.float64)
        R = np.zeros((T, A), dtype=np.float64)

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
