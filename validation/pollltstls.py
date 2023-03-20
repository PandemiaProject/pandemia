
import logging
import numpy as np
import os
import datetime
import matplotlib.pyplot as plt
import matplotlib
import csv

from collections import defaultdict
from tqdm import tqdm

def get_data(directory, label):

    dates_to_data = defaultdict(list)
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        if os.path.isfile(f):
            length_label = len(label)
            if filename[0:length_label] == label:
                with open(f, newline='') as csvfile:
                    next(csvfile)
                    data = csv.reader(csvfile, delimiter=',')
                    for row in data:
                        date = datetime.datetime.strptime(str(row[1]), '%d/%m/%Y')
                        total = float(row[-1])
                        dates_to_data[date].append(total)

    return dates_to_data

def smooth(array):
    kernel_size = 14
    kernel = np.ones(kernel_size) / kernel_size
    array = np.convolve(array, kernel, mode='same')
    return array

def get_days(dates_to_data):

    day = 0
    days = []
    day_to_date = {}
    start_date = min(dates_to_data)
    end_date = max(dates_to_data)
    date = start_date
    while date <= end_date:
        days.append(day)
        day_to_date[day] = date
        date += datetime.timedelta(days=1)
        day += 1
    number_of_days = len(days)

    return days, day_to_date, number_of_days

def get_arrays(dates_to_data):

    averages = []
    start_date = min(dates_to_data)
    end_date = max(dates_to_data)
    date = start_date
    while date <= end_date:
        averages.append(dates_to_data[date])
        date += datetime.timedelta(days=1)

    averages = np.array(averages).astype(np.float64)

    y_mean = np.mean(averages, axis=1, dtype=np.float64)
    y_median = np.median(averages, axis=1)
    y_025 = np.percentile(averages, 2.5, axis=1)
    y_25 = np.percentile(averages, 25, axis=1)
    y_75 = np.percentile(averages, 75, axis=1)
    y_975 = np.percentile(averages, 97.5, axis=1)

    y_mean = smooth(y_mean)
    y_median = smooth(y_median)
    y_025 = smooth(y_025)
    y_25 = smooth(y_25)
    y_75 = smooth(y_75)
    y_975 = smooth(y_975)

    return y_mean, y_median, y_025, y_25, y_75, y_975

plt.figure(figsize=(12, 6))
font = {'size' : 10}
plt.rc('font', **font)

fname = "deaths_simulation.png"
label = "death_counts"
colour = 'black'
directory = "C:/tmp/"

dates_to_data = get_data(directory, label)
days, day_to_date, number_of_days = get_days(dates_to_data)
y_mean, y_median, y_025, y_25, y_75, y_975 = get_arrays(dates_to_data)

total = np.sum(y_mean)

plot_label = "Total Deaths: " + str(f'{int(total):,}')

if (y_025 is not None) and (y_975 is not None):
    plt.fill_between(range(number_of_days), y_025, y_975, color=colour, alpha=0.1)
if (y_25 is not None) and (y_75 is not None):
    plt.fill_between(range(number_of_days), y_25, y_75, color=colour, alpha=0.25)
# if (y_median is not None):
#     plt.plot(range(number_of_days), y_median, colour, linestyle='dashed', label="Median")
plt.plot(range(number_of_days), y_mean, colour, label=plot_label)

############################

historical_data_fp = "C:/Users/jt511/Documents/GitHub/pandemia/Scenarios/Homogeneous/data/historical_data/who_data.csv"

historical_daily_deaths_dict = {day_to_date[day]: 0 for day in days}
if historical_data_fp is not None:
    with open(historical_data_fp, newline='') as csvfile:
        next(csvfile)
        deaths_data = csv.reader(csvfile, delimiter=',')
        for row in deaths_data:
            date = datetime.datetime.strptime(str(row[0]), '%d/%m/%Y')
            historical_daily_deaths_dict[date] = float(row[2])
hist_daily_deaths = []
for day in range(number_of_days):
    hist_daily_deaths.append(historical_daily_deaths_dict[day_to_date[day]])
hist_daily_deaths = np.array(hist_daily_deaths).astype(np.float64)
hist_daily_deaths = smooth(hist_daily_deaths)

plt.plot(range(number_of_days), hist_daily_deaths, colour, linestyle='dashed', label="Deaths - Reported")

############################

# plt.xlabel('Date')
increment = int(min(number_of_days // 20, number_of_days))
if increment > 0:
    plt.xticks(ticks=[i*increment for i in range((number_of_days // increment))],
        labels=[day_to_date[days[i*increment]].strftime('%d/%m/%Y%Z') for i in range((number_of_days // increment))])
    plt.xticks(rotation=90)

plt.legend(loc='upper right')
plt.gca().ticklabel_format(axis='y', style='plain')

plt.gca().yaxis.set_major_formatter(plt.matplotlib.ticker.StrMethodFormatter('{x:,.0f}'))

plt.xlim([0, number_of_days])
# plt.ylim([0, ylim])
plt.gca().set_ylim(bottom=0)
plt.grid(False)
fig = plt.gcf()

# dirname = os.path.dirname(fname)
# if dirname != '':
#     os.makedirs(dirname, exist_ok=True)
# fig.savefig(fname, bbox_inches='tight')

plt.show()

# plt.close()
