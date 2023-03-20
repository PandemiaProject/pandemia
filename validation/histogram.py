
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from datetime import datetime
import csv

def get_dates(filename):

    df = pd.read_csv(filename)

    dates = [datetime.strptime(str(date), '%d/%m/%Y') for date in list(df['iso8601'])]

    return dates

def get_date_of_first_infection_simulated(filename, dates):

    df = pd.read_csv(filename)

    df = df.drop(columns=['day', 'Total for Strain: 0', 'Total for Strain: 1', 'Total for Strain: 2', 'Total for Strain: 3', 'Total :'])

    iso2s_to_data = defaultdict(list)
    for col in df.columns[1:]:
        iso2 = col[0:2]
        if len(iso2s_to_data[iso2]) == 0:
            iso2s_to_data[iso2] = list(df[col])
        else:
            iso2s_to_data[iso2] = [x + y for x, y in zip(iso2s_to_data[iso2], list(df[col]))]

    dates = [datetime.strptime(str(date), '%d/%m/%Y') for date in list(df['iso8601'])]

    date_of_first_infection = {}
    for iso2 in iso2s_to_data:
        for day in range(len(dates)):
            if iso2s_to_data[iso2][day] > 0:
                date_of_first_infection[iso2] = dates[day]
                break
    
    return date_of_first_infection

def get_iso_and_name_dicts():

    iso2s_to_names = {}
    iso2s_to_iso3s = {}
    all_countries = []
    with open('C:/Users/jt511/Documents/Github/Pandemia/Scenarios/Heterogeneous/data/country_data.csv', newline='') as csvfile:
        next(csvfile)
        data = csv.reader(csvfile, delimiter=',')
        for row in data:
            name = row[0]
            iso2 = row[1]
            iso3 = row[2]
            iso2s_to_names[iso2] = name
            iso2s_to_iso3s[iso2] = iso3
            all_countries.append(name)

    return all_countries, iso2s_to_names, iso2s_to_iso3s

def get_date_of_first_infection_reported():

    date_of_first_infection_reported = {}
    iso3s_to_dates = defaultdict(list)
    with open('C:/Users/jt511/Documents/Github/Pandemia/Scenarios/Validation/data/owid-covid-data.csv', newline='') as csvfile:
        next(csvfile)
        data = csv.reader(csvfile, delimiter=',')
        for row in data:
            iso3 = row[0]
            date = datetime.strptime(str(row[3]), '%Y-%m-%d')
            new_cases = str(row[5])
            if len(new_cases) > 0:
                if float(new_cases) > 0:
                    iso3s_to_dates[iso3].append(date)
    for iso3 in iso3s_to_dates:
        min_date = min(iso3s_to_dates[iso3])
        date_of_first_infection_reported[iso3] = min_date # str(min_date.strftime('%d/%m/%Y'))
    # print(date_of_first_infection_reported)

    return date_of_first_infection_reported 

def calculate_differences(hist_data, date_of_first_infection_simulated, date_of_first_infection_reported, iso2s_to_names, iso2s_to_iso3s):

    for iso2 in date_of_first_infection_simulated:
        if iso2 not in ['KP']:
            if iso2s_to_iso3s[iso2] in date_of_first_infection_reported:
                date_1 = date_of_first_infection_simulated[iso2]
                date_2 = date_of_first_infection_reported[iso2s_to_iso3s[iso2]]
                # if date_1 > date_2:
                # print(iso2s_to_names[iso2], "Simulated: " + str(date_1.strftime('%d/%m/%Y')), "Reported: " + str(date_2.strftime('%d/%m/%Y')))
                # print(iso2s_to_names[iso2], "Difference: " + str((date_2 - date_1).days))
                hist_data[iso2].append(int((date_2 - date_1).days))
    
    return hist_data

def print_dates(dates, date_of_first_infection_simulated, all_countries, iso2s_to_names):

    countries_hit = []
    for date in dates:
        countries = []
        for iso2 in date_of_first_infection_simulated:
            if date_of_first_infection_simulated[iso2] == date:
                countries.append(iso2s_to_names[iso2])
                countries_hit.append(iso2s_to_names[iso2])
        print("Date: " + str(date.strftime('%d/%m/%Y')))
        print(countries)

    countries_not_hit = []
    for country in all_countries:
        if country not in countries_hit:
            countries_not_hit.append(country)
    print("Countries without infections: ")
    print(countries_not_hit)

all_countries, iso2s_to_names, iso2s_to_iso3s = get_iso_and_name_dicts()
date_of_first_infection_reported = get_date_of_first_infection_reported()

hist_data = defaultdict(list)

filename = 'C:/tmp/strain_counts_validation_1.csv'
dates = get_dates(filename)
date_of_first_infection_simulated = get_date_of_first_infection_simulated(filename, dates)
hist_data = calculate_differences(hist_data, date_of_first_infection_simulated, date_of_first_infection_reported, iso2s_to_names, iso2s_to_iso3s)

print_dates(dates, date_of_first_infection_simulated, all_countries, iso2s_to_names)

filename = 'C:/tmp/strain_counts_validation_2.csv'
dates = get_dates(filename)
date_of_first_infection_simulated = get_date_of_first_infection_simulated(filename, dates)
hist_data = calculate_differences(hist_data, date_of_first_infection_simulated, date_of_first_infection_reported, iso2s_to_names, iso2s_to_iso3s)

# filename = 'C:/tmp/strain_counts_validation_3.csv'
# dates = get_dates(filename)
# date_of_first_infection_simulated = get_date_of_first_infection_simulated(filename, dates)
# hist_data = calculate_differences(hist_data, date_of_first_infection_simulated, date_of_first_infection_reported, iso2s_to_names, iso2s_to_iso3s)

# filename = 'C:/tmp/strain_counts_validation_4.csv'
# dates = get_dates(filename)
# date_of_first_infection_simulated = get_date_of_first_infection_simulated(filename, dates)
# hist_data = calculate_differences(hist_data, date_of_first_infection_simulated, date_of_first_infection_reported, iso2s_to_names, iso2s_to_iso3s)

# filename = 'C:/tmp/strain_counts_validation_5.csv'
# dates = get_dates(filename)
# date_of_first_infection_simulated = get_date_of_first_infection_simulated(filename, dates)
# hist_data = calculate_differences(hist_data, date_of_first_infection_simulated, date_of_first_infection_reported, iso2s_to_names, iso2s_to_iso3s)

# filename = 'C:/tmp/strain_counts_validation_6.csv'
# dates = get_dates(filename)
# date_of_first_infection_simulated = get_date_of_first_infection_simulated(filename, dates)
# hist_data = calculate_differences(hist_data, date_of_first_infection_simulated, date_of_first_infection_reported, iso2s_to_names, iso2s_to_iso3s)

# filename = 'C:/tmp/strain_counts_validation_7.csv'
# dates = get_dates(filename)
# date_of_first_infection_simulated = get_date_of_first_infection_simulated(filename, dates)
# hist_data = calculate_differences(hist_data, date_of_first_infection_simulated, date_of_first_infection_reported, iso2s_to_names, iso2s_to_iso3s)

# filename = 'C:/tmp/strain_counts_validation_8.csv'
# dates = get_dates(filename)
# date_of_first_infection_simulated = get_date_of_first_infection_simulated(filename, dates)
# hist_data = calculate_differences(hist_data, date_of_first_infection_simulated, date_of_first_infection_reported, iso2s_to_names, iso2s_to_iso3s)

# filename = 'C:/tmp/strain_counts_validation_9.csv'
# dates = get_dates(filename)
# date_of_first_infection_simulated = get_date_of_first_infection_simulated(filename, dates)
# hist_data = calculate_differences(hist_data, date_of_first_infection_simulated, date_of_first_infection_reported, iso2s_to_names, iso2s_to_iso3s)

# filename = 'C:/tmp/strain_counts_validation_10.csv'
# dates = get_dates(filename)
# date_of_first_infection_simulated = get_date_of_first_infection_simulated(filename, dates)
# hist_data = calculate_differences(hist_data, date_of_first_infection_simulated, date_of_first_infection_reported, iso2s_to_names, iso2s_to_iso3s)

averages = []
for iso2 in hist_data:
    average = sum(hist_data[iso2]) / len(hist_data[iso2])
    averages.append(average)

averages = [a for a in averages if a > -100]

print(averages)
print(sum(averages) / len(averages))

plt.hist(averages, bins=40)
plt.show()
