"""Default policy model, specifiying a predefined set of policy updates, on regional transmission and
to travel between regions"""

from collections import defaultdict
import logging
import pandas
import numpy as np
import csv
import datetime

from ..policy_maker_model import PolicyMakerModel

log = logging.getLogger("validation_policy_maker_model")
#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class ValidationPolicyMakerModel(PolicyMakerModel):
    """Default model of a policy maker, specifiying a predefined set of policy interventions loaded from
    data files, affecting transmission in regions and travel between regions"""

    def __init__(self, config, scale_factor, clock, number_of_regions,
                 number_of_vaccines, age_groups):
        """Initialize component"""
        super().__init__(config, scale_factor)

        self.scale_factor = scale_factor

        self.age_groups = age_groups
        self.number_of_vaccines = number_of_vaccines
        self.number_of_vaccination_age_groups = len(age_groups)

        self.simulation_length_days = clock.simulation_length_days
        self.epoch = clock.epoch
        self.number_of_regions = number_of_regions

        self.stay_at_home_data_fp = config['stay_at_home_data_fp']
        self.travel_data_fp       = config['travel_data_fp']
        self.vaccination_data_fp  = config['vaccination_data_fp']

        self.max_transmission_control = config['max_transmission_control']
        self.max_travel_control       = config['max_travel_control']

        self.stay_at_home = defaultdict(dict)
        self.travel_dict = defaultdict(dict)
        self.raw_vaccination_data = defaultdict(dict)

    def new_policy(self, policy):
        """Set new policy"""

        self.day_to_date = {}
        for day in range(self.simulation_length_days):
            self.day_to_date[day] =\
                (self.epoch + datetime.timedelta(days=day))

        self.transmission_control_policy =\
            np.ones((self.simulation_length_days, self.number_of_regions), dtype=np.float64)
        self.border_closure_policy =\
            np.ones((self.simulation_length_days, self.number_of_regions), dtype=np.float64)
        self.vaccination_policy =\
            np.full((self.simulation_length_days, self.number_of_regions,
                     self.number_of_vaccination_age_groups,
                     self.number_of_vaccines), 0, dtype=np.int64)

        with open(self.stay_at_home_data_fp, newline='') as csvfile:
            next(csvfile)
            stay_at_home_data = csv.reader(csvfile, delimiter=',')
            for row in stay_at_home_data:
                date = datetime.datetime.strptime(str(row[2]), '%Y-%m-%d')
                iso3 = str(row[1])
                self.stay_at_home[date][iso3] = float(float(str(row[3])) / 3)

        with open(self.travel_data_fp, newline='') as csvfile:
            next(csvfile)
            travel_data = csv.reader(csvfile, delimiter=',')
            for row in travel_data:
                date = datetime.datetime.strptime(str(row[3]), '%Y-%m-%d')
                iso3 = str(row[1])
                self.travel_dict[date][iso3] = float(int(str(row[4])) / 4)

        with open(self.vaccination_data_fp, newline='') as csvfile:
            next(csvfile)
            data = csv.reader(csvfile, delimiter=',')
            for row in data:
                iso3 = str(row[0])
                date = datetime.datetime.strptime(str(row[3]), '%Y-%m-%d')
                total_first_doses = 0
                if row[35] != "":
                    total_first_doses = int(float(row[35]))
                total_third_doses = 0
                if row[37] != "":
                    total_third_doses = int(float(row[37]))
                self.raw_vaccination_data[iso3][date] = [total_first_doses, total_third_doses]

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        # Transmission control and border closure
        for day in range(self.simulation_length_days):
            date = self.day_to_date[day]
            if date in self.stay_at_home:
                if vector_region.other_name in self.stay_at_home[date]:
                    self.transmission_control_policy[day][vector_region.id] =\
                        1 - (self.max_transmission_control *\
                             self.stay_at_home[date][vector_region.other_name])
            else:
                self.transmission_control_policy[day][vector_region.id] = 1.0
            if date in self.travel_dict:
                if vector_region.other_name in self.travel_dict[date]:
                    self.border_closure_policy[day][vector_region.id] =\
                        1 - (self.max_travel_control *\
                             self.travel_dict[date][vector_region.other_name])
            else:
                self.border_closure_policy[day][vector_region.id] = 1.0

        # Vaccination
        vaccination_data = defaultdict(list)
        data_dict = self.raw_vaccination_data
        iso3 = vector_region.other_name
        if iso3 not in data_dict:
            vaccination_data[self.epoch] = [0, 0]
            min_date = self.epoch
            max_date = self.epoch
        else:
            min_date = min(data_dict[iso3])
            max_date = max(data_dict[iso3])
            date = min_date
            while date <= max_date:
                if date in data_dict[iso3]:
                    yesterday = date - datetime.timedelta(days=1)
                    if yesterday in data_dict[iso3]:
                        new_first_doses =\
                            int(data_dict[iso3][date][0] - data_dict[iso3][yesterday][0])
                        new_third_doses =\
                            int(data_dict[iso3][date][1] - data_dict[iso3][yesterday][1])
                        if new_first_doses >= 0 and new_third_doses >= 0:
                            new = [new_first_doses, new_third_doses]
                        else:
                            new = [0, 0]
                    else:
                        new = [0, 0]
                    row = new
                else:
                    row = [0, 0]
                vaccination_data[date] = row
                date += datetime.timedelta(days=1)

        # Fill missing values
        start_date = self.epoch
        end_date = start_date + datetime.timedelta(days=self.simulation_length_days + 1)
        dates = pandas.date_range(min_date,max_date,freq='d').tolist()
        required_dates = pandas.date_range(start_date,end_date,freq='d').tolist()
        missing_dates = [date for date in required_dates if date not in dates]
        for date in missing_dates:
            vaccination_data[date] = [0,0]

        # Get number of agents in each age group
        num_by_age_groups = [0, 0, 0]
        for n in range(vector_region.number_of_agents):
            if vector_region.age[n] < 18:
                num_by_age_groups[0] += 1
            if (vector_region.age[n] >= 18) and (vector_region.age[n] < 65):                        # TODO do age groups properly...
                num_by_age_groups[1] += 1
            if vector_region.age[n] >= 65:
                num_by_age_groups[2] += 1

        # Construct vaccinations arrays
        date = start_date
        age_group_first_doses = 2
        age_group_third_doses = 2
        num_first_doses_by_age_group = [0,0,0]
        num_third_doses_by_age_group = [0,0,0]
        for day in range(self.simulation_length_days):
            first_doses = int(self.scale_factor * vaccination_data[date][0])
            third_doses = int(self.scale_factor * vaccination_data[date][1])
            self.vaccination_policy[day][vector_region.id][age_group_first_doses][0] = first_doses
            self.vaccination_policy[day][vector_region.id][age_group_third_doses][1] = third_doses
            num_first_doses_by_age_group[age_group_first_doses] += first_doses
            num_third_doses_by_age_group[age_group_third_doses] += third_doses
            if num_first_doses_by_age_group[age_group_first_doses] >\
                0.8 * num_by_age_groups[age_group_first_doses]:
                age_group_first_doses -= 1
                age_group_first_doses = max(age_group_first_doses, 0)
            if num_third_doses_by_age_group[age_group_third_doses] >\
                0.8 * num_by_age_groups[age_group_third_doses]:
                age_group_third_doses -= 1
                age_group_third_doses = max(age_group_third_doses, 0)
            date += datetime.timedelta(days=1)

    def initial_conditions(self, vector_region):
        """Initial policy"""

        vector_region.num_to_vaccinate =\
            self.vaccination_policy[0][vector_region.id]

        pass

    def dynamics(self, vector_region, day):
        """Changes to policy"""

        vector_region.current_region_transmission_multiplier *=\
            self.transmission_control_policy[day][vector_region.id]

        vector_region.current_border_closure_multiplier =\
            self.border_closure_policy[day][vector_region.id]

        vector_region.num_to_vaccinate =\
            self.vaccination_policy[day][vector_region.id]
