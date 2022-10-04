"""Default input model, specifiying a predefined set of policy updates, on regional transmission and
to travel between regions"""

from collections import defaultdict
import logging
import numpy as np
import csv
import datetime

from pandemia.components.input_model import InputModel

log = logging.getLogger("validation_input_model")
#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class ValidationInputModel(InputModel):
    """Default model of input, specifiying a predefined set of policy interventions loaded from
    data files, affecting transmission in regions and travel between regions"""

    def __init__(self, config, scale_factor, clock, number_of_regions,
                 number_of_vaccines, age_groups):
        """Initialize component"""
        super().__init__(config, scale_factor)

        self.scale_factor = scale_factor

        self.simulation_length_days = clock.simulation_length_days
        self.epoch = clock.epoch
        self.number_of_regions = number_of_regions

        self.stay_at_home_data_fp = config['stay_at_home_data_fp']
        self.travel_data_fp       = config['travel_data_fp']

        self.max_transmission_control = config['max_transmission_control']
        self.max_travel_control       = config['max_travel_control']

        self.stay_at_home = defaultdict(dict)
        self.travel_dict = defaultdict(dict)

    def new_input(self, input_arrays):
        """Set new input"""

        self.day_to_date = {}
        for day in range(self.simulation_length_days):
            self.day_to_date[day] =\
                (self.epoch + datetime.timedelta(days=day)).strftime('%d/%m/%Y%Z')

        self.transmission_control_input =\
            np.ones((self.simulation_length_days, self.number_of_regions), dtype=float)
        self.border_closure_input =\
            np.ones((self.simulation_length_days, self.number_of_regions), dtype=float)

        with open(self.stay_at_home_data_fp, newline='') as csvfile:
            next(csvfile)
            stay_at_home_data = csv.reader(csvfile, delimiter=',')
            for row in stay_at_home_data:
                date = datetime.datetime.strptime(str(row[2]), '%Y-%m-%d').strftime('%d/%m/%Y%Z')
                iso3 = str(row[1])
                self.stay_at_home[date][iso3] = float(float(str(row[3])) / 3)

        with open(self.travel_data_fp, newline='') as csvfile:
            next(csvfile)
            travel_data = csv.reader(csvfile, delimiter=',')
            for row in travel_data:
                date = datetime.datetime.strptime(str(row[3]), '%Y-%m-%d').strftime('%d/%m/%Y%Z')
                iso3 = str(row[1])
                self.travel_dict[date][iso3] = float(int(str(row[4])) / 4)

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        for day in range(self.simulation_length_days):
            date = self.day_to_date[day]
            if date in self.stay_at_home:
                if vector_region.other_name in self.stay_at_home[date]:
                    self.transmission_control_input[day][vector_region.id] =\
                        1 - (self.max_transmission_control *\
                             self.stay_at_home[date][vector_region.other_name])
            else:
                self.transmission_control_input[day][vector_region.id] = 1.0
            if date in self.travel_dict:
                if vector_region.other_name in self.travel_dict[date]:
                    self.border_closure_input[day][vector_region.id] =\
                        1 - (self.max_travel_control *\
                             self.travel_dict[date][vector_region.other_name])
            else:
                self.border_closure_input[day][vector_region.id] = 1.0

    def initial_conditions(self, vector_region):
        """Initial input"""

        pass

    def dynamics(self, vector_region, day):
        """Changes to input"""

        vector_region.current_region_transmission_multiplier *=\
            self.transmission_control_input[day][vector_region.id]

        vector_region.current_border_closure_multiplier =\
            self.border_closure_input[day][vector_region.id]
