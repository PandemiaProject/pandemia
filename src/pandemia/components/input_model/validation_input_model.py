"""Simple input model"""

from collections import defaultdict
import logging
import numpy as np
import csv
import datetime
from dateutil.parser import parse
from pytz import country_names

from pandemia.components.input_model import InputModel

log = logging.getLogger("validation_input_model")
#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class ValidationInputModel(InputModel):
    """Simple model of input"""

    def __init__(self, config, scale_factor, clock, number_of_regions,
                 number_of_vaccines, age_groups, enable_ctypes):
        """Initialize component"""
        super().__init__(config, scale_factor)

        self.enable_ctypes = enable_ctypes

        self.scale_factor = scale_factor

        self.simulation_length_days = clock.simulation_length_days
        self.epoch = clock.epoch
        self.number_of_regions = number_of_regions

        self.containment_data_fp = config['containment_data_fp']
        self.travel_data_fp      = config['travel_data_fp']

        self.max_transmission_control = config['max_transmission_control']
        self.max_travel_control       = config['max_travel_control']

        self.containment_dict = defaultdict(dict)
        self.travel_dict = defaultdict(dict)

    def new_input(self, input_arrays):
        """Set new input"""

        self.day_to_date = {}
        for day in range(self.simulation_length_days):
            self.day_to_date[day] =\
                (self.epoch + datetime.timedelta(days=day)).strftime('%m/%d/%Y%Z')

        self.transmission_control_input =\
            np.ones((self.simulation_length_days, self.number_of_regions), dtype=float)
        self.border_closure_input =\
            np.ones((self.simulation_length_days, self.number_of_regions), dtype=float)

        with open(self.containment_data_fp, newline='') as csvfile:
            next(csvfile)
            containment_data = csv.reader(csvfile, delimiter=',')
            for row in containment_data:
                date = parse(str(row[3]), dayfirst=True).strftime('%m/%d/%Y%Z')
                iso2 = str(row[2])
                self.containment_dict[date][iso2] = float(float(str(row[4])) / 100)

        with open(self.travel_data_fp, newline='') as csvfile:
            next(csvfile)
            travel_data = csv.reader(csvfile, delimiter=',')
            for row in travel_data:
                date = parse(str(row[3]), dayfirst=True).strftime('%m/%d/%Y%Z')
                iso2 = str(row[2])
                self.travel_dict[date][iso2] = float(int(str(row[4])) / 4)

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        for day in range(self.simulation_length_days):
            date = self.day_to_date[day]
            if date in self.containment_dict:
                if vector_region.name in self.containment_dict[date]:
                    self.transmission_control_input[day][vector_region.id] =\
                        1 - (self.max_transmission_control *\
                             self.containment_dict[date][vector_region.name])
            else:
                self.transmission_control_input[day][vector_region.id] = 1.0
            if date in self.travel_dict:
                if vector_region.name in self.travel_dict[date]:
                    self.border_closure_input[day][vector_region.id] =\
                        1 - (self.max_travel_control *\
                             self.travel_dict[date][vector_region.name])
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
