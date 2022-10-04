"""Default hospitalization and death model"""

import logging
import numpy as np

from ctypes import c_void_p, c_int, c_double, cdll

from pandemia.components.hospitalization_and_death_model import HospitalizationAndDeathModel

log = logging.getLogger("default_hospitalization_and_death_model")

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class DefaultHospitalizationAndDeathModel(HospitalizationAndDeathModel):
    """Default model of agent hospitalization and death"""

    def __init__(self, config):
        """Initialize component"""
        super().__init__(config)

        lib = cdll.LoadLibrary("./src/pandemia/components/hospitalization_and_death_model/"
                                "default_hospitalization_and_death_model_functions.dll")

        self.dynamics_hospitalization_and_death = lib.dynamics_hospitalization_and_death

        self.dynamics_hospitalization_and_death.restype = None

        self.hospital_threshold     = config['hospital_threshold']
        self.hospital_location_type = config['hospital_location_type']
        self.cemetery_location_type = config['cemetery_location_type']

    def vectorize_component(self, vector_region):
        """Initializes numpy arrays associated to this component"""

        number_of_agents = vector_region.number_of_agents
        number_of_locations = vector_region.number_of_locations

        hospitals = []
        cemeteries = []
        for m in range(number_of_locations):
            if vector_region.location_typ_strings[m] == self.hospital_location_type:
                hospitals.append(m)
            if vector_region.location_typ_strings[m] == self.cemetery_location_type:
                cemeteries.append(m)

        vector_region.number_of_hospitals = len(hospitals)
        vector_region.number_of_cemeteries = len(cemeteries)
        vector_region.hospitals = np.array(hospitals).astype(int)
        vector_region.cemeteries = np.array(cemeteries).astype(int)
        vector_region.in_hospital = np.zeros((number_of_agents), dtype=int)
        vector_region.in_cemetery = np.zeros((number_of_agents), dtype=int)

    def initial_conditions(self, vector_region):
        """Initial hospitalization and death"""

        for n in range(vector_region.number_of_agents):
            vector_region.in_hospital[n] = 0
            vector_region.in_cemetery[n] = 0

    def dynamics(self, vector_region):
        """Changes agent locations"""

        self.dynamics_hospitalization_and_death(
            c_int(vector_region.number_of_agents),
            c_int(vector_region.number_of_hospitals),
            c_int(vector_region.number_of_cemeteries),
            c_int(vector_region.id),
            c_double(self.hospital_threshold),
            c_void_p(vector_region.current_location.ctypes.data),
            c_void_p(vector_region.current_region.ctypes.data),
            c_void_p(vector_region.current_disease.ctypes.data),
            c_void_p(vector_region.requesting_location_update.ctypes.data),
            c_void_p(vector_region.requested_location_update.ctypes.data),
            c_void_p(vector_region.requesting_facemask_update.ctypes.data),
            c_void_p(vector_region.requested_facemask_update.ctypes.data),
            c_void_p(vector_region.in_hospital.ctypes.data),
            c_void_p(vector_region.hospitals.ctypes.data),
            c_void_p(vector_region.in_cemetery.ctypes.data),
            c_void_p(vector_region.cemeteries.ctypes.data),
            c_void_p(vector_region.random_state.ctypes.data)
        )
