
import numpy as np
from typing import Union
import pickle

from .region import Region, VectorRegion

class World:
    """A world consiting of regions.
    
    For example, a country consisting of administrative divisions, or the whole world consisting of
    countries.

    Parameters:
    -----------
    scale_factor : float
        Attached to each world is a scale factor, indicating that scaling may have taken place in
        the creation of the world. For example, a world with 1 million agents representing a real
        world with 100 million agents would have the scale factor 0.01. This is stored as an
        attribute since the scale factor is needed to scale quantities during input and output.

    Attributes:
    -----------
    regions : list[Region]
        The list of regions in this world.
    number_of_regions : int
        How many regions are in this world.
    scale_factor : float
        Attached to each world is a scale factor, indicating that scaling may have taken place in
        the creation of the world. For example, a world with 1 million agents representing
        a real world with 100 million agents would have the scale factor 0.01. This is stored as an
        attribute since the scale factor is needed to scale quantities during input and output.
    travel_matrix : np.ndarray
        An array of integers of dimension number_of_regions x number_of_regions.
    """

    def __init__(self, scale_factor):

        self.regions: list[Region] = []
        self.number_of_regions: int = None
        self.scale_factor: float = scale_factor
        self.travel_matrix: Union[None, np.ndarray] = None

    def vectorize_world(self):
        """Converts an object of type World to an object of type VectorWorld. A VectorWorld contains
        a list of objects of type VectorRegion, as opposed to objects of type Region, as in a World.

        Returns:
        --------
        new_vector_world : VectorWorld
            A vector representation of the world.
        """

        new_vector_world = VectorWorld(self.scale_factor)

        # Vectorize regions
        for region in self.regions:
            new_vector_region =\
                region.vectorize_region()
            new_vector_world.vector_regions.append(new_vector_region)

        new_vector_world.number_of_regions = self.number_of_regions
        new_vector_world.travel_matrix = self.travel_matrix

        return new_vector_world

class VectorWorld:
    """A world consiting of vectorized regions.
    
    For example, a country consisting of administrative divisions, or the whole world consisting of
    countries.

    Parameters:
    -----------
    scale_factor : float
        Attached to each world is a scale factor, indicating that scaling may have taken place in
        the creation of the world. For example, a world with 1 million agents representing a real
        world with 100 million agents would have the scale factor 0.01. This is stored as an
        attribute since the scale factor is needed to scale quantities during input and output.

    Attributes:
    -----------
    vector_regions : list[VectorRegion]
        The list of regions in this vectorized world.
    number_of_regions : int
        How many vectorized regions are in this vectorized world.
    scale_factor : float
        Attached to each world is a scale factor, indicating that scaling may have taken place in
        the creation of the world. For example, a world with 1 million agents representing
        a real world with 100 million agents would have the scale factor 0.01. This is stored as an
        attribute since the scale factor is needed to scale quantities during input and output.
    travel_matrix : np.ndarray
        An array of integers of dimension number_of_regions x number_of_regions.
    """

    def __init__(self, scale_factor):

        self.vector_regions: list[VectorRegion] = []
        self.number_of_regions = None
        self.scale_factor = scale_factor
        self.travel_matrix = None

    def to_file(self, output_filename: str):
        """Write an object to disk at the filename given.

        Parameters:
            output_filename (str):The filename to write to.  Files get overwritten
                                  by default.

        Returns:
            None
        """

        with open(output_filename, 'wb') as fout:
            pickle.dump(self, fout, protocol=pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def from_file(input_filename: str):
        """Read an object from disk from the filename given.

        Parameters:
            input_filename (str):The filename to read from.

        Returns:
            obj(Object):The python object read from disk
        """

        with open(input_filename, 'rb') as fin:
            payload = pickle.load(fin)

        return payload
