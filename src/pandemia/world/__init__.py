
from .region import Region, VectorRegion

class World:
    """Represents a world consiting of regions, for example a country consisting of administrative
    division, or the whole world consisting of countries.

    Parameters:
    ----------
    scale_factor : float
        Attached to each world is a scale factor, indicating that scaling may have taken place in
        the creation of the world (`float`). For example, a world with 1 million agents representing
        a real world with 100 million agents would have the scale factor 0.01. This is stored as an
        attribute since the scale factor is needed to scale quantities during input and output.
    """

    regions = None
    """The list of regions in this world (`list[Region]`).
    """

    number_of_regions = None
    """How many regions are in this world (`int`).
    """

    scale_factor = None
    """Attached to each world is a scale factor, indicating that scaling may have taken place in the
    creation of the world (`float`). For example, a world with 1 million agents representing a real
    world with 100 million agents would have the scale factor 0.01. This is stored as an attribute
    since the scale factor is needed to scale quantities during input and output.
    """

    travel_matrix = None
    """An array of integers of dimension number_of_regions x number_of_regions
    (`Union[None, np.ndarray]`).
    """

    def __init__(self, scale_factor):

        self.regions: list[Region] = []
        self.number_of_regions: int = None
        self.scale_factor: float = scale_factor
        self.travel_matrix: Union[None, np.ndarray] = None

    def vectorize_world(self):
        """Converts an object of type World to an object of type VectorWorld. A VectorWorld contains
        a list of objects of type VectorRegion, as opposed to objects of type Region, as in a World.

        Returns
        -------
        new_vector_world : `VectorWorld`
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
    """Represents a world consiting of vectorized regions.

    Parameters:
    ----------
    scale_factor : float
        Attached to each world is a scale factor, indicating that scaling may have taken place in
        the creation of the world (`float`). For example, a world with 1 million agents representing
        a real world with 100 million agents would have the scale factor 0.01. This is stored as an
        attribute since the scale factor is needed to scale quantities during input and output.
    """

    regions = None
    """The list of vectorized regions in this vectorized world (`list[VectorRegion]`).
    """

    number_of_regions = None
    """How many regions are in this world (`int`).
    """

    scale_factor = None
    """Attached to each world is a scale factor, indicating that scaling may have taken place in the
    creation of the world (`float`). For example, a world with 1 million agents representing a real
    world with 100 million agents would have the scale factor 0.01. This is stored as an attribute
    since the scale factor is needed to scale quantities during input and output.
    """

    travel_matrix = None
    """An array of integers of dimension number_of_regions x number_of_regions
    (`Union[None, np.ndarray]`).
    """

    def __init__(self, scale_factor):

        self.vector_regions: list[VectorRegion] = []
        self.number_of_regions = None
        self.scale_factor = scale_factor
        self.travel_matrix = None
