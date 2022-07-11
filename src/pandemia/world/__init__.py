"""Represents a world consisting of regions"""

from pandemia.world.region import Region, VectorRegion

class VectorWorld:
    """Represents a world consisting of regions using numpy arrays"""

    def __init__(self, scale_factor):

        self.vector_regions: list[VectorRegion] = []
        self.number_of_regions = None
        self.scale_factor = scale_factor
        self.travel_matrix = None

class World:
    """Represents a world consisting of regions"""

    def __init__(self, scale_factor):

        self.regions: list[Region] = []
        self.number_of_regions = None
        self.scale_factor = scale_factor
        self.travel_matrix = None

    def vectorize_world(self):

        new_vector_world = VectorWorld(self.scale_factor)

        # Vectorize regions
        for region in self.regions:
            new_vector_region =\
                region.vectorize_region()
            new_vector_world.vector_regions.append(new_vector_region)

        new_vector_world.number_of_regions = self.number_of_regions
        new_vector_world.travel_matrix = self.travel_matrix

        return new_vector_world
