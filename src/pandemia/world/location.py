"""Represents locations with a type and coordinates in space."""

# Allows classes to return their own type, e.g. from_file below
from __future__ import annotations

import uuid

from math import sqrt
from pyproj import Transformer

LocationTuple = tuple[float, float]

class Location:

    def __init__(self, typ: str, coord: LocationTuple):
        """Represents a location.

        Parameters:
          typ (str): The type of location, as a string, for example House, Restaurant etc
          coord (tuple): 2-tuple with x, y grid coordinates
        """

        self.uuid = uuid.uuid4().hex
        self.typ = typ
        self.coord = coord

    def distance_euclidean(self, other: Location) -> float:
        """Return the distance between the two locations."""

        return sqrt(((self.coord[0]-other.coord[0])**2)\
                  + ((self.coord[1]-other.coord[1])**2))

# pylint: disable=invalid-name
def ETRS89_to_WGS84(coord: LocationTuple) -> LocationTuple:
    """Convert from ETRS89 grid format to WGS84 lat, lon format"""

    # 4326 is the EPSG identifier of WGS84, 3035 is the EPSG identifier of ETRS89
    return Transformer.from_crs('epsg:3035', 'epsg:4326').transform(coord[0], coord[1])

def WGS84_to_ETRS89(coord: LocationTuple) -> LocationTuple:
    """Convert from WGS84 lat, lon format ETRS89 to grid format"""

    # 4326 is the EPSG identifier of WGS84, 3035 is the EPSG identifier of ETRS89
    return Transformer.from_crs('epsg:4326', 'epsg:3035').transform(coord[0], coord[1])
