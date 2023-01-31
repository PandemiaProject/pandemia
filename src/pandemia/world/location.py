
# Allows classes to return their own type
from __future__ import annotations

import uuid
from math import sqrt
from pyproj import Transformer

LocationTuple = tuple[float, float]

class Location:
    """A location, for example an area of land or a building.

    During a simulation, agents move between locations according to their daily and weekly routines.

    Parameters:
    -----------
    typ : str
        The type of location, as a string. For example "House" or "Restaurant".
    coord : tuple
        A 2-tuple of floats, representing x, y coordinates.

    Attributes:
    -----------
    uuid : str
        A universally unique identifier for this location.
    typ : str
        The type of location, as a string. For example "House" or "Restaurant".
    coord : tuple
        A 2-tuple of floats, representing x, y coordinates.
    """

    def __init__(self, typ: str, coord: LocationTuple):

        self.uuid: str = uuid.uuid4().hex
        self.typ: str = typ
        self.coord: LocationTuple = coord

    def distance_euclidean(self, other: Location) -> float:
        """Calculate the Euclidean distance between two locations.

        Parameters:
        -----------
        other : `Location`
            The other location.

        Returns:
        --------
        distance : `float`
            The Euclidean distance to the other location.
        """

        distance = sqrt(((self.coord[0]-other.coord[0])**2) + ((self.coord[1]-other.coord[1])**2))

        return distance

# pylint: disable=invalid-name
def ETRS89_to_WGS84(coord: LocationTuple) -> LocationTuple:
    """Convert from ETRS89 grid format to WGS84 lat, lon format.

    Parameters:
    -----------
    coord : `LocationTuple`
        A 2-tuple of floats.

    Returns:
    --------
    new_coord : `LocationTuple`
        If coord represents the coordinates of a location in ETRS89 format, then new_coord
        represents the coordinates of that location in WGS84 format.
    """

    # 4326 is the EPSG identifier of WGS84, 3035 is the EPSG identifier of ETRS89
    new_coord = Transformer.from_crs('epsg:3035', 'epsg:4326').transform(coord[0], coord[1])

    return new_coord

def WGS84_to_ETRS89(coord: LocationTuple) -> LocationTuple:
    """Convert from WGS84 lat, lon format to ETRS89 grid format.

    Parameters:
    -----------
    coord : `LocationTuple`
        A 2-tuple of floats.

    Returns:
    --------
    new_coord : `LocationTuple`
        If coord represents the coordinates of a location in WGS84 format, then new_coord
        represents the coordinates of that location in ETRS89 format.
    """

    # 4326 is the EPSG identifier of WGS84, 3035 is the EPSG identifier of ETRS89
    new_coord = Transformer.from_crs('epsg:4326', 'epsg:3035').transform(coord[0], coord[1])
    
    return new_coord
