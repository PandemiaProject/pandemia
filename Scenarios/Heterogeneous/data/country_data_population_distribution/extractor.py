"""Extracts national population density grid data from the GPWv4 data set:

Center for International Earth Science Information Network - CIESIN - Columbia University. 2018.
Gridded Population of the World, Version 4 (GPWv4): Population Density, Revision 11. Palisades,
New York: NASA Socioeconomic Data and Applications Center (SEDAC). https://doi.org/10.7927/H49C6VHW.
Accessed 31 OCTOBER 2022.

"""

import numpy as np
from tqdm import tqdm

# Load iso3 to value lookup
iso3_to_value_dict = {}
data = open("raw/gpw_v4_national_identifier_grid_rev11_lookup.txt", "r")
data.readline()
for row in data:
    split_row = row.split("\t")
    value = int(split_row[0])
    iso3 = str(split_row[1])
    iso3_to_value_dict[iso3] = value

# Load population density grid and national identifier grid
national_identifier_grid =\
    np.genfromtxt("raw/gpw_v4_national_identifier_grid_rev11_2pt5_min.asc", skip_header=6)
population_density_grid =\
    np.genfromtxt("raw/gpw_v4_population_density_rev11_2020_2pt5_min.asc", skip_header=6)
ncols = 8640
nrows = 4320
xllcorner = -180
yllcorner = -90
cellsize = 0.041666666666667
NODATA_value = -9999

def create_asc(iso3, values):
    """Extracts a density map for the country of the given iso3 and value"""

    density_grid = np.zeros_like(population_density_grid)
    id_grid = np.zeros_like(national_identifier_grid)

    # Zero grid outside the desired country or countries
    for value in values:

        new_density_grid = np.copy(population_density_grid)
        new_density_grid[national_identifier_grid != value] = 0
        new_density_grid[new_density_grid == NODATA_value] = 0
        density_grid = density_grid + new_density_grid

        new_id_grid = np.copy(national_identifier_grid)
        new_id_grid[national_identifier_grid != value] = 0
        new_id_grid[national_identifier_grid == value] = 1
        id_grid = id_grid + new_id_grid

    nonzero = np.nonzero(id_grid)

    if nonzero[0].size != 0:

        # Extract subgrid tightly containing the desired country
        x_min = np.min(nonzero[1])
        x_max = np.max(nonzero[1])
        y_min = np.min(nonzero[0])
        y_max = np.max(nonzero[0])
        width = x_max - x_min
        height = y_max - y_min
        new_ncols = width + 1
        new_nrows = height + 1
        density_grid = density_grid[y_min: y_min + new_nrows, x_min: x_min + new_ncols]
        id_grid = id_grid[y_min: y_min + new_nrows, x_min: x_min + new_ncols]

        # Pad arrays so that lower left corner has integer coordinates
        error_x = (xllcorner + (x_min * cellsize)) % 1
        error_y = (yllcorner + ((nrows - (y_max + 1)) * cellsize)) % 1
        pad_x = int(error_x / cellsize)
        pad_y = int(error_y / cellsize)
        density_grid = np.pad(density_grid, ((0, pad_y), (pad_x, 0)), 'constant',
                              constant_values=((0, 0), (0, 0)))
        id_grid = np.pad(id_grid, ((0, pad_y), (pad_x, 0)), 'constant',
                         constant_values=((0, 0), (0, 0)))
        id_grid = id_grid.astype(int)

        # Record the coordinates of the lower left corner
        grid_xllcorner = int(xllcorner + ((x_min - pad_x) * cellsize))
        grid_yllcorner = int(yllcorner + ((nrows - (y_max + pad_y + 1)) * cellsize))

        # Save the density grid as an asc file with required header data
        fmt = " ".join(["%1.9g"] * (density_grid.shape[1]))
        with open('processed/density_grids/' + iso3.lower() + '.asc', 'wb') as f:
            f.write(b'ncols ' + str(new_ncols + pad_x).encode() + b'\n')
            f.write(b'nrows ' + str(new_nrows + pad_y).encode() + b'\n')
            f.write(b'xllcorner ' + str(grid_xllcorner).encode() + b'\n')
            f.write(b'yllcorner ' + str(grid_yllcorner).encode() + b'\n')
            f.write(b'cellsize ' + str(cellsize).encode() + b'\n')
            f.write(b'NODATA_value ' + str(NODATA_value).encode() + b'\n')
            np.savetxt(f, density_grid, fmt=fmt, delimiter=" ")

        # Save the ids grid as an asc file with required header data
        fmt = " ".join(["%d"] * (id_grid.shape[1]))
        with open('processed/id_grids/' + iso3.lower() + '.asc', 'wb') as f:
            f.write(b'ncols ' + str(new_ncols + pad_x).encode() + b'\n')
            f.write(b'nrows ' + str(new_nrows + pad_y).encode() + b'\n')
            f.write(b'xllcorner ' + str(grid_xllcorner).encode() + b'\n')
            f.write(b'yllcorner ' + str(grid_yllcorner).encode() + b'\n')
            f.write(b'cellsize ' + str(cellsize).encode() + b'\n')
            f.write(b'NODATA_value ' + str(NODATA_value).encode() + b'\n')
            np.savetxt(f, id_grid, fmt=fmt, delimiter=" ")

values = [iso3_to_value_dict['SRB'], iso3_to_value_dict['MNE']]
create_asc('SCG', values)

for iso3 in tqdm(iso3_to_value_dict):
    create_asc(iso3, [iso3_to_value_dict[iso3]])
