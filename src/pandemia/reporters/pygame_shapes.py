"""Reporters that output a pygame plot"""

import os
os.environ['PYGAME_HIDE_SUPPORT_PROMPT'] = "hide"

from collections import defaultdict
import matplotlib.pyplot as plt
from pandemia.reporters import Reporter
import random
import pygame
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class PygameShapes(Reporter):
    """Reporter that returns a plot at the end of the simulation. This reporter uses the shapefiles
    of each region to plot the number of people infected each day (daily prevalence) in each region.
    The style of the plot is that of a heat map, with region colour intensity corresponding to the
    prevalence for that region."""

    def __init__(self, telemetry_bus, config):
        super().__init__(telemetry_bus)

        self.font_size          = config['font_size']
        self.max_norm_prev      = config['max_norm_prev']
        self.refresh_rate       = config['refresh_rate']
        self.font_size          = config['font_size']
        self.display_width      = config['display_width']
        self.display_height     = config['display_height']
        self.points_per_polygon = config['points_per_polygon']
        self.fullscreen         = config['fullscreen']
        cmap                    = config['cmap']

        self.subscribe("vector_world.data", self.initialize)
        self.subscribe("strain_counts.update", self.update)
        self.subscribe("simulation.end", self.stop_sim)

        self.ids_to_population_size  = defaultdict(float)
        self.ids_to_multipolygons    = defaultdict(MultiPolygon)
        self.ids_to_list_of_points   = defaultdict(list)
        self.ids_to_colours          = defaultdict(str)

        self.scalar_cmap = self.get_scalar_cmap(cmap, 0, 1)

        self.number_of_strains = None
        self.number_of_regions = None

        self.scale_factor = None

        self.done = False

    def initialize(self, world, number_of_strains):
        """Intializes reporter"""

        pygame.init()

        if self.fullscreen:
            infoObject = pygame.display.Info()
            self.display_width = infoObject.current_w
            self.display_height = infoObject.current_h

        self.get_polygons(world)

        self.number_of_strains = number_of_strains
        self.number_of_regions = world.number_of_regions

        if self.fullscreen:
            self.screen = pygame.display.set_mode((0, 0), pygame.FULLSCREEN)
        else:
            self.screen = pygame.display.set_mode((self.display_width,
                                                   self.display_height + self.font_size))

        self.pygame_clock = pygame.time.Clock()

    def update(self, clock, strain_counts):
        """Draws polygonal plots of regions"""

        if not self.done:

            prevalence = {}
            infected = 0
            for r in range(self.number_of_regions):
                normalized_prevalence = 0
                for s in range(self.number_of_strains):
                    normalized_prevalence += strain_counts[(r * self.number_of_strains) + s]
                infected += normalized_prevalence
                normalized_prevalence = normalized_prevalence / self.ids_to_population_size[r]
                prevalence[r] = min(normalized_prevalence, self.max_norm_prev) / self.max_norm_prev

            # Refresh screen
            self.screen.fill("white")

            # Colour regions by prevalence
            for r in range(self.number_of_regions):
                colour = self.scalar_cmap.to_rgba(prevalence[r], bytes = True)
                if self.ids_to_list_of_points[r] is not None:
                    for points in self.ids_to_list_of_points[r]:
                        pygame.draw.polygon(self.screen, colour, points, width=0)

            # Print day and infections
            BLACK = (0, 0, 0)
            FONT = pygame.font.Font("freesansbold.ttf", self.font_size)
            surface = pygame.display.get_surface()

            text_surf = FONT.render("Infected: " + str(int(infected)), True, BLACK)
            text_rect = text_surf.get_rect(topleft=(0, surface.get_height() - (2 * self.font_size)))
            self.screen.blit(text_surf, text_rect)

            text_surf = FONT.render("Date: " + str(clock.iso8601()), True, BLACK)
            text_rect = text_surf.get_rect(topleft=(0, surface.get_height() - self.font_size))
            self.screen.blit(text_surf, text_rect)

            # Update display
            pygame.display.update()

            self.pygame_clock.tick(self.refresh_rate)

            for event in pygame.event.get():
                if event.type == pygame.MOUSEBUTTONDOWN:
                    if event.button == 1:
                        pygame.quit()
                        self.done = True
                if event.type == pygame.QUIT:
                    pygame.quit()
                    self.done = True

    def get_polygons(self, world):
        """Assigns points and polygons to each region if coordinate data exists"""

        self.scale_factor = world.scale_factor

        # Get minimum and maxiumum coordinates
        min_x_by_region = []
        min_y_by_region = []
        max_x_by_region = []
        max_y_by_region = []
        for vector_region in world.vector_regions:
            if vector_region.coordinates is not None:
                coordinates = vector_region.coordinates
                # Transform coordinates and create new polygon
                for part in coordinates:
                    points = part[0]
                    x = [i[0] for i in points[:]]
                    y = [i[1] for i in points[:]]
                    min_x_by_region.append(min(x))
                    min_y_by_region.append(min(y))
                    max_x_by_region.append(max(x))
                    max_y_by_region.append(max(y))

        if len(min_x_by_region) > 0:
            min_x = min(min_x_by_region)
        else:
            min_x = 0
        if len(min_y_by_region) > 0:
            min_y = min(min_y_by_region)
        else:
            min_y = 0
        if len(max_x_by_region) > 0:
            max_x = max(max_x_by_region)
        else:
            max_x = 0
        if len(max_y_by_region) > 0:
            max_y = max(max_y_by_region)
        else:
            max_y = 0

        # Assign polygonal shape to regions for rendering using coordinate data
        for vector_region in world.vector_regions:
            self.ids_to_population_size[vector_region.id] =\
                vector_region.number_of_agents * (1 / self.scale_factor)
            if vector_region.coordinates is not None:
                self.ids_to_colours[vector_region.id] = "white"
                coordinates = vector_region.coordinates
                # Transform coordinates and create new polygon
                list_of_points = []
                for part in coordinates:
                    points = part[0]
                    # Transform coordinates to display size
                    x = [int(self.display_width * (i[0] - min_x) / (max_x - min_x))
                         for i in points[:]]
                    y = [int(self.display_height * (1 - ((i[1] - min_y) / (max_y - min_y))))
                         for i in points[:]]
                    # Put coordinates into a list of tuples
                    new_points = list(zip(x,y))
                    # Subsample
                    new_points = [new_points[i] for i in
                                  sorted(random.sample(range(len(new_points)),
                                  min(len(new_points), self.points_per_polygon)))]
                    # Add to list of points
                    list_of_points.append(new_points)
                # Update region with shape data
                self.ids_to_multipolygons[vector_region.id] =\
                    MultiPolygon([Polygon(points) for points in list_of_points])
                self.ids_to_list_of_points[vector_region.id] = list_of_points

    def stop_sim(self):
        """Called when the simulation ends."""

        pygame.quit()

    def get_scalar_cmap(self, cmap, min_val, max_val):
        """Constucts scalar colour map to represent disease prevalence when rendering"""

        cm = plt.get_cmap(cmap)
        cNorm  = colors.Normalize(vmin=min_val, vmax=max_val)

        return cmx.ScalarMappable(norm=cNorm, cmap=cm)
