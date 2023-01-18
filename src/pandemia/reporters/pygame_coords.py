"""Reporters that output a pygame plot"""

import os
os.environ['PYGAME_HIDE_SUPPORT_PROMPT'] = "hide"

import numpy as np
import matplotlib.pyplot as plt
from . import Reporter
import pygame
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

#pylint: disable=unused-argument
#pylint: disable=attribute-defined-outside-init
class PygameCoords(Reporter):
    """Reporter that returns a plot at the end of the simulation. This reporter uses the coordinates
    of each location to plot the number of people infected each day (daily prevalence) within the
    geographical area corresponding to each pixel in a rectangular plot of the regions. The style of
    the plot is that of a heat map, with pixel colour intensity corresponding to the prevalence for
    that pixel."""

    def __init__(self, telemetry_bus, config):
        super().__init__(telemetry_bus)

        self.font_size            = config['font_size']
        self.plot_prevalence      = config['plot_prevalence']
        self.max_prev_to_plot     = config['max_prev_to_plot']
        self.plot_num_infected    = config['plot_num_infected']
        self.max_num_inf_to_plot  = config['max_num_inf_to_plot']
        self.refresh_rate         = config['refresh_rate']
        self.font_size            = config['font_size']
        self.display_width        = config['display_width']
        self.display_height       = config['display_height']
        self.cellsize_config  = config['cellsize']
        self.fullscreen           = config['fullscreen']
        cmap                      = config['cmap']

        self.subscribe("vector_world.data", self.initialize)
        self.subscribe("data.update", self.update)
        self.subscribe("simulation.end", self.stop_sim)

        self.scalar_cmap = self.get_scalar_cmap(cmap, 0, 1)

        self.scale_factor = None

        self.render_dict = {}
        self.rectangle_width = None
        self.rectangle_height = None

        self.done = False

    def initialize(self, world, number_of_strains):
        """Intializes reporter"""

        pygame.init()

        if self.fullscreen:
            infoObject = pygame.display.Info()
            self.display_width = infoObject.current_w
            self.display_height = infoObject.current_h

        self.scale_factor = world.scale_factor

        # Simplify colour map for faster rendering
        self.prevalence_to_colour = {}
        for prevalence in range(11):
            self.prevalence_to_colour[prevalence] =\
                self.scalar_cmap.to_rgba((prevalence + 1) / 11, bytes = True)

        # Store number of agents in each region
        self.number_of_agents =\
            {vr.id: vr.number_of_agents for vr in world.vector_regions}

        # Determine minimum and maximum spatial coordinates that need to be plotted
        min_x = min([np.min(vr.location_x_coords) for vr in world.vector_regions])
        min_y = min([np.min(vr.location_y_coords) for vr in world.vector_regions])
        max_x = max([np.max(vr.location_x_coords) for vr in world.vector_regions])
        max_y = max([np.max(vr.location_y_coords) for vr in world.vector_regions])

        spatial_width = max_x - min_x
        spatial_height = max_y - min_y

        # Preserving aspect ratio, determine width and height of area to plot
        r_old = spatial_width / spatial_height
        r_new = self.display_width / self.display_height

        if (r_old > r_new):
            self.width = int(self.display_width)
            self.height = int(self.width / r_old)
        else:
            self.height = int(self.display_height)
            self.width = int(self.height * r_old)

        if self.cellsize_config is not None:
            self.cellsize = int(self.width / (spatial_width / self.cellsize_config)) + 1
            self.square_height = int(self.height / (spatial_height / self.cellsize_config)) + 1
        else:
            self.cellsize = 1
            self.square_height = 1

        # Determine the mapping of spatial coordinates to display coordinates on screen
        self.locations_to_pixels_by_region = {}
        for vector_region in world.vector_regions:
            self.locations_to_pixels_by_region[vector_region.id] =\
                np.zeros((vector_region.number_of_locations, ), dtype=np.int64)
        max_x_display = 0
        max_y_display = 0
        self.region_coordinates = []
        for vector_region in world.vector_regions:
            for m in range(vector_region.number_of_locations):
                x_spatial = vector_region.location_x_coords[m]
                y_spatial = vector_region.location_y_coords[m]
                x_display =\
                    int((self.width - 1) * (x_spatial - min_x) / (max_x - min_x))
                y_display =\
                    int((self.height - 1) * (1 - (y_spatial - min_y) / (max_y - min_y)))
                self.locations_to_pixels_by_region[vector_region.id][m] =\
                    (x_display * self.height) + y_display
                if x_display > max_x_display:
                    max_x_display = x_display
                if y_display > max_y_display:
                    max_y_display = y_display
            if vector_region.region_coordinates is not None:
                for coord in vector_region.region_coordinates:
                    x_spatial = coord[0]
                    y_spatial = coord[1]
                    x_display =\
                        int((self.width - 1) * (x_spatial - min_x) / (max_x - min_x))
                    y_display =\
                        int((self.height - 1) * (1 - (y_spatial - min_y) / (max_y - min_y)))
                    self.region_coordinates.append((x_display, y_display))
                    if x_display > max_x_display:
                        max_x_display = x_display
                    if y_display > max_y_display:
                        max_y_display = y_display

        # Calculate offsets to allow the image to be displayed in the center of the screen
        self.x_offset = int((self.display_width - max_x_display) / 2)
        self.y_offset = int((self.display_height - max_y_display) / 2)

        if self.fullscreen:
            self.screen = pygame.display.set_mode((0, 0), pygame.FULLSCREEN)
        else:
            self.screen = pygame.display.set_mode((self.display_width,
                                                   self.display_height + self.font_size))

        # Create background surface
        self.screen.fill("white")
        if len(self.region_coordinates) > 0:
            for coord in self.region_coordinates:
                colour = self.prevalence_to_colour[0]
                pygame.draw.rect(self.screen, colour, pygame.Rect(coord[0] + self.x_offset,
                                 coord[1] + self.y_offset, self.cellsize,
                                 self.square_height))
        self.background_screen = self.screen.copy()

        self.pygame_clock = pygame.time.Clock()

    def update(self, clock, current_location, current_infectiousness):
        """Draws polygonal plots of regions"""

        if not self.done:

            num_infected = np.zeros((self.width * self.height), dtype=np.float64)
            num_total    = np.zeros((self.width * self.height), dtype=np.float64)

            infected = 0
            for id in current_location:
                for n in range(self.number_of_agents[id]):
                    m = current_location[id][n]
                    pixel = self.locations_to_pixels_by_region[id][m]
                    num_total[pixel] += 1
                    if current_infectiousness[id][n] > 0:
                        num_infected[pixel] += 1
                        infected += 1

            empty_pixels = num_total == 0
            num_total += empty_pixels.astype(int)
            if self.plot_prevalence and not self.plot_num_infected:
                prevalence = np.minimum(np.divide(num_infected, num_total),
                                        self.max_prev_to_plot) / self.max_prev_to_plot
            if self.plot_num_infected and not self.plot_prevalence:
                prevalence = np.minimum(num_infected, self.max_num_inf_to_plot) /\
                                        self.max_num_inf_to_plot
            prevalence = (prevalence * 10).astype(int)

            # Refresh screen
            self.screen.blit(self.background_screen, (0,0))

            # Colour squares
            for x in range(self.width):
                for y in range(self.height):
                    if not empty_pixels[(x * self.height) + y]:
                        colour = self.prevalence_to_colour[prevalence[(x * self.height) + y]]
                        pygame.draw.rect(self.screen, colour, pygame.Rect(x + self.x_offset,
                                         y + self.y_offset, self.cellsize, self.square_height))

            infected = int((1 / self.scale_factor) * infected)

            # Print day and infections
            FONT = pygame.font.Font("freesansbold.ttf", self.font_size)
            surface = pygame.display.get_surface()

            text_surf = FONT.render("Infected: " + str(int(infected)), True, "black")
            text_rect = text_surf.get_rect(topleft=(0, surface.get_height() - (2 * self.font_size)))
            self.screen.blit(text_surf, text_rect)

            text_surf = FONT.render("Date: " + str(clock.iso8601()), True, "black")
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

    def stop_sim(self):
        """Called when the simulation ends."""

        pygame.quit()

    def get_scalar_cmap(self, cmap, min_val, max_val):
        """Constucts scalar colour map to represent disease prevalence when rendering"""

        cm = plt.get_cmap(cmap)
        cNorm  = colors.Normalize(vmin=min_val, vmax=max_val)

        return cmx.ScalarMappable(norm=cNorm, cmap=cm)
