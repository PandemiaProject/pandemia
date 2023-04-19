"""Reporters that output to the terminal"""

from tqdm import tqdm

from . import Reporter

#pylint: disable=unused-argument
class TimeReporter(Reporter):
    """Uses TQDM to plot a progress bar"""

    def __init__(self, telemetry_bus, config):

        super().__init__(telemetry_bus)

        self.subscribe('sim.time', self.event)

        self.tqdm = None

    def event(self, clock):
        """Plots progress bar with tick and clock time"""

        if self.tqdm is None:
            self.tqdm = tqdm(total=clock.simulation_length_days - 1)

        self.tqdm.n = clock.day
        self.tqdm.set_description(f"{clock.iso8601()}")

        if clock.day == clock.simulation_length_days - 1:
            self.tqdm.close()
