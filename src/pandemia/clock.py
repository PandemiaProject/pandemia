"""Used to represent time in the model"""

import logging
from datetime import datetime, timedelta
from typing import Union

from dateutil.parser import parse

log = logging.getLogger("clock")

class Clock:
    """Tracks time"""
    # They're appropriate in this case.
    # pylint: disable=too-many-instance-attributes

    def __init__(self, tick_length_s: int, simulation_length_days: int,
                 epoch: Union[datetime, str]=datetime.now()):
        """Create a new clock.

        The clock is an iterator that counts forward in time by one day, ignoring timezone changes
        and other complexities. Each day is subdivided into ticks.

        Parameters:
            tick_length_s (int):Number of wall-clock seconds per simulation tick
            epoch (str or datetime):A datetime representing the starting point for the simulation
            simulation_length_days (int):Number of days the simulation will run for
        """

        # Check that day length (in seconds) is divisible by tick length
        if 86400 % tick_length_s != 0:
            raise ValueError("Day length must be divisible by tick length")

        if isinstance(epoch, str):
            parsed_epoch = parse(epoch)
            if parsed_epoch is None:
                raise ValueError(f"Failed to parse epoch: {epoch}")
            self.epoch = parsed_epoch
        else:
            self.epoch = epoch

        self.simulation_length_days = simulation_length_days

        self.tick_length_s   = tick_length_s
        self.ticks_in_second = 1          / self.tick_length_s
        self.ticks_in_minute = 60         / self.tick_length_s
        self.ticks_in_hour   = 3600       / self.tick_length_s
        self.ticks_in_day    = int(86400  / self.tick_length_s)
        self.ticks_in_week   = int(604800 / self.tick_length_s)

        self.epoch_week_offset = int(self.epoch.weekday() * self.ticks_in_day \
                                 + self.epoch.hour        * self.ticks_in_hour \
                                 + self.epoch.minute      * self.ticks_in_minute \
                                 + self.epoch.second      * self.ticks_in_second)

        self.day = 0
        self.started = False
        self.reset()

        log.debug("New clock created starting %s, tick_length=%i, simulation_days=%i",
                  self.epoch, tick_length_s, simulation_length_days)

    def reset(self) -> None:
        """Reset the clock to the start"""
        log.debug("Resetting clock at day=%i", self.day)
        self.day = 0
        self.started = False

    def __iter__(self):
        self.reset()
        return self

    def __next__(self):

        if self.started:
            self.day += 1
        else:
            self.started = True

        if self.day >= self.simulation_length_days:
            raise StopIteration()

        return self.day

    def __len__(self):
        return self.simulation_length_days

    def iso8601(self) -> str:
        """Return ISO 8601 time as a string"""

        return self.now().strftime('%d/%m/%Y%Z')

    def ticks_through_week(self) -> int:
        """Returns the number of whole ticks through the week this is"""
        return int((self.epoch_week_offset + (self.day * self.ticks_in_day)) % self.ticks_in_week)

    def now(self) -> datetime:
        """Return a datetime.datetime showing the clock time"""
        return self.epoch + timedelta(days=self.day)
