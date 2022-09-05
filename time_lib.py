from typing import Optional
from dataclasses import dataclass

import numpy as np

from numpy import float64 as f64

@dataclass
class PCS_Time:
    """a time library to help me convert between different time bases and display time in a more meaninful way"""
    _years: Optional[f64] = 0
    _weeks: Optional[f64] = 0
    _days: Optional[f64] = 0
    _hours: Optional[f64] = 0
    _minutes: Optional[f64] = 0
    _seconds: Optional[f64] = 0

    def __post_init__(self):
        # compute how many total seconds we have from all contributions
        seconds_per_minute = 60
        seconds_per_hour = np.multiply(seconds_per_minute, 60)
        seconds_per_day = np.multiply(seconds_per_hour, 24)
        seconds_per_week = np.multiply(seconds_per_day, 7)
        seconds_per_year = 31556952

        seconds_from_years = self._years * seconds_per_year
        seconds_from_weeks = self._weeks * seconds_per_week
        seconds_from_days = self._days * seconds_per_day
        seconds_from_hours = self._hours * seconds_per_hour
        seconds_from_minutes = self._minutes * seconds_per_minute
        
        self.seconds = (self._seconds + seconds_from_minutes + seconds_from_hours + seconds_from_days + seconds_from_weeks + seconds_from_years)

        # update each field to represent to total time
        self.minutes = np.divide(self.seconds, seconds_per_minute)
        self.hours = np.divide(self.seconds, seconds_per_hour)
        self.days = np.divide(self.seconds, seconds_per_day)
        self.weeks = np.divide(self.seconds, seconds_per_week)
        self.years = np.divide(self.seconds, seconds_per_year)

    def __str__(self) -> str:
        if self._years > 1:
            time_left = self.years
            (years, time_left) = divmod(time_left, 1)
            time_left *= 52 # weeks/year
            (weeks, time_left) = divmod(time_left, 1)
            time_left *= 7 # days/week
            (days, time_left) = divmod(time_left, 1)
            time_left *= 24 # hours/day
            (hours, time_left) = divmod(time_left, 1)
            time_left *= 60 # minutes/hour
            (minutes, time_left) = divmod(time_left, 1)
            time_left *= 60 # seconds/minute
            (seconds, time_left) = divmod(time_left, 1)
            return "{:.0f} years, {:.0f} weeks, {:.0f} days, {:.0f} hours, {:.0f} minutes, {:.0f} seconds".format(years, weeks, days, hours, minutes, seconds)
        elif self._weeks > 1:
            time_left = self.weeks
            (weeks, time_left) = divmod(time_left, 1)
            time_left *= 7 # days/week
            (days, time_left) = divmod(time_left, 1)
            time_left *= 24 # hours/day
            (hours, time_left) = divmod(time_left, 1)
            time_left *= 60 # minutes/hour
            (minutes, time_left) = divmod(time_left, 1)
            time_left *= 60 # seconds/minute
            (seconds, time_left) = divmod(time_left, 1)
            return "{:.0f} weeks, {:.0f} days, {:.0f} hours, {:.0f} minutes, {:.0f} seconds".format(weeks, days, hours, minutes, seconds)
        elif self._days > 1:
            time_left = self.days
            (days, time_left) = divmod(time_left, 1)
            time_left *= 24 # hours/day
            (hours, time_left) = divmod(time_left, 1)
            time_left *= 60 # minutes/hour
            (minutes, time_left) = divmod(time_left, 1)
            time_left *= 60 # seconds/minute
            (seconds, time_left) = divmod(time_left, 1)
            return "{:.0f} days, {:.0f} hours, {:.0f} minutes, {:.0f} seconds".format(days, hours, minutes, seconds)
        elif self._hours > 1:
            time_left = self.hours
            (hours, time_left) = divmod(time_left, 1)
            time_left *= 60 # minutes/hour
            (minutes, time_left) = divmod(time_left, 1)
            time_left *= 60 # seconds/minute
            (seconds, time_left) = divmod(time_left, 1)
            return "{:.0f} hours, {:.0f} minutes, {:.0f} seconds".format(hours, minutes, seconds)
        elif self._minutes > 1:
            time_left = self.minutes
            (minutes, time_left) = divmod(time_left, 1)
            time_left *= 60 # seconds/minute
            (seconds, time_left) = divmod(time_left, 1)
            return "{:.0f} minutes, {:.0f} seconds".format(minutes, seconds)
        else:
            return "{:.0f} seconds".format(self.seconds)


    @property
    def seconds(self) ->f64:
        return self._seconds
    
    @seconds.setter
    def seconds(self, s: f64) -> None:
        self._seconds = s

    @property
    def minutes(self) ->f64:
        return self._minutes
    
    @minutes.setter
    def minutes(self, m: f64) -> None:
        self._minutes = m

    @property
    def hours(self) ->f64:
        return self._hours
    
    @hours.setter
    def hours(self, h: f64) -> None:
        self._hours = h

    @property
    def days(self) ->f64:
        return self._days
    
    @days.setter
    def days(self, d: f64) -> None:
        self._days = d

    @property
    def weeks(self) ->f64:
        return self._weeks
    
    @weeks.setter
    def weeks(self, w: f64) -> None:
        self._weeks = w

    @property
    def years(self) ->f64:
        return self._years
    
    @years.setter
    def years(self, y: f64) -> None:
        self._years = y


if __name__ == "__main__":

    some_time = PCS_Time(_minutes=10, _seconds=30)

    print("Some time is {:.3f} seconds".format(some_time.seconds))
    print("Which is {:.3f} minutes".format(some_time.minutes))
    print("Which is {:.6f} hours".format(some_time.hours))
    print("Which is {:.6f} days".format(some_time.days))
    print("Which is {:.6f} weeks".format(some_time.weeks))
    print("Which is {:.6f} years".format(some_time.years))


