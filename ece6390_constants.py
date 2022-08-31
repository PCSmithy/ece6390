from typing import Optional
from dataclasses import dataclass

from numpy import float64 as f64

from time_lib import PCS_Time

@dataclass
class Quantity:
    _value: Optional[f64] = 0.0
    unit: Optional[str] = "-"

    def __call__(self) -> f64:
        return self.value

    def __str__(self) -> str:
        return "{:.3e} {}".format(self._value, self.unit)

    def __add__(self, other) -> None:
        result = None
        if isinstance(other, (f64, int)):
            result = self.value + other
        else:
            try:
                result = self.value + other.value
            except ArithmeticError as e:
                print("Error: {}, Unable to add {} and {}".format(e, self, other))

        return result

    def __sub__(self, other) -> None:
        result = None
        if isinstance(other, (f64, int)):
            result = self.value - other
        else:
            try:
                result = self.value - other.value
            except ArithmeticError as e:
                print("Error: {}, Unable to multiply {} and {}".format(e, self, other))

        return result

    def __mul__(self, other) -> None:
        result = None
        if isinstance(other, (f64, int)):
            result = self.value * other
        else:
            try:
                result = self.value * other.value
            except ArithmeticError as e:
                print("Error: {}, Unable to multiply {} and {}".format(e, self, other))

        return result

    def __truediv__(self, other) -> None:
        result = None
        if isinstance(other, (f64, int)):
            result = self.value / other
        else:
            try:
                result = self.value / other.value
            except ArithmeticError as e:
                print("Error: {}, Unable to divide {} by {}".format(e, self, other))

        return result

    @property
    def value(self) -> None:
        return self._value

    @value.setter
    def value(self, s: f64) -> None:
        self._value = s


Earth_radius = Quantity(6380, "km")
Earth_mass = Quantity(5.974e24, "kg")
Gravitational_constant = Quantity(6.672e-20, "km^3/(kg*s^2)")
Sidereal_day = PCS_Time(_hours=23, _minutes=56, _seconds=4.09)
Keplers_constant_on_Earth = Quantity(3.986004418e5, "km^3/s^2")

if __name__ == "__main__":

    print("Mp = {}, value: {}, unit: {}".format(Earth_mass, Earth_mass.value, Earth_mass.unit))
    print("R_earth = {}".format(Earth_radius))
    print("G = {}".format(Gravitational_constant))
    print("Sidereal Day = {} seconds".format(Sidereal_day.seconds))


