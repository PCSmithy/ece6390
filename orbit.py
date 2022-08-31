from typing import Optional

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from numpy import float64 as f64
from numpy import pi
from numpy import sqrt
from numpy import float_power as fpow

from ece6390_constants import Earth_mass as Me
from ece6390_constants import Gravitational_constant as Gravity
from ece6390_constants import Sidereal_day
from ece6390_constants import Quantity

from time_lib import PCS_Time

class EllipticalOrbit(object):
    T = Quantity(unit="seconds")
    R = Quantity(unit="km")
    G = Quantity()
    M = Quantity()
    V = Quantity(unit="km/s")
    e = Quantity(unit="")
    a = Quantity(unit="km")
    b = Quantity(unit="km")

    def __init__(self, Grav=Gravity, Mp=Me, T_seconds: Optional[Quantity] = None, apoapsis: Optional[Quantity] = None, periapsis: Optional[Quantity] = None) -> None:
        if T_seconds is None and apoapsis is None and periapsis is None:
            raise AttributeError("Cannot compute orbit parameters without either the Period T, or the apoapsis / periapsis")

        self.G = Grav
        self.M = Mp

        if T_seconds is None:
            alpha = periapsis / apoapsis
            self.e.value = (1 - alpha)/(1 + alpha)

            if periapsis.unit not in ["m", "km"] or apoapsis.unit not in ["m", "km"]:
                raise AttributeError("Incorrect Radius R units")
            
            if apoapsis.unit == "m":
                apoapsis.value /= 1000
                apoapsis.unit = "km"

            if periapsis.unit == "m":
                periapsis.value /= 1000
                periapsis.unit = "km"
                
            self.a.value = (apoapsis + periapsis) / 2
            self.T.value = ( 2*pi*fpow(self.a(), (3/2)) ) / sqrt(self.G * self.M)

        else: # R is none
            if T_seconds.unit != self.T.unit:
                raise AttributeError("Incorrect Period T units")
            else:
                self.T = T_seconds
                self.R.value = fpow( (self.G*self.M*fpow(self.T(), 2)) / (4*fpow(pi,2)), (1/3)) 

    @property
    def semi_major_axis(self):
        return self.a

    @property
    def eccentricity(self):
        return self.e

    @property
    def period(self):
        return self.T

    @property
    def radius(self):
        return self.R

    @property
    def V_max(self):
        return Quantity(f64(sqrt(self.G*self.M / (self.a*(1-self.e())))), "km/s")
    
    @property
    def V_min(self):
        return Quantity(f64(sqrt(self.G*self.M / (self.a*(1+self.e())))), "km/s")

class CircularOrbit(EllipticalOrbit):
    def __init__(self, Grav=Gravity, Mp=Me, T_seconds: Optional[Quantity] = None, R: Optional[Quantity] = None) -> None:
        if T_seconds is None and R is None:
            raise AttributeError("Cannot compute circular orbit without either the Period T, or the Radius R")
        super().__init__(Grav, Mp, T_seconds, apoapsis=R, periapsis=R)

        self.V.value = sqrt(self.G*self.M / self.R())

    @property
    def angular_velocity(self):
        return Quantity(self.V / self.R, "rad/s")

    @property
    def velocity(self):
        return self.V



if __name__ == "__main__":

    T = Quantity(Sidereal_day.seconds, "seconds")
    geostationary_orbit = CircularOrbit(T_seconds=T)
    print("Geostationary Orbit around Earth has a Period of {} and a Radius of {}, with a velocity of {}".format(geostationary_orbit.period, geostationary_orbit.radius, geostationary_orbit.angular_velocity))

    # R2 = Quantity(4.216417e4, "km")
    R2 = Quantity(20e3, "km")
    geostat_orbit = CircularOrbit(R=R2)
    print("Geostationary Orbit2 around Earth has a Period of {} and a Radius of {}, with a velocity of {}".format(PCS_Time(_seconds=geostat_orbit.period()), geostat_orbit.radius, geostat_orbit.velocity))


