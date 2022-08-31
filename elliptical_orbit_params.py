import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from numpy import float64 as f64
from numpy import pi
from numpy import sqrt
from numpy import float_power as fpow

from ece6390_constants import Earth_mass as Me
from ece6390_constants import Earth_radius as R_earth
from ece6390_constants import Gravitational_constant as Gravity
from ece6390_constants import Sidereal_day
from ece6390_constants import Quantity

from time_lib import PCS_Time

from orbit import EllipticalOrbit

if __name__ == "__main__":
    
    apogee = Quantity(f64(20200 + R_earth()), "km")
    perigee = Quantity(f64(1000 + R_earth()), "km")

    eo = EllipticalOrbit(apoapsis=apogee, periapsis=perigee)

    period_hours = PCS_Time(_seconds=(eo.period()/2))

    print("Orbit Period: {}".format(period_hours))
    print("Orbit Eccentricity: {}".format(eo.eccentricity))
    print("Orbit Semi-major axis: {}".format(eo.a))
    print("Orbit Velocity at apogee: {}".format(eo.V_min))
    print("Orbit Velocity at perigee: {}".format(eo.V_max))

    # print("Time to do a Holmann Transfer (Period / 2): {} {}".format(period_hours/2, period_hours.unit))

    # apogee = Quantity(f64(20200), "km")
    # perigee = Quantity(f64(1000), "km")
    
    # # eccentricity
    # alpha = (perigee + R_earth) / (apogee + R_earth)
    # e = (1 - alpha) / (1 + alpha)

    # print("Eccentricity: {:.3f}".format(e))

    # # semi-major axis
    # a = Quantity(f64((apogee + perigee + (2*R_earth())) / 2), "km")
    # print("Semi-major axis: {}".format(a))

    # # period
    # T_s = (2*pi*fpow(a(), (3/2))) / (sqrt(Gravity*Me))
    # print("Period: {:.5f} seconds".format(T_s))

    # # orbital speed
    # V_theta_m_per_s = sqrt(Gravity * Me / a())
    # print("Speed: {:.5f}km/s".format(V_theta_m_per_s/1000))


    # T_squared = (4*fpow(pi, 2)*fpow(a, 3)) / (Gravity*Me)
    # T_s = sqrt((4*fpow(pi,2)*fpow(a,3))/(Gravity*Me))

    # print("Period: {:.3f} seconds".format(T_s))


