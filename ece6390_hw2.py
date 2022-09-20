import matplotlib.pyplot as plt
import numpy as np

from numpy import float64 as f64
from numpy import pi

from ece6390_constants import Atlanta_latitude as ES_LAT
from ece6390_constants import Atlanta_longitude as ES_LON
from ece6390_constants import Earth_angular_veolcity as wEarth
from ece6390_constants import Sidereal_day
from ece6390_constants import Earth_mass as Mearth
from ece6390_constants import Sun_mass as Msun
from ece6390_constants import Earth_radius as R_earth
from ece6390_constants import Quantity, Sun_radius

from orbit import CircularOrbit, NumericalOrbitSimulator

from time_lib import PCS_Time

if __name__ == "__main__":

    # calculate some initial conditions for a GEO orbit
    T = Quantity(Sidereal_day.seconds, "seconds")
    geostationary_orbit = CircularOrbit(T_seconds=T)
    print("Geostationary Orbit around Earth has a Period of {} and a Radius of {}, with a velocity of {}".format(geostationary_orbit.period, geostationary_orbit.radius, geostationary_orbit.angular_velocity))


    Pa_to_test = [3, 4, 5, 6, 7, 10, 15, 25, 50]
    azimuth_range = np.zeros(len(Pa_to_test))
    azimuth_deviation = np.zeros(len(Pa_to_test))

    for i, Pa in enumerate(Pa_to_test):

        cross_sectional_area = Quantity(10, "m^2")
        mass = Quantity(Pa*cross_sectional_area(), "kg")

        t0 = f64(0)
        tf = PCS_Time(_years=1).seconds
        delta_t_1 = f64(10)          # seconds
        n = int((tf-t0)/delta_t_1)  # number of samples
        t = np.linspace(t0, tf, n)
        OS1 = NumericalOrbitSimulator(R_init            =geostationary_orbit.radius, 
                                        theta_init      =ES_LON, 
                                        V_r_init        =Quantity(f64(0), "km/s"), 
                                        V_theta_init    =geostationary_orbit.velocity, 
                                        delta_t         =Quantity(delta_t_1, "s"),
                                        n               =n,
                                        As              =cross_sectional_area,
                                        Ms              =mass)

        OS1.compute_look_angles()

        azimuth_range[i] = np.rad2deg(OS1.look_angle_azimuth.max() - OS1.look_angle_azimuth.min())
        azimuth_deviation[i] = azimuth_range[i] / 2
        print("i: {}, Density: {}, Initial Azimuth: {:.4f}, Azimuth Deviation: {:.4f}".format(i, Pa, 
                                                                                                np.rad2deg(OS1.look_angle_azimuth[0]), 
                                                                                                azimuth_deviation[i]))

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_title("Patrick Smith, ECE6390 HW2\nAzimuth look angle deviation vs area density")
    ax.plot(Pa_to_test, azimuth_deviation, ".-")
    ax.grid(True)
    ax.set_xlabel("Area Density, kg/m^2")
    ax.set_ylabel("Azimuth Deviation, degrees")

    plt.show()

