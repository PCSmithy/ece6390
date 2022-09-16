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

from orbit import CircularOrbit, OrbitSimulatorWithLookAngles

from time_lib import PCS_Time

if __name__ == "__main__":

    # calculate some initial conditions for a GEO orbit
    T = Quantity(Sidereal_day.seconds, "seconds")
    geostationary_orbit = CircularOrbit(T_seconds=T)
    print("Geostationary Orbit around Earth has a Period of {} and a Radius of {}, with a velocity of {}".format(geostationary_orbit.period, geostationary_orbit.radius, geostationary_orbit.angular_velocity))


    t0 = f64(0)
    tf = PCS_Time(_weeks=4).seconds
    delta_t_1 = f64(5)          # seconds
    n = int((tf-t0)/delta_t_1)  # number of samples
    t = np.linspace(t0, tf, n)
    OS1 = OrbitSimulatorWithLookAngles(R_init        =geostationary_orbit.radius, 
                                        theta_init    =ES_LON, 
                                        V_r_init      =Quantity(f64(0), "km/s"), 
                                        V_theta_init  =Quantity(geostationary_orbit.velocity()*1.0, "km/s"), 
                                        delta_t       =Quantity(delta_t_1, "s"),
                                        n             =n,
                                        ssp_lat_init  =Quantity(0, "deg"),
                                        es_lat        =ES_LAT,
                                        es_lon        =ES_LON)



    # plots for part 1-3
    earth_r = np.full((360,1), R_earth())
    earth_theta = np.linspace(0, 2*pi, 360)

    fig = plt.figure()
    ax = fig.add_subplot(projection="polar", facecolor="lightgoldenrodyellow")
    ax.set_title("Patrick Smith, ECE6390 HW2\nTimestep: {}".format(PCS_Time(_seconds=delta_t_1)))

    ax.plot(earth_theta, earth_r, lw=1, color="tab:green", label='Earth, Radius 6390 km')

    ax.plot(OS1.theta, OS1.r, lw=2, color="tab:orange", label='GEO orbit, Eccentricity: {:.3f}, Period: {}'.format(OS1.eccentricity, OS1.period))
    ax.tick_params(grid_color="palegoldenrod")

    ax.set_rmax(150e3)
    ax.set_rlabel_position(40)  # Move radial labels away from plotted line
    angle = np.deg2rad(67.5)
    ax.legend(loc="lower left",
            bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2))
    ax.grid(True)


    fig2 = plt.figure()
    ax2 = fig2.add_subplot()
    ax2.set_title("Patrick Smith, ECE6390 HW1\nTimestep: {}".format(PCS_Time(_seconds=delta_t_1)))

        
    ax2.plot(t, np.rad2deg(OS1.rotation_of_earth_in_fixed_frame), lw=2, color="tab:blue", label='rotation of earth, '.format())
    ax2.plot(t, np.rad2deg(OS1.theta), lw=2, color="tab:green", label='rotation of sat, '.format())
    ax2.plot(t, np.rad2deg(OS1.ssp_lon), lw=2, color="tab:orange", label='Sat longitude, '.format())
    ax2.plot(t, np.rad2deg(OS1.look_angle_azimuth), lw=2, color="tab:red", label='Look Angle Azimuth, '.format())
    ax2.plot(t, np.rad2deg(OS1.look_angle_elevation), lw=2, color="tab:pink", label='Look Angle Elevation, '.format())

    angle = np.deg2rad(67.5)
    ax2.legend(loc="lower left",
            bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2))
    ax2.grid(True)

    fig3 = plt.figure()
    ax3 = fig3.add_subplot(projection="polar", facecolor="lightgoldenrodyellow")
    ax3.set_title("Patrick Smith, ECE6390 HW1\nTimestep: {}".format(PCS_Time(_seconds=delta_t_1)))

    ax3.plot(OS1.look_angle_azimuth, OS1.look_angle_elevation, lw=2, color="tab:red", label='Look Angles, '.format())

    angle = np.deg2rad(67.5)
    ax3.legend(loc="lower left",
            bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2))
    ax3.grid(True)



    plt.show()
