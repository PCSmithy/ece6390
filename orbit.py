from typing import Optional

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from numpy import float64 as f64
from numpy import pi, sqrt, square, cos, sin, arctan2, arctan
from numpy import float_power as fpow

from ece6390_constants import Atlanta_latitude as ES_LAT
from ece6390_constants import Atlanta_longitude as ES_LON
from ece6390_constants import Earth_angular_veolcity as wEarth
from ece6390_constants import Earth_mass as Me
from ece6390_constants import Earth_radius as R_earth
from ece6390_constants import Sun_mass as Msun
from ece6390_constants import Solar_pressure as Fsun
from ece6390_constants import Gravitational_constant as Gravity
from ece6390_constants import Sidereal_day
from ece6390_constants import Quantity

from look_angles import LookAngles
from time_lib import PCS_Time

class EllipticalOrbit(object):
    T = Quantity(unit="seconds")    # period
    R = Quantity(unit="km")         # radius
    G = Quantity()                  # gravitational constant
    M = Quantity()                  # mass of body to orbit
    V = Quantity(unit="km/s")       # velocity
    e = Quantity(unit="")           # eccentricity
    a = Quantity(unit="km")         # semi-major axis
    b = Quantity(unit="km")         # semi-minor axis

    def __init__(self, Grav=Gravity, Mp=Me, T_seconds: Optional[Quantity] = None, apoapsis: Optional[Quantity] = None, periapsis: Optional[Quantity] = None) -> None:
        if T_seconds is None and apoapsis is None and periapsis is None:
            raise AttributeError("Cannot compute orbit parameters without either the Period T, or the apoapsis / periapsis")

        self.G = Grav
        self.M = Mp

        if T_seconds is None:
            alpha = periapsis / apoapsis
            self.e.value = (1 - alpha)/(1 + alpha)

            if periapsis.unit not in ["m", "km", "au"] or apoapsis.unit not in ["m", "km", "au"]:
                raise AttributeError("Incorrect Radius R units")
                
            self.a.value = (apoapsis.km + periapsis.km) / 2
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
        return self.e.value

    @property
    def period(self):
        return PCS_Time(_seconds=self.T.value)

    @property
    def radius(self):
        return self.R

    @property
    def V_max(self):
        return Quantity(f64(sqrt(self.G*self.M*(1+self.e()) / (self.a*(1-self.e())))), "km/s")
    
    @property
    def V_min(self):
        return Quantity(f64(sqrt(self.G*self.M*(1-self.e()) / (self.a*(1+self.e())))), "km/s")


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


class NumericalOrbitSimulator(object):
    delta_t = Quantity(unit="s")
    G = Quantity()
    M = Quantity()
    T = None

    # orbit params
    r           = None
    delta_r     = None
    r_dot       = None
    theta       = None
    delta_theta = None
    theta_dot   = None

    # solar pressure params
    alpha = None
    As = None
    Ms = None

    # look angle params
    rotation_of_earth_in_fixed_frame = None
    ssp_lat = None
    ssp_lon = None
    look_angle_azimuth = None
    look_angle_elevation = None
    es_lat = Quantity(0, "deg")
    es_lon = Quantity(0, "deg")

    def __init__(self,  R_init: Quantity, 
                        theta_init: Quantity, 
                        V_r_init: Quantity,  
                        V_theta_init: Quantity, 
                        delta_t: Quantity, 
                        n: int, 
                        As: Quantity=Quantity(0, "m^2"),
                        Ms: Quantity=Quantity(10, "kg"),
                        alpha=0.5,
                        Grav=Gravity, 
                        Mp=Me, ) -> None:

        if delta_t.unit != self.delta_t.unit:
            raise AttributeError("Wrong units on delta_t arguement, must be (s)")
        self.delta_t = delta_t
        self.G = Grav
        self.M = Mp

        self.alpha = alpha
        self.As = As
        self.Ms = Ms

        self.rotation_of_earth_in_fixed_frame = f64(np.zeros(n))
        self.r           = f64(np.zeros(n))
        self.delta_r     = f64(np.zeros(n))
        self.r_dot       = f64(np.zeros(n))
        self.theta       = f64(np.zeros(n))
        self.delta_theta = f64(np.zeros(n))
        self.theta_dot   = f64(np.zeros(n))

        # integrate from t0 to tf
        for i in range(n):
            # initialize arrays
            if i == 0:
                self.r[i] = R_init.km
                self.r_dot[i] = V_r_init()
                self.theta[i] = theta_init.rad
                self.theta_dot[i] = V_theta_init() / self.r[i]

                self.delta_r[i] = self.r_dot[i]*self.delta_t()
                self.r[i+1] = self.r[i] + self.delta_r[i]

                self.delta_theta[i] = self.theta_dot[i]*self.delta_t()
                self.theta[i+1] = self.theta[i] + self.delta_theta[i]

                self.delta_r[i+1] = self.delta_r_n_plus_one(self.r[i], self.theta[i], self.delta_r[i], self.delta_theta[i])
                self.delta_theta[i+1] = self.delta_theta_n_plus_one(self.r[i], self.theta[i], self.delta_r[i], self.delta_theta[i])
                continue

            # exit out at the end to avoid array size mismatches
            if i == n-1:
                break

            self.rotation_of_earth_in_fixed_frame[i+1] = (self.rotation_of_earth_in_fixed_frame[i] + wEarth * delta_t) % (2*np.pi)

            self.r[i+1] = self.r[i] + self.delta_r[i]
            self.theta[i+1] = (self.theta[i] + self.delta_theta[i]) % (2*pi)

            self.delta_r[i+1] = self.delta_r_n_plus_one(self.r[i], self.theta[i], self.delta_r[i], self.delta_theta[i])
            self.delta_theta[i+1] = self.delta_theta_n_plus_one(self.r[i], self.theta[i], self.delta_r[i], self.delta_theta[i])

            # detect when we've reached 1 full period
            # if (self.theta[i] > (1.9*pi)) and (self.theta[i+1] < (0.1*pi)):
            #     self.T = PCS_Time(_seconds=i*self.delta_t())
            #     print("1 Period took {} samples, or {}".format(i, PCS_Time(_seconds=i*self.delta_t())))

    def delta_r_n_plus_one(self, r_n, theta_n, delta_r_n, delta_theta_n):
        dt = self.delta_t()
        r = r_n + (1/2)*delta_r_n
        phi = theta_n + (1/2)*delta_theta_n
        mu = self.M*self.G
        As = self.As()
        Ms = self.Ms()

        a = r*square(delta_theta_n)
        b = mu*square(dt)/square(r)
        c = Fsun()*self.alpha*As*square(dt)*cos(phi)/Ms

        return (delta_r_n + a - b + c)

    def delta_theta_n_plus_one(self, r_n, theta_n, delta_r_n, delta_theta_n):
        dt = self.delta_t()
        r = r_n + (1/2)*delta_r_n
        phi = theta_n + (1/2)*delta_theta_n
        As = self.As()
        Ms = self.Ms()

        a = 2*delta_r_n*delta_theta_n
        b = Fsun()*self.alpha*As*square(dt)*sin(phi)/Ms

        return (delta_theta_n - (a + b)/r)

    def compute_look_angles(self, es_lat: Quantity = ES_LAT, es_lon: Quantity = ES_LON) -> None:
        n = len(self.r)
        self.ssp_lat = f64(np.zeros(n))
        self.ssp_lon = f64(np.zeros(n))
        self.look_angle_azimuth = f64(np.zeros(n))
        self.look_angle_elevation = f64(np.zeros(n))
        self.es_lat = es_lat
        self.es_lon = es_lon

        for i in range(n):
            # hard-coded orbit on the equator for now...
            self.ssp_lat[i] = 0.0

            # the longitude of the satellite is the difference in the rotation of the earth in a fixed frame and rotation of the satellite in the fixed frame
            # 0-2pi rad
            self.ssp_lon[i] = (self.theta[i] - self.rotation_of_earth_in_fixed_frame[i]) if (self.theta[i] - self.rotation_of_earth_in_fixed_frame[i]) > 0 else (2*np.pi + self.theta[i] - self.rotation_of_earth_in_fixed_frame[i])

            LA = LookAngles(self.ssp_lat[i], self.ssp_lon[i], self.r[i])

            self.look_angle_azimuth[i] = LA.azimuth()
            self.look_angle_elevation[i] = LA.elevation()

    @property
    def period(self):
        return self.T

    @property
    def eccentricity(self):
        return (self.apogee - self.perigee) / (self.apogee + self.perigee)

    @property
    def apogee(self):
        return Quantity(max(self.r), "km")

    @property
    def perigee(self):
        return Quantity(min(self.r), "km")



if __name__ == "__main__":

    T = Quantity(Sidereal_day.seconds, "seconds")
    geostationary_orbit = CircularOrbit(T_seconds=T)
    print("Geostationary Orbit around Earth has a Period of {} and a Radius of {}, with a velocity of {}".format(geostationary_orbit.period, geostationary_orbit.radius, geostationary_orbit.angular_velocity))

    # R2 = Quantity(4.216417e4, "km")
    R2 = Quantity(20e3, "km")
    geostat_orbit = CircularOrbit(R=R2)
    print("Geostationary Orbit2 around Earth has a Period of {} and a Radius of {}, with a velocity of {}".format(geostat_orbit.period, geostat_orbit.radius, geostat_orbit.velocity))

    Halley_perihelion = Quantity(f64(0.59278), "au")
    Halley_aphelion = Quantity(f64(35.14), "au")
    hc = EllipticalOrbit(Mp=Msun, apoapsis=Halley_aphelion, periapsis=Halley_perihelion)
    print("Halley's comet period: {},\nEccentricity: {:.3f},\nSemi-major axis: {},\nV_max: {},\nV_min: {},\naphilion: {},\nperihelion: {}\n".format(hc.period, hc.eccentricity, hc.a, hc.V_max, hc.V_min, Halley_aphelion, Halley_perihelion))

    cross_sectional_area = Quantity(30, "m^2")
    mass = Quantity(1000, "kg")

    area_density = mass / cross_sectional_area

    t0 = f64(0)
    tf = PCS_Time(_years=1).seconds
    delta_t_1 = f64(10)          # seconds
    n = int((tf-t0)/delta_t_1)  # number of samples
    t = np.linspace(t0, tf, n)

    OS1 = NumericalOrbitSimulator(R_init        =geostationary_orbit.radius, 
                                theta_init    =ES_LON, 
                                V_r_init      =Quantity(f64(0), "km/s"), 
                                V_theta_init  =geostationary_orbit.velocity, 
                                delta_t       =Quantity(delta_t_1, "s"),
                                n             =n,
                                As              =cross_sectional_area,
                                Ms              =mass)

    OS1.compute_look_angles()

    earth_r = np.full((360,1), R_earth())
    earth_theta = np.linspace(0, 2*pi, 360)

    fig = plt.figure()
    ax = fig.add_subplot(projection="polar", facecolor="lightgoldenrodyellow")
    ax.set_title("Patrick Smith, ECE6390 HW2\nTimestep: {}".format(PCS_Time(_seconds=delta_t_1)))

    ax.plot(earth_theta, earth_r, lw=1, color="tab:green", label='Earth, Radius 6390 km')

    ax.plot(OS1.theta, OS1.r, lw=2, color="tab:orange", label='GEO orbit, Eccentricity: {:.3f}, Period: {}'.format(OS1.eccentricity, OS1.period))
    ax.tick_params(grid_color="palegoldenrod")

#     ax.set_rmax(150e3)
    ax.set_rlabel_position(40)  # Move radial labels away from plotted line
    angle = np.deg2rad(67.5)
    ax.legend(loc="lower left",
            bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2))
    ax.grid(True)

    fig2 = plt.figure()
    ax2 = fig2.add_subplot()
    ax2.set_title("Patrick Smith, ECE6390 HW1\nTimestep: {}\nArea Density: {:.3f}".format(PCS_Time(_seconds=delta_t_1), area_density))

    # ax2.plot(t, np.rad2deg(OS1.rotation_of_earth_in_fixed_frame), lw=2, color="tab:blue", label='rotation of earth, '.format())
    # ax2.plot(t, np.rad2deg(OS1.theta), lw=2, color="tab:green", label='rotation of sat, '.format())
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
    ax3.set_theta_zero_location("S")
    ax3.set_rlabel_position(80)  # Move radial labels away from plotted line
    ax3.plot(OS1.look_angle_azimuth, np.rad2deg(OS1.look_angle_elevation), lw=2, color="tab:red", label='Look Angles, '.format())
    ax3.set_rmax(90)

    angle = np.deg2rad(67.5)
    ax3.legend(loc="lower left",
            bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2))
    ax3.grid(True)

    plt.show()

    