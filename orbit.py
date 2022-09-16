from typing import Optional

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from numpy import float64 as f64
from numpy import pi
from numpy import sqrt
from numpy import float_power as fpow

from ece6390_constants import Atlanta_latitude as ES_LAT
from ece6390_constants import Atlanta_longitude as ES_LON
from ece6390_constants import Earth_angular_veolcity as wEarth
from ece6390_constants import Earth_mass as Me
from ece6390_constants import Earth_radius as R_earth
from ece6390_constants import Sun_mass as Msun
from ece6390_constants import Gravitational_constant as Gravity
from ece6390_constants import Sidereal_day
from ece6390_constants import Quantity

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

    # arrays to hold the state data r, r_dot, theta, and theta_dot
    r           = None
    delta_r     = None
    r_dot       = None
    theta       = None
    delta_theta = None
    theta_dot   = None

    def __init__(self, R_init: Quantity, 
                        theta_init: Quantity, 
                        V_r_init: Quantity,  
                        V_theta_init: Quantity, 
                        delta_t: Quantity, 
                        n: int, 
                        Grav=Gravity, 
                        Mp=Me, ) -> None:

        if delta_t.unit != self.delta_t.unit:
            raise AttributeError("Wrong units on delta_t arguement, must be (s)")
        self.delta_t = delta_t
        self.G = Grav
        self.M = Mp

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

                self.delta_r[i+1] = self.delta_r_n_plus_one(self.delta_r[i], self.r[i], self.delta_theta[i])
                self.delta_theta[i+1] = self.delta_theta_n_plus_one(self.delta_theta[i], self.delta_r[i], self.r[i])
                continue

            # exit out at the end to avoid array size mismatches
            if i == n-1:
                break

            self.r[i+1] = self.r[i] + self.delta_r[i]
            self.theta[i+1] = (self.theta[i] + self.delta_theta[i]) % (2*pi)

            self.delta_r[i+1] = self.delta_r_n_plus_one(self.delta_r[i], self.r[i], self.delta_theta[i])
            self.delta_theta[i+1] = self.delta_theta_n_plus_one(self.delta_theta[i], self.delta_r[i], self.r[i])

            # detect when we've reached 1 full period
            if (self.theta[i] > (1.9*pi)) and (self.theta[i+1] < (0.1*pi)):
                self.T = PCS_Time(_seconds=i*self.delta_t())
                print("1 Period took {} samples, or {}".format(i, PCS_Time(_seconds=i*self.delta_t())))

    def delta_r_n_plus_one(self, delta_r_n, r_n, delta_theta_n):
        a = (r_n+(1/2)*delta_r_n)*fpow(delta_theta_n, 2)
        b = self.G*self.M*fpow(self.delta_t(), 2)/fpow((r_n+(0.5*delta_r_n)), 2)

        return (delta_r_n + (a - b))

    def delta_theta_n_plus_one(self, delta_theta_n, delta_r_n, r_n):
        return (delta_theta_n - (2*delta_r_n*delta_theta_n)/(r_n+(1/2)*delta_r_n))

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


class OrbitSimulatorWithLookAngles(NumericalOrbitSimulator):
    rotation_of_earth_in_fixed_frame = None
    ssp_lat = None
    ssp_lon = None
    look_angle_azimuth = None
    look_angle_elevation = None
    es_lat = Quantity(0, "deg")
    es_lon = Quantity(0, "deg")

    def __init__(self, R_init: Quantity, theta_init: Quantity, V_r_init: Quantity, V_theta_init: Quantity, delta_t: Quantity, n: int, ssp_lat_init: Quantity, es_lat: Quantity, es_lon: Quantity, Grav=Gravity, Mp=Me) -> None:
        super().__init__(R_init, theta_init, V_r_init, V_theta_init, delta_t, n, Grav, Mp)

        if es_lat.unit != "deg" or es_lon.unit != "deg":
            raise AttributeError("Wrong units on Earth State coordinates, must be (deg)")
        
        self.es_lon = es_lon
        self.es_lat = es_lat

        self.rotation_of_earth_in_fixed_frame = f64(np.zeros(n))
        self.ssp_lat                = f64(np.zeros(n))
        self.ssp_lon                = f64(np.zeros(n))
        self.look_angle_azimuth     = f64(np.zeros(n))
        self.look_angle_elevation   = f64(np.zeros(n))

        for i in range(n):
            # hard-coded orbit on the equator for now...
            self.ssp_lat[i] = np.deg2rad(ES_LAT())

            # the longitude of the satellite is the difference in the rotation of the earth in a fixed frame and rotation of the satellite in the fixed frame
            # 0-2pi rad
            self.ssp_lon[i] = (self.rotation_of_earth_in_fixed_frame[i] - self.theta[i]) if (self.rotation_of_earth_in_fixed_frame[i] - self.theta[i]) > 0 else (2*np.pi + self.rotation_of_earth_in_fixed_frame[i] - self.theta[i])

            le = self.es_lon()
            Le = self.es_lat()

            ls = self.ssp_lon[i]
            Ls = self.ssp_lat[i]

            # compute the azimuth
            gamma = np.arccos(np.sin(Ls)*np.sin(Le)+np.cos(Ls)*np.cos(Le)*np.cos(ls-le))
            alpha = np.arcsin(np.sin(abs(le-ls)*np.cos(Ls)/np.sin(gamma)))

            self.look_angle_azimuth[i] = (np.pi + alpha) if (self.ssp_lon[i] < np.deg2rad(self.es_lon())) else (np.pi - alpha)

            # compute elevation
            re = R_earth()
            rs = self.r[i]

            self.look_angle_elevation[i] = np.arccos(np.sin(gamma)/np.sqrt(1+np.square(re/rs)-2*(re/rs)*np.cos(gamma)))

            if i == (n-1):
                break
            self.rotation_of_earth_in_fixed_frame[i+1] = (self.rotation_of_earth_in_fixed_frame[i] + wEarth * delta_t) % (2*np.pi)



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


    t0 = f64(0)
    tf = PCS_Time(_weeks=1).seconds
    delta_t_1 = f64(5)          # seconds
    n = int((tf-t0)/delta_t_1)  # number of samples

    OS1 = OrbitSimulatorWithLookAngles(R_init        =geostationary_orbit.radius, 
                                        theta_init    =ES_LON, 
                                        V_r_init      =Quantity(f64(0), "km/s"), 
                                        V_theta_init  =geostationary_orbit.velocity, 
                                        delta_t       =Quantity(delta_t_1, "s"),
                                        n             =n,
                                        ssp_lat_init  =Quantity(0, "deg"),
                                        es_lat        =ES_LAT,
                                        es_lon        =ES_LON)