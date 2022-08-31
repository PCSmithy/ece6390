from typing import Optional

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

class OrbitSimulator(object):
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

    def __init__(self, R_init: Quantity, theta_init: Quantity, V_r_init: Quantity,  V_theta_init: Quantity, delta_t: Quantity, n: int, Grav=Gravity, Mp=Me, ) -> None:
        if delta_t.unit != self.delta_t.unit:
            raise AttributeError("Wrong units on delta_t arguement, must be (s)")
        self.delta_t = delta_t
        self.G = Grav
        self.M = Mp

        self.r           = np.zeros(n)
        self.delta_r     = np.zeros(n)
        self.r_dot       = np.zeros(n)
        self.theta       = np.zeros(n)
        self.delta_theta = np.zeros(n)
        self.theta_dot   = np.zeros(n)

        # integrate from t0 to tf
        for i in range(n):
            # initialize arrays
            if i == 0:
                self.r[i] = R_init()
                self.r_dot[i] = V_r_init()
                self.theta[i] = theta_init()
                self.theta_dot[i] = V_theta_init() / self.r[i]

                self.delta_r[i] = self.r_dot[i]*timestep
                self.r[i+1] = self.r[i] + self.delta_r[i]

                self.delta_theta[i] = self.theta_dot[i]*timestep
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

            if (self.theta[i] > (1.9*pi)) and (self.theta[i+1] < (0.1*pi)):
                self.T = PCS_Time(_seconds=i*timestep)
                print("1 Period took {} samples, or {}".format(i, PCS_Time(_seconds=i*timestep)))

    def delta_r_n_plus_one(self, delta_r_n, r_n, delta_theta_n):
        return (delta_r_n + ((r_n+(1/2)*delta_r_n)*fpow(delta_theta_n, 2) - self.G*self.M*fpow(self.delta_t(), 2)/fpow(r_n, 2)))

    def delta_theta_n_plus_one(self, delta_theta_n, delta_r_n, r_n):
        return (delta_theta_n - (2*delta_r_n*delta_theta_n)/(r_n+(1/2)*delta_r_n))

    @property
    def period(self):
        return self.T

    @property
    def apogee(self):
        return Quantity(max(self.r), "km")

    @property
    def perigee(self):
        return Quantity(min(self.r), "km")


if __name__ == "__main__":

    t0 = f64(0)
    tf = f64(PCS_Time(_hours=12.5).seconds)
    timestep = f64(5)
    n = int((tf-t0)/timestep)

    t = np.linspace(t0, tf, n)

    # initial values

    OS1 = OrbitSimulator(R_init=Quantity(f64(20e3), "km"), 
                        theta_init=Quantity(f64(0), "rad"), 
                        V_r_init=Quantity(f64(0), "km/s"), 
                        V_theta_init=Quantity(f64(5), "km/s"), 
                        delta_t=Quantity(timestep, "s"),
                        n=n)


    tf = f64(PCS_Time(_days=2, _hours=1, _minutes=5).seconds)
    n = int((tf-t0)/timestep)
    OS2 = OrbitSimulator(R_init=Quantity(f64(20e3), "km"), 
                        theta_init=Quantity(f64(0), "rad"), 
                        V_r_init=Quantity(f64(-3), "km/s"), 
                        V_theta_init=Quantity(f64(5), "km/s"), 
                        delta_t=Quantity(timestep, "s"),
                        n=n)


    tf = f64(PCS_Time(_weeks=0.5).seconds)
    n = int((tf-t0)/timestep)
    OS3 = OrbitSimulator(R_init=Quantity(f64(20e3), "km"), 
                        theta_init=Quantity(f64(0), "rad"), 
                        V_r_init=Quantity(f64(-6), "km/s"), 
                        V_theta_init=Quantity(f64(5), "km/s"), 
                        delta_t=Quantity(timestep, "s"),
                        n=n)


    earth_r = np.full((360,1), R_earth())
    earth_theta = np.linspace(0, 2*pi, 360)

    fig = plt.figure()
    ax = fig.add_subplot(projection="polar", facecolor="lightgoldenrodyellow")
    ax.set_title("Patrick Smith, ECE6390 HW1")

    ax.plot(earth_theta, earth_r, lw=1, color="tab:green", label='Earth, Radius 6390 km')

    ax.plot(OS1.theta, OS1.r, lw=2, color="tab:orange", label='Part A, Apogee {}, Period {}'.format(OS1.apogee, OS1.T))
    ax.plot(OS2.theta, OS2.r, ls="--", color="tab:blue", label='Part B, Apogee {}, Period {}'.format(OS2.apogee, OS2.T))
    ax.plot(OS3.theta, OS3.r, ls="-.", color="tab:red", label='Part C, Escaped')
    ax.tick_params(grid_color="palegoldenrod")

    ax.set_rmax(150e3)
    ax.set_rlabel_position(40)  # Move radial labels away from plotted line
    angle = np.deg2rad(67.5)
    ax.legend(loc="lower left",
            bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2))
    ax.grid(True)

    plt.show()
