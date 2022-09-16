import matplotlib.pyplot as plt
import numpy as np

from numpy import float64 as f64
from numpy import pi

from ece6390_constants import Sun_mass as Msun
from ece6390_constants import Earth_radius as R_earth
from ece6390_constants import Quantity, Sun_radius

from orbit import EllipticalOrbit, NumericalOrbitSimulator

from time_lib import PCS_Time

if __name__ == "__main__":

    t0 = f64(0)
    tf = PCS_Time(_hours=12.5).seconds
    delta_t_1 = f64(5)          # seconds
    n = int((tf-t0)/delta_t_1)  # number of samples

    # part a
    OS1 = NumericalOrbitSimulator(R_init        =Quantity(f64(20e3), "km"), 
                                  theta_init    =Quantity(f64(0), "rad"), 
                                  V_r_init      =Quantity(f64(0), "km/s"), 
                                  V_theta_init  =Quantity(f64(5), "km/s"), 
                                  delta_t       =Quantity(delta_t_1, "s"),
                                  n             =n)

    # part b
    # tf = PCS_Time(_days=2, _hours=1, _minutes=5).seconds
    tf = PCS_Time(_days=2, _hours=2).seconds
    n = int((tf-t0)/delta_t_1)
    OS2 = NumericalOrbitSimulator(R_init        =Quantity(f64(20e3), "km"), 
                                  theta_init    =Quantity(f64(0), "rad"), 
                                  V_r_init      =Quantity(f64(-3), "km/s"), 
                                  V_theta_init  =Quantity(f64(5), "km/s"), 
                                  delta_t       =Quantity(delta_t_1, "s"),
                                  n             =n)

    # part c
    tf = PCS_Time(_weeks=0.5).seconds
    n = int((tf-t0)/delta_t_1)
    OS3 = NumericalOrbitSimulator(R_init        =Quantity(f64(20e3), "km"), 
                                  theta_init    =Quantity(f64(0), "rad"), 
                                  V_r_init      =Quantity(f64(-6), "km/s"), 
                                  V_theta_init  =Quantity(f64(5), "km/s"), 
                                  delta_t       =Quantity(delta_t_1, "s"),
                                  n             =n)

    # part d - Halley's Comet
    tf = PCS_Time(_years=78).seconds
    delta_t_2 = 1000 # seconds
    n = int((tf-t0)/delta_t_2)
    Halley_perihelion   = Quantity(f64(0.59278), "au") # from wikipedia
    Halley_aphelion     = Quantity(f64(35.14), "au")

    # compute V_min / V_max / Period of the Halley's Comet orbit to give us our initial conditions for the numerical simulation
    hc = EllipticalOrbit(Mp=Msun, apoapsis=Halley_aphelion, periapsis=Halley_perihelion)

    Halley = NumericalOrbitSimulator(R_init         =Halley_perihelion, 
                                     theta_init     =Quantity(0, "rad"), 
                                     V_r_init       =Quantity(0, "km/s"),    
                                     V_theta_init   =Quantity(hc.V_max, "km/s"), 
                                     delta_t        =Quantity(delta_t_2, "s"), 
                                     Mp             =Msun, 
                                     n              =n)


    # plots for part 1-3
    earth_r = np.full((360,1), R_earth())
    earth_theta = np.linspace(0, 2*pi, 360)

    fig = plt.figure()
    ax = fig.add_subplot(projection="polar", facecolor="lightgoldenrodyellow")
    ax.set_title("Patrick Smith, ECE6390 HW1\nTimestep: {}".format(PCS_Time(_seconds=delta_t_1)))

    ax.plot(earth_theta, earth_r, lw=1, color="tab:green", label='Earth, Radius 6390 km')

    ax.plot(OS1.theta, OS1.r, lw=2, color="tab:orange", label='Part A, Eccentricity: {:.3f}, Period: {}'.format(OS1.eccentricity, OS1.period))
    ax.plot(OS2.theta, OS2.r, ls="--", color="tab:blue", label='Part B, Eccentricity: {:.3f}, Period: {}'.format(OS2.eccentricity, OS2.period))
    ax.plot(OS3.theta, OS3.r, ls="-.", color="tab:red", label='Part C, Escaped')
    ax.tick_params(grid_color="palegoldenrod")

    ax.set_rmax(150e3)
    ax.set_rlabel_position(40)  # Move radial labels away from plotted line
    angle = np.deg2rad(67.5)
    ax.legend(loc="lower left",
            bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2))
    ax.grid(True)

    # plot for part 4
    sun_r = np.full((360, 1), Sun_radius())
    sun_theta = np.linspace(0, 2*pi, 360)

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(projection="polar", facecolor="lightgoldenrodyellow")
    ax2.set_title("Patrick Smith, ECE6390 HW1, Halley's Comet\nEccentricity: {:.3f}, Period: {}\nTimestep: {} (or {} seconds)".format(Halley.eccentricity, Halley.period, PCS_Time(_seconds=delta_t_2), delta_t_2))

    ax2.plot(sun_theta, sun_r, lw=3, color="tab:orange", label='Sun, Radius {}'.format(Sun_radius))
    ax2.plot(Halley.theta, Halley.r, ls="-.", color="tab:blue", label="Halley's Comet")

    ax2.tick_params(grid_color="palegoldenrod")

    ax2.set_rlabel_position(40)  # Move radial labels away from plotted line
    angle = np.deg2rad(67.5)
    ax2.legend(loc="lower left",
            bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2))
    ax2.grid(True)

    plt.show()
