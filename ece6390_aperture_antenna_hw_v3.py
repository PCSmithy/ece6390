import matplotlib.pyplot as plt
import numpy as np

from numpy import pi, cos, sin, power, exp, log10

from scipy.special import j0
import scipy.integrate as integrate

from prettytable import PrettyTable

j = complex(0,1)
eta = 377

def E_diamond(rho):
    return E_0*power(cos((pi*rho)/(2*R)), n)

def E_diamond_obs(r, theta):
    return (j*k*cos(theta)/r)*exp(-j*k*r)*(integrate.quad(lambda rho: E_diamond(rho)*rho*j0(k*rho*sin(theta)), 0, R)[0])

def S_av(r, theta):
    return (1/(2*eta))*power(np.real(E_diamond_obs(r, theta)), 2)

def Pt():
    return 2*pi*(r**2)*(integrate.quad(lambda theta: np.real(S_av(r, theta)*sin(theta)), 0, pi/2)[0])

colors = ["red", "blue", "green", "orange"]

if __name__ == "__main__":
    E_0 = 1                             # arbitrary
    wavelength = 0.1                    # 10cm
    k = 2*pi/wavelength                 # wavenumber

    lower = 0
    upper = pi/2
    N = 1000

    theta_list = np.linspace(lower, upper, N)
    delta_theta = (upper - lower) / N

    pt = PrettyTable(["R scale", "n", "Peak Gain (dB)", "SLL (dB)", "3dB Beamwidth (deg)"])
    pt.custom_format["Peak Gain (dB)"] = lambda f, v: f"{v:.3f}"
    pt.custom_format["SLL (dB)"] = lambda f, v: f"{v:.3f}"
    pt.custom_format["3dB Beamwidth (deg)"] = lambda f, v: f"{v:.3f}"

    for scale in [5, 10]:
        # compute dish radius
        R = scale*wavelength

        # ensure we are well into the far field 
        r = (power((2*R), 2) / wavelength) * 2

        fig = plt.figure()
        ax = fig.add_subplot(projection="polar", facecolor="lightgoldenrodyellow")

        for n in range(4):
            E_theta = np.zeros(N, dtype=np.complex_)
            Sav_theta = np.zeros(N)
            D_theta = np.zeros(N)
            D_theta_log = np.zeros(N)

            # compute the total power radiated
            total_power = Pt()

            # compute the directivity at every theta
            for i, theta in enumerate(theta_list):
                Sav_theta[i] = S_av(r, theta)
                D_theta[i] = Sav_theta[i] / (total_power/(4*pi*(r**2)))

            # convert the directivity to the log scale
            D_theta_log = 10*log10(D_theta)

            theta_pi = np.linspace(-pi/2, pi/2, N*2)
            D_theta_pi = np.zeros(2*N)
            d1_theta = 0
            d1_theta_prev = 0
            side_lobe_gain = 0

            peak_gain = D_theta_log[0]
            peak_gain_minus_3db = peak_gain - 3
            _3dB_beamwidth = 0

            # inspect the directivity across theta to find key radiation parameters,
            # mirror the other half of the data into [-pi/2 to 0] for plotting
            for i, _ in enumerate(D_theta_pi):
                if i >= N: # [0, pi/2]
                    gain = D_theta_log[-1*(N-i)]
                    D_theta_pi[i] = gain

                    # find the inflection point of the first side lobe
                    d1_theta_prev = d1_theta
                    d1_theta = (D_theta_pi[i] - D_theta_pi[i-1]) / delta_theta
                    
                    if (d1_theta < 0) and (d1_theta_prev > 0) and (side_lobe_gain == 0):
                        side_lobe_gain = gain
                        SLL = side_lobe_gain - peak_gain

                    # find the 3dB beam width
                    if (gain < peak_gain_minus_3db) and (_3dB_beamwidth == 0):
                        _3dB_beamwidth = (theta_pi[i]*2)*180/pi

                else: # [-pi/2, 0]
                    D_theta_pi[i] = D_theta_log[N-i-1]

            label = f"R: {scale}X, n: {n}"
            pt.add_row([scale, n, peak_gain, SLL, _3dB_beamwidth])

            ax.plot(theta_pi, D_theta_pi, lw=1, color=f"tab:{colors[n]}", label=label)

        ax.set_title("R = {}*wavelengths".format(scale))
        ax.legend()
        ax.set_rmin(0)
        ax.set_rmax(39)
        angle = 90
        ax.set_thetamin(-angle)
        ax.set_thetamax(angle)
        ax.grid(True)
        ax.set_theta_zero_location("N")

    print(pt)
    plt.show()    


