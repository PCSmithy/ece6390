
from ece6390_constants import Quantity
from ece6390_constants import Atlanta_longitude
from ece6390_constants import Atlanta_latitude
from ece6390_constants import Earth_radius

from numpy import arccos, arcsin, sin, cos, tan, pi, sqrt, square, abs

class LookAngles(object):
    elevation = Quantity(0.0, "rad")
    azimuth = Quantity(0.0, "rad")
    
    def __init__(self, ssp_lat, ssp_lon, sat_r: Quantity, es_lat=Atlanta_latitude, es_lon=Atlanta_longitude, es_r: Quantity = Earth_radius) -> None:
        # convert everything to radians
        # 0.0
        Ls = ssp_lat if not isinstance(ssp_lat, Quantity) else ssp_lat.rad
        
        # -1.833 / 4.45
        ls = ssp_lon if not isinstance(ssp_lon, Quantity) else ssp_lon.rad

        # 0.589
        Le = es_lat if not isinstance(es_lat, Quantity) else es_lat.rad
        
        # -1.473
        le = es_lon if not isinstance(es_lon, Quantity) else es_lon.rad

        re = es_r if not isinstance(es_r, Quantity) else es_r.km
        rs = sat_r if not isinstance(sat_r, Quantity) else sat_r.km

        gamma = arccos(sin(Ls)*sin(Le)+cos(Ls)*cos(Le)*cos(ls-le))

        self.elevation.value = arccos(sin(gamma)/sqrt(1+square(re/rs)-2*(re/rs)*cos(gamma)))

        alpha = arcsin(sin(abs(le-ls))*cos(Ls)/sin(gamma))

        ssp_is_north_of_es = tan(Le)*cos(ls-le) < tan(Ls)
        ssp_is_east_of_es = ls > le

        if ssp_is_north_of_es:
            if ssp_is_east_of_es:
                self.azimuth.value = alpha
            else:
                self.azimuth.value = 2*pi-alpha
        else:
            if ssp_is_east_of_es:
                self.azimuth.value = pi - alpha
            else:
                self.azimuth.value = pi + alpha



if __name__ == "__main__":
    LA1 = LookAngles(Quantity(0, "deg"), Quantity(-105, "deg"), Quantity(42000, "km"))

    print("Azimuth: {}\tElevation: {}".format(LA1.azimuth, LA1.elevation))
