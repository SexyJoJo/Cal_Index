import pandas as pd
import metpy.calc as mpcalc
from metpy.units import units

def test_CAPE(prs, tem, dewp):
    def get_tv(tem):
        pass
    tem = list(map(get_tv, tem))
    p = pd.Series(prs).values * units.hPa
    T = pd.Series(tem).values * units.degC
    Td = pd.Series(dewp).values * units.degC
    prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degC')
    cape, cin = mpcalc.cape_cin(p, T, Td, prof, which_el='most_cape')
    return float(str(cape).split()[0])
