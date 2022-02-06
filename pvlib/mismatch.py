
import pvlib
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


# TODO: how do you rescale CEC parameters to the cell level?
# TODO: some kind of interpolation/cache layer on top of pvlib SDM
# TODO: make voltage consistent across all currents for an I-V curve?

def _balanced_linspace(start, stop, num, pivot=0):
    left = np.linspace(start, pivot, num//2, endpoint=False)
    right = np.linspace(pivot, stop, num//2)
    return np.concatenate([left, right])


class IVCurve:

    def __init__(self, i, v):
        self.i = i
        self.v = v

    def get_voltage(self, i):
        f_interp = interp1d(np.flipud(self.i), np.flipud(self.v),
                            kind='linear', fill_value='extrapolate')
        return f_interp(i)

    def get_current(self, v):
        f_interp = interp1d(self.v, self.i,
                            kind='linear', fill_value='extrapolate')
        return f_interp(v)

    def combine(self, other, connection):
        if connection == 'series':
            isc_max = np.maximum(np.max(self.i[self.v > 0], axis=0),
                                 np.max(other.i[other.v > 0], axis=0))
            i = _balanced_linspace(np.minimum(np.min(self.i, axis=0), np.min(other.i, axis=0)),
                                   np.maximum(np.max(self.i, axis=0), np.max(other.i, axis=0)),
                                   1000, pivot=isc_max)
            return IVCurve(i, self.get_voltage(i) + other.get_voltage(i))
        if connection == 'parallel':
            v = _balanced_linspace(np.minimum(np.min(self.v, axis=0), np.min(other.v, axis=0)),
                                   np.maximum(np.max(self.v, axis=0), np.max(other.v, axis=0)),
                                   num=1000, pivot=0)
            return IVCurve(self.get_current(v) + other.get_current(v), v)

        raise ValueError("connection must be one of 'series' or 'parallel', "
                         f"got '{connection}'")

    @property
    def p_mp(self):
        # TODO: delete this function?
        p = self.i * self.v
        return np.max(p)

    def plot(self):
        # TODO: delete this function
        plt.plot(self.v, self.i)
        Isc = np.max(self.i[self.v > 0])
        plt.ylim(0, Isc*1.1)


class Circuit:

    def __init__(self, conditions, connection, diode_voltage=None):
        self.conditions = conditions
        self.connection = connection
        self.diode_voltage = diode_voltage

    def apply(self, function):
        curves = []
        for condition in self.conditions:
            if isinstance(condition, dict):
                curves.append(function(**condition))
            else:
                curves.append(condition.apply(function))

        overall_curve = curves[0]
        for curve in curves[1:]:
            overall_curve = overall_curve.combine(curve, self.connection)

        Vd = self.diode_voltage
        if Vd is not None:
            overall_curve.v[overall_curve.v < -Vd] = -Vd

        return overall_curve

# %%

def make_cec_cell(params):
    """
    hacky function to use module-level CEC parameters to generate cell-level
    I-V curves.  See the multiplication and division by 72 below.
    """
    
    def get_curve(irradiance, temperature):
        
        sde_args = pvlib.pvsystem.calcparams_cec(
            np.array(irradiance),
            np.array(temperature),
            alpha_sc=params['alpha_sc'],
            a_ref=params['a_ref'],
            I_L_ref=params['I_L_ref'],
            I_o_ref=params['I_o_ref'],
            R_sh_ref=params['R_sh_ref'],
            R_s=params['R_s'],
            Adjust=params['Adjust'],
            EgRef=1.121,
            dEgdT=-0.0002677
        )
        kwargs = {
            'breakdown_factor': 2e-3,
            'breakdown_exp': 3,
            'breakdown_voltage': -15*72,
        }
        v_oc = pvlib.singlediode.bishop88_v_from_i(
            0.0, *sde_args, **kwargs
        )
        vd = np.linspace(0.99*kwargs['breakdown_voltage'], 1.01*v_oc, 1000)

        ivcurve_i, ivcurve_v, _ = pvlib.singlediode.bishop88(vd, *sde_args, **kwargs)
        return IVCurve(ivcurve_i, ivcurve_v/72)

    return get_curve

# %%

cecmods = pvlib.pvsystem.retrieve_sam('cecmod')
params = cecmods['Canadian_Solar_Inc__CS3W_400P']

# %% shadow across bottom of module (portrait)

function = make_cec_cell(params)

submodules = [
    Circuit([
        {'irradiance': 1000 if i > 0 else 200, 'temperature': 25} for i in range(24)
    ], connection='series', diode_voltage=0.5)
    for j in range(3)
]

module = Circuit(submodules, connection='series')

curve = module.apply(function)
curve.plot()

# %% uneven soiling

function = make_cec_cell(params)

submodules = [
    Circuit([
        {'irradiance': 1000 - j*100, 'temperature': 25} for i in range(24)
    ], connection='series', diode_voltage=0.5)
    for j in range(3)
]

module = Circuit(submodules, connection='series')

curve = module.apply(function)
curve.plot()


# %% time series

function = make_cec_cell(params)

times = pd.date_range('2019-01-01 10:00', '2019-01-01 15:00', freq='h', tz='Etc/GMT+5')

submodules = [
    Circuit([
        {'irradiance': pd.Series(1000 - j*100, index=times),
         'temperature': pd.Series(25, index=times)
         } for i in range(24)
    ], connection='series', diode_voltage=0.5)
    for j in range(3)
]

module = Circuit(submodules, connection='series')

curve = module.apply(function)
curve.plot()
