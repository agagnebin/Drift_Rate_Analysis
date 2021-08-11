from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
import numpy as np
#import dateutil
from astropy.convolution import convolve, Gaussian1DKernel
from astropy import constants as Const
from argparse import Namespace


def to_Angle(val, unit='degree'):
    if isinstance(val, (u.quantity.Quantity, Angle)):
        return val
    return Angle(val, unit)


#def to_Time(times, toffset=None):
 #   """Convert to Time."""

def to_Time(times, toffset = None):
    """Convert to Time."""
    if dates is None or dates is False:
        return None
    if isinstance(dates, list) or isinstance(times, (list, ndarray)):
        if isinstance(times, ndarray):
            times = list(times)
        if not isinstance(times, list):
            times = [times] * len(dates)
        if not isinstance(dates, list):
            dates = [dates] * len(times)
        if len(dates) != len(times):
            raise ValueError("dates/times list lengths must match.")
        return_Time = []
        if len(dates) > 1000:
            print("Converting {} time entries - could take a moment.".format(len(dates)))
        for _date, _time in zip(dates, times):
            return_Time.append(get_astropytime(_date, _time))
        return Time(return_Time)
    if isinstance(dates, str):
        if dates.lower() == 'none':
            return None
        if dates == '<':
            return Time('2000-01-01', scale='utc')
        if dates == '>':
            return Time.now() + TimeDelta(1000, format='jd')
        if dates.lower() == 'now' or dates.lower() == 'current':
            return Time.now()
    if isinstance(dates, Time):
        return_Time = dates
    elif isinstance(dates, datetime):
        return_Time = Time(dates, format='datetime')
    else:
        try:
            dates = float(dates)
            if dates > 1000000000.0:
                return_Time = Time(dates, format='gps')
            elif dates > 2400000.0 and dates < 2500000.0:
                return_Time = Time(dates, format='jd')
            else:
                raise ValueError(f'Invalid format:  date as a number should be gps time '
                                 f'or julian date, not {dates}.')
        except ValueError:
            dates = dates.replace('/', '-')
            try:
                return_Time = Time(dates, scale='utc')
            except ValueError:
                raise ValueError(
                    f'Invalid format:  YYYY[/-]M[/-]D [HH:MM:SS], not {dates}')
    # add on times
    if times is None or abs(times) < 1E-6:
        return return_Time
    try:
        times = float(times)
        return return_Time + TimeDelta(times * 3600.0, format='sec')
    except ValueError:
        pass
    sign_of_times = 1.0
    if times[0] == '-':
        sign_of_times = -1.0
        times = times[1:]
    add_time = 0.0
    for i, d in enumerate(times.split(':')):
        add_time += (float(d)) * 3600.0 / (60.0**i)
    add_time *= sign_of_times
    return return_Time + TimeDelta(add_time, format='sec')

def to_separation(c1, c2):
    if not isinstance(c1, SkyCoord):
        c1 = SkyCoord(to_Angle(c1[0]), to_Angle(c1[1]))
    if not isinstance(c2, SkyCoord):
        c2 = SkyCoord(to_Angle(c2[0]), to_Angle(c2[1]))
    return c1.separation(c2)


class BaseEphem:
    param = ['times', 'ra', 'dec', 'az', 'el', 'x', 'y', 'z', 'D', 'dt',
             'radot', 'decdot', 'azdot', 'eldot', 'xdot', 'ydot', 'zdot',
             'Ddot', 'Ddotdot']

    def __init__(self):
        """
        Provides an init'd base class for ephemerides per:
            times :  Time
            dt :  np.array
            ra, dec, az, el:  Angle
            radot, decdot, azdot, eldot:  deg/sec (Quantity)
            x, y, z:  m (Quantity)
            xdot, ydot, zdot, Ddot:  m/s (Quantity)
            Ddotdot: m/s/s (Quantity)
        """
        self.initall()
        self._E = None  # Class archive for interp
        self.c0 = Const.c  # Get used constants handy
        self.au = Const.au
        self.ly = 1.0 * u.lyr.to('m')

    def initall(self, **kwargs):
        for par in self.param:
            if par in kwargs.keys():
                setattr(self, par, kwargs[par])
            else:
                setattr(self, par, [])

    def elapsed(self, unit='hour'):
        elpsd = (self.times.jd - self.times[0].jd)
        if unit == 'hour':
            elpsd *= 24.0
        elif unit.startswith('m'):
            elpsd *= (24.0 * 60.0)
        elif unit.startswith('s'):
            elpsd *= (24.0 * 3600.0)
        return elpsd

    def to_Time(self, times='times', toffset=None):
        """
        Convert the array to Time and set as self.times.

        If str, then convert class 'times' variable.
        """
        if isinstance(times, str) and times == 'times':
            times = self.times
        self.times = to_Time(times, toffset)

    def to_Angle(self, ang, angle_unit=u.degree, angv=None):
        """
        Convert attribute ang to class array.

        """
        if angv is None:
            angv = getattr(self, ang)
        setattr(self, ang, to_Angle(angv, angle_unit))

    def interp(self, par, times):
        """
        Interpolate attribute par onto times.

        par : str
            string of the attribute to use
        times : Time
            times to be interpolated onto.
        """
        if self._E is None:
            self.reset_base()
        clpar = getattr(self._E, f"{par}")
        if not len(clpar):
            return clpar
        try:
            parunit = clpar.unit
        except AttributeError:
            return np.interp(times.jd, self._E.times.jd, clpar)
        return u.Quantity(np.interp(times.jd, self._E.times.jd, clpar), parunit)

    def reset_base(self):
        if self._E is None:
            self._E = Namespace()
        for par in self.param:
            setattr(self._E, par, getattr(self, par))

    def at(self, times=None):
        """
        Put all data at "times".  None to self._E data.
        """
        if self._E is None:  # Set ephem archive data
            self.reset_base()
        if times is None:
            for par in self.param:
                setattr(self, par, getattr(self._E, par))
        else:
            self.to_Time(times)
            for par in self.param:
                if par != 'times':
                    setattr(self, par, self.interp(par, times))

    def calc_dt(self):
        """
        Compute self.dt in sec as ndarray Quantity.
        """
        self.dt = [0.0]  # in seconds
        for i, _t in enumerate(self.times[1:]):
            self.dt.append((_t - self.times[i]).value * 3600.0 * 24.0)
        self.dt[0] = self.dt[1]
        self.dt = np.array(self.dt) * u.second

    def smooth_array(self, arr, smooth):
        try:
            anit = arr.unit
        except AttributeError:
            anit = None
        if smooth:
            if anit is not None:
                return u.Quantity(convolve(arr, Gaussian1DKernel(smooth), boundary='extend'), anit)
            else:
                return convolve(arr, Gaussian1DKernel(smooth), boundary='extend')
        else:
            return arr

    def vis(self, arr, val=0.0, horizon=0.0):
        """
        Get filter for visible above horizon.  Those below the horizon are set to val.

        Usage is to call e.g. visible_doppler = self.vis(self.doppler)
        """
        varr = np.array(arr)
        varr[np.where(self.el < horizon)] = val
        return varr

    def dbydt(self, par, smooth=None, unwrap=False):
        if '.' in par:
            ns = getattr(self, par.split('.')[0])
            par = par.split('.')[1]
        else:
            ns = self
        if self.dt is None or len(self.dt) != len(getattr(ns, par)):
            self.calc_dt()
        deriv = f"{par}dot"
        setattr(ns, deriv, [0.0])
        if unwrap:
            _param = np.unwrap(getattr(ns, par))
        else:
            _param = getattr(ns, par)
        _param = self.smooth_array(_param, smooth)
        for i, _pp in enumerate(_param[1:]):
            getattr(ns, deriv).append((_pp - _param[i])/self.dt[i+1])
        getattr(ns, deriv)[0] = getattr(ns, deriv)[1]
        setattr(ns, deriv, u.Quantity(getattr(ns, deriv)))
