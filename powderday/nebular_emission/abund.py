from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as InterpUS
from powderday.nebular_emission.cloudy_tools import sym_to_name

"""
------------------------------------------------------------------------------------------
From cloudyfsps written by Nell Byler.
(Source https://github.com/nell-byler/cloudyfsps/blob/master/cloudyfsps/nebAbundTools.py 
retrieved in October 2019)
------------------------------------------------------------------------------------------
"""


def getNebAbunds(set_name, logZ, dust=True, re_z=False, **kwargs):
    """
    neb_abund.get_abunds(set_name, logZ, dust=True, re_z=False)
    set_name must be 'dopita', 'newdopita', 'cl01' or 'yeh'
    """
    allowed_names = ['dopita', 'newdopita', 'cl01', 'yeh',
                     'varyNO', 'gutkin', 'UVbyler', 'varyCO']
    if set_name in allowed_names:
        return eval('{}({}, dust={}, re_z={})'.format(set_name, logZ, dust, re_z))
    else:
        raise IOError(allowed_names)


class abundSet(object):
    def __init__(self, set_name, logZ):
        """
        overarching class for abundance sets.
        abundSet('dopita', 0.0)
        """
        self.logZ = logZ
        self.abund_0 = load_abund(set_name)
        self.depl = load_depl(set_name)
        self.calcSpecial()
        self.calcFinal()
        self.inputStrings()

    def calcSpecial(self):
        return

    def calcFinal(self):
        return

    def inputStrings(self):
        self.solarstr = 'abundances {} {}'.format(self.solar, self.grains)
        elem_strs = []
        names = sym_to_name()
        for key in self.abund_0.keys():
            elm = names[key]
            abund = self.__getattribute__(key)
            # if hasattr(self, 're_z'):
            #    if key != 'He':
            #        abund -= self.re_z
            outstr = 'element abundance {0} {1:.2f} log'.format(elm, abund)
            elem_strs.append(outstr)
        self.__setattr__('elem_strs', elem_strs)
        return


class dopita(abundSet):
    solar = 'old solar 84'

    def __init__(self, logZ, dust=True, re_z=False):
        """
        Dopita+2001: old solar abundances = 0.019
        ISM grains
        """
        if dust:
            self.grains = 'no grains\ngrains ISM'
        else:
            self.grains = 'no grains'
        if re_z:
            self.re_z = logZ
        else:
            self.re_z = 0.0
        abundSet.__init__(self, 'dopita', logZ)

    def calcSpecial(self):
        """
        piece-wise function for nitrogen abund (step-function)
        functional form for helium
        """

        def calc_N(logZ):
            if logZ <= -0.63:
                return -4.57 + logZ
            else:
                return -3.94 + (2.0 * logZ)

        def calc_He(logZ):
            return np.log10(0.08096 + (0.02618 * (10.0 ** logZ)))

        self.__setattr__('He', calc_He(self.logZ))
        self.__setattr__('N', calc_N(self.logZ) + self.depl['N'])
        return

    def calcFinal(self):
        """
        apply depletions and scale with logZ
        """
        [self.__setattr__(key, val + self.logZ + self.depl[key])
         for key, val in self.abund_0.items() if not hasattr(self, key)]
        return


class newdopita(abundSet):
    solar = 'GASS10'

    def __init__(self, logZ, dust=True, re_z=False):
        """
        Abundances from Dopita (2013)
            Solar Abundances from Grevasse 2010 - z= 0.013
            includes smooth polynomial for N/O, C/O relationship
            functional form for He(z)
            new depletion factors
            ISM grains
        """
        if dust:
            self.grains = 'no grains\ngrains ISM'
        else:
            self.grains = 'no grains'
        self.re_z = re_z
        abundSet.__init__(self, 'newdopita', logZ)

    def calcSpecial(self):
        def calc_He(logZ):
            return np.log10(0.0737 + (0.024 * (10.0 ** logZ)))

        def calc_CNO(logZ):
            oxy = np.array([7.39, 7.50, 7.69, 7.99, 8.17,
                            8.39, 8.69, 8.80, 8.99, 9.17, 9.39])
            nit = np.array([-6.61, -6.47, -6.23, -5.79, -5.51,
                            -5.14, -4.60, -4.40, -4.04, -3.67, -3.17])
            car = np.array([-5.58, -5.44, -5.20, -4.76, -4.48,
                            -4.11, -3.57, -3.37, -3.01, -2.64, -2.14])
            O = self.abund_0['O'] + logZ
            C = float(InterpUS(oxy, car, k=1)(O + 12.0))
            N = float(InterpUS(oxy, nit, k=1)(O + 12.0))
            return C, N, O

        self.__setattr__('He', calc_He(self.logZ))
        C, N, O = calc_CNO(self.logZ)
        [self.__setattr__(key, val + self.depl[key])
         for key, val in zip(['C', 'N', 'O'], [C, N, O])]
        return

    def calcFinal(self):
        [self.__setattr__(key, val + self.logZ + self.depl[key])
         for key, val in self.abund_0.items() if not hasattr(self, key)]
        return


class UVbyler(abundSet):
    solar = 'GASS10'

    def __init__(self, logZ, dust=True, re_z=False):
        """
        Abundances from Dopita (2013)
            Solar Abundances from Grevasse 2010 - z= 0.013
            New fit for N/O, C/O relationship
            functional form for He(z)
            new depletion factors
            ISM grains
        """
        if dust:
            self.grains = 'no grains\ngrains ISM'
        else:
            self.grains = 'no grains'
        self.re_z = re_z
        abundSet.__init__(self, 'UVbyler', logZ)

    def calcSpecial(self):
        def calc_He(logZ):
            return np.log10(0.0737 + (0.024 * (10.0 ** logZ)))

        def calc_CNO(logZ):
            O = self.abund_0['O'] + logZ
            # C = np.log10((1.0*10.**O)*(10.**-1.1 + 10.**(2.96 + O)))
            C = np.log10((10. ** O) * (10. ** -0.7 + 10. ** (4.8 + 1.45 * O)))
            # N = np.log10((1.0*10.**O)*(10.**-1.8 + 10.**(2.2 + O)))
            # N = np.log10((10.**O)*(10.**-1.5 + 10.**(2.5 + 1.2*O)))
            N = np.log10((1.0 * 10. ** O) * (10. ** -1.55 + 10. ** (2.3 + 1.1 * O)))
            # N  = -4.81 + logZ if logZ <= -0.3 else -4.51 + 2.0*logZ
            return C, N, O

        self.__setattr__('He', calc_He(self.logZ))
        C, N, O = calc_CNO(self.logZ)
        [self.__setattr__(key, val + self.depl[key])
         for key, val in zip(['C', 'N', 'O'], [C, N, O])]
        return

    def calcFinal(self):
        [self.__setattr__(key, val + self.logZ + self.depl[key])
         for key, val in self.abund_0.items() if not hasattr(self, key)]
        return


class gutkin(abundSet):
    solar = 'GASS10'

    def __init__(self, logZ, dust=True, re_z=False):
        """
        Gutkin+2016
            PARSEC metallicity (Bressan+2012)
            based on Grevesse+Sauvel (1998) and Caffau+2011
        """
        if dust:
            self.grains = 'no grains\ngrains ISM'
        else:
            self.grains = 'no grains'
        self.re_z = re_z
        abundSet.__init__(self, 'gutkin', logZ)

    def calcSpecial(self):
        def calc_He(logZ):
            Z = (10. ** logZ) * 0.01524
            Y = 0.2485 + 1.7756 * Z
            X = 1. - Y - Z
            return np.log10(Y / X / 4.)

        def calc_CNO(logZ):
            O = self.abund_0['O'] + logZ
            N = np.log10((0.41 * 10. ** O) * (10. ** -1.6 + 10. ** (2.33 + O)))
            C = self.abund_0['C'] + logZ
            return C, N, O

        self.__setattr__('He', calc_He(self.logZ))
        C, N, O = calc_CNO(self.logZ)
        [self.__setattr__(key, val)
         for key, val in zip(['C', 'N', 'O'], [C, N, O])]
        return

    def calcFinal(self):
        [self.__setattr__(key, val)
         for key, val in self.abund_0.items() if not hasattr(self, key)]
        return


def load_abund(set_name):
    if set_name == 'dopita':
        adict = dict(He=-1.01,
                     C=-3.44,
                     N=-3.95,
                     O=-3.07,
                     Ne=-3.91,
                     Mg=-4.42,
                     Si=-4.45,
                     S=-4.79,
                     Ar=-5.44,
                     Ca=-5.64,
                     Fe=-4.33,
                     F=-7.52,
                     Na=-5.69,
                     Al=-5.53,
                     P=-6.43,
                     Cl=-6.73,
                     K=-6.87,
                     Ti=-6.96,
                     Cr=-6.32,
                     Mn=-6.47,
                     Co=-7.08,
                     Ni=-5.75,
                     Cu=-7.73,
                     Zn=-7.34)
    elif set_name == 'newdopita':
        adict = dict(He=-1.01,
                     C=-3.57,
                     N=-4.60,
                     O=-3.31,
                     Ne=-4.07,
                     Na=-5.75,
                     Mg=-4.40,
                     Al=-5.55,
                     Si=-4.49,
                     S=-4.86,
                     Cl=-6.63,
                     Ar=-5.60,
                     Ca=-5.66,
                     Fe=-4.50,
                     Ni=-5.78,
                     F=-7.44,
                     P=-6.59,
                     K=-6.97,
                     Cr=-6.36,
                     Ti=-7.05,
                     Mn=-6.57,
                     Co=-7.01,
                     Cu=-7.81,
                     Zn=-7.44)
    elif set_name == 'UVbyler':
        adict = dict(He=-1.01,
                     C=-3.57,
                     N=-4.17,
                     O=-3.31,
                     Ne=-4.07,
                     Na=-5.75,
                     Mg=-4.40,
                     Al=-5.55,
                     Si=-4.49,
                     S=-4.86,
                     Cl=-6.63,
                     Ar=-5.60,
                     Ca=-5.66,
                     Fe=-4.50,
                     Ni=-5.78,
                     F=-7.44,
                     P=-6.59,
                     K=-6.97,
                     Cr=-6.36,
                     Ti=-7.05,
                     Mn=-6.57,
                     Co=-7.01,
                     Cu=-7.81,
                     Zn=-7.44)
    elif set_name == 'gutkin':
        adict = dict(He=-1.01,
                     C=-3.53,
                     N=-4.32,
                     O=-3.17,
                     F=-7.47,
                     Ne=-4.01,
                     Na=-5.70,
                     Mg=-4.45,
                     Al=-5.56,
                     Si=-4.48,
                     P=-6.57,
                     S=-4.87,
                     Cl=-6.53,
                     Ar=-5.63,
                     K=-6.92,
                     Ca=-5.67,
                     Sc=-8.86,
                     Ti=-7.01,
                     V=-8.03,
                     Cr=-6.36,
                     Mn=-6.64,
                     Fe=-4.51,
                     Co=-7.11,
                     Ni=-5.78,
                     Cu=-7.82,
                     Zn=-7.43)
    return adict


def load_depl(set_name):
    if set_name == 'dopita':
        ddict = dict(C=-0.30,
                     N=-0.22,
                     O=-0.22,
                     Ne=0.0,
                     Mg=-0.70,
                     Si=-1.0,
                     S=0.0,
                     Ar=0.0,
                     Ca=-2.52,
                     Fe=-2.0,
                     F=0.0,
                     Na=0.0,
                     Al=0.0,
                     P=0.0,
                     Cl=0.0,
                     K=0.0,
                     Ti=0.0,
                     Cr=0.0,
                     Mn=0.0,
                     Co=0.0,
                     Ni=0.0,
                     Cu=0.0,
                     Zn=0.0)
    elif set_name == 'newdopita':
        ddict = dict(He=0.00,
                     C=-0.30,
                     N=-0.05,
                     O=-0.07,
                     Ne=0.00,
                     Na=-1.00,
                     Mg=-1.08,
                     Al=-1.39,
                     Si=-0.81,
                     S=0.00,
                     Cl=-1.00,
                     Ar=0.00,
                     Ca=-2.52,
                     Fe=-1.31,
                     Ni=-2.00,
                     F=0.0,
                     P=0.0,
                     K=0.0,
                     Cr=0.0,
                     Ti=0.0,
                     Mn=0.0,
                     Co=0.0,
                     Cu=0.0,
                     Zn=0.0)
    elif set_name == 'UVbyler':
        ddict = dict(He=0.00,
                     C=-0.30,
                     N=-0.05,
                     O=-0.07,
                     Ne=0.00,
                     Na=-1.00,
                     Mg=-1.08,
                     Al=-1.39,
                     Si=-0.81,
                     S=0.00,
                     Cl=-1.00,
                     Ar=0.00,
                     Ca=-2.52,
                     Fe=-1.31,
                     Ni=-2.00,
                     F=0.0,
                     P=0.0,
                     K=0.0,
                     Cr=0.0,
                     Ti=0.0,
                     Mn=0.0,
                     Co=0.0,
                     Cu=0.0,
                     Zn=0.0)
    elif set_name == 'gutkin':
        ddict = dict(He=0.00,
                     Li=-0.8,
                     C=-0.30,
                     O=-0.15,
                     Na=-0.60,
                     Mg=-0.70,
                     Al=-1.70,
                     Si=-1.00,
                     Cl=-0.30,
                     Ca=-2.52,
                     Fe=-2.00,
                     Ni=-1.40)
    return ddict
