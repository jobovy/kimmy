# OneZone.py: simple one-zone chemical evolution models
from functools import wraps
import numpy
import hashlib
from astropy import units as u
def _recalc_model(method):
    @wraps(method)
    def wrapper(*args,**kwargs):
        new_model_hash= args[0]._model_hash()
        if new_model_hash != args[0]._current_model_hash:
            args[0]._update_timescales()
            args[0]._calc_equilibrium()
            args[0]._current_model_hash= new_model_hash
        new_solar_hash= args[0]._solar_hash()
        if new_solar_hash != args[0]._current_solar_hash:
            args[0]._calc_solar()
            args[0]._current_solar_hash= new_solar_hash
        return method(*args,**kwargs)
    return wrapper   
_defaults= {'eta':2.5,
            'tau_SFE':   1.*u.Gyr,
            'tau_SFH':   6.*u.Gyr,
            'tau_Ia':    1.5*u.Gyr,
            'min_dt_Ia': 0.15*u.Gyr,
            'sfh':       'exp',
            'mCC_O':     0.015,
            'mCC_Fe':    0.0012,
            'mIa_O':     0.,
            'mIa_Fe':    0.0017,
            'r':         0.4,
            'solar_O':   8.69,
            'solar_Fe':  7.47}
class OneZone(object):
    """OneZone: simple one-zone chemical evolution models"""
    def __init__(self,**kwargs):
        """
        NAME:
           __init__
        PURPOSE:
           Setup a simple OneZone chemical-evolution object using the analytical formalism from Weinberg, Andrews, & Freudenberg (2017)
        INPUT:
           The following can be updated on-the-fly:
              eta= (2.5) outflow mass loading factor as a fraction of SFR (float)
              tau_SFE= (1 Gyr) star-formation efficiency time scale (Quantity with units of time)
              tau_SFH= (6 Gyr) star-formation exponential decay time scale
              sfh= ('exp') use an exponential ['exp'] or linear-exponential ['linexp'] star-formation history
              tau_Ia= (1.5 Gyr) SNe Ia exponential decay time scale
              min_dt_Ia= (150 Myr) minimum time delay for SNe Ia
              mCC_O= (0.015) mass fraction of oxygen returned by core-collapse SNe (mass of O / stellar mass formed)
              mCC_Fe= (0.0012) mass fraction of iron returned by core-collapse SNe (mass of O / stellar mass formed)
              mIa_O= (0.) mass fraction of oxygen returned by SNe Ia (mass of O / stellar mass formed)
              mIa_Fe= (0.0017) mass fraction of iron returned by SNe Ia (mass of O / stellar mass formed)
              r= (0.4) mass recycling parameter (core-collapse SNe + AGB returns): amount of mass returned at abundances of star at birth
              solar_O= (8.69) solar oxygen number density on the x_O = 12 + log10(X_O/H) scale
              solar_Fe=7.47 solar iron number density on the x_O = 12 + log10(X_O/H) scale
        OUTPUT:
           instance
        HISTORY:
           2018-07-09 - Written - Bovy (UofT)
        """
        self._initialize_params(**kwargs)
        # Setup hash for storing models
        self._current_model_hash= None
        self._current_solar_hash= None
        return None

    def _initialize_params(self,**kwargs):
        self._init_params= {} # To store initial
        for key in _defaults.keys():
            self.__dict__[key]= kwargs.get(key,_defaults[key])
            self._init_params[key]= kwargs.get(key,_defaults[key])
        return None

    def initial(self):
        for key in _defaults.keys():
            self.__dict__[key]= self._init_params[key]
        return None

    def default(self):
        for key in _defaults.keys():
            self.__dict__[key]= _defaults[key]
        return None

    def __str__(self):
        out= ''
        for key in sorted(_defaults.keys()):
            out+= '{0:<10}:\t{1}\n'.format(key,self.__dict__[key])
        return out[:-1]

    def _calc_solar(self):
        self._logZO_solar= -2.25+self.solar_O-8.69
        self._logZFe_solar= -2.93+self.solar_Fe-7.47
        return None

    # Equilibrium and model parameters
    def _update_timescales(self):
        # Update all relevant timescales for the model based on the current
        # model parameters
        self._tau_dep= self.tau_SFE/(1.+self.eta-self.r)
        self._tau_dep_SFH= 1./(1./self._tau_dep-1./self.tau_SFH)
        self._tau_dep_Ia= 1./(1./self._tau_dep-1./self.tau_Ia)
        self._tau_Ia_SFH= 1./(1./self.tau_Ia-1./self.tau_SFH)
        return None

    def _calc_equilibrium(self):
        self._ZO_CC_eq= self.mCC_O*self._tau_dep_SFH/self.tau_SFE
        self._ZO_Ia_eq= self.mIa_O*self._tau_dep_SFH/self.tau_SFE\
            *self._tau_Ia_SFH/self.tau_Ia\
            *numpy.exp(self.min_dt_Ia/self.tau_SFH)
        self._ZFe_CC_eq= self.mCC_Fe*self._tau_dep_SFH/self.tau_SFE
        self._ZFe_Ia_eq= self.mIa_Fe*self._tau_dep_SFH/self.tau_SFE\
            *self._tau_Ia_SFH/self.tau_Ia\
            *numpy.exp(self.min_dt_Ia/self.tau_SFH)
        return None

    # Time evolution equations
    def _evol_CC(self,t):
        if self.sfh.lower() == 'exp':
            return (1.-numpy.exp(-t/self._tau_dep_SFH))
        else:
            return (1.-self._tau_dep_SFH/t
                       *(1.-numpy.exp(-t/self._tau_dep_SFH)))
    
    def _evol_Ia(self,t):
        # Ia contribution
        dt= t-self.min_dt_Ia
        idx= dt > 0.
        out= numpy.zeros(t.shape)
        if self.sfh.lower() == 'exp':
            out[idx]+= \
                (1.-numpy.exp(-dt[idx]/self._tau_dep_SFH)
                                -self._tau_dep_Ia/self._tau_dep_SFH
                                *(numpy.exp(-dt[idx]/self._tau_Ia_SFH)
                                  -numpy.exp(-dt[idx]/self._tau_dep_SFH)))\
                                  .to(u.dimensionless_unscaled).value
        else:
            out[idx]+= \
                (self._tau_Ia_SFH/t[idx]\
                *(dt[idx]/self._tau_Ia_SFH+self._tau_dep_Ia/self._tau_dep_SFH
                                         *numpy.exp(-dt[idx]/self._tau_Ia_SFH)
                  +(1.+self._tau_dep_SFH/self._tau_Ia_SFH
                    -self._tau_dep_Ia/self._tau_dep_SFH)\
                      *numpy.exp(-dt[idx]/self._tau_dep_SFH)
                  -(1.+self._tau_dep_SFH/self._tau_Ia_SFH)))\
                  .to(u.dimensionless_unscaled).value
        return out

    # Abundances
    @_recalc_model
    def O_H(self,t):
        # CCSNe contribution
        ZO_t= self._ZO_CC_eq*self._evol_CC(t)
        # Ia contribution
        ZO_t+= self._ZO_Ia_eq*self._evol_Ia(t)
        # DO WE NEED TO ADD HYDROGEN EVOLUTION AS WELL? SMALL EFFECT?
        return numpy.log10(ZO_t)-self._logZO_solar

    @_recalc_model
    def Fe_H(self,t):
        # CCSNe contribution
        ZFe_t= self._ZFe_CC_eq*self._evol_CC(t)
        # Ia contribution
        ZFe_t+= self._ZFe_Ia_eq*self._evol_Ia(t)
        # DO WE NEED TO ADD HYDROGEN EVOLUTION AS WELL? SMALL EFFECT?
        return numpy.log10(ZFe_t)-self._logZFe_solar

    def O_Fe(self,t):
        return self.O_H(t)-self.Fe_H(t)

    def _model_hash(self):
        return hashlib.md5(numpy.array([self.eta,
                                        self.tau_SFE.to(u.Gyr).value,
                                        self.tau_SFH.to(u.Gyr).value,
                                        self.tau_Ia.to(u.Gyr).value,
                                        self.min_dt_Ia.to(u.Gyr).value,
                                        self.mCC_O,
                                        self.mCC_Fe,
                                        self.mIa_O,
                                        self.mIa_Fe,
                                        self.r])).hexdigest()

    def _solar_hash(self):
        return hashlib.md5(numpy.array([self.solar_O,
                                        self.solar_Fe])).hexdigest()

