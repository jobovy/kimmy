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
        return method(*args,**kwargs)
    return wrapper   
class OneZone(object):
    """OneZone: simple one-zone chemical evolution models"""
    def __init__(self,
                 eta=2.5,tau_SFE=1.*u.Gyr,tau_SFH=6.*u.Gyr,
                 tau_Ia=1.5*u.Gyr,min_dt_Ia=0.15*u.Gyr,
                 solar_O=8.69,solar_Fe=7.47,
                 mCC_O=0.015,mCC_Fe=0.0012,mIa_O=0.,mIa_Fe=0.0017,
                 r=0.4):
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
              tau_Ia= (1.5 Gyr) SNe Ia exponential decay time scale
              min_dt_Ia= (150 Myr) minimum time delay for SNe Ia
              mCC_O= (0.015) mass fraction of oxygen returned by core-collapse SNe (mass of O / stellar mass formed)
              mCC_Fe= (0.0012) mass fraction of iron returned by core-collapse SNe (mass of O / stellar mass formed)
              mIa_O= (0.) mass fraction of oxygen returned by SNe Ia (mass of O / stellar mass formed)
              mIa_Fe= (0.0017) mass fraction of iron returned by SNe Ia (mass of O / stellar mass formed)
              r= (0.4) mass recycling parameter (core-collapse SNe + AGB returns): amount of mass returned at abundances of star at birth
           The next parameters are fixed for the instance:
              solar_O= (8.69) solar oxygen number density on the x_O = 12 + log10(X_O/H) scale
              solar_Fe=7.47 solar iron number density on the x_O = 12 + log10(X_O/H) scale
        OUTPUT:
           instance
        HISTORY:
           2018-07-09 - Written - Bovy (UofT)
        """
        # Store everything internally
        self.eta= eta
        self.tau_SFE= tau_SFE
        self.tau_SFH= tau_SFH
        self.tau_Ia= tau_Ia
        self.min_dt_Ia= min_dt_Ia
        self.mCC_O= mCC_O
        self.mCC_Fe= mCC_Fe
        self.mIa_O= mIa_O
        self.mIa_Fe= mIa_Fe
        self.r= r
        # Set solar
        self._solar_O= solar_O
        self._solar_Fe= solar_Fe
        self._calc_solar()
        # Setup hash
        self._current_model_hash= None
        return None

    def _calc_solar(self):
        self._logZO_solar= -2.25+self._solar_O-8.69
        self._logZFe_solar= -2.93+self._solar_Fe-7.47
        return None

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

    # Abundances
    @_recalc_model
    def O_H(self,t):
        # CCSNe contribution
        ZO_t= self._ZO_CC_eq*(1.-numpy.exp(-t/self._tau_dep_SFH))
        # Ia contribution
        dt= t-self.min_dt_Ia
        idx= dt > 0.
        ZO_t[idx]+= self._ZO_Ia_eq*(1.-numpy.exp(-dt[idx]/self._tau_dep_SFH)
                                    -self._tau_dep_Ia/self._tau_dep_SFH
                                    *(numpy.exp(-dt[idx]/self._tau_Ia_SFH)
                                      -numpy.exp(-dt[idx]/self._tau_dep_SFH)))
        # DO WE NEED TO ADD HYDROGEN EVOLUTION AS WELL? SMALL EFFECT?
        return numpy.log10(ZO_t)-self._logZO_solar

    @_recalc_model
    def Fe_H(self,t):
        self._update_timescales()
        self._calc_equilibrium()
        # CCSNe contribution
        ZFe_t= self._ZFe_CC_eq*(1.-numpy.exp(-t/self._tau_dep_SFH))
        # Ia contribution
        dt= t-self.min_dt_Ia
        idx= dt > 0.
        ZFe_t[idx]+= self._ZFe_Ia_eq*(1.-numpy.exp(-dt[idx]/self._tau_dep_SFH)
                                      -self._tau_dep_Ia/self._tau_dep_SFH
                                      *(numpy.exp(-dt[idx]/self._tau_Ia_SFH)
                                       -numpy.exp(-dt[idx]/self._tau_dep_SFH)))
        # DO WE NEED TO ADD HYDROGEN EVOLUTION AS WELL? SMALL EFFECT?
        return numpy.log10(ZFe_t)-self._logZFe_solar

    @_recalc_model
    def O_Fe(self,t):
        return self.O_H(t)-self.Fe_H(t)

