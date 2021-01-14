#! /usr/bin/python
#-*- coding: utf-8 -*-

import math
## Last updated: August 17, 2015


def md_pwl():
    """ Compute the post-seismic deformation/correction "d" using parametric
      model PWL (Piece-Wise Linear Function)

      Returns
      ----------------------------------------------------------------------
      float
            post-seismic correction in mm (at dtq years after the earthquake)

  """
    return 0e0


def md_log(dtq, a1, t1):
    """ Compute the post-seismic deformation/correction "d" using parametric
      model Logarithmic Function

      Parameters
      ----------------------------------------------------------------------
      dtq : float 
            time difference (t-t_Earthquake) in decimal year (but see note below)
      a1: float
          amplitude 1 of the parametric model, in mm
      t1: float
          relaxation time 1, in decimal years

      Returns
      ----------------------------------------------------------------------
      float
            post-seismic correction in mm (at dtq years after the earthquake)
      
      Notes
      ----------------------------------------------------------------------
      Time unit is decimal year. It is advised to compute "dtq" by:
      (MJD - MJD_Earthquake)/365.25 where MJD is the modified julian day.

  """
    return a1 * math.log(1e0 + dtq / t1)


def md_exp(dtq, a1, t1):
    """ Compute the post-seismic deformation/correction "d" using parametric
      model Exponential Function

      Parameters
      ----------------------------------------------------------------------
      dtq : float 
            time difference (t-t_Earthquake) in decimal year (but see note below)
      a1: float
          amplitude 1 of the parametric model, in mm
      t1: float
          relaxation time 1, in decimal years

      Returns
      ----------------------------------------------------------------------
      float
            post-seismic correction in mm (at dtq years after the earthquake)
      
      Notes
      ----------------------------------------------------------------------
      Time unit is decimal year. It is advised to compute "dtq" by:
      (MJD - MJD_Earthquake)/365.25 where MJD is the modified julian day.

  """
    te1 = dtq / t1
    return a1 * (1e0 - math.exp(-te1))


def md_logexp(dtq, a1, t1, a2, t2):
    """ Compute the post-seismic deformation/correction "d" using parametric
      model Logarithmic and Exponential Function

      Parameters
      ----------------------------------------------------------------------
      dtq : float 
            time difference (t-t_Earthquake) in decimal year (but see note below)
      a1: float
          amplitude 1 of the parametric model, in mm
      a2: float
          amplitude 2 of the parametric model, in mm
      t1: float
          relaxation time 1, in decimal years
      t2: float
          relaxtaion time 2, in decimal years

      Returns
      ----------------------------------------------------------------------
      float
            post-seismic correction in mm (at dtq years after the earthquake)
      
      Notes
      ----------------------------------------------------------------------
      Time unit is decimal year. It is advised to compute "dtq" by:
      (MJD - MJD_Earthquake)/365.25 where MJD is the modified julian day.

  """
    te2 = dtq / t2
    return a1 * math.log(1e0 + dtq / t1) + a2 * (1e0 - math.exp(-te2))


def md_expexp(dtq, a1, t1, a2, t2):
    """ Compute the post-seismic deformation/correction "d" using parametric
      model Two Exponential Functions

      Parameters
      ----------------------------------------------------------------------
      dtq : float 
            time difference (t-t_Earthquake) in decimal year (but see note below)
      a1: float
          amplitude 1 of the parametric model, in mm
      a2: float
          amplitude 2 of the parametric model, in mm
      t1: float
          relaxation time 1, in decimal years
      t2: float
          relaxtaion time 2, in decimal years

      Returns
      ----------------------------------------------------------------------
      float
            post-seismic correction in mm (at dtq years after the earthquake)
      
      Notes
      ----------------------------------------------------------------------
      Time unit is decimal year. It is advised to compute "dtq" by:
      (MJD - MJD_Earthquake)/365.25 where MJD is the modified julian day.

  """
    te1 = dtq / t1
    te2 = dtq / t2
    return a1 * (1e0 - math.exp(-te1)) + a2 * (1e0 - math.exp(-te2))


def parametric(model='pwl', *args):
    """ Compute the post-seismic deformation/correction "d" using the 
      parametric model specified by the (input) variable 'model'.

      Parameters
      ----------------------------------------------------------------------
      model: string or int
          can be any of the following:
          - 'pwl' or 0 to denote a Piece-Wise Linear Model
          - 'log' or 1 to denote a Logarithmic Model
          - 'exp' or 2 to denote an Exponential Model
          - 'logexp' or 3 to denote a Logarithmic + Exponential Model
          - 'expexp' or 4 to denote a Two_part Exponential Model
      dtq : float 
            time difference (t-t_Earthquake) in decimal year (but see note below)
            Not needed in the case where model is 0.
      a1: float
          amplitude 1 of the parametric model, in mm
          Not needed in the case where model is 0.
      t1: float
          relaxation time 1, in decimal years
          Not needed in the case where model is 0.
      a2: float
          amplitude 2 of the parametric model, in mm
          Not needed in the case where model is 0, 1 or 2
      t2: float
          relaxtaion time 2, in decimal years
          Not needed in the case where model is 0, 1 or 2

      Returns
      ----------------------------------------------------------------------
      float
            post-seismic correction in mm (at dtq years after the earthquake)
      
      Notes
      ----------------------------------------------------------------------
      Time unit is decimal year. It is advised to compute "dtq" by:
      (MJD - MJD_Earthquake)/365.25 where MJD is the modified julian day.

  """
    model_dict = {
        'pwl': md_pwl,
        'log': md_log,
        'exp': md_exp,
        'logexp': md_logexp,
        'expexp': md_expexp
    }

    int_model_dict = {0: 'pwl', 1: 'log', 2: 'exp', 3: 'logexp', 4: 'expexp'}
    try:
        model + 1
        return (model_dict[int_model_dict[model]](
            *args) if model > 0 else model_dict[int_model_dict[model]]())
    except:
        return (model_dict[model](
            *args) if model != 'pwl' else model_dict[model]())
