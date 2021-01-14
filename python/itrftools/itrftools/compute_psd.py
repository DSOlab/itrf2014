#! /usr/bin/python
#-*- coding: utf-8 -*-

from __future__ import print_function
import sys, datetime, math
import numpy as np
from itrftools import parametric


def time_str2dt(time_str):
    """ Resolve a datetime string of type: YY:DDD:SSSSS to a python datetime
      instance.

      Parameters:
      -------------------
      time_str : string
                 A datetime string of type YY:DDD:SSSSS, e.g. '09:280:80331'

      Returns
      -------------------
      datetime
                 A datetime instance
  """
    cyr, doy, isec = map(int, time_str.split(':'))
    yr = cyr + 2000 if datetime.datetime.now().year > cyr + 2000 else cyr + 1900
    dt = (datetime.datetime(yr, 1, 1) + datetime.timedelta(doy - 1) +
          datetime.timedelta(seconds=isec))
    return dt


def get_psd_model(line, cmp):
    """ This function extracts the PSD model and parameters off from a line 
      of a ITRF-psd-*.dat file.
      Such lines have the form:
      +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ANTC  A 41713S001 10:058:23656 E 3 -192.03  0.5969  -72.74  0.0799     GPS
                                     N 3   61.57  2.1357   26.26  0.2294
                                     U 4  157.62  3.3132   25.61  0.1854
      +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      The function disregards everything before column 32 and after column 72
      (hence the technique is not read). We are only interested in the
      component, model type (int) and model parameters.

      Parameters:
      ----------------
      line: string
            A line of a PSD station record, i.e. a line describing the
            PSD model for the North, East or Up component.
      cmp: character
           The component to be read. The function will assert that cmp is the
           same as the component read at column 32.

      Returns:
      ----------------
      tuple
          A 2-element tuple where the first element is an int [0,4] describing
          the model type and the second element is a list containing the model
          parameters; this second element has variable size depending on the
          model type.

  """
    assert (line[32] == cmp)
    model = int(line[34])
    assert (model >= 0 and model < 5)
    params = map(float, line[36:72].split())
    return model, params


def get_next_psd(fin, line):
    """ This function extracts the PSD model and parameters for a station, off
      from an ITRF-psd-*.dat file. It uses the function get_psd_model() to
      extract the relevant fields per component.
      
      Parameters:
      ----------------
      line: string
            The first line of a PSD station record, i.e. a line describing the
            PSD model for the East component.
      fin: input file stream
            The (opened) file from which we are reading from. This is normaly
            a file of type ITRF2014-psd-*.dat; The function will assert that
            the first line passed in is for the East component and after
            resolving it, will try to read and resolve two more lines for the
            North and Up components respectively.

      Returns:
      ----------------
      tuple
          The elements of the returned tuple are:
          [0] : (string) Station name (4-char id)
          [1] : (string) Station domes number
          [2] : (datetime.datetime) Datatime (instance) of the earthquake
          [3] : (int) Model id for the PSD of the East component
          [4] : (list of floats) Parameters for the PSD of the East component
          [5] : (int) Model id for the PSD of the North component
          [6] : (list of floats) Parameters for the PSD of the North component
          [7] : (int) Model id for the PSD of the Up component
          [8] : (list of floats) Parameters for the PSD of the Up component

  """
    sta_name = line[1:5]
    char_id = line[7]
    domes = line[9:18]
    dtime = time_str2dt(line[19:31])
    modele, parame = get_psd_model(line, 'E')
    line = fin.readline()
    modeln, paramn = get_psd_model(line, 'N')
    line = fin.readline()
    modelu, paramu = get_psd_model(line, 'U')
    return sta_name, domes, dtime, modele, parame, modeln, paramn, modelu, paramu


def compute_psd(psd_file, t=datetime.datetime.now(), station=None, domes=None):
    """ Given an (ITRF) .PSD file, aka a file containing ITRF-like post seismic
      deformation parametrs, compute the PSD correction per [e,n,u] component
      for a given station at a given time t. The station can be described by
      either passing in its 4-char id, or its DOMES number.
      If the station (or domes) is not found in the file, then the tuple
      <name_give, domes_given, 0, 0, 0> is returned.

      Parameters:
      -----------
      psd_file: string
          The PSD (.dat) file to extract the PSD corrections from. (see e.g.
          ftp://itrf.ign.fr/pub/itrf/itrf2014/ITRF2014-psd-gnss.dat
      t: datetime.datetime
          The time we want the PSD at.
      station: string
          The name of the station (4-char id)
      domes:
          The domes number.

      Returns:
      --------
      tuple (of size 5)
          First two elements are the station name and domes as written in the
          PSD file if the station is found; else they are the values provided.
          The next 3 tuple elements are the (total) PSD in [e,n,u] components
          respectively in mm.

      Note:
      -----
      To specify the station, you must use either the 'station' parameter, or
      the 'domes' parameter, but **not** both.
  """
    de = dn = du = 0e0
    num_of_psd = 0
    sta_def = ' ' * 4 if not station else station
    dms_def = ' ' * 9 if not domes else domes
    found = False
    assert (not station or not domes)
    if station:
        station = station.upper()
    with open(psd_file) as fin:
        line = fin.readline()
        while line:
            sta, dms, et, me, pe, mn, pn, mu, pu = get_next_psd(fin, line)
            if (not station and domes == dms) or (not domes and station == sta):
                dt = t - et
                dyr = (dt.days + dt.seconds / 86400e0) / 365.25e0
                de += parametric.parametric(me, dyr, *pe)
                dn += parametric.parametric(mn, dyr, *pn)
                du += parametric.parametric(mu, dyr, *pu)
                num_of_psd += 1
                found = True
                sta_def = sta
                dms_def = dms
            line = fin.readline()
        # print 'Number of individuals PSDs applied: {}'.format(num_of_psd)
        return sta_def, dms_def, de, dn, du


def xyz2llh(x, y, z, a=6378137e0, f=0.003352810681183637418):
    """ Cartesian to ellispoidal coordinates, based on
      Transformation from Cartesian to geodetic coordinates accelerated by 
      Halley's method, J. Geodesy (2006), 79(12): 689-693
  """
    # Functions of ellipsoid parameters.
    aeps2 = a * a * 1e-32
    e2 = (2.0e0 - f) * f
    e4t = e2 * e2 * 1.5e0
    ep2 = 1.0e0 - e2
    ep = math.sqrt(ep2)
    aep = a * ep
    # Compute Coefficients of (Modified) Quartic Equation
    # Remark: Coefficients are rescaled by dividing by 'a'
    # Compute distance from polar axis squared.
    p2 = x * x + y * y
    # Compute longitude lambda.
    if p2:
        lon = math.atan2(y, x)
    else:
        lon = .0e0
    # Ensure that Z-coordinate is unsigned.
    absz = abs(z)
    if p2 > aeps2:  # Continue unless at the poles
        # Compute distance from polar axis.
        p = math.sqrt(p2)
        # Normalize.
        s0 = absz / a
        pn = p / a
        zp = ep * s0
        # Prepare Newton correction factors.
        c0 = ep * pn
        c02 = c0 * c0
        c03 = c02 * c0
        s02 = s0 * s0
        s03 = s02 * s0
        a02 = c02 + s02
        a0 = math.sqrt(a02)
        a03 = a02 * a0
        d0 = zp * a03 + e2 * s03
        f0 = pn * a03 - e2 * c03
        # Prepare Halley correction factor.
        b0 = e4t * s02 * c02 * pn * (a0 - ep)
        s1 = d0 * f0 - b0 * s0
        cp = ep * (f0 * f0 - b0 * c0)
        # Evaluate latitude and height.
        phi = math.atan(s1 / cp)
        s12 = s1 * s1
        cp2 = cp * cp
        h = (p * cp + absz * s1 -
             a * math.sqrt(ep2 * s12 + cp2)) / math.sqrt(s12 + cp2)
    else:  # // Special case: pole.
        phi = math.pi / 2e0
        h = absz - aep
    # Restore sign of latitude.
    if z < 0.e0:
        phi = -phi
    return phi, lon, h


def enu2xyz(e, n, u, x, y, z):
    """ Transform a [e,n,u] vector (i.e.  local East, North, Up coordinates)
      to cartesian [X,Y,Z].
      Reference: http://www.navipedia.net/index.php/Transformations_between_ECEF_and_ENU_coordinates
      Parameters:
      -----------
      Returns:
      -----------
      list (of floats witsh size = 3)
          the list is: [x, y, z]
  """
    lat, lon, hgt = xyz2llh(x, y, z)
    sl = np.sin(lon)
    cl = np.cos(lon)
    sf = np.sin(lat)
    cf = np.cos(lat)
    R = np.matrix([[-sl, -cl * sf, cl * cf], [cl, -sl * sf, sl * cf],
                   [0e0, cf, sf]])
    enu = np.matrix([[e], [n], [u]])
    return [item for sublist in (R * enu).tolist() for item in sublist]


## Example usage
##if __name__ == "__main__":
##    de, dn, du = compute_psd('ITRF2014-psd-gnss.dat', t=datetime.datetime.now(), station='ANKR')
##    print('PSD correction in [e,n,u] = [{}, {}, {}]'.format(de, dn, du))
##    de, dn, du = compute_psd('ITRF2014-psd-gnss.dat', t=datetime.datetime.now(), domes='20805M002')
##    print('PSD correction in [e,n,u] = [{}, {}, {}]'.format(de, dn, du))
##    de, dn, du = compute_psd('ITRF2014-psd-gnss.dat', t=datetime.datetime.now(), domes='11401M001')
##    print('PSD correction in [e,n,u] = [{}, {}, {}]'.format(de, dn, du))
##    de, dn, du = compute_psd('ITRF2014-psd-gnss.dat', t=datetime.datetime.now(), station='COCO')
##    print('PSD correction in [e,n,u] = [{}, {}, {}]'.format(de, dn, du))
##    de, dn, du = compute_psd('ITRF2014-psd-gnss.dat', t=datetime.datetime.now(), station='PDEL')
##    print('PSD correction in [e,n,u] = [{}, {}, {}]'.format(de, dn, du))
