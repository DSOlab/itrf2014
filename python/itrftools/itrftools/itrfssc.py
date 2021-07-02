#! /usr/bin/python
#-*- coding: utf-8 -*-

from __future__ import print_function
import sys, datetime
from itrftools.compute_psd import time_str2dt


def read_header(fstream):
    """ Read the header of a .SSC file (i.e. a reference frame info file).
        This function will read the first line and resolve the reference frame
        name and reference epoch. Next, it will read and skip the next 6 lines.
        An example .SSC file, can be found at 
        http://itrf.ign.fr/ITRF_solutions/2008/doc/ITRF2008_GNSS.SSC.txt.

        Parameters:
        -----------
        fstream: An input file stream
                 The (already opened) .SSC file.

        Returns:
        ---------
        tuple
            A two-element tuple, where the first element is the reference frame
            name and the second is the reference epoch, as a datetime.datetime
            instance.
    """
    line = fstream.readline()
    l = line.split()
    assert (len(l) == 8)
    assert (l[1] == 'STATION' and l[2] == 'POSITIONS' and l[3] == 'AT' and
            l[4] == 'EPOCH' and l[6] == 'AND' and l[7] == 'VELOCITIES')
    refframe = l[0]
    refepoch = int(float(l[5]))
    assert (refepoch == float(l[5]))
    for i in range(6):
        fstream.readline()
    return refframe, datetime.datetime(year=refepoch,
                                       month=1,
                                       day=1,
                                       hour=0,
                                       minute=0,
                                       second=0)


def read_next_record(line, fstream):
    """ Read a record line (for one station) off from a SSC file. Actually, the
        first line of the record is passed in as the 'line' parameter, while the
        second is read (from this function) using the input file stream fstream.
        The function will read and resolve all fields, and return them as a
        dictionary. If the record has a validity interval (i.e. the fields
        'DATA_START' and 'DATA_END' are not empty), they will be used; otherwise
        they are set to datetime.min and datetime.max respectively.

        Parameters:
        -----------
        line: string
            The first line of a .SSC record for a station (e.g.
        (ITRF2008)    10002M006 Grasse (OCA)     GNSS GRAS  4581690.900   556114.837  4389360.793 0.001 0.001 0.001  2 03:113:00000 04:295:43200
        (ITRF2014)    10002M006 Grasse (OCA)     GNSS GRAS  4581690.8267   556114.9242  4389360.8453 0.0006 0.0006 0.0006  1 00:000:00000 96:277:00000
            see also http://itrf.ign.fr/ITRF_solutions/2008/doc/ITRF2008_GNSS.SSC.txt)
        fstream: Input file stream
            An input SSC file stream, where the read position is exactly after
            the line passed in the function (i.e. next line to be read, is the
            second line of the station record).

        Returns:
        --------
        dictionary
            A dictionary where the keys are the fields read off from the records
    """
    domes = line[0:9]
    name = line[10:27]
    tqn = line[27:31]
    cid = line[32:36]
    tstart = datetime.datetime.min
    tstop = datetime.datetime.max
    try:  # start and stop times included
        stop_col = line[36:].index(':') - 2 + 36
        dstrs = line[stop_col:].split()
        assert (len(dstrs) == 2)
        if dstrs[0].strip() != '00:000:00000':
            tstart = time_str2dt(dstrs[0])
        if dstrs[1].strip() != '00:000:00000':
            tstop = time_str2dt(dstrs[1])
        stop_col -= 4
    except:
        stop_col = len(line) - 1
    x, y, z, sx, sy, sz = map(float, line[36:stop_col].split())
    line = fstream.readline()
    l = line.split()
    assert (l[0] == domes and len(l) == 7)
    vx, vy, vz, svx, svy, svz = map(float, l[1:])
    return {
        'domes': domes,
        'name': name,
        'tqn': tqn,
        'id': cid,
        'start': tstart,
        'stop': tstop,
        'x': x,
        'y': y,
        'z': z,
        'vx': vx,
        'vy': vy,
        'vz': vz,
        'sx': sx,
        'sy': sy,
        'sz': sz,
        'svx': svx,
        'svy': svy,
        'svz': svz
    }


def extrapolate(dtq, x0, vx):
    """ Given a linear model described by x0 and Vx (constant coef. and velocity)
        at a reference epoch t0, compute the y value at time ti. The parameter
        dtq is the datetime difference in (fractional) years, from ti to t0.
        Parameters:
        -----------
        dtq: float
             Delta years from ti to t0 (i.e. time difference from reference epoch
             to the time we want the computation at).
        x0: float
            constant term of the linear model
        vx: float
            Velocity of the linear mode (must be annual)
        Returns:
        ---------
        float
            Value of the model at delta time dtq (i.e. x0+vx(ti-t0) = x0+vx*dtq).
    """
    return x0 + vx * dtq


def itrf_extrapolate(ssc_file,
                     t0,
                     t=datetime.datetime.now(),
                     station=[],
                     domes=[]):
    """ Given an (ITRF) .SSC file, compute the coordinates of a station list
        at the given epoch (t).The stations can be described by either passing
        in their 4-char id, or their DOMES numbers.

        Parameters:
        -----------
        ssc_file: string
            The .SSC (.txt) file to extract the coordinates from. (see e.g.
            http://itrf.ign.fr/ITRF_solutions/2014/doc/ITRF2014_GNSS.SSC.txt)
        t: datetime.datetime
            The time we want the coordinates at.
        t0: datetime.datetime
            The reference epoch
        station: list of strings
            The names of the stations (4-char id)
        domes: list of strings
            The domes numbers.

        Returns:
        --------
        a list of dictionaries.
            Each entry of the returned list, is a dictionary of type:
            {'station':, 'domes': , 'x':, 'y':, 'z':}, where x, y and z are 
            the extrapolated coordinates in m.
            Stations which were not matched are not included in the list.

    """
    results = []
    if station:
        station = [x.upper() for x in station]
    with open(ssc_file) as fin:
        frame, refepoch = read_header(fin)
        line = fin.readline()
        while line:
            dic = read_next_record(line, fin)
            if dic['domes'] in domes or dic['id'] in station:
                if t >= dic['start'] and t < dic['stop']:
                    dt = t - t0
                    dyr = (dt.days + dt.seconds / 86400e0) / 365.25e0
                    print('[DEBUG] Found crd for station {:} in SSC: {:15.5f} {:15.5f} {:15.5f}'.format(dic['id']+ ' ' +dic['domes'], dic['x'], dic['y'], dic['z']))
                    print('                               Velocity : {:15.5f} {:15.5f} {:15.5f} delta_year={:20.15f}'.format(dic['vx'], dic['vy'], dic['vz'], dyr))
                    x = extrapolate(dyr, dic['x'], dic['vx'])
                    y = extrapolate(dyr, dic['y'], dic['vy'])
                    z = extrapolate(dyr, dic['z'], dic['vz'])
                    results.append({
                        'station': dic['id'],
                        'domes': dic['domes'],
                        'x': x,
                        'y': y,
                        'z': z
                    })
            line = fin.readline()
        return results


# example usage
#if __name__ == "__main__":
#    x,y,z = itrf_extrapolate('ITRF2008_GNSS.SSC.txt', t=datetime.datetime.now(), station='ANKR')
#    print('Coordinates X={}, Y={}, Z={}'.format(x,y,z))
#    x,y,z = itrf_extrapolate('ITRF2014_GNSS.SSC.txt', t=datetime.datetime.now(), station='ANKR')
#    print('Coordinates X={}, Y={}, Z={}'.format(x,y,z))
