# itrftools

itrftools is a collection of Python and C++ modules/programs to assist the 
processing and manipulation of ITRF14-related files; ease computation and 
extrapolation of station coordinates.

Specifically, the program can be used to perform:

* Coordinate extrapolation to a reference epoch using a SSC file, optionaly
      including a Post Seismic Deformation file (for ITRF2014)

* Post Seismic Deformation (PSD) computation (per component)

Stations to be considered can be identified either by their 4-char ID or via
their DOMES. See the script's help message for more info.

To make use of the package and/or script, you will need the relevant files
(e.g. PSD and SSC files). All of these are publicly available at the 
[ITRF website](http://itrf.ensg.ign.fr/) and are **NOT** included within the 
package.

The package contains two implementations, one in Python and one in C++, which 
can be installed individually; both give identical results and work exactly 
the same way. Users only need to install one of the two implementations.

## Python Version

The package requires [numpy](http://www.numpy.org/) and should work both for
Python 2.x and Python 3.x. No other requirement exists. The package is made up 
by a module (which needs to be installed) and a program under `/bin/itrftool`, 
which gets automatically installed (during the package installation process).

To install the package and the executable, run (in the folder `python/itrftools`)
`$> python setup.py install`
After that, you should have an executable named `itrftools` installed

## C++ Version

To install the C++ version you need to first install the 
[ggdatetime](https://github.com/xanthospap/ggdatetime.git) library. Once you 
have that, just configure and make, aka go into the `cpp` directory and type:
```
$> ./configure
$> make
```


# Usage

## Compute PSD Values
Compute PSD values for a given date and a list of stations; note that the stations 
can be specified either by name (aka their 4-char id) or by DOMES number (or both).
Here we compute the PSD values for stations with id's NRMD, COCO and TONG and stations with 
DOME's 97401M003, 50902M001 and 49971M001 for day of year 150 of year 2020 (that is
2020-05-29). Note that the DOMES 50902M001 and the id TONG correspond to the same station, 
hence only one record is written for it. Information on computed the PSD values are
extracted from the file `../data/ITRF2014-psd-gnss.dat`. Stations with PSD values of 0e0 
do not have corresponding records in the input PSD information file.
```
$> itrftool -s NRMD COCO TONG -m 97401M003 50902M001 49971M001 -p ../data/ITRF2014-psd-gnss.dat -y 2020 -d 150 --psd-only
NAME   DOMES   East(mm) North(mm) Up(mm)        EPOCH
---- --------- -------- -------- -------- ------------------
     97401M003     0.00     0.00     0.00 2020-05-29 00:00:00
COCO 50127M001    14.51    24.50     0.00 2020-05-29 00:00:00
NRMD               0.00     0.00     0.00 2020-05-29 00:00:00
TONG 50902M001    43.78   -13.65     0.00 2020-05-29 00:00:00
```
## Extrapolate ITRF2014 Coordinates
Extrapolate coordinates in ITRF2014 at epoch 2020-05-29 for the given stations (some 
specified by id others by domes). Note that the DOMES 50902M001 and the id TONG correspond to the same station,
hence only one record is written for it. Also, there is no record for a station with 
domes number 49971M001, hence no result is printed.
```
$> src/itrftool -s NRMD COCO TONG -m 97401M003 50902M001 49971M001  -y 2020 -d 150 -c ../data/ITRF2014_GNSS.SSC.txt -p ../data/ITRF2014-psd-gnss.dat

Reference Frame: ITRF2014, Reference Epoch: 2010-01-01 00:00:00
NAME   DOMES         X(m)           Y(m)            Z(m)        EPOCH
---- --------- --------------- --------------- --------------- ------------------
COCO 50127M001    -741951.09602   6190961.71574  -1337767.36193 2020-05-29 00:00:00
NRMD 92701M005   -5743538.11585   1380503.86427  -2397895.98837 2020-05-29 00:00:00
REUN 97401M003    3364098.92612   4907944.67286  -2293466.68314 2020-05-29 00:00:00
TONG 50902M001   -5930303.53647   -500148.80597  -2286366.30075 2020-05-29 00:00:00
```

> Minor format changes may be exhibeted between the C++ and the Python implementation; 
> e.g. Python results are not sorted (alphabeticaly)

# todo:

 - [ ] Compute standard deviation values for the extrapolated station coordinates.
 - [ ] install C++ binary (e.g. in /usr/local/bin)
 - [ ] configure Makefile.am for production
 - [ ] check clang++

# bug and comments

Please send any bugs, feedback, suggestions, comments, etc ..... to
xanthos@mail.ntua.gr or dganastasiou@gmail.com
