#ifndef __ITRF_CRD_TOOLS_HPP__
#define __ITRF_CRD_TOOLS_HPP__

#include <algorithm>
#include <cassert>
#include <fstream>
#include <functional>
#include <string>
#include <vector>
#ifdef DEBUG
#include <iostream>
#endif
#include "ggdatetime/datetime_write.hpp"
#include "ggdatetime/dtcalendar.hpp"

namespace ngpt {

/// A structure to hold a station and its coordinates.
/// It holds the station's name (4char-id + ' ' + domes number), plus the
/// 3 cartesian components (x, y, z). This is clearly a dummy class, to ease
/// the use of such simple collections.
struct sta_crd {
  double x, y, z;   ///< Coordinates in m, [x,y,z] components
  std::string site; ///< NAME+' '+DOMES = 9+1+4 chars

  /// Constructor.
  explicit sta_crd(const std::string &s, double xc, double yc,
                   double zc) noexcept
      : x{xc}, y{yc}, z{zc}, site{s} {};
};

namespace itrf_details {

/// A structure to hold a (full, two-line) SSC record for a station.
template <typename S> struct ssc_record {
  std::string site;       ///< NAME+' '+DOMES = 9+1+4 chars
  ngpt::datetime<S> from, ///< Validity interval, from ...
      to;                 ///< Vlidity interval, to ...
  double x, y, z,         ///< Coordinates in m
      vx, vy, vz,         ///< Velocities in m/y
      sx, sy, sz,         ///< Coordinate sigmas
      svx, svy, svz;      ///< Velocity sigmas
};

/// A simple class to hold PSD records.
template <typename S> struct psd_record {
  std::string site;          ///< NAME + ' ' + DOMES
  ngpt::datetime<S> teq;     ///< Time of earthquake
  int emdn, nmdn, umdn;      ///< Model numbers for e, n and u components
  double ea1, et1, ea2, et2; ///< a1,t1,a2,t2 parameters for the east component
  double na1, nt1, na2,
      nt2;                   ///< a1,t1,a2,t2 parameters for the north component
  double ua1, ut1, ua2, ut2; ///< a1,t1,a2,t2 parameters for the up component
};

/// Compare the first 4 chars of two strings. This function is mean to implement
/// (in the context it is used) a station 4-char-id comparisson.
inline int compare_sta_id(const std::string &str1, const std::string &str2) {
  return str1.compare(0, 4, str2, 0, 4);
}

/// Compare chars in the range [5, 9) of two strings. This function is mean to
/// implement (in the context it is used) a station DOMES number comparisson.
inline int compare_sta_domes(const std::string &str1, const std::string &str2) {
  return str1.compare(5, 9, str2, 5, 9);
}

/// Compute the post-seismic deformation/correction using parametric models:
///  - PWL (Piece-Wise Linear Function)
///  - Logarithmic Function
///  - Exponential Function
///  - Logarithmic + Exponential
///  - Two Exponential Functions
/// This is a translation of the fortran module parametric.f, found at
/// ftp://itrf.ign.fr/pub/itrf/itrf2014/parametric.f
///
/// @parameter[in] model The model to use to compute the PSD; valid are ints
///                      in range [0,4]
/// @parameter[in] dtq   Time difference (t - t_Earthquake) in decimal year. It
///                      It is advised to compute "dtq" by
///                      (MJD - MJD_Earthquake)/365.25 where MJD is the modified
///                      julian day.
/// @parameter[in] a1    Amplitude 1 of the parametric model, if modn = 1 or 2
///                      (or 3 or 4, if a2 & t2 are supplied) in mm
/// @parameter[in] a2    Amplitude 2 of the parametric model, if modn = 3 or 4
///                      in mm
/// @parameter[in] t1    Relaxation time 1, if modn = 1 or 2 (or 3 or 4, if a2
///                      & t2 are supplied) in decimal years
/// @parameter[in] t2    Relaxation time 2, if modn = 3 or 4 in decimal years
/// @return              The Post Seismic Deformation correction in mm.
///
double parametric(int model, double dtq = 0e0, double a1 = 0e0, double t1 = 0e0,
                  double a2 = 0e0, double t2 = 0e0) noexcept;

/// Read PSD model number and respective parameters off from a
/// ITRF2014-psd-*.dat file. This function will take (as input) a line from a
/// PSD .dat file, and resolve the fields recorded (hence we are really
/// interested at line columns >= 34). Depending on the model (number), the
/// corresponding number of coefficients will be collected (e.g. if the model
/// is 0, then no parameters will be read, but if the model is 1, 2 paramters
/// will be read).
///
/// @param[in]  line     The line (of a PSD .dat file) to be resolved.
/// @param[out] model_nr The model number collected from the line (int in
///                      range [0,4])
/// @param[out] a1       The a1 parameter of the model in mm; read if
///                      model_nr > 0.
/// @param[out] t1       The t1 parameter of the model in fractional years;
///                      read if model_nr > 0.
/// @param[out] a2       The a2 parameter of the model in mm; read if
///                      model_nr > 3.
/// @param[out] t2       The t2 parameter of the model in fractional years;
///                      read if model_nr > 3.
/// @return              The number of parameters read/collected. Anything
///                      other that an integer in the range [0, 4] denotes
///                      an error.
int read_psd_parameters(const std::string &line, int &model_nr, double &a1,
                        double &t1, double &a2, double &t2);

/// Read a station PSD record off from a PSD .dat file. This function will take
/// in a file stream (actually an open PSD .dat file) and try to read the PSD
/// record for a station, that is the record at which the file's get pointer is
/// set at. It will try to read 3 lines (one per component) and resolve all
/// respective fields.
/// PSD .dat files use a strict format, and they are available here:
/// http://itrf.ensg.ign.fr/ITRF_solutions/2014/ITRF2014_files.php
///
/// @param[in]  psd_stream The input file stream (an open PSD .dat file).
/// @param[out] rec        An instance of type psd_record<S>, where the resolved
///                        parameters/fields will be stored. Note that depending
///                        on the PSD model, only a subset of the instance will
///                        have valid values; e.g. if the model for the east
///                        component is 2 (i.e. rec.emdn = 2), then only the
///                        members rec.ea1 and rec.et1 will have valid values;
///                        the values of members rec.ea2 and rec.et2 will not
///                        be changed (kept as in input) and hence will have no
///                        meaning.
/// @return An integer value denoting the function status. A value of 0,
/// signifies
///         that the function did everything correctly and all three lines were
///         resolved successefuly. Else (i.e. in case of error), an integer
///         other thatn 0 is returned.
///
/// @warning The function will not check the validity of the resolved date (e.g.
///          will not check if year, day of year or seconds are valid numbers
///          and that the datetime instance formed is correct).
template <typename S>
int read_next_record_psd(std::ifstream &psd_stream, psd_record<S> &rec) {
  constexpr int max_chars{256};
  std::string line;
  line.reserve(max_chars);
  std::size_t idx;

  if (std::getline(psd_stream, line)) {
    rec.site = line.substr(1, 5);
    rec.site += line.substr(9, 9);
    int iyr, idoy;
    long isec;
    iyr = std::stoi(line.substr(19, 5), &idx);
    assert(idx == 2);
    iyr += (iyr > 70) ? (1900) : (2000);
    idoy = std::stoi(line.substr(22, 5), &idx);
    assert(idx == 3);
    isec = std::stoi(line.substr(26, 6), &idx);
    isec *= S::template sec_factor<long>(); // cast seconds to whatever S is
    ngpt::datetime<S> tmp{ngpt::year{iyr}, ngpt::day_of_year{idoy}, S{isec}};
    rec.teq = tmp;
    assert(line[32] == 'E');
    if (read_psd_parameters(line, rec.emdn, rec.ea1, rec.et1, rec.ea2,
                            rec.et2) < 0)
      return -1;
  } else {
    return 1;
  }

  if (std::getline(psd_stream, line)) {
    assert(line[32] == 'N');
    if (read_psd_parameters(line, rec.nmdn, rec.na1, rec.nt1, rec.na2,
                            rec.nt2) < 0)
      return -1;
  } else {
    return 1;
  }

  if (std::getline(psd_stream, line)) {
    assert(line[32] == 'U');
    if (read_psd_parameters(line, rec.umdn, rec.ua1, rec.ut1, rec.ua2,
                            rec.ut2) < 0)
      return -1;
  } else {
    return 1;
  }

  return 0;
}

/// Function to read the header off from a SSC-type file.
///
/// Given an SSC-type input file stream (ascii), this function will read the
/// first line and extract information (i.e. reference frame name and
/// reference epoch). Then it will read (and skip) the next 6 lines, so that
/// the input stream is at a position of reading record lines (i.e. next line
/// to be read from the stream will be the first record line).
/// The function will automatically go to the top of the file to read the
/// header.
///
/// @param[in] ssc_stream An SSC-type input file stream (aka std::ifstream)
/// @param[out] ref_frame The name of the reference frame as extracted from
///                       the file header (atually the first word of the first
///                       line).
/// @return               The reference epoch as float, i.e. as fractional year.
///                       or -1 to signal that somethin went wrong while
///                       resolving the fields of the first line.
/// @note  Always check that the returned year is greater than 0. If not, then
///        the header was not read properly.
///
float read_ssc_header(std::ifstream &ssc_stream, std::string &ref_frame);

/// Read a station record froma an SSC files (stream).
///
/// This function is used to read a (full, two-line) station record off from a
/// SSC file. It takes a great deal of prequations to check the validity of
/// the lines read. Two lines are sequentially read, and their fields are
/// resolved to fill in an ssc_record instance. The function will check for
/// the fields of validity interval (from, to dates); if they do not exist or
/// all their fields (i.e. year, day of month and seconds) are equal to 0, the
/// respective ssc_record entries are set to datetime<S>::min() and
/// datetime<S>::max(); otherwise they are extracted and filled in the
/// ssc_record. Note however that the validity of the dates is not explicitely
/// checked (i.e. if doy's are in the range 0-365/366). Note also that the
/// stream must be in a position that the next line to be read is the first of a
/// two-line ssc station record.
///
/// @param[in] ssc_stream The input SSC file stream; it should be in a position
///                       so that the next line to be read is the first of a
///                       two-line record.
/// @param[in] record     An instance of type ssc_record, where the resolved
///                       lines/fields are stored.
/// @return               On success, the function will return 0; else it will
///                       return some other int.
///
/// @note  Before calling this function (on an open SSC files), you should have
///        alredy called the read_ssc_header function (so that the stream is
///        set to the right first position).
template <typename S>
int read_next_record(std::ifstream &ssc_stream, ssc_record<S> &record) {
  constexpr int max_chars{256};
  std::string line;
  line.reserve(max_chars);

  // first line has domes, site_id, position info and validity interval
  if (std::getline(ssc_stream, line)) {
    record.site = line.substr(32, 5);  // 4-char id
    record.site += line.substr(0, 10); // domes
    std::size_t pos{36}, idx;
    record.x = std::stod(line.substr(pos, 20), &idx); // x
    pos += idx;
    record.y = std::stod(line.substr(pos, 20), &idx); // y
    pos += idx;
    record.z = std::stod(line.substr(pos, 20), &idx); // z
    pos += idx;
    record.sx = std::stod(line.substr(pos, 20), &idx); // sx
    pos += idx;
    record.sy = std::stod(line.substr(pos, 20), &idx); // sy
    pos += idx;
    record.sz = std::stod(line.substr(pos, 20), &idx); // sz
    pos += idx;
    record.from = ngpt::datetime<S>::min(); // preset from to min date
    record.to = ngpt::datetime<S>::max();   // preset to to max date
    if ((pos = line.find_first_of(':', pos)) != std::string::npos) {
      int iyr, idoy;
      long isec;
      // Resolve 'from' date string ...
      pos -= 2;
      iyr = std::stoi(line.substr(pos, 2), &idx);
      assert(idx == 2);
      pos += idx + 1;
      idoy = std::stoi(line.substr(pos, 3), &idx);
      assert(idx == 3);
      pos += idx + 1;
      isec = std::stol(line.substr(pos, 6), &idx); // seconds
      isec *= S::template sec_factor<long>(); // cast seconds to whatever S is
      if ((iyr + idoy + isec) != 0) {
        iyr += (iyr > 70) ? (1900) : (2000);
        ngpt::datetime<S> tmp{ngpt::year{iyr}, ngpt::day_of_year{idoy},
                              S{isec}};
        record.from = tmp;
      }
      // Resolve 'to' date sting ...
      pos += idx;
      pos = line.find_first_of(':', pos) - 2;
      iyr = std::stoi(line.substr(pos, 2), &idx);
      pos += idx + 1;
      idoy = std::stoi(line.substr(pos, 3), &idx);
      assert(idx == 3);
      pos += idx + 1;
      isec = std::stol(line.substr(pos));
      isec *= S::template sec_factor<long>();
      if ((iyr + idoy + isec) != 0) {
        iyr += (iyr > 70) ? (1900) : (2000);
        ngpt::datetime<S> tmp{ngpt::year{iyr}, ngpt::day_of_year{idoy},
                              S{isec}};
        record.to = tmp;
      }
    }
  } else {
    return 1;
  }

  // second line has velocity info
  if (std::getline(ssc_stream, line)) {
    assert(!line.compare(0, 9, record.site, 5, 9));
    std::size_t pos{36}, idx;
    record.vx = std::stod(line.substr(pos, 20), &idx); // vx
    pos += idx;
    record.vy = std::stod(line.substr(pos, 20), &idx); // vy
    pos += idx;
    record.vz = std::stod(line.substr(pos, 20), &idx); // vz
    pos += idx;
    record.svx = std::stod(line.substr(pos, 20), &idx); // svx
    pos += idx;
    record.svy = std::stod(line.substr(pos, 20), &idx); // svy
    pos += idx;
    record.svz = std::stod(line.substr(pos, 20), &idx); // svz
  } else {
    return 1;
  }

  return 0;
}

} // namespace itrf_details

template <typename S>
int ssc_extrapolate(std::ifstream &fin,
                    const std::vector<std::string> &stations,
                    const ngpt::datetime<S> &t, const ngpt::datetime<S> &t0,
                    std::vector<sta_crd> &results, bool use_domes = false) {
  std::function<int(const std::string &, const std::string &)> cmp =
      itrf_details::compare_sta_id;
  if (use_domes)
    cmp = itrf_details::compare_sta_domes;
  ngpt::datetime_interval<S> dt{ngpt::delta_date(t, t0)};
  double dyr = dt.as_mjd() / 365.25;

  results.clear();
  results.reserve(stations.size());

  itrf_details::ssc_record<S> record;
  std::vector<std::string> sta{stations};
  //  warning! if domes are provided, then the strings in the stations vector
  //  are
  //+ e.g. '97401M003', '92701M003', etc ..... Now, we need to transform these
  // to
  //+ ID+' '+DOMES, i.e. add 5* whitespaces at the begining
  if (use_domes)
    std::transform(sta.begin(), sta.end(), sta.begin(),
                   [](std::string &s) { return ("     " + s); });

  auto it = sta.begin();
  std::string site;
  while (!itrf_details::read_next_record<S>(fin, record) && sta.size()) {
    site = record.site;
    if ((it = std::find_if(sta.begin(), sta.end(), [=](const std::string &str) {
           return !cmp(site, str);
         })) != sta.end()) {
      if (t >= record.from && t < record.to) {
        sta.erase(it);
        auto x = record.x + (record.vx * dyr);
        auto y = record.y + (record.vy * dyr);
        auto z = record.z + (record.vz * dyr);
        results.emplace_back(site, x, y, z);
      }
    }
  }
  return results.size(); // number of stations actually found
}

template <typename S>
int compute_psd(const char *psd_file, const std::vector<std::string> &stations,
                const ngpt::datetime<S> &t, std::vector<sta_crd> &results,
                bool use_domes = false) {
  results.clear();
  results.reserve(stations.size());

  std::ifstream fin(psd_file);
  if (!fin.is_open())
    return -1;

  std::vector<std::string> sta(stations);
  //  warning! if domes are provided, then the strings in the stations vector
  //  are
  //+ e.g. '97401M003', '92701M003', etc ..... Now, we need to transform these
  // to
  //+ ID+' '+DOMES, i.e. add 5* whitespaces at the begining
  if (use_domes)
    std::transform(sta.begin(), sta.end(), sta.begin(),
                   [](std::string &s) { return ("     " + s); });

  std::function<int(const std::string &, const std::string &)> cmp =
      itrf_details::compare_sta_id;
  if (use_domes)
    cmp = itrf_details::compare_sta_domes;

  itrf_details::psd_record<S> rec;
  auto send = sta.end();
  auto rit = results.begin();
  std::string site;
  double dyr;
  while (!itrf_details::read_next_record_psd<S>(fin, rec)) {
    site = rec.site;
    // is the station read of interest?
    if (auto it =
            std::find_if(sta.begin(), sta.end(),
                         [=](const auto &str) { return !cmp(site, str); });
        it != send) {
      site = rec.site;
      // do we already have it in the results vector?
      if (rit = std::find_if(results.begin(), results.end(),
                             [&site = std::as_const(rec.site)](
                                 const sta_crd &a) { return site == a.site; });
          rit != results.end()) {
        ;
      } else {
        results.emplace_back(site, 0e0, 0e0, 0e0);
        rit = results.end() - 1;
      }
      // compute/append PSD
      if (t >= rec.teq) {
        ngpt::datetime_interval<S> dt(ngpt::delta_date(t, rec.teq));
        dyr = dt.as_mjd() / 365.25;
        rit->x += itrf_details::parametric(rec.emdn, dyr, rec.ea1, rec.et1,
                                           rec.ea2, rec.et2);
        rit->y += itrf_details::parametric(rec.nmdn, dyr, rec.na1, rec.nt1,
                                           rec.na2, rec.nt2);
        rit->z += itrf_details::parametric(rec.umdn, dyr, rec.ua1, rec.ut1,
                                           rec.ua2, rec.ut2);
      }
    }
  }
  return results.size(); // number of stations actually found
}

} // namespace ngpt

#endif
