#ifndef __ITRF_CRD_TOOLS_HPP__
#define __ITRF_CRD_TOOLS_HPP__

#include "ggdatetime/datetime_write.hpp"
#include "ggdatetime/dtcalendar.hpp"
#include <cstdint>
#include <fstream>
#include <vector>

namespace itrf {

/// @class StationId A struct to hold a statio id.
/// A station id has two (string) fields: a 4-character name and
/// a DOMES number (9 chars).
/// Members are c-strings; at construction last character of the strings is
/// '\0'. Do not change that.
struct StationId {
  char mname[5] = {'\0'};
  char mdomes[10] = {'\0'};
  const char *name() const noexcept { return mname; }
  const char *domes() const noexcept {
    return mdomes;
  }
}; //StationId

/// A structure to hold a station and its coordinates.
/// It holds the station's name (as StationId), plus the
/// 3 cartesian components (x, y, z). This is clearly a 
/// wrapper class to ease such simple collections.
struct sta_crd {
  StationId staid;
  double x, y, z; ///< Coordinates in m, [x,y,z] components
  const char *name() const noexcept { return staid.name(); }
  const char *domes() const noexcept { return staid.domes(); }
}; // sta_crd

/// @enum psd_model
/// Possible ITRF-defined Post Seismic Deformation Models
enum class psd_model : uint_fast8_t {
  pwl = 0, // piece-wise linear
  logarithmic,
  exponential,
  logexp, // logarithmic plus exponential
  twoexp  // two exponential functions
}; // psd_model

namespace itrf_details {

/// @brief Cast an integer to a psd_model enum
psd_model int2model(int model_nr);

/// @enum stationComparissonPolicy
/// An enum class to enable different comparisson options when 
/// comparing stationId's.
enum class stationComparissonPolicy : uint_fast8_t {
  use_name_id = 0, // only name (aka 4char id)
  use_domes,       // only domes number
  use_full_name    // name plus domes
}; // stationComparissonPolicy

/// A structure to hold a (full, two-line) SSC record for a station.
struct ssc_record {
  StationId staid;
  ngpt::datetime<ngpt::seconds> from, ///< Validity interval, from ...
      to;                             ///< Vlidity interval, to ...
  double x, y, z,                     ///< Coordinates in m
      vx, vy, vz,                     ///< Velocities in m/y
      sx, sy, sz,                     ///< Coordinate sigmas
      svx, svy, svz;                  ///< Velocity sigmas
  const char *name() const noexcept { return staid.name(); }
  const char *domes() const noexcept { return staid.domes(); }
}; // ssc_record

/// A simple class to hold PSD records.
struct psd_record {
  StationId staid;
  psd_model emdn, nmdn, umdn; ///< Model numbers for e, n and u components
  ngpt::datetime<ngpt::seconds> teq; ///< Time of earthquake
  double ea1, et1, ea2,
      et2; ///< a1,t1,a2,t2 parameters for the east component
  double na1, nt1, na2,
      nt2;                   ///< a1,t1,a2,t2 parameters for the north component
  double ua1, ut1, ua2, ut2; ///< a1,t1,a2,t2 parameters for the up component
  const char *name() const noexcept { return staid.name(); }
  const char *domes() const noexcept { return staid.domes(); }
}; // psd_record

template <typename T1, typename T2>
inline int compare_stations(const T1 &s1, const T2 &s2,
                            stationComparissonPolicy policy) noexcept {
  switch (policy) {
  case stationComparissonPolicy::use_name_id:
    return std::strncmp(s1.name(), s2.name(), 4);
  case stationComparissonPolicy::use_domes:
    return std::strncmp(s1.domes(), s2.domes(), 9);
  case stationComparissonPolicy::use_full_name:
    return std::strncmp(s1.name(), s2.name(), 4) &&
           std::strncmp(s1.domes(), s2.domes(), 9);
  }
  return -1;
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
/// @parameter[in] dtq   Time difference (t - t_Earthquake) in decimal year.
/// It
///                      It is advised to compute "dtq" by
///                      (MJD - MJD_Earthquake)/365.25 where MJD is the
///                      modified julian day.
/// @parameter[in] a1    Amplitude 1 of the parametric model, if modn = 1 or 2
///                      (or 3 or 4, if a2 & t2 are supplied) in mm
/// @parameter[in] a2    Amplitude 2 of the parametric model, if modn = 3 or 4
///                      in mm
/// @parameter[in] t1    Relaxation time 1, if modn = 1 or 2 (or 3 or 4, if a2
///                      & t2 are supplied) in decimal years
/// @parameter[in] t2    Relaxation time 2, if modn = 3 or 4 in decimal years
/// @return              The Post Seismic Deformation correction in mm.
///
double parametric(psd_model model, double dtq = 0e0, double a1 = 0e0,
                  double t1 = 0e0, double a2 = 0e0, double t2 = 0e0) noexcept;

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
int read_psd_parameters(const char *line, psd_model &model_nr, double &a1,
                        double &t1, double &a2, double &t2) noexcept;

/// Read a station PSD record off from a PSD .dat file. This function will
/// take in a file stream (actually an open PSD .dat file) and try to read the
/// PSD record for a station, that is the record at which the file's get
/// pointer is set at. It will try to read 3 lines (one per component) and
/// resolve all respective fields. PSD .dat files use a strict format, and
/// they are available here:
/// http://itrf.ensg.ign.fr/ITRF_solutions/2014/ITRF2014_files.php
///
/// @param[in]  psd_stream The input file stream (an open PSD .dat file).
/// @param[out] rec        An instance of type psd_record<S>, where the
/// resolved
///                        parameters/fields will be stored. Note that
///                        depending on the PSD model, only a subset of the
///                        instance will have valid values; e.g. if the model
///                        for the east component is 2 (i.e. rec.emdn = 2),
///                        then only the members rec.ea1 and rec.et1 will have
///                        valid values; the values of members rec.ea2 and
///                        rec.et2 will not be changed (kept as in input) and
///                        hence will have no meaning.
/// @return An integer value denoting the function status. A value of 0,
/// signifies
///         that the function did everything correctly and all three lines
///         were resolved successefuly. Else (i.e. in case of error), an
///         integer other thatn 0 is returned.
///
/// @warning The function will not check the validity of the resolved date
/// (e.g.
///          will not check if year, day of year or seconds are valid numbers
///          and that the datetime instance formed is correct).
int read_next_record_psd(std::ifstream &psd_stream, psd_record &rec) noexcept;

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
/// @return               The reference epoch as float, i.e. as fractional
/// year.
///                       or -1 to signal that somethin went wrong while
///                       resolving the fields of the first line.
/// @note  Always check that the returned year is greater than 0. If not, then
///        the header was not read properly.
///
double read_ssc_header(std::ifstream &ssc_stream, char *ref_frame) noexcept;

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
/// stream must be in a position that the next line to be read is the first of
/// a two-line ssc station record.
///
/// @param[in] ssc_stream The input SSC file stream; it should be in a
/// position
///                       so that the next line to be read is the first of a
///                       two-line record.
/// @param[in] record     An instance of type ssc_record, where the resolved
///                       lines/fields are stored.
/// @return               On success, the function will return 0; else it will
///                       return some other int.
///
/// @note  Before calling this function (on an open SSC files), you should
/// have
///        alredy called the read_ssc_header function (so that the stream is
///        set to the right first position).
int read_next_record(std::ifstream &ssc_stream, ssc_record &record) noexcept;

} // namespace itrf_details

int ssc_extrapolate(std::ifstream &fin, const std::vector<StationId> &stations,
                    const ngpt::datetime<ngpt::seconds> &t,
                    const ngpt::datetime<ngpt::seconds> &t0,
                    std::vector<sta_crd> &results,
                    itrf_details::stationComparissonPolicy policy) noexcept;

int compute_psd(const char *psd_file, const std::vector<StationId> &stations,
                const ngpt::datetime<ngpt::seconds> &t,
                std::vector<sta_crd> &results,
                itrf_details::stationComparissonPolicy policy) noexcept;

} // namespace itrf

#endif
