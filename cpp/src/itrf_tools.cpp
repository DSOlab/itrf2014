#include "itrf_tools.hpp"
#include <stdexcept>
#include <algorithm>

/// Compute the post-seismic deformation/correction using a parametric model.
double itrf::itrf_details::parametric(itrf::psd_model model, double dtq,
                                      double a1, double t1, double a2,
                                      double t2) noexcept {
  double d{0e0}, te1, te2;
  switch (model) {
  case psd_model::pwl: // PWL (Piece-Wise Linear Function)
    d = 0e0;
    break;
  case psd_model::logarithmic: // Logarithmic Function
    d = a1 * std::log(1e0 + dtq / t1);
    break;
  case psd_model::exponential: // Exponential Function
    te1 = dtq / t1;
    d = a1 * (1e0 - std::exp(-te1));
    break;
  case psd_model::logexp: // Logarithmic + Exponential
    te2 = dtq / t2;
    d = a1 * std::log(1e0 + dtq / t1) + a2 * (1e0 - std::exp(-te2));
    break;
  case psd_model::twoexp: // Two Exponential Functions
    te1 = dtq / t1;
    te2 = dtq / t2;
    d = a1 * (1e0 - std::exp(-te1)) + a2 * (1e0 - std::exp(-te2));
  }

  return d;
}

itrf::psd_model itrf::itrf_details::int2model(int mnr) {
  switch (mnr) {
  case 0:
    return psd_model::pwl;
  case 1:
    return psd_model::logarithmic;
  case 2:
    return psd_model::exponential;
  case 3:
    return psd_model::logexp;
  case 4:
    return psd_model::twoexp;
  default:
    throw std::runtime_error("[ERROR] Failed to transform int to psd model!");
  }
}

int itrf::compute_psd(
    const char *psd_file, const std::vector<itrf::StationId> &stations,
    const ngpt::datetime<ngpt::seconds> &t, std::vector<itrf::sta_crd> &results,
    itrf::itrf_details::stationComparissonPolicy policy) noexcept {
  results.clear();
  results.reserve(stations.size());

  std::ifstream fin(psd_file);
  if (!fin.is_open())
    return -1;

  itrf_details::psd_record rec;
  double dyr;

  auto j = results.end();
  while (!itrf_details::read_next_record_psd(fin, rec)) {
    if (auto it = std::find_if(stations.cbegin(), stations.cend(),
                               [&](const StationId& staid) {
                                 return !itrf_details::compare_stations(
                                     rec, staid, policy);
                               });
        it != stations.cend()) {
      // do we arleady have the station?
      if (j = std::find_if(results.begin(), results.end(),
                           [&](const sta_crd &crd) {
                             return !itrf_details::compare_stations(rec, crd,
                                                                    policy);
                           });
          j == results.end()) {
        results.emplace_back(sta_crd{});
        j = results.begin() + results.size() - 1;
        std::strncmp(j->staid.name, rec.name, 4);
        std::strncmp(j->staid.domes, rec.domes, 9);
      }
      // compute/append PSD
      if (t >= rec.teq) {
        ngpt::datetime_interval<ngpt::seconds> dt(ngpt::delta_date(t, rec.teq));
        dyr = dt.as_mjd() / 365.25e0;
        j->x += itrf_details::parametric(rec.emdn, dyr, rec.ea1, rec.et1,
                                         rec.ea2, rec.et2);
        j->y += itrf_details::parametric(rec.nmdn, dyr, rec.na1, rec.nt1,
                                         rec.na2, rec.nt2);
        j->z += itrf_details::parametric(rec.umdn, dyr, rec.ua1, rec.ut1,
                                         rec.ua2, rec.ut2);
      }
    }
  }
  return results.size(); // number of stations actually found
}

int itrf::ssc_extrapolate(std::ifstream &fin,
                    const std::vector<itrf::StationId> &stations,
                    const ngpt::datetime<ngpt::seconds> &t, const ngpt::datetime<ngpt::seconds> &t0,
                    std::vector<itrf::sta_crd> &results, itrf::itrf_details::stationComparissonPolicy policy) noexcept {
  
  ngpt::datetime_interval<ngpt::seconds> dt{ngpt::delta_date(t, t0)};
  double dyr = dt.as_mjd() / 365.25e0;

  results.clear();
  results.reserve(stations.size());

  itrf_details::ssc_record record;
  auto it = stations.begin();
  int stations_found =0;

  while (!itrf_details::read_next_record(fin, record) && stations_found<stations.size()) {
    if ((it = std::find_if(stations.begin(), stations.end(), [=](const StationId& &str) {
           return !itrf_details::compare_stations<itrf_details::ssc_record, itrf::StationId>(record, str, policy);
         })) != stations.end()) {
      if (t >= record.from && t < record.to) {
        results.emplace_back(sta_crd{});
        auto j = results.end() - 1;
        j->x = record.x + (record.vx * dyr);
        j->y = record.y + (record.vy * dyr);
        j->z = record.z + (record.vz * dyr);
        ++stations_found;
      }
    }
  }
  return results.size(); // number of stations actually found
}