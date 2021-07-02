#include "itrf_tools.hpp"
#include "sinex.hpp"
#include <algorithm>
#include <stdexcept>

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


/// Find (and return a const_iterator to) a sinex::SolutionEstimate in the 
/// est_vec vector, to match the given ref_est.
auto find_matching_psd_parameter(const std::vector<sinex::SolutionEstimate>& est_vec,
  const sinex::SolutionEstimate& ref_est) noexcept {
    return std::find_if(est_vec.cbegin(), est_vec.cend(), [&](const sinex::SolutionEstimate& est) {
      return !std::strncmp(ref_est.m_param_type, est.m_param_type, 6)
        && !std::strncmp(ref_est.m_site_code, est.m_site_code, 4)
        && !std::strncmp(ref_est.m_point_code, est.m_point_code, 2)
        && (ref_est.m_epoch == est.m_epoch);}
        );
}

void filter_stations(const std::vector<itrf::StationId>& refsta, const std::vector<sinex::SiteId>& snxsta, 
  std::vector<sinex::SiteId>& ressta, itrf::itrf_details::stationComparissonPolicy policy) noexcept {
  
  ressta.reserve(refsta.size());

  switch (policy) {
    
    case (itrf::itrf_details::stationComparissonPolicy::use_name_id):
      for (const auto& ssta : snxsta) {
        auto it = std::find_if(refsta.cbegin(), refsta.cend(), [&](const itrf::StationId& ista){ return !std::strncmp(ista.mname, ssta.m_site_code, 4);});
        if (it != refsta.cend()) {
          ressta.emplace_back(ssta);

        }
      }
      return;
    
    case (itrf::itrf_details::stationComparissonPolicy::use_domes):
      for (const auto& ssta : snxsta) {
        auto it = std::find_if(refsta.cbegin(), refsta.cend(), [&](const itrf::StationId& ista){ return !std::strncmp(ista.mdomes, ssta.m_domes, 9);});
        if (it != refsta.cend()) {
          ressta.emplace_back(ssta);
        }
      }
      return;
    
    case (itrf::itrf_details::stationComparissonPolicy::use_full_name):
      for (const auto& ssta : snxsta) {
        auto it = std::find_if(refsta.cbegin(), refsta.cend(),[&](const itrf::StationId& ista){ 
          return !std::strncmp(ista.mname, ssta.m_site_code, 4) && !std::strncmp(ista.mdomes, ssta.m_domes, 9);});
        if (it != refsta.cend()) {
          ressta.emplace_back(ssta);
        }
      }
      return;
  }

  return;
}

auto find_corresponding_psd_param(const sinex::SolutionEstimate& psd, const std::vector<sinex::SolutionEstimate>& psd_vec) noexcept {
      char cor_type[7] = {'\0'};
 
      if (!std::strncmp(psd.m_param_type, "AEXP", 4)) {
        std::memcpy(cor_type, "TEXP", 4*sizeof(char));
        std::memcpy(cor_type+4, psd.m_param_type+4, 2);
      } else if (!std::strncmp(psd.m_param_type, "ALOG", 4)) {
        std::memcpy(cor_type, "TLOG", 4*sizeof(char));
        std::memcpy(cor_type+4, psd.m_param_type+4, 2);
      } else if (!std::strncmp(psd.m_param_type, "TEXP", 4) || !std::strncmp(psd.m_param_type, "TLOG", 4)) {
            std::memset(cor_type, '\0', sizeof(char)*7);
      } else {
        fprintf(stderr, "[ERROR] Unknow parameter type in PSD Sinex \"%s\"\n", psd.m_param_type);
        return psd_vec.cend();
      }

      if (*cor_type) {
        sinex::SolutionEstimate tpsd(psd);
        std::memcpy(tpsd.m_param_type, cor_type, 6);
        return find_matching_psd_parameter(psd_vec, tpsd);
      }

      return psd_vec.cend();
}


int itrf::compute_psd_from_sinex(
    const char *sinex_file, const std::vector<itrf::StationId> &stations,
    const ngpt::datetime<ngpt::seconds> &t, std::vector<itrf::sta_crd> &results,
    itrf::itrf_details::stationComparissonPolicy policy) noexcept {
  
  results.clear();
  results.reserve(stations.size());

  Sinex snx(sinex_file);

  std::vector<sinex::SiteId> sites_in_snx;
  sites_in_snx.reserve(200); /* approximately nr of sites in Sinex */
  if (snx.parse_block_site_id(sites_in_snx)) {
    fprintf(stderr, "[ERROR] Failed to collect sites from input SINEX file\n");
    return 1;
  }

  // filter collected station to only include the given ones
  std::vector<sinex::SiteId> sites_vec;
  filter_stations(stations, sites_in_snx, sites_vec, policy);
  /*printf("[DEBUG] These are the filtered solutions:\n");
  for (const auto& s : sites_vec) printf("-- \"%s %s\"\n", s.m_site_code, s.m_domes);*/

  // collect PSD parameters for each station in Sinex
  std::vector<sinex::SolutionEstimate> psd_vec;
  psd_vec.reserve(sites_vec.size()*8); // approximate hint for size */
  if (snx.parse_block_solution_estimate(psd_vec, sites_vec)) {
    fprintf(stderr, "[ERROR] Failed to collect solution estimates from SINEX\n");
    return 1;
  }

  // iterate through collected PSD's and compute total PSD per station. We are
  // interested in parameters of type "AEXP_?" and "ALOG_?", if and only if they
  // correspond to a station in sites_vec
  for (auto psd_it = psd_vec.cbegin(); psd_it != psd_vec.cend(); ++psd_it) {
    // printf("[DEBUG] iterating over PSD type \"%s\"\n", psd_it->m_param_type);
    if ( !std::strncmp(psd_it->m_param_type, "AEXP_", 5) || !std::strncmp(psd_it->m_param_type, "ALOG_", 5) ) {
      if (auto psd_site_it = std::find_if(
        sites_vec.cbegin(), sites_vec.cend(), [&](const sinex::SiteId &staid) {
            return !std::strncmp(psd_it->m_site_code, staid.m_site_code, 4) 
              && !std::strncmp(psd_it->m_point_code, staid.m_point_code, 2);
          }); psd_site_it != sites_vec.cend()) {
            // ok, so this psd referes to the site psd_site_it and we need to compute PSD
            // find the corresponding element in the resulting vector where we need to add
            // the PSD values
            auto res_it = std::find_if(results.begin(), results.end(), [&](const itrf::sta_crd& sta) {
              return !std::strcmp(sta.staid.name(), psd_site_it->m_site_code) 
                && !std::strcmp(sta.staid.domes(), psd_site_it->m_domes); 
              });
            if (res_it == results.end()) {
              results.emplace_back(itrf::sta_crd{});
              res_it = results.end() - 1;
              std::memcpy(res_it->staid.mname, psd_site_it->m_site_code, 4*sizeof(char));
              std::memcpy(res_it->staid.mdomes, psd_site_it->m_domes, 9*sizeof(char));
            }
            // good .... now we need the corresponding time parameter for the PSD we have
            if (t >= psd_it->m_epoch) {
              ngpt::datetime_interval<ngpt::seconds> dt(ngpt::delta_date(t, psd_it->m_epoch));
              double dyr = dt.as_mjd() / 365.25e0;
                auto t_param = find_corresponding_psd_param(*psd_it, psd_vec);
                if (t_param == psd_vec.cend()) {
                  fprintf(stderr, "[ERROR] Failed to find corresponding parameter for \"%s\"\n", psd_it->m_param_type);
                  return 1;
                }
                double* cmp = &(res_it->x);
                if (psd_it->m_param_type[5] == 'E') cmp = &(res_it->y);
                else if (psd_it->m_param_type[5] == 'U') cmp = &(res_it->z);
                psd_model model = psd_model::logarithmic;
                if (std::strncmp(psd_it->m_param_type+1, "EXP", 3))
                  model = psd_model::exponential;
                // printf("--> units \"%s\" and \"%s\"\n", psd_it->m_units, t_param->m_units);
                // printf("[DEBUG] adding psd to station %s name %s, with estimate a=%10.5f, t=%10.3f\n", res_it->name(), psd_it->m_param_type, psd_it->m_estimate, t_param->m_estimate);
                assert( !std::strcmp(psd_it->m_units, "m   ") && !std::strcmp(t_param->m_units, "y   ") );
                *cmp += itrf_details::parametric(model, dyr, psd_it->m_estimate * 1e3, t_param->m_estimate, 0e0, 0e0);
            }
          }
      }
    }

  return 0;
}

int itrf::compute_psd(
    const char *psd_file, const std::vector<itrf::StationId> &stations,
    const ngpt::datetime<ngpt::seconds> &t, std::vector<itrf::sta_crd> &results,
    itrf::itrf_details::stationComparissonPolicy policy) noexcept {

  results.clear();
  results.reserve(stations.size());

  std::ifstream fin(psd_file);
  if (!fin.is_open()) {
    fprintf(stderr, "[ERROR] Failed to open file \'%s\'\n", psd_file);
    return 1;
  }

  itrf_details::psd_record rec;
  double dyr;

  auto j = results.end();
  while (!itrf_details::read_next_record_psd(fin, rec)) {
    auto it = std::find_if(
        stations.cbegin(), stations.cend(), [&](const StationId &staid) {
          return !itrf_details::compare_stations(rec, staid, policy);
        });
    if (it != stations.cend()) {
      j = std::find_if(results.begin(), results.end(), [&](const sta_crd &crd) {
        return !itrf_details::compare_stations(rec, crd, policy);
      });
      if (j == results.end()) {
        results.emplace_back(sta_crd{});
        j = results.begin() + results.size() - 1;
        std::memcpy(j->staid.mname, rec.name(), 4);
        std::memcpy(j->staid.mdomes, rec.domes(), 9);
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
  return 0;
}

int itrf::ssc_extrapolate(
    std::ifstream &fin, const std::vector<itrf::StationId> &stations,
    const ngpt::datetime<ngpt::seconds> &t,
    const ngpt::datetime<ngpt::seconds> &t0,
    std::vector<itrf::sta_crd> &results,
    itrf::itrf_details::stationComparissonPolicy policy) noexcept {

  ngpt::datetime_interval<ngpt::seconds> dt{ngpt::delta_date(t, t0)};
  double dyr = dt.as_mjd() / 365.25e0;

  results.clear();
  results.reserve(stations.size());

  itrf_details::ssc_record record;
  auto it = stations.begin();
  int stations_found = 0;

  while (!itrf_details::read_next_record(fin, record) &&
         stations_found < (int)stations.size()) {
    if ((it = std::find_if(
             stations.begin(), stations.end(), [=](const StationId &str) {
               return !itrf_details::compare_stations(record, str, policy);
             })) != stations.end()) {
      if (t >= record.from && t < record.to) {
        results.emplace_back(sta_crd{});
        auto j = results.end() - 1;
        std::memcpy(j->staid.mname, record.name(), 4);
        std::memcpy(j->staid.mdomes, record.domes(), 9);
        j->x = record.x + (record.vx * dyr);
        j->y = record.y + (record.vy * dyr);
        j->z = record.z + (record.vz * dyr);
        ++stations_found;
        // printf("[DEBUG] Found crd for station %s %s in SSC: %15.5f %15.5f %15.5f\n", j->staid.mname, j->staid.mdomes, record.x, record.y, record.z);
        // printf("[DEBUG]                          Velocity : %15.5f %15.5f %15.5f delta-years=%20.15f\n", record.vx, record.vy, record.vz, dyr);
        // printf("[DEBUG] That is %20.9f + %20.6f * %20.15f = %20.5f\n", record.x, record.vx, dyr, j->x);
        printf("[DEBUG] Station %s %s, X=%20.5f\n", j->staid.mname, j->staid.mdomes, j->x);
      }
    }
  }
  return results.size(); // number of stations actually found
}
