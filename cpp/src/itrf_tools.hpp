#ifndef __ITRF_CRD_TOOLS_HPP__
#define __ITRF_CRD_TOOLS_HPP__

#include "itrf_details.hpp"
#include <vector>

namespace itrf {

using itrf_details::psd_model;
using itrf_details::StationId;

/// A structure to hold a station and its coordinates.
/// It holds the station's name (as StationId), plus the
/// 3 cartesian components (x, y, z). This is clearly a 
/// wrapper class to ease such simple collections.
struct sta_crd {
  StationId staid;
  double x{0e0}, y{0e0}, z{0e0}; ///< Coordinates in m, [x,y,z] components
  const char *name() const noexcept { return staid.name(); }
  const char *domes() const noexcept { return staid.domes(); }
  bool is_same(const sta_crd& other) const noexcept {
    return !std::strncmp(staid.mname, other.staid.mname, 4) && !std::strncmp(staid.mdomes, other.staid.mdomes, 9);
  }
}; // sta_crd


int ssc_extrapolate(std::ifstream &fin, const std::vector<StationId> &stations,
                    const ngpt::datetime<ngpt::seconds> &t,
                    const ngpt::datetime<ngpt::seconds> &t0,
                    std::vector<sta_crd> &results,
                    itrf_details::stationComparissonPolicy policy) noexcept;

int compute_psd(const char *psd_file, const std::vector<StationId> &stations,
                const ngpt::datetime<ngpt::seconds> &t,
                std::vector<sta_crd> &results,
                itrf_details::stationComparissonPolicy policy) noexcept;

int compute_psd_from_sinex(const char *sinex_file, const std::vector<StationId> &stations,
                const ngpt::datetime<ngpt::seconds> &t,
                std::vector<sta_crd> &results,
                itrf_details::stationComparissonPolicy policy) noexcept;

} // namespace itrf

#endif
