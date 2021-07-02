#include "ggeodesy/car2ell.hpp"
#include "ggeodesy/geodesy.hpp"
#include "itrf_tools.hpp"
#include "ggdatetime/datetime_write.hpp"
#include <algorithm>
#include <array>
#include <cstring>
#include <map>

using itrf::StationId;

struct psd_delta {
  StationId staid;
  std::array<double, 6> dr = {0e0}; // de, dn, du, dx, dy, dz
  psd_delta(const StationId &str, const std::array<double, 6> &a) noexcept
      : staid(str), dr(a){};
};

void help_message() noexcept {
  std::cout << "usage: itrftool [-h] [-s [STATION_LIST [STATION_LIST ...]]] "
               "[-m [DOMES_LIST [DOMES_LIST ...]]] [-c SSC_FILE] [-p PSD_FILE] "
               "-y YEAR -d DOY [--psd-only]\n";
  std::cout
      << "\nExtrapolate coordinates from a SSC file, optionaly including a "
         "PSD file\n";
  std::cout << "\noptional arguments:";
  std::cout << "-h, --help            show this help message and exit\n";
  std::cout << "-s [STATION_LIST [STATION_LIST ...]], --stations "
               "[STATION_LIST [STATION_LIST ...]]\n";
  std::cout
      << "                      A whitespace seperated list of stations "
         "to compute coordinates for. The given names are checked against the "
         "\"4-char ID\" in the input files.\n";
  std::cout << "-m [DOMES_LIST [DOMES_LIST ...]], --domes [DOMES_LIST "
               "[DOMES_LIST ...]]\n";
  std::cout << "                      A whitespace seperated list of station "
               "domes to compute coordinates for. The given deomes are checked "
               "against "
               "the \"DOMES\" field in the input files.\n";
  std::cout << "-c SSC_FILE, --ssc SSC_FILE\n";
  std::cout << "                      A SSC ascci file to extract coordinates "
               "and velocities from. These files are normaly accessible at: "
               "http://itrf.ign.fr/ITRF_solutions/\n";
  std::cout << "-x SINEX_FILE, --sinex SINEX_FILE\n";
  std::cout << "                      A SINEX psd file to extract psd parameters "
               "from. These files are normaly accessible at: "
               "http://itrf.ign.fr/ITRF_solutions/\n";
  std::cout << "-p PSD_FILE, --psd PSD_FILE\n";
  std::cout
      << "                      A PSD ascci file to extract "
         "Post-Seismic-Deformation models andparameters from. These files are "
         "normaly accessible at: "
         "ftp://itrf.ign.fr/pub/itrf/itrf2014/ITRF2014-psd-gnss.dat\n";
  std::cout
      << "-y YEAR, --year YEAR  The year to extrapolate coordinates at.\n";
  std::cout << "-d DOY, --doy DOY     The day of year to extrapolate "
               "coordinates at.\n";
  std::cout << "--psd-only            If this switch is on, then only the PSD "
               "corrections are computed (per component); no extrapolation of "
               "coordinates "
               "is performed\n";
  std::cout << "\nNational Technical University of Athens";
  std::cout << "\nDionysos Satellite Observatory - 2020\n";
  return;
}

int parse_cmd(int argc, char *argv[],
              std::map<char, std::vector<std::string>> &cmd_map) {
  int dummy = 0;
  // setup default positional arguments
  cmd_map['n'] = std::vector<std::string>{"0"}; //< psd only

  for (int i = 1; i < argc;) {
    if (++dummy > 1000)
      return 50; /* just to be safe ... */
    if (!std::strcmp(argv[i], "-s") || !std::strcmp(argv[i], "--stations")) {
      if (argc <= i + 1)
        return 1;
      ++i;
      cmd_map['s'] = std::vector<std::string>();
      while (i < argc) {
        if (argv[i][0] == '-')
          break;
        cmd_map['s'].emplace_back(std::string(argv[i]));
        ++i;
      }
    } else if (!std::strcmp(argv[i], "-m") ||
               !std::strcmp(argv[i], "--domes")) {
      if (argc <= i + 1)
        return 1;
      ++i;
      cmd_map['m'] = std::vector<std::string>();
      while (i < argc) {
        if (argv[i][0] == '-')
          break;
        cmd_map['m'].emplace_back(std::string(argv[i]));
        ++i;
      }
    } else if (!std::strcmp(argv[i], "-c") || !std::strcmp(argv[i], "--ssc")) {
      if (argc <= i + 1)
        return 1;
      cmd_map['c'] = std::vector<std::string>{std::string(argv[i + 1])};
      i += 2;
    } else if (!std::strcmp(argv[i], "-p") || !std::strcmp(argv[i], "--psd")) {
      if (argc <= i + 1)
        return 1;
      cmd_map['p'] = std::vector<std::string>{std::string(argv[i + 1])};
      i += 2;
    } else if (!std::strcmp(argv[i], "-x") || !std::strcmp(argv[i], "--sinex")) {
      if (argc <= i + 1)
        return 1;
      cmd_map['x'] = std::vector<std::string>{std::string(argv[i + 1])};
      i += 2;
    } else if (!std::strcmp(argv[i], "-y") || !std::strcmp(argv[i], "--year")) {
      if (argc <= i + 1)
        return 1;
      cmd_map['y'] = std::vector<std::string>{std::string(argv[i + 1])};
      i += 2;
    } else if (!std::strcmp(argv[i], "-d") || !std::strcmp(argv[i], "--doy")) {
      if (argc <= i + 1)
        return 1;
      cmd_map['d'] = std::vector<std::string>{std::string(argv[i + 1])};
      i += 2;
    } else if (!std::strcmp(argv[i], "--psd-only")) {
      cmd_map['n'][0] = "1";
      i += 1;
    } else if (!std::strcmp(argv[i], "-h") || !std::strcmp(argv[i], "--help")) {
      cmd_map['h'] = std::vector<std::string>{"0"};
      i += 1;
    } else {
      std::cerr << "\n[WARNING] Invalid command line argument \"" << argv[i]
                << "\". Skipping";
      ++i;
    }
  }

  // make sure the user didn't mess up the args ...
  auto mend = cmd_map.end();
  auto i1 = cmd_map.find('n');
  if (i1->second[0] == "1" && ( cmd_map.find('p') == mend && cmd_map.find('x')==mend)) {
    std::cerr << "[ERROR] If you need the PSD values, you need to supply a "
                 "PSD or SINEX file!\n";
    return 1;
  }
  if ((i1->second[0] == "0" && cmd_map.find('c') == mend) && cmd_map.find('h') == mend) {
    std::cerr << "[ERROR] You need to supply an SSC file for coordinate "
                 "extrapolation\n";
    return 1;
  }
  if (cmd_map.find('p') != mend && cmd_map.find('x')!=mend ) {
    std::cerr << "[ERROR] Only provide one of PSD/SINEX files\n";
    return 1;
  }
  if ((cmd_map.find('y') == mend || cmd_map.find('d') == mend) && cmd_map.find('h')==mend) {
    std::cerr << "[ERROR] Need to provide a year and a day_of_year\n";
    return 1;
  }
  return 0;
}

std::vector<StationId>
str2stationId(const std::vector<std::string> &vstr,
              itrf::itrf_details::stationComparissonPolicy policy) noexcept {
  std::vector<StationId> res;
  if (policy == itrf::itrf_details::stationComparissonPolicy::use_name_id) {
    for (const auto &str : vstr) {
      StationId n;
      std::strncpy(n.mname, str.c_str(), 4);
      res.push_back(n);
    }
  } else if (policy ==
             itrf::itrf_details::stationComparissonPolicy::use_domes) {
    for (const auto &str : vstr) {
      StationId n;
      std::strncpy(n.mdomes, str.c_str(), 9);
      res.push_back(n);
    }
  } else {
    for (const auto &str : vstr) {
      StationId n;
      std::strncpy(n.mname, str.c_str(), 4);
      const char *c = str.c_str();
      std::strncpy(n.mdomes, c + 5, 9);
      res.push_back(n);
    }
  }
  return res;
}

std::vector<itrf::sta_crd> merge_sort_unique(std::vector<itrf::sta_crd> &v1,
                                             std::vector<itrf::sta_crd> &v2) {
  using itrf::sta_crd;
  // concatenate vectors to v1
  v1.insert(v1.end(), std::make_move_iterator(v2.begin()),
            std::make_move_iterator(v2.end()));
  // sort based on station name
  std::sort(v1.begin(), v1.end(), [](const sta_crd &a, const sta_crd &b) {
    return std::strncmp(a.staid.name(), b.staid.name(), 4);
  });
  // move duplicates to end .....
  auto last =
      std::unique(v1.begin(), v1.end(), [](const sta_crd &a, const sta_crd &b) {
        return !(std::strncmp(a.staid.domes(), b.staid.domes(), 9));
      });
  // and delete them
  v1.erase(last, v1.end());
  return v1;
}
  
int compute_psd__(const char* psd_fn, const char* snx_fn, const std::vector<std::string>& sites, ngpt::datetime<ngpt::seconds>& t, std::vector<itrf::sta_crd>& results,
  itrf::itrf_details::stationComparissonPolicy policy) {

    if (psd_fn)
    
      return itrf::compute_psd(
          psd_fn,
          str2stationId(
              sites,
              policy),
          t, results, policy);
    
    return itrf::compute_psd_from_sinex(
          snx_fn,
          str2stationId(
              sites,
              policy),
          t, results, policy);
  }


int main(int argc, char *argv[]) {
  using sec = ngpt::seconds;

  std::map<char, std::vector<std::string>> cmd_map;
  if (parse_cmd(argc, argv, cmd_map))
    return 10;

  // if help message requested, print help and exit
  if (cmd_map.find('h') != cmd_map.end()) {
    help_message();
    return 0;
  }

  // epoch t to extrapolate at
  auto mend = cmd_map.end();
  auto it = cmd_map.find('y');
  int year = std::stoi(it->second[0]);
  it = cmd_map.find('d');
  int doy = std::stoi(it->second[0]);
  ngpt::datetime<sec> t(ngpt::year{year}, ngpt::day_of_year{doy}, sec{0});

  // vectors to store results
  std::vector<itrf::sta_crd> res1, res2, psd_results;

  // Compute PSD if a sinex or PSD file is given
  if (cmd_map.find('p') != mend) {
    std::string psd_file = cmd_map['p'][0];
    if (auto its = cmd_map.find('s'); its != mend) {
      if (compute_psd__(psd_file.c_str(), nullptr, its->second, t, res1, itrf::itrf_details::stationComparissonPolicy::use_name_id)) return 11;
    }
    if (auto its = cmd_map.find('m'); its != mend) {
      if (compute_psd__(psd_file.c_str(), nullptr, its->second, t, res2, itrf::itrf_details::stationComparissonPolicy::use_domes)) return 11;
    }
    psd_results = merge_sort_unique(res1, res2);
    printf("[DEBUG] PSD file found; computed deformation for %3d sites\n", (int)psd_results.size());
  } else if (cmd_map.find('x') != mend) {
    std::string snx_file = cmd_map['x'][0];
    if (auto its = cmd_map.find('s'); its != mend) {
      if (compute_psd__(nullptr, snx_file.c_str(), its->second, t, res1, itrf::itrf_details::stationComparissonPolicy::use_name_id)) return 12;
    }
    if (auto its = cmd_map.find('m'); its != mend) {
      if (compute_psd__(nullptr, snx_file.c_str(), its->second, t, res2, itrf::itrf_details::stationComparissonPolicy::use_domes)) return 12;
    }
    psd_results = merge_sort_unique(res1, res2);
    printf("[DEBUG] SINEX file found; computed deformation for %3d sites\n", (int)psd_results.size());
  }

  // easy case: We have a PSD file but no SSC; Only compute PSD in [e,n,u]
  it = cmd_map.find('n');
  if (auto tit = cmd_map.find('c'); it->second[0] == "1" && tit == mend) {
    printf("\nNAME   DOMES   East(mm) North(mm) Up(mm)        EPOCH");
    printf("\n---- --------- -------- -------- -------- ------------------");
    for (const auto &i : psd_results) {
      printf("\n%s %s %+8.2f %+8.2f %+8.2f %s", i.staid.name(), i.staid.domes(),
             i.x, i.y, i.z, ngpt::strftime_ymd_hms(t).c_str());
    }
    std::cout << "\n";
    return 0;
  }

  // open the SSC file and extract reference epoch and frame
  char ref_frame[25];
  double refyear;
  std::ifstream fin(cmd_map['c'][0]);
  if ((refyear = itrf::itrf_details::read_ssc_header(fin, ref_frame)) < 0e0) {
    std::cerr << "\n[ERROR] Failed reading SSC header for \"" << cmd_map['c'][0]
              << "\"";
    return -1;
  }

  // reference epoch (from float) as datetime instance
  int t0_yr = (int)refyear;
  assert((double)t0_yr - refyear == 0e0);
  ngpt::datetime<sec> t0(ngpt::year(t0_yr), ngpt::day_of_year(1), sec(0));

  // extrapolate coordinates to epoch t using the ssc file
  std::vector<itrf::sta_crd> ssc_results;
  if (!res1.empty()) res1.clear();
  if (!res2.empty()) res2.clear();
  if (auto its = cmd_map.find('s'); its != mend) {
    itrf::ssc_extrapolate(
        fin,
        str2stationId(
            its->second,
            itrf::itrf_details::stationComparissonPolicy::use_name_id),
        t, t0, res1, itrf::itrf_details::stationComparissonPolicy::use_name_id);
  }
  if (auto its = cmd_map.find('m'); its != mend) {
    itrf::ssc_extrapolate(
        fin,
        str2stationId(its->second,
                      itrf::itrf_details::stationComparissonPolicy::use_domes),
        t, t0, res2, itrf::itrf_details::stationComparissonPolicy::use_domes);
  }
  ssc_results = merge_sort_unique(res1, res2);
  // printf("[DEBUG] Coordinates extrapolated for %2d sites\n", (int)results.size());

  // add PSD (after convertion to cartesian) to the xtrapolated coordinates
  std::vector<psd_delta> psd_info; // for reporting later on ...
  if (psd_results.size()) {
    for (auto &ssc : ssc_results) {
      auto psd = std::find_if(psd_results.cbegin(), psd_results.cend(), [&](const itrf::sta_crd& psd_neu){return psd_neu.is_same(ssc);});
      if (psd != psd_results.cend()) {
          double lat, lon, hgt, dx, dy, dz;
          ngpt::car2ell<ngpt::ellipsoid::grs80>(ssc.x, ssc.y, ssc.z, lat, lon, hgt);
          ngpt::top2car(psd->y * 1e-3, psd->x * 1e-3, psd->z * 1e-3, lat, lon,
                        dx, dy, dz);
          ssc.x += dx;
          ssc.y += dy;
          ssc.z += dz;
          if (cmd_map['n'][0] == "1") {
            psd_info.emplace_back(
                ssc.staid, std::array<double, 6>{psd->x, psd->y, psd->z,
                                                 dx * 1e3, dy * 1e3, dz * 1e3});
          }
      }
    }   
  }

  std::cout << "\nReference Frame: " << ref_frame
            << ", Reference Epoch: " << ngpt::strftime_ymd_hms(t0);
  if (cmd_map['n'][0] != "1") {
    printf("\nNAME   DOMES         X(m)           Y(m)            Z(m)        "
           "EPOCH");
    printf("\n---- --------- --------------- --------------- --------------- "
           "------------------");
    for (const auto &i : ssc_results) {
      printf("\n%s %s %15.5f %15.5f %15.5f %s", i.staid.name(), i.staid.domes(),
             i.x, i.y, i.z, ngpt::strftime_ymd_hms(t).c_str());
    }
  } else {
    printf(
        "\nNAME   DOMES   East(mm) North(mm) Up(mm)   X(mm)    Y(mm)     Z(mm) "
        "     EPOCH");
    printf("\n---- --------- -------- -------- -------- -------- -------- "
           "-------- ------------------");
    for (const auto &i : psd_info) {
      printf("\n%s %s", i.staid.name(), i.staid.domes());
      for (int k = 0; k < 6; ++k)
        printf("%+8.2f ", i.dr[k]);
      ngpt::strftime_ymd_hms(t).c_str();
    }
  }

  std::cout << "\n";
  return 0;
}