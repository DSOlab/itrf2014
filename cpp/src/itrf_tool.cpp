#include <array>
#include <cstring>
#include <map>

#include "ggeodesy/car2ell.hpp"
#include "ggeodesy/geodesy.hpp"
#include "itrf_tools.hpp"

struct psd_delta {
  std::string site;
  std::array<double, 6> dr = {0e0}; // de, dn, du, dx, dy, dz
  psd_delta(std::string str, std::array<double, 6> &&a) noexcept
      : site(std::move(str)), dr(std::move(a)){};
};

template<int N>
struct charArray { 
  char ar[N] = {'\0'}; 
  explicit charArray(const char* str) noexcept {
    std::memcpy(ar, str, N-1);
  }
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
      cmd_map['s'] = std::vector<charArray<5>>();
      while (i < argc) {
        if (argv[i][0] == '-')
          break;
        cmd_map['s'].emplace_back(argv[i]);
        ++i;
      }
    } else if (!std::strcmp(argv[i], "-m") ||
               !std::strcmp(argv[i], "--domes")) {
      if (argc <= i + 1)
        return 1;
      ++i;
      cmd_map['m'] = std::vector<charArray<10>>();
      while (i < argc) {
        if (argv[i][0] == '-')
          break;
        cmd_map['m'].emplace_back(argv[i]);
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
  if (i1->second[0] == "1" && cmd_map.find('p') == mend) {
    std::cerr << "[ERROR] If you need the PSD values, you need to supply a "
                 "PSD file!\n";
    return 1;
  }
  if (i1->second[0] == "0" && cmd_map.find('c') == mend) {
    std::cerr << "[ERROR] You need to supply an SSC file for coordinate "
                 "extrapolation\n";
    return 1;
  }
  if (cmd_map.find('y') == mend || cmd_map.find('d') == mend) {
    std::cerr << "[ERROR] Need to provide a year and a day_of_year\n";
    return 1;
  }
  return 0;
}

/*
std::vector<itrf::sta_crd> merge_sort_unique(std::vector<itrf::sta_crd> &v1,
                                             std::vector<itrf::sta_crd> &v2) {
  using itrf::sta_crd;
  // concatenate vectors to v1
  v1.insert(v1.end(), std::make_move_iterator(v2.begin()),
            std::make_move_iterator(v2.end()));
  // sort based on station name
  std::sort(v1.begin(), v1.end(),
            [](const sta_crd &a, const sta_crd &b) { return a.site < b.site; });
  // move duplicates to end .....
  auto last =
      std::unique(v1.begin(), v1.end(), [](const sta_crd &a, const sta_crd &b) {
        return a.site == b.site;
      });
  // and delete them
  v1.erase(last, v1.end());
  return std::move(v1);
}
*/

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
  std::vector<itrf::sta_crd> res1, res2;

  // easy case: We have a PSD file but no SSC; Only compute PSD in [e,n,u]
  it = cmd_map.find('n');
  if (auto tit = cmd_map.find('c'); it->second[0] == "1" && tit == mend) {
    std::string psd_file = cmd_map['p'][0];
    if (auto its = cmd_map.find('s'); its != mend) {
      ngpt::compute_psd(psd_file.c_str(), its->second, t, res1, false);
    }
    if (auto its = cmd_map.find('m'); its != mend) {
      ngpt::compute_psd(psd_file.c_str(), its->second, t, res2, true);
    }
    auto results = merge_sort_unique(res1, res2);
    printf("\nNAME   DOMES   East(mm) North(mm) Up(mm)        EPOCH");
    printf("\n---- --------- -------- -------- -------- ------------------");
    for (const auto &i : results) {
      printf("\n%s %+8.2f %+8.2f %+8.2f %s", i.site.c_str(), i.x, i.y, i.z,
             ngpt::strftime_ymd_hms(t).c_str());
    }
    std::cout << "\n";
    return 0;
  }

  // open the SSC file and extract reference epoch and frame
  std::string ref_frame;
  float refyear;
  std::ifstream fin(cmd_map['c'][0]);
  if ((refyear = ngpt::itrf_details::read_ssc_header(fin, ref_frame)) < 0e0) {
    std::cerr << "\n[ERROR] Failed reading SSC header for \"" << cmd_map['c'][0]
              << "\"";
    return -1;
  }

  // reference epoch (from float) as datetime instance
  int t0_yr = (int)refyear;
  assert((float)t0_yr - refyear == 0e0);
  ngpt::datetime<mlsec> t0(ngpt::year(t0_yr), ngpt::day_of_year(1), mlsec(0));

  // extrapolate coordinates to epoch t using the ssc file
  if (auto its = cmd_map.find('s'); its != mend) {
    ngpt::ssc_extrapolate(fin, its->second, t, t0, res1, false);
  }
  if (auto its = cmd_map.find('m'); its != mend) {
    ngpt::ssc_extrapolate(fin, its->second, t, t0, res2, true);
  }
  auto results = merge_sort_unique(res1, res2);

  // do we need to add the PSD's ?
  std::vector<psd_delta> psd_info;
  if (it = cmd_map.find('p'); it != mend) {
    std::vector<ngpt::sta_crd> pres;
    std::string psd_file = it->second[0];
    // compute the psd values for the station-id vector  ....
    if (auto its = cmd_map.find('s'); its != mend) {
      ngpt::compute_psd(psd_file.c_str(), its->second, t, pres, false);
      // add results to the original results vector (psd are in [n,e,u] and mm)
      for (auto &rec : results) {
        if (auto rit =
                std::find_if(pres.begin(), pres.end(),
                             [&ref = std::as_const(rec)](const auto &a) {
                               return !(a.site.compare(0, 4, ref.site, 0, 4));
                             });
            rit != pres.end()) {
          double lat, lon, hgt, dx, dy, dz;
          ngpt::car2ell<ngpt::ellipsoid::grs80>(rec.x, rec.y, rec.z, lat, lon,
                                                hgt);
          ngpt::top2car(rit->y * 1e-3, rit->x * 1e-3, rit->z * 1e-3, lat, lon,
                        dx, dy, dz);
          rec.x += dx;
          rec.y += dy;
          rec.z += dz;
          if (cmd_map['n'][0] == "1") {
            psd_info.emplace_back(
                rec.site, std::array<double, 6>{rit->x, rit->y, rit->z,
                                                dx * 1e3, dy * 1e3, dz * 1e3});
          }
        } /* else {
           std::cout<<"\n[DEBUG] No PSD record for station \""<<rec.site<<"\"";
         }*/
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
    for (const auto &i : results) {
      printf("\n%s %15.5f %15.5f %15.5f %s", i.site.c_str(), i.x, i.y, i.z,
             ngpt::strftime_ymd_hms(t).c_str());
    }
  } else {
    printf(
        "\nNAME   DOMES   East(mm) North(mm) Up(mm)   X(mm)    Y(mm)     Z(mm) "
        "     EPOCH");
    printf("\n---- --------- -------- -------- -------- -------- -------- "
           "-------- ------------------");
    for (const auto &i : psd_info) {
      printf("\n%s", i.site.c_str());
      for (int k = 0; k < 6; ++k)
        printf("%+8.2f ", i.dr[k]);
      ngpt::strftime_ymd_hms(t).c_str();
    }
  }

  std::cout << "\n";
  return 0;

  /*
  using mlsec = ngpt::milliseconds;

  const char* ssc_file = "ITRF2008_GNSS.SSC.txt";
  std::string reff;
  float       reft;

  std::ifstream fin (ssc_file);
  reft = ngpt::itrf_details::read_ssc_header(fin, reff);
  std::cout<<"\nFrame is \""<<reff<<"\", time is "<<reft<<"\n";

  int t0_yr = (int)reft;
  assert((float)t0_yr-reft == 0e0);
  ngpt::datetime<mlsec> t0 {ngpt::year{t0_yr}, ngpt::day_of_year{1}, mlsec{0}};

  std::vector<std::string> stations;
  stations.emplace_back("NRMD 92701M005");
  stations.emplace_back("COCO");
  stations.emplace_back("REUN 97401M003");
  stations.emplace_back("AZRY 49971M001");
  stations.emplace_back("ANKR");

  // ngpt::datetime<mlsec> t {ngpt::year{2017}, ngpt::day_of_year{143},
  ngpt::hours{12},
  // ngpt::minutes{40}, mlsec{365554}};
  ngpt::datetime<mlsec> t{ngpt::year{2017}, ngpt::day_of_year{143}, mlsec{0}};

  std::vector<ngpt::sta_crd> sta_crd;

  ngpt::ssc_extrapolate(fin, stations, t, t0, sta_crd);
  printf("\nNAME   DOMES         X(m)           Y(m)            Z(m) EPOCH");
  printf("\n---- --------- --------------- --------------- ---------------
  ------------------"); for (const auto& i : sta_crd) { printf("\n%s %15.5f
  %15.5f %15.5f %s", i.site.c_str(), i.x, i.y, i.z,
  ngpt::strftime_ymd_hms(t).c_str());
  }

  ngpt::compute_psd("ITRF2014-psd-gnss.dat", stations, t, sta_crd);
  printf("\nNAME   DOMES         X(mm)          Y(mm)           Z(mm) EPOCH");
  printf("\n---- --------- --------------- --------------- ---------------
  ------------------"); for (const auto& i : sta_crd) { printf("\n%s %15.5f
  %15.5f %15.5f %s", i.site.c_str(), i.x, i.y, i.z,
  ngpt::strftime_ymd_hms(t).c_str());
  }

  std::cout<<"\n";
  return 0;
  */
}
