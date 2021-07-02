#include "itrf_details.hpp"
#include <cassert>
#include <cstdio>

/// Read a station PSD record off from a PSD .dat file.
int itrf::itrf_details::read_psd_parameters(const char *line, psd_model &model,
                                            double &a1, double &t1, double &a2,
                                            double &t2) noexcept {
  int model_nr = line[34] - '0';
  if (model_nr < 0 || model_nr > 4) {
    fprintf(stderr, "[ERROR] Invalid PSD model at line:\n%s", line);
    return -1;
  }

  model = int2model(model_nr);

  const char *start = line + 35;
  char *end;
  switch (model_nr) {
  case 0:
    return 0;
  case 1:
    // same as case 2
    [[fallthrough]];
  case 2:
    a1 = std::strtod(start, &end);
    if (start == end)
      return -1;
    start = end;
    a2 = std::strtod(start, &end);
    if (start == end)
      return -1;
    return 2;
  case 3:
    // same as case 4
    [[fallthrough]];
  case 4:
    a1 = std::strtod(start, &end);
    if (start == end)
      return -1;
    start = end;
    t1 = std::strtod(start, &end);
    if (start == end)
      return -1;
    start = end;
    a2 = std::strtod(start, &end);
    if (start == end)
      return -1;
    start = end;
    t2 = std::strtod(start, &end);
    if (start == end)
      return -1;
    return 4;
  }
  return -1;
}

int itrf::itrf_details::read_next_record_psd(std::ifstream &psd_stream,
                                             psd_record &rec) noexcept {
  constexpr int max_chars{256};
  char line[max_chars];
  char *start, *end;

  if (psd_stream.getline(line, max_chars)) {
    std::memcpy(rec.staid.mname, line + 1, 4);
    std::memcpy(rec.staid.mdomes, line + 9, 9);
    long idate[3];
    start = line + 19;
    for (int i = 0; i < 3; i++) {
      idate[i] = std::strtol(start, &end, 10);
      if (start == end)
        return 1;
      start = end + 1;
    }
    idate[0] += (idate[0] > 70) ? (1900) : (2000);
    ngpt::datetime<ngpt::seconds> tmp{ngpt::year(idate[0]),
                                      ngpt::day_of_year(idate[1]),
                                      ngpt::seconds(idate[2])};
    rec.teq = tmp;
    assert(line[32] == 'E');
    if (read_psd_parameters(line, rec.emdn, rec.ea1, rec.et1, rec.ea2,
                            rec.et2) < 0)
      return -1;
  } else {
    if (psd_stream.eof())
      return -1;
    fprintf(stderr, "[ERROR] Failed reading line from PSD file (#1).\n");
    return 1;
  }

  if (psd_stream.getline(line, max_chars)) {
    assert(line[32] == 'N');
    if (read_psd_parameters(line, rec.nmdn, rec.na1, rec.nt1, rec.na2,
                            rec.nt2) < 0)
      return -1;
  } else {
    return 1;
  }

  if (psd_stream.getline(line, max_chars)) {
    assert(line[32] == 'U');
    if (read_psd_parameters(line, rec.umdn, rec.ua1, rec.ut1, rec.ua2,
                            rec.ut2) < 0)
      return -1;
  } else {
    fprintf(stderr, "[ERROR] Failed reading line from PSD file (#2).\n");
    return 1;
  }

  return 0;
}
