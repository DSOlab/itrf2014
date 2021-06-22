#include "itrf_tools.hpp"
#include <cassert>
#include <cstdio>
#include <fstream>

void next_nonwp(const char *start, char **end) noexcept {
  while (*start && *start != ' ')
    ++start;
  *end = (char *)start;
}

/// Function to read the header off from a SSC-type file.
double itrf::itrf_details::read_ssc_header(std::ifstream &ssc_stream,
                                           char *ref_frame) noexcept {
  constexpr std::size_t max_chars = 256;
  char line[max_chars];

  ssc_stream.seekg(0, std::ios::beg);
  if (!ssc_stream.getline(line, max_chars)) {
    fprintf(stderr, "[ERROR] Failed reading SSC file.\n");
    return 1;
  }

  // get the reference frame, which is the frst word in the line
  char *start = line, *end;
  while (*start && *start == ' ')
    ++start;
  if (!*start) {
    fprintf(stderr, "[ERROR] Failed resolving first line in SSC file.\n");
    return 1;
  }
  std::size_t rf_len = 0;
  end = start;
  while (*end && *end != ' ') {
    ++rf_len;
    ++end;
  }
  if (!*end) {
    fprintf(stderr, "[ERROR] Failed resolving first line in SSC file.\n");
    return 1;
  }
  std::memcpy(ref_frame, start, rf_len);
  ref_frame[rf_len] = '\0';

  // the header is actually pretty standard .... check the middle part
  const char *middle_part = "STATION POSITIONS AT EPOCH";
  std::size_t middle_part_len = std::strlen(middle_part);
  if (std::strncmp(middle_part, end + 1, middle_part_len)) {
    fprintf(stderr, "[ERROR] Expected string \"%s\" but did not find it!\n",
            middle_part);
    return 1;
  }
  start = end + middle_part_len + 1;

  // get the reference epoch
  double ref_epoch = std::strtod(start, &end);
  if (!ref_epoch || start == end) {
    fprintf(stderr, "[ERROR] Failed to resolve reference epoch.\n");
    return 1;
  }

  // check the last part
  const char *last_part = "AND VELOCITIES";
  std::size_t last_part_len = std::strlen(last_part);
  if (std::strncmp(last_part, end + 1, last_part_len)) {
    fprintf(stderr, "[ERROR] Expected string \"%s\" but did not find it!\n",
            last_part);
    return 1;
  }

  // read a bunch of no-info lines ....
  for (int i = 0; i < 6; i++)
    ssc_stream.getline(line, max_chars);

  // retun reference epoch as float
  return ref_epoch;
}

int itrf::itrf_details::read_next_record(
    std::ifstream &ssc_stream,
    itrf::itrf_details::ssc_record &record) noexcept {
  constexpr int max_chars{256};
  char line[max_chars];
  double dvec[3];

  // first line has domes, site_id, position info and validity interval
  if (ssc_stream.getline(line, max_chars)) {
    std::memcpy(record.staid.mname, line + 32, sizeof(char) * 4); // 4-char id
    std::memcpy(record.staid.mdomes, line, sizeof(char) * 9);     // domes

    char *end, *start;
    start = line + 36;

    // resolve coordinates
    for (int i = 0; i < 3; i++) {
      dvec[i] = std::strtod(start, &end);
      assert(dvec[i] && start != end);
      start = end;
    }
    record.x = dvec[0];
    record.y = dvec[1];
    record.z = dvec[2];

    // resolve std. deviations
    for (int i = 0; i < 3; i++) {
      dvec[i] = std::strtod(start, &end);
      assert(start != end);
      start = end;
    }
    record.sx = dvec[0];
    record.sy = dvec[1];
    record.sz = dvec[2];

    record.from =
        ngpt::datetime<ngpt::seconds>::min();         // preset from to min date
    record.to = ngpt::datetime<ngpt::seconds>::max(); // preset to to max date

    // find first occurance of ':' if any
    start = line;
    while (*start && *start != ':')
      ++start;

    if (*start) {
      long idate[3];
      assert(*start == ':');
      start -= 2;

      // Resolve 'from' date string ...
      for (int i = 0; i < 3; i++) {
        idate[i] = std::strtol(start, &end, 10);
        assert(start != end);
        start = end + 1;
      }
      if ((idate[0] + idate[1] + idate[2]) != 0) {
        idate[0] += (idate[0] > 70) ? (1900) : (2000);
        ngpt::datetime<ngpt::seconds> tmp{ngpt::year(idate[0]),
                                          ngpt::day_of_year(idate[1]),
                                          ngpt::seconds(idate[2])};
        record.from = tmp;
      }

      // Resolve 'to' date sting ...
      while (*start && *start != ':')
        ++start;
      if (!*start || *start != ':') {
        fprintf(stderr, "[ERROR] Failed to resolve endind date in line\n/%s",
                line);
        return 1;
      }
      start -= 2;
      for (int i = 0; i < 3; i++) {
        idate[i] = std::strtol(start, &end, 10);
        assert(start != end);
        start = end + 1;
      }
      if ((idate[0] + idate[1] + idate[2]) != 0) {
        idate[0] += (idate[0] > 70) ? (1900) : (2000);
        ngpt::datetime<ngpt::seconds> tmp{ngpt::year(idate[0]),
                                          ngpt::day_of_year(idate[1]),
                                          ngpt::seconds(idate[2])};
        record.to = tmp;
      }
    }
  } else {
    if (ssc_stream.eof())
      return -1;
    fprintf(stderr, "[ERROR] Failed to read line from SSC file.\n");
    return 1;
  }

  // second line has velocity info
  if (ssc_stream.getline(line, max_chars)) {
    if (std::strncmp(line, record.domes(), 9)) {
      fprintf(stderr, "[ERROR] Expected velocity info line, found:\n%s", line);
      return 2;
    }
    char *end, *start;
    start = line + 36;
    for (int i = 0; i < 3; i++) {
      dvec[i] = std::strtod(start, &end);
      assert(/*dvec[i] &&*/ start != end);
      start = end;
    }
    record.vx = dvec[0];
    record.vy = dvec[1];
    record.vz = dvec[2];
    for (int i = 0; i < 3; i++) {
      dvec[i] = std::strtod(start, &end);
      assert(dvec[i] && start != end);
      start = end;
    }
    record.svx = dvec[0];
    record.svy = dvec[1];
    record.svz = dvec[2];
  } else {
    fprintf(stderr, "[ERROR] Failed to read line from SSC file.\n");
    return 2;
  }

  return 0;
}