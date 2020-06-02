#include "itrf_tools.hpp"

/// Compute the post-seismic deformation/correction using a parametric model.
double
ngpt::itrf_details::parametric(int model, double dtq, double a1,
    double t1, double a2, double t2)
noexcept
{
  double d {0e0}, te1, te2;
  switch (model) {
    case 0: // PWL (Piece-Wise Linear Function)
      d = 0e0;
      break;
    case 1: // Logarithmic Function
      d = a1*std::log(1e0+dtq/t1);
      break;
    case 2: // Exponential Function
      te1 = dtq/t1;
      d   =  a1*(1e0-std::exp(-te1));
      break;
    case 3: // Logarithmic + Exponential
      te2 = dtq/t2;
      d = a1*std::log(1e0+dtq/t1) + a2*(1e0-std::exp(-te2));
      break;
    case 4: // Two Exponential Functions
      te1 = dtq/t1;
      te2 = dtq/t2;
      d   = a1*(1e0-std::exp(-te1)) + a2*(1e0-std::exp(-te2));
  }

  return d;
}

/// Read a station PSD record off from a PSD .dat file.
int
ngpt::itrf_details::read_psd_parameters(const std::string& line, int& model_nr,
    double& a1, double& t1, double& a2, double& t2)
{
  model_nr = line[34] - '0';
  if (model_nr<0 || model_nr>4) return -1;

  std::size_t pos {35}, idx;
  switch (model_nr) {
    case 0:
      return 0;
    case 1:
      // same as case 2
    case 2:
      a1 = std::stod(line.substr(pos, 10), &idx);   // a1
      pos += idx;
      t1 = std::stod(line.substr(pos, 10), &idx);   // t1
      return 2;
    case 3:
      // same as case 4
    case 4:
      a1 = std::stod(line.substr(pos, 10), &idx);   // a1
      pos += idx;
      t1 = std::stod(line.substr(pos, 10), &idx);   // t1
      pos += idx;
      a2 = std::stod(line.substr(pos, 10), &idx);   // a1
      pos += idx;
      t2 = std::stod(line.substr(pos, 10), &idx);   // t1
      return 4;
  }
  return -1;
}

/// Function to read the header off from a SSC-type file.
float
ngpt::itrf_details::read_ssc_header(std::ifstream& ssc_stream,
    std::string& ref_frame)
{
  using pos_t =  std::string::size_type;

  constexpr pos_t max_chars {256};
  const std::string middle_part {"STATION POSITIONS AT EPOCH"},
                      last_part {"AND VELOCITIES"};
  const pos_t mdp_sz {middle_part.size()},
              ltp_sz {last_part.size()};
  const char whitesp = ' ';
  auto npos = std::string::npos;
  std::string line;
  line.reserve(max_chars);

  ssc_stream.seekg(0, std::ios::beg);
  std::getline(ssc_stream, line);

  // get the reference frame, which is the frst word in the line
  pos_t length = line.size();
  pos_t pos1   = line.find_first_not_of(whitesp);
  pos_t pos2   = line.find_first_of(whitesp, pos1);
  if (!((pos1!=npos && pos2!=npos) && (pos2>pos1))) return -1e0;
  ref_frame    = line.substr(pos1, pos2-pos1);

  // the header is actually pretty standard .... check the middle part
  if (!(length>(pos2+mdp_sz)) ||
      line.compare(pos2+1, mdp_sz, middle_part)) return -1e0;

  // get the reference epoch
  pos1 = pos2 + mdp_sz + 1;
  pos2 = line.find_first_of(whitesp, pos1+1);
  std::string ref_epoch {line.substr(pos1, pos2-pos1)};

  // check the last part
  if (!(length>=(pos2 + ltp_sz)) ||
     line.compare(pos2+1, ltp_sz, last_part)) return -1e0;

  // read a bunch of no-info lines ....
  for (int i=0; i<6; i++) std::getline(ssc_stream, line);

  // retun reference epoch as float
  return std::stof(ref_epoch);
}
