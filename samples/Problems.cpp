#include <fstream>
#include <iostream>
#include <samples/Problems.hpp>

#include <iomanip>
#include <sstream>

bool Problem::load(const std::string &file) {
  bool ok = true;
  std::ifstream f;
  f.open(file);

  if (ok && !f.good())
    ok = false;
  if (ok && !f.is_open())
    ok = false;

  if (ok) {
    f >> nCameras;
    f >> nPoints;
    int ne;
    f >> ne;
    edges.resize(ne);
    for (int i = 0; i < edges.size(); i++) {
      f >> edges[i].first;
      f >> edges[i].second;
    }
  }
  f.close();

  return ok;
}

std::string Problem::toString() const {
  std::ostringstream s;
  s << "nCameras " << nCameras << ", nPoints " << nPoints << ", nEdges " << edges.size() << std::endl;
  return s.str();
}

// TODO(sc) add the txt file as resource file and open it, here was too slow to compile
const Problem ProblemCollections::ladybug_49_7776 = Problem{49, 7776, {}};
