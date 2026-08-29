// Minimalkandidat 3: ein einziger Meniskus-Loesungsschritt.
#include <cstdio>
#include "es/meniscus.hpp"
using namespace es;
int main() {
  OpenCapillaryParams cp; cp.r_bore=1e-5; cp.r_outer=2e-5; cp.shank_length=4e-4; cp.z_rim=0;
  MeniscusParams mp; mp.r_contact=1e-5; mp.gamma=0.0452; mp.delta_p=0; mp.n_nodes=61;
  mp.max_outer=6; mp.tol=1e-3;
  MeniscusSolver s(make_capillary_open(cp), mp);
  MeniscusSolution m = s.solve_at_height(3e-6);
  std::printf("U = %.3f V, conv=%d\n", m.voltage, m.converged ? 1 : 0);
  return 0;
}
