#include <cmath>
#include "nr.h"
using namespace std;

DP NR::dist_si_sio2(Vec_IO_DP &len, Vec_IO_DP &hlen, Mat_IO_DP &x, 
		    const int a, const int b)
{
  static DP dx, dy, dz, dab;
  dx = fabs(x[a][0] - x[b][0]);
  if (dx > hlen[0])
    dx = len[0] - dx;
  dy = fabs(x[a][1] - x[b][1]);
  if (dy > hlen[1])
    dy = len[1] - dy;
  dz = fabs(x[a][2] - x[b][2]);
  if (dz > hlen[2])
    dz = len[2] - dz;
  dab = sqrt(dx*dx + dy*dy + dz*dz);
  return dab;
}
