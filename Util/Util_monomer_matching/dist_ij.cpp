#include <cmath>
#include "nr.h"
using namespace std;

extern int numC;
extern int numN;

void NR::dist_ij(Mat_IO_DP &x, Mat_IO_DP &len, Mat_IO_DP &hlen, 
		 const int a, const int b, DP &r_ij, Vec_IO_DP &v_ij)
{
  int i, j, ix, iy, iz, iter;
  DP d_ab;
  Vec_IO_INT d_cell(numC);
  Vec_IO_DP dx(numC), d_ij(numC), d_ij_hold(numC);

  for (j = 0; j < numC; j++) {
    d_ij[j] = x[b][j] - x[a][j];
  }
  
  d_ij_hold = d_ij;

  iter = 0;
  for (ix = -1; ix < 2; ix++) {
    d_cell[0] = ix;
    for (iy = -1; iy < 2; iy++) {
      d_cell[1] = iy;
      for (iz = -1; iz < 2; iz++) {
	d_cell[2] = iz;

	d_ij = d_ij_hold;
	for (i = 0; i < numC; i++) {
	  for (j = 0; j < numC; j++) {
	    d_ij[j] += d_cell[i]*len[i][j];
	  }
	}
	d_ab = sqrt(d_ij[0]*d_ij[0] + d_ij[1]*d_ij[1] + d_ij[2]*d_ij[2]);
	if (iter == 0) {
	  r_ij = d_ab;
	  v_ij = d_ij;
	}
	else {
	  if (d_ab < r_ij) {
	    r_ij = d_ab;
	    v_ij = d_ij;
	  }
	}
	iter++;
      }
    }
  }
}
