/************************************************************
 * November 29, 2005
 * Name: Sillinger Weber for Alloy, Neighbor Part  
 * File name: swNeighbor.cpp
 *
 * Developers: Sangheon Lee and Gyeong S. Hwang
 * 
 *
 *
 *************************************************************/

#include <cmath>
#include "nr.h"
using namespace std;

extern int nemax, namax;
extern int na, nnmax;
extern int npmax;
extern DP cutoff1;
extern int nsamplemax;

void NR::swNeighbor_mode4(Vec_IO_DP &len, Vec_IO_DP &hlen, Mat_IO_DP &x, 
			  Mat_IO_INT &np, Vec_IO_INT &cntnp)
{
  int i, j, k, l, ii, jj, kk, ll;  
  //  int nnmax1;
  DP dx, dy, dz, dij;
  //  Mat_IO_DP r(na, nnmax);

  np = -1;

  //r = 0.0;

  for (i = 0; i < namax; i++) 
    cntnp[i] = 0;

  for (i = 0; i < idsize2; i++) {
    ii = id[i];
    for (j = 0; j < na; j++) {
      jj = j;
      if (jj != ii) {
	dx = x[ii][0] - x[jj][0];
	if (dx >= -hlen[0] && dx < hlen[0])
	  ;
	else if (dx < -hlen[0]) {
	  dx += len[0];
	}
	else {
	  dx -= len[0];  
	}
	dy = x[ii][1] - x[jj][1];
	if (dy >= -hlen[1] && dy < hlen[1])
	  ;
	else if (dy < -hlen[1]) {
	  dy += len[1];
	}
	else {
	  dy -= len[1];
	}
	dz = x[ii][2] - x[jj][2];
	if (dz >= -hlen[2] && dz < hlen[2])
	  ;
	else if (dz < -hlen[2]) {
	  dz += len[2];
	}
	else {
	  dz -= len[2];
	}
      
	dij = sqrt(dx*dx + dy*dy + dz*dz);

	if (dij < cutoff1) {
	  np[ii][cntnp[ii]] = jj;
	  //	  np[jj][cntnp[jj]] = ii;
	  cntnp[ii]++;
	  //	  cntnp[jj]++;
	}
      }
    }
  }
}

