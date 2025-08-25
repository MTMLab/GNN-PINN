#ifndef _NR_H_
#define _NR_H_
#include <fstream>
#include <complex>
#include "nrutil.h"
#include "nrtypes.h"

using namespace std;

namespace NR 
{

  void dist_ij(Mat_IO_DP &x, Mat_IO_DP &len, Mat_IO_DP &hlen, 
	       const int a, const int b, DP &r_ij, Vec_IO_DP &v_ij);

  void neighbor_list(Mat_IO_DP &x, Mat_IO_DP &len, Mat_IO_DP &hlen, Vec_IO_INT &elem, Mat_IO_INT &nn);

}
#endif /* _NR_H_ */







