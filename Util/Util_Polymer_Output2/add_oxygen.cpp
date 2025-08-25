#include "nr.h"
#include <cmath>
#include <iostream>
using namespace std;
using std::cout;
using std::cin;
using std::endl;

#include <iomanip>
#include <cstdlib>

void NR::add_oxygen(Vec_IO_DP &len, Vec_IO_DP &hlen, 
		    Mat_IO_DP &x, Mat_IO_DP &x1,  Mat_IO_INT &nn)
{ 
  int i, j, k, l;
  int count;
  DP b;
  const int n = x.nrows();
  const int n1 = nn.nrows();
  const int m1 = nn.ncols();
  const int n2 = x1.nrows(); // Number of Si-Si bonds or number of Oxygen atoms
  Mat_IO_INT bond(n, n);
  count = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < 4; j++) {
      if (bond[i][nn[i][j]] == 0) {
	for (k = 0; k < 3; k++) {
	  b = fabs(x[nn[i][j]][k]-x[i][k]);
	  if (b > hlen[k]) 
	    x1[count][k] = (x[i][k] + x[nn[i][j]][k] + len[k]) / 2.0;
	  else 
	    x1[count][k] = (x[i][k] + x[nn[i][j]][k]) / 2.0;
	}
	count++;
	bond[i][nn[i][j]] = bond[nn[i][j]][i] = 1;
      }
    }
  }  
} 



