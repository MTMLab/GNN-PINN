#include "nr.h"
#include <cmath>
#include <iostream>
using namespace std;
using std::cout;
using std::cin;
using std::endl;

#include <iomanip>
#include <cstdlib>

void NR::neighbor_sinw(Vec_IO_DP &len, Vec_IO_DP &hlen,
 		       const int nTotal, const int nSiTotal, 
		       const int nOTotal, Mat_IO_DP &x, 
		       Mat_IO_INT &nn)
{ 
  
  int i, j, k, k1;
  const int m = 3;
  const int num1 = 4;
  const int num2 = 2;
  Vec_IO_INT count_si(nSiTotal), count_o(nOTotal);
  Vec_IO_INT count_si_hold(nSiTotal), count_o_hold(nOTotal);
  Mat_IO_INT nn_hold(nTotal, num1);
  for (i = 0; i < nSiTotal; i++) {
    count_si[i] = 0;
    count_si_hold[i] = 0;
  }
  for (i = 0; i < nOTotal; i++) {
     count_o[i] = 0;
     count_o_hold[i] = 0;
  }
  for (i = 0; i < nTotal; i++) {
    for (j = 0; j < num1; j++) {
      nn[i][j] = -1;      
    }
  }

  DP d_ss;
  DP dij;
  int si1, si2, o1, o2;
  int flag;

  for (i = nSiTotal; i < nTotal; i++) {
    cout << "i = " << i << endl;
    nn_hold = nn;
    count_si_hold = count_si;
    count_o_hold = count_o;
    d_ss = 0.0;
    while (count_o[i-nSiTotal] < num2) {
      count_si = count_si_hold;
      count_o = count_o_hold;
      nn = nn_hold;
      for (j = 0; j < nSiTotal; j++) {
	flag = 0;
	if (count_si[j] < num1 && count_o[i-nSiTotal] < num2) {
	  dij = dist_si_sio2(len, hlen, x, i, j);
	  if (dij < d_ss) {
	    for (k = 0; k < count_si[j]; k++) {
	      if (nn[j][k] == i) 
		flag = 1;
	      /*
	      else {
		if (nn[j][k] < nSiTotal) {
		  for (k1 = 0; k1 < count_si[nn[j][k]]; k1++) {
		    if (nn[nn[j][k]][k1] == i) {
		      flag = 1;
		    }
		  }
		}
		else {
		  for (k1 = 0; k1 < count_o[nn[j][k]-nSiTotal]; k1++) {
		    if (nn[nn[j][k]][k1] == i) {
		      flag = 1;
		    }
		  }
		}
	      }
	      */

	    }
	    if (flag == 0) {
	      nn[i][count_o[i-nSiTotal]] = j;
	      nn[j][count_si[j]] = i;
	      count_o[i-nSiTotal]++;
	      count_si[j]++;
	    }
	  }      
	}
      }
      d_ss += 0.3;
    }
  }


 

  for (i = 0; i < nSiTotal; i++) {
    cout << "i = " << i << endl;
    nn_hold = nn;
    count_si_hold = count_si;
    count_o_hold = count_o;
    d_ss = 0.0;
    while (count_si[i] < num1) {
      count_si = count_si_hold;
      count_o = count_o_hold;
      nn = nn_hold;
      for (j = 0; j < nSiTotal; j++) {
	flag = 0;
	dij = dist_si_sio2(len, hlen, x, i, j);
    
	if (i != j) {
	  if (count_si[j] < num1 && count_si[i] < num1) {
	    if (dij < d_ss) {
	      for (k = 0; k < count_si[j]; k++) {
		if (nn[j][k] == i) {
		  flag = 1;
		} 
	      
		else {
		  if (nn[j][k] < nSiTotal) {
		    for (k1 = 0; k1 < count_si[nn[j][k]]; k1++) {
		      if (nn[nn[j][k]][k1] == i) {
			flag = 1;
		      }
		    }
		  }
		  else {
		    for (k1 = 0; k1 < count_o[nn[j][k]-nSiTotal]; k1++) {
		      if (nn[nn[j][k]][k1] == i) {
			flag = 1;
		      }
		    }
		  }
		}
	       
	      }
	      
	      if (flag == 0) {
		nn[i][count_si[i]] = j;
		count_si[i]++;
		
		nn[j][count_si[j]] = i;		      
		count_si[j]++;		    
	      }
	    }
	  }
	}	
      }
      d_ss += 0.3;
    }
  } 
} 
