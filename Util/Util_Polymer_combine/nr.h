#ifndef _NR_H_
#define _NR_H_
#include <fstream>
#include <complex>
#include "nrutil.h"
#include "nrtypes.h"

using namespace std;

namespace NR 
{

  DP brent(const DP ax, const DP bx, const DP cx,
	   DP f(const DP), const DP tol, DP &xmin);

  DP f1dim(const DP x);

  void frprmn(Vec_IO_DP &len, Vec_IO_DP &hlen, Vec_IO_INT &id, Mat_IO_DP &x,
	      Mat_IO_INT &nn, Vec_IO_INT &cnt, Vec_IO_INT &elem,
	      Mat_IO_INT &np, Vec_IO_INT &cntnp, int &penalty,
	      const DP ftol, int &iter, 
	      DP &fret, int &MAX1, int &MAX2,  
	      DP func(Vec_IO_DP &, Vec_IO_DP &, Vec_IO_INT &, Mat_IO_DP &,  
		      Mat_IO_INT &, Vec_IO_INT &, Vec_IO_INT &, 
		      Mat_IO_INT &, Vec_IO_INT &, int &), 
	      void dfunc(Vec_IO_DP &, Vec_IO_DP &, Vec_IO_INT &, 
			 Mat_IO_DP &, Mat_IO_INT &, Vec_IO_INT &, 
			 Vec_IO_INT &, Mat_IO_INT &, Vec_IO_INT &, 
			 int &, Mat_IO_DP &));

  void ktForce(Vec_IO_DP &len, Vec_IO_DP &hlen, Vec_IO_INT &idcom,  
	       Mat_IO_DP &x, Mat_IO_INT &nn, Vec_IO_INT &cnt, 
	       Vec_IO_INT &elem, Mat_IO_INT &np, Vec_IO_INT &cntnp, 
	       int &penalty, Mat_IO_DP &g);

  DP ktPotential(Vec_IO_DP &len, Vec_IO_DP &hlen, Vec_IO_INT &idcom,  
		 Mat_IO_DP &x, Mat_IO_INT &nn, Vec_IO_INT &cnt, 
		 Vec_IO_INT &elem, Mat_IO_INT &np, Vec_IO_INT &cntnp, 
		 int &penalty );

  void linmin(Mat_IO_DP &x,  Mat_IO_INT &nn, Vec_IO_INT &cnt, 
	      Mat_IO_INT &np, Vec_IO_INT &cntnp, Mat_IO_DP &xi, DP &fret,
	      DP func(Vec_IO_DP &, Vec_IO_DP &, Vec_IO_INT &, Mat_IO_DP &, 
		      Mat_IO_INT &, Vec_IO_INT &, Vec_IO_INT &, 
		      Mat_IO_INT &, Vec_IO_INT &, int &));

  void mnbrak(DP &ax, DP &bx, DP &cx, DP &fa, DP &fb, DP &fc, 
	      DP func(const DP));

  void swForce(Vec_IO_DP &len, Vec_IO_DP &hlen, Mat_IO_DP &x, 
	       Mat_IO_INT &nn, Vec_IO_INT &cnt, Vec_IO_INT &elem, 
	       Mat_IO_DP &g);

  void swNeighbor(Vec_IO_DP &len, Vec_IO_DP &hlen, Vec_IO_INT &cntnp, 
		  Mat_IO_DP &x, Mat_IO_INT &np);


  DP swPotential(Vec_IO_DP &len, Vec_IO_DP &hlen, Mat_IO_DP &x, 
		 Mat_IO_INT &nn, Vec_IO_INT &cnt, Vec_IO_INT &elem);
}
#endif /* _NR_H_ */







