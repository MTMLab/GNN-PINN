/************************************************************
 *
 *************************************************************/

#include <cmath>
#include "nr.h"
using namespace std;

//
// Global variables.
//
extern int nemax, namax;
extern int na, nnmax;
extern int aSi;
extern int aGe;
extern int aO;
extern int aH;
extern DP E_Si_Si, E_Si_O, E_Ge_Ge, E_Ge_O, E_Si_Ge;

//extern int nQCmax;
//extern int nQC;
//extern Vec_IO_DP E_QC, R_QC;
//
// parameters
//

extern Mat3D_IO_DP R0;
extern Mat3D_IO_DP kb_L;
extern Mat3D_IO_DP kb_R;
extern Mat3D_IO_DP kb_N_L;
extern Mat3D_IO_DP kb_N_R;
extern Mat3D_IO_DP cos0;
extern Mat3D_IO_DP cos0_BALANCE;
extern DP cos0_O_SI_O_SILANONE;
extern Mat3D_IO_DP ka_L;
extern Mat3D_IO_DP ka_R;
extern DP k_ijk_O_Si_O_SILANONE_L;
extern DP k_ijk_O_Si_O_SILANONE_R;
extern Mat3D_IO_DP ka_N_L;
extern Mat3D_IO_DP ka_N_R;
extern DP N_FACTOR_O_Si_O_SILANONE_L;
extern DP N_FACTOR_O_Si_O_SILANONE_R;
extern Mat3D_IO_DP sub;
extern DP gamma2, p_dist;
extern Mat_IO_DP pen; 
extern Vec_IO_DP coord_pen;
extern Vec_IO_DP ring_penalty;


DP NR::ktPotential(Vec_IO_DP &len, Vec_IO_DP &hlen, Mat_IO_DP &x, 
		   Mat_IO_INT &nn, Vec_IO_INT &cnt, Vec_IO_INT &elem, 
		   Mat_IO_INT &np, Vec_IO_INT &cntnp, Vec_IO_INT &id, 
		   int &nid, int &penalty )
{
 int i, j, k, l, ii, iii, jj, kk, ll, i1, i2, i3, i4, k_id, l_id,kkk,kkkkk, n_dihedral, bb_id, ki, bbb_id;

  int nSi, nO, nGe;
  int cnt0, cnt1, cnt2;
  int flag;
  int charge1;
  DP H1, H2, H3, H4, H5;
  DP dx, dy, dz, dxji, dxjk, dyji, dyjk, dzji, dzjk;
  DP rij, rjk, yij, yjk, rji, yji;
  DP cosijk;
  DP E_strain, E_sub, E_pen, E_sys, E_tor;
  int nSiSi, nSiO, nGeGe, nGeO, nSiGe;
  DP dp1; 
  DP coeff_hold;

  Mat_IO_DP r(na, nnmax);
  Mat_IO_DP dxmat(na, nnmax), dymat(na, nnmax), dzmat(na, nnmax);
  Vec_IO_INT oxi(namax), oxi2(namax);
  Mat_IO_INT mat1(namax, nnmax);
  Vec_IO_INT idat(namax);
  
    DP vb1x, vb1y, vb1z, vb2x, vb2y, vb2z, vb2xm, vb2ym, vb2zm, vb3x, vb3y, vb3z, ax, ay, az, bx, by, bz;
  DP rasq, rbsq, rgsq, rg, rginv, ra2inv, rb2inv, rabinv, c, s;

  Vec_IO_DP ax1(3), ay1(3), az1(3), vb3x1(3), vb3y1(3), vb3z1(3), bx1(3), by1(3), bz1(3), cs(9), vb1x1(3),vb1y1(3), vb1z1(3);
  Vec_IO_INT a_id(3), b_id(3);
  Mat_IO_INT cs_id(9,2);
  DP cs_min, c1,c2,c3,c4, s1, s2, s3, s4, theta1, theta2, theta3, theta4;
  int cs_id_min1, cs_id_min2;
  Mat_IO_INT dihedrallist(80000,4);
  Mat_IO_INT nn_dihedral(80000,2);

 
  //
  // Silicon oxidation states
  //
  nSi = nGe = 0;
  oxi2 = 0;
  for (i = 0; i < na; i++) {
    ii = i;
    if (elem[ii] == aSi) {
      nSi++;
      charge1 = 0;
      for (j = 0; j < cnt[ii]; j++) {
	if (elem[nn[ii][j]] == aO) {
	  charge1++;
	}
	if (elem[nn[ii][j]] == aGe) {
	  oxi2[ii]++;
	}
      }
    }
    if (elem[ii] == aGe) {
      nGe++;
      charge1 = 0;
      for (j = 0; j < cnt[ii]; j++) {
	if (elem[nn[ii][j]] == aO) {
	  charge1++;
	}
	if (elem[nn[ii][j]] == aSi) {
	  oxi2[ii]++;
	}
      }      
    }
    oxi[ii] = charge1;    
  }
  //
  // identify atoms
  // 

  
  nSiSi = nSiO = nGeGe = nGeO = nSiGe = 0;
  for (i = 0; i < nid; i++) {
    ii = id[i];
    for (j = 0; j < cnt[ii]; j++) {
      if (elem[ii] == aSi && elem[nn[ii][j]] == aSi) {
	nSiSi++;
      }
      else if (elem[ii] == aSi && elem[nn[ii][j]] == aO 
	       || elem[ii] == aO && elem[nn[ii][j]] == aSi) {
	nSiO++;
      }
      else if (elem[ii] == aGe && elem[nn[ii][j]] == aGe) {
	nGeGe++;
      }
      else if (elem[ii] == aGe && elem[nn[ii][j]] == aO 
	       || elem[ii] == aO && elem[nn[ii][j]] == aGe) {
	nGeO++;
      }	
      else if (elem[ii] == aSi && elem[nn[ii][j]] == aGe 
	       || elem[ii] == aGe && elem[nn[ii][j]] == aSi) {
	nSiGe++;
      }	
    }
  }


  for (i = 0; i < nid; i++) {
    cnt1 = cnt[id[i]];
    ii = id[i];
    //    cout << cnt1 << endl;
    for (j = 0; j < cnt1; j++) {
      jj = nn[ii][j];

      dx = x[ii][0] - x[jj][0];
      if (dx >= -hlen[0] && dx < hlen[0])
	;
      else if (dx < -hlen[0]) {
	dx += len[0];
      }
      else {
	dx -= len[0];  
      }
      dxmat[ii][j] = dx;
      
      dy = x[ii][1] - x[jj][1];
      if (dy >= -hlen[1] && dy < hlen[1])
	;
      else if (dy < -hlen[1]) {
	dy += len[1];
      }
      else {
	dy -= len[1];
      }
      dymat[ii][j] = dy;

      dz = x[ii][2] - x[jj][2];
      if (dz >= -hlen[2] && dz < hlen[2])
	;
      else if (dz < -hlen[2]) {
	dz += len[2];
      }
      else {
	dz -= len[2];
      }
      dzmat[ii][j] = dz;

      rij = sqrt(dx*dx + dy*dy + dz*dz);	
      r[ii][j] = rij;      

      if ((elem[ii] == aSi && elem[jj] == aSi)
	  ||(elem[ii] == aGe && elem[jj] == aGe)) {
	if (oxi[ii] == 0 && oxi[jj] == 0)
	  mat1[ii][j] = 0;
	else if ((oxi[ii] == 0 && oxi[jj] == 1) 
		 || (oxi[ii] == 1 && oxi[jj] == 0)) 
	  mat1[ii][j] = 1;
	else if ((oxi[ii] == 0 && oxi[jj] == 2) 
		 || (oxi[ii] == 2 && oxi[jj] == 0))
	  mat1[ii][j] = 2;
	else if ((oxi[ii] == 0 && oxi[jj] == 3) 
		 || (oxi[ii] == 3 && oxi[jj] == 0))
	  mat1[ii][j] = 3;
	else if (oxi[ii] == 1 && oxi[jj] == 1)
	  mat1[ii][j] = 4;
	else if ((oxi[ii] == 1 && oxi[jj] == 2) 
		 || (oxi[ii] == 2 && oxi[jj] == 1))
	  mat1[ii][j] = 5;
	else if ((oxi[ii] == 1 && oxi[jj] == 3) 
		 || (oxi[ii] == 3 && oxi[jj] == 1))
	  mat1[ii][j] = 6;
	else if (oxi[ii] == 2 && oxi[jj] == 2)
	  mat1[ii][j] = 7;	
	else if ((oxi[ii] == 2 && oxi[jj] == 3) 
		 || (oxi[ii] == 3 && oxi[jj] == 2))
	  mat1[ii][j] = 8;
	else if (oxi[ii] == 3 && oxi[jj] == 3)
	  mat1[ii][j] = 9;
      }
      else if((elem[ii] == aSi && elem[jj] == aGe)
	      ||(elem[ii] == aGe && elem[jj] == aSi)) {
	if (oxi[ii] == 0 && oxi[jj] == 0)
	  mat1[ii][j] = 0;
	else if (oxi[ii] == 0 && oxi[jj] == 1)
	  mat1[ii][j] = 1;
	else if (oxi[ii] == 0 && oxi[jj] == 2)
	  mat1[ii][j] = 2;
	else if (oxi[ii] == 0 && oxi[jj] == 3)
	  mat1[ii][j] = 3;
	else if (oxi[ii] == 1 && oxi[jj] == 0)
	  mat1[ii][j] = 4;
	else if (oxi[ii] == 1 && oxi[jj] == 1)
	  mat1[ii][j] = 5;
	else if (oxi[ii] == 1 && oxi[jj] == 2)
	  mat1[ii][j] = 6;
	else if (oxi[ii] == 1 && oxi[jj] == 3)
	  mat1[ii][j] = 7;
	else if (oxi[ii] == 2 && oxi[jj] == 0)
	  mat1[ii][j] = 8;
	else if (oxi[ii] == 2 && oxi[jj] == 1)
	  mat1[ii][j] = 9;
	else if (oxi[ii] == 2 && oxi[jj] == 2)
	  mat1[ii][j] = 10;
	else if (oxi[ii] == 2 && oxi[jj] == 3)
	  mat1[ii][j] = 11;
	else if (oxi[ii] == 3 && oxi[jj] == 0)
	  mat1[ii][j] = 12;
	else if (oxi[ii] == 3 && oxi[jj] == 1)
	  mat1[ii][j] = 13;
	else if (oxi[ii] == 3 && oxi[jj] == 2)
	  mat1[ii][j] = 14;
	else if (oxi[ii] == 3 && oxi[jj] == 3)
	  mat1[ii][j] = 15;
      }

      else if((elem[ii] == aSi || elem[ii] == aGe) 
	      && (elem[jj] == aO || elem[jj] == aO)) {
	if (oxi[ii] == 1) 
	  mat1[ii][j] = 1;
	else if (oxi[ii] == 2) 
	  mat1[ii][j] = 2;
	else if (oxi[ii] == 3) 
	  mat1[ii][j] = 3;
	else if (oxi[ii] == 4) 
	  mat1[ii][j] = 4;
      }

      else if((elem[ii] == aO || elem[ii] == aO) 
	      && (elem[jj] == aSi || elem[jj] == aGe)) {
	if (oxi[jj] == 1) 
	  mat1[ii][j] = 1;
	else if (oxi[jj] == 2) 
	  mat1[ii][j] = 2;
	else if (oxi[jj] == 3) 
	  mat1[ii][j] = 3;
	else if (oxi[jj] == 4) 
	  mat1[ii][j] = 4;
      }

      else {
	mat1[ii][j] = 0;
      }
    }
  }

  E_sys = 0.0;

  //
  // Suboxide penalty energy.
  //
  H4 = 0.0;
  for (i = 0; i < nid; i++) {
    ii = id[i];
    if (elem[ii] == aSi || elem[ii] == aGe) {
      if (cnt[ii] == 4) {
	H4 += sub[elem[ii]][oxi2[ii]][oxi[ii]];
      }
      else if (cnt[ii] == 0) {
	H4 += coord_pen[0];
      }
      else if (cnt[ii] == 1) {
	H4 += coord_pen[1];
      }
      else if (cnt[ii] == 2) {
	H4 += coord_pen[2];
      }
      else if (cnt[ii] == 3) {
	H4 += coord_pen[3];
      }      
    } 
  }


  //
  // One body interaction (not considered at current stage)
  //
  H1 = 0.0;

  //
  // Two body interaction and suboxide penalty energy.
  //

  H2 = 0.0;

  for (i = 0; i < nid; i++) {
    cnt1 = cnt[id[i]];

    ii = id[i];

    for (j = 0; j < cnt1; j++) {
      jj = nn[ii][j];
      rij = r[ii][j];

      /*
      cout << elem[ii] << "   " 
	   << elem[jj] << "   "
	   << mat1[ii][j] << "   " 
	   << R0[elem[ii]][elem[jj]][mat1[ii][j]] << endl;
      */

      if (rij < R0[elem[ii]][elem[jj]][mat1[ii][j]]) {
	H2 += 0.25*kb_L[elem[ii]][elem[jj]][mat1[ii][j]]
	  *pow(R0[elem[ii]][elem[jj]][mat1[ii][j]] - rij,kb_N_L[elem[ii]][elem[jj]][mat1[ii][j]]);

	/*
	cout << elem[ii] << "   " 
	     << elem[jj] << "   "
	     << mat1[ii][j] << "   " 
	     << kb_N_L[elem[ii]][elem[jj]][mat1[ii][j]] << endl;
	*/	
      }	
      else {
 	H2 += 0.25*kb_R[elem[ii]][elem[jj]][mat1[ii][j]]
	  *pow(rij - R0[elem[ii]][elem[jj]][mat1[ii][j]],kb_N_R[elem[ii]][elem[jj]][mat1[ii][j]]);
	/*
	cout << elem[ii] << "   " 
	     << elem[jj] << "   "
	     << mat1[ii][j] << "   " 
	     << kb_N_R[elem[ii]][elem[jj]][mat1[ii][j]] << endl;
	*/
      }
    }
  }

  //
  // Three body interaction.
  // 
  H3 = 0.0;

  for (j = 0; j < nid; j++) {
    jj = id[j];
    cnt1 = cnt[jj];
    for (i = 0; i < cnt1; i++) {
      ii = nn[jj][i];
      for (k = i+1; k < cnt1; k++) {
	kk = nn[jj][k];
	//
	// Distance between atom jj and atom ii.
	//	
	rji = r[jj][i];
	//
	// Distance between atom jj and atom kk.
	//
	rjk = r[jj][k];

	//
	// cosine therm of three atom ii, jj and  kk with jj vertex. 
	// 
	
	// dxji, dxjk;
	dxji = -dxmat[jj][i];
	dxjk = -dxmat[jj][k];
	
	// dyji, dyjk;
	dyji = -dymat[jj][i];
	dyjk = -dymat[jj][k];
	
	// dzji, dzjk;
	dzji = -dzmat[jj][i];
	dzjk = -dzmat[jj][k];
	cosijk = (dxji*dxjk + dyji*dyjk + dzji*dzjk) / (rji*rjk);

	/*
	if ((elem[jj] == aSi || elem[jj] == aGe) && cnt[jj] == 3) {
	  if (cnt[ii] != 1 && cnt[kk] != 1) {   
	    if (cosijk < cos0_O_SI_O_SILANONE) {
	      H3 += 0.5*k_ijk_O_Si_O_SILANONE_L
		*pow(cos0_O_SI_O_SILANONE - cosijk, 
		     N_FACTOR_O_Si_O_SILANONE_L);	  	    
	    }
	    else {
	      H3 += 0.5*k_ijk_O_Si_O_SILANONE_R
		*pow(cosijk - cos0_O_SI_O_SILANONE, 
		     N_FACTOR_O_Si_O_SILANONE_R);
	    }
	  }
	}
	*/
	//	else {
	if (cosijk < cos0_BALANCE[elem[ii]][elem[jj]][elem[kk]]) {
	  if (cosijk < cos0[elem[ii]][elem[jj]][elem[kk]]) {
	    H3 += 0.5*ka_L[elem[ii]][elem[jj]][elem[kk]]
	      *pow(cos0[elem[ii]][elem[jj]][elem[kk]] - cosijk, 
		   ka_N_L[elem[ii]][elem[jj]][elem[kk]]);	  
	  }
	  else {
	    H3 += 0.5*ka_L[elem[ii]][elem[jj]][elem[kk]]
	      *pow(cosijk - cos0[elem[ii]][elem[jj]][elem[kk]], 
		   ka_N_L[elem[ii]][elem[jj]][elem[kk]]);	  
	  }
	}
	else {
	  if (cosijk < cos0[elem[ii]][elem[jj]][elem[kk]]) {	  
	    H3 += 0.5*ka_R[elem[ii]][elem[jj]][elem[kk]]
	      *pow(cos0[elem[ii]][elem[jj]][elem[kk]] - cosijk, 
		   ka_N_R[elem[ii]][elem[jj]][elem[kk]]);	  	    
	  }
	  else {
	    H3 += 0.5*ka_R[elem[ii]][elem[jj]][elem[kk]]
	      *pow(cosijk - cos0[elem[ii]][elem[jj]][elem[kk]], 
		   ka_N_R[elem[ii]][elem[jj]][elem[kk]]);	  	    
	  }
	}
	//	}
      }       
    }
  }

  //  cout << "H2: " << H2 << endl;
  //  cout << "H3: " << H3 << endl;
  //  cout << "H4: " << H4 << endl;

  E_strain = H2 + H3;

  E_sub = H4;
  H5 = 0.0;
 // if (penalty == 0) 
    E_sys = E_strain + E_sub; 
 /*
  else if (penalty == 1) {
    H5 = 0.0; 
    for (i = 0; i < nid; i++) {
      ii = id[i];
      cnt2 = cntnp[ii];

      for (j = 0; j < cnt2; j++) {
	jj = np[ii][j];
	flag = 0;
	cnt0 = cnt[ii];
	for (k = 0; k < cnt0; k++) {

	  kk = nn[ii][k];
	  if (jj == kk)
	    flag = 1;
	  else {
	    cnt1 = cnt[kk];
	    for (l = 0; l < cnt1; l++) {

	      ll = nn[kk][l];
	      if (jj == ll)
		flag = 1;
	    }
	  }
	}
	
	if (flag == 0) {
	  
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
	  rij = sqrt(dx*dx + dy*dy + dz*dz);	

	  if (rij < pen[elem[ii]][elem[jj]]) {

	    H5 += 0.5*gamma2*pow(pen[elem[ii]][elem[jj]]-rij,3);
	  }
	  else {
	    ;
	  }
	} 	
      }
    }
    E_pen = H5;
    
    E_sys = E_strain + E_pen + E_sub ;

    //    cout << E_pen << endl;
  }


    */
    E_tor = 0;
    
    // torsion
  //  cout << "tor" << endl;
 

 
 kkkkk = 0;
 dihedrallist = 0;
 n_dihedral = 0;
   for (i = 0; i < nid; i++) {
    ii = id[i];
    if ( elem[ii] == aGe ) { 
 //if (elem[nn[ii][0]] == aGe && elem[nn[ii][1]] == aGe && elem[nn[ii][2]] == aGe && elem[nn[ii][3]] == aGe) {
       for (j = 0; j < 4; j++) {        
           jj = nn[ii][j];
     if ( elem[jj] == aGe ) { 
 if ( (elem[nn[jj][0]] == aGe && elem[nn[jj][1]] == aGe && elem[nn[jj][2]] == aGe && elem[nn[jj][3]] == aGe ) || (elem[nn[jj][0]] == aSi && elem[nn[jj][1]] == aGe && elem[nn[jj][2]] == aGe && elem[nn[jj][3]] == aGe ) || (elem[nn[jj][0]] == aGe && elem[nn[jj][1]] == aSi && elem[nn[jj][2]] == aGe && elem[nn[jj][3]] == aGe ) || (elem[nn[jj][0]] == aGe && elem[nn[jj][1]] == aGe && elem[nn[jj][2]] == aSi && elem[nn[jj][3]] == aGe ) || (elem[nn[jj][0]] == aGe && elem[nn[jj][1]] == aGe && elem[nn[jj][2]] == aGe && elem[nn[jj][3]] == aSi ) ) {
 //   if (elem[nn[jj][0]] == aGe && elem[nn[jj][1]] == aGe && elem[nn[jj][2]] == aGe && elem[nn[jj][3]] == aGe) {   
//	if ( ii > jj ) continue;
	//	n_dihedral = n_dihedral + 1;
        nn_dihedral[n_dihedral][0] = ii;
        nn_dihedral[n_dihedral][1] = jj;
		    n_dihedral = n_dihedral + 1;
		 //   	    cout << ii << " : " << jj << endl;
      }
      }
    }
   }
}
//cout << "AA1" << endl;
   for ( iii = 0; iii < n_dihedral; iii++ ) {

     ii = nn_dihedral[iii][0];
     jj = nn_dihedral[iii][1];   
                    i2 = ii;
                    i3 = jj;
                    
                    for ( k = 0; k < 4; k++ ) {
                         if ( nn[ii][k] == jj ) {
                              k_id = k;
                              }
                              }   
                              
                    for ( ki = 0; ki < 4; ki++ ) {
                         if ( nn[jj][ki] == ii ) {
                              l_id = ki;
                              }
                              }               
		    k = 0;                                                                   
                            for ( l = 0; l < 4; l++ ) {
                         if ( nn[ii][l] == jj ) continue;
                           i1 = nn[ii][l];
        
    vb1x1[k] = x[i1][0] - x[i2][0];
    vb1y1[k] = x[i1][1] - x[i2][1];
    vb1z1[k] = x[i1][2] - x[i2][2];
         if (vb1x1[k] >= -hlen[0] && vb1x1[k] < hlen[0])	
       ;
      else if (vb1x1[k] < -hlen[0]) {
	vb1x1[k] += len[0];
      }
      else {
	vb1x1[k] -= len[0];  
      }

      if (vb1y1[k]  >= -hlen[1] && vb1y1[k]  < hlen[1])
	;
      else if (vb1y1[k]  < -hlen[1]) {
	vb1y1[k]  += len[1];
      }
      else {
	vb1y1[k]  -= len[1];
      }
           if (vb1z1[k] >= -hlen[2] && vb1z1[k] < hlen[2])
	;
      else if (vb1z1[k] < -hlen[2]) {
	vb1z1[k] += len[2];
      }
      else {
	vb1z1[k] -= len[2];
      }

    // 2nd bond

    vb2x = x[i3][0] - x[i2][0];
    vb2y = x[i3][1] - x[i2][1];
    vb2z = x[i3][2] - x[i2][2];
    if (vb2x >= -hlen[0] && vb2x < hlen[0])	
       ;
      else if (vb2x < -hlen[0]) {
	vb2x += len[0];
      }
      else {
	vb2x -= len[0];  
      }

      if (vb2y  >= -hlen[1] && vb2y  < hlen[1])
	;
      else if (vb2y  < -hlen[1]) {
	vb2y  += len[1];
      }
      else {
	vb2y  -= len[1];
      }
           if (vb2z >= -hlen[2] && vb2z < hlen[2])
	;
      else if (vb2z < -hlen[2]) {
	vb2z += len[2];
      }
      else {
	vb2z -= len[2];
      }

    vb2xm = -vb2x;
    vb2ym = -vb2y;
    vb2zm = -vb2z;

    // 3rd bond
    
    // c,s calculation

    ax1[k] = vb1y1[k]*vb2zm - vb1z1[k]*vb2ym;
    ay1[k] = vb1z1[k]*vb2xm - vb1x1[k]*vb2zm;
    az1[k] = vb1x1[k]*vb2ym - vb1y1[k]*vb2xm;
     
    a_id[k] = i1;
    k = k + 1;
			    } 

			    k = 0;
                l = 0;
    vb3x = 0;
    vb3y = 0;
    vb3z = 0;
      for ( l = 0; l < 4; l++ ) {
            if ( nn[jj][l] == ii ) continue;
            i4 = nn[jj][l];   
            
            
                vb2x = x[i3][0] - x[i2][0];
    vb2y = x[i3][1] - x[i2][1];
    vb2z = x[i3][2] - x[i2][2];
    if (vb2x >= -hlen[0] && vb2x < hlen[0])	
       ;
      else if (vb2x < -hlen[0]) {
	vb2x += len[0];
      }
      else {
	vb2x -= len[0];  
      }

      if (vb2y  >= -hlen[1] && vb2y  < hlen[1])
	;
      else if (vb2y  < -hlen[1]) {
	vb2y  += len[1];
      }
      else {
	vb2y  -= len[1];
      }
           if (vb2z >= -hlen[2] && vb2z < hlen[2])
	;
      else if (vb2z < -hlen[2]) {
	vb2z += len[2];
      }
      else {
	vb2z -= len[2];
      }

    vb2xm = -vb2x;
    vb2ym = -vb2y;
    vb2zm = -vb2z;
                        
            
     vb3x1[k] = x[i4][0] - x[i3][0];
    vb3y1[k] = x[i4][1] - x[i3][1];
    vb3z1[k] = x[i4][2] - x[i3][2];
    if (vb3x1[k] >= -hlen[0] && vb3x1[k] < hlen[0])	
       ;
      else if (vb3x1[k] < -hlen[0]) {
	vb3x1[k] += len[0];
      }
      else {
	vb3x1[k] -= len[0];  
      }

      if (vb3y1[k]  >= -hlen[1] && vb3y1[k]  < hlen[1])
	;
      else if (vb3y1[k]  < -hlen[1]) {
	vb3y1[k]  += len[1];
      }
      else {
	vb3y1[k]  -= len[1];
      }
           if (vb3z1[k] >= -hlen[2] && vb3z1[k] < hlen[2])
	;
      else if (vb3z1[k] < -hlen[2]) {
	vb3z1[k] += len[2];
      }
      else {
	vb3z1[k] -= len[2];
      }   
          
    
    bx1[k] = vb3y1[k]*vb2zm - vb3z1[k]*vb2ym;
    by1[k] = vb3z1[k]*vb2xm - vb3x1[k]*vb2zm;
    bz1[k] = vb3x1[k]*vb2ym - vb3y1[k]*vb2xm;
    b_id[k] = i4;
    k = k + 1;
      }

//  cout << "CCCC" << endl;


kkk = 0;
    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            
   // abx = ay[i]*bz[j] - az[i]*by[j];
  //  aby = az[i]*bx[j] - ax[i]*bz[j];
  //  abz = ax[i]*by[j] - ay[i]*bx[j];
    
 //   torsion = (abx*vb2xm + aby*vb2ym + abz*vb2zm);    
    
    rasq = ax1[i]*ax1[i] + ay1[i]*ay1[i] + az1[i]*az1[i];
    rbsq = bx1[j]*bx1[j] + by1[j]*by1[j] + bz1[j]*bz1[j];
    rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
    rg = sqrt(rgsq);
    
    rginv = ra2inv = rb2inv = 0.0;
    rginv = 1.0/rg;
    ra2inv = 1.0/rasq;
    rb2inv = 1.0/rbsq;
    rabinv = sqrt(ra2inv*rb2inv);

    c = (ax1[i]*bx1[j] + ay1[i]*by1[j] + az1[i]*bz1[j])*rabinv;
    s = rg*rabinv*(ax1[i]*vb3x1[j] + ay1[i]*vb3y1[j] + az1[i]*vb3z1[j]);
    cs[kkk] = c;
    cs_id[kkk][0] = a_id[i];
    cs_id[kkk][1] = b_id[j];
    kkk = kkk + 1;    
	}
    } // l;
    // cout << "CCCC" << endl;
   cs_min = cs[0];
   cs_id_min1 = cs_id[0][0];
   cs_id_min2 = cs_id[0][1];
   bb_id = 0;
     if ( cs[1] < cs_min ) {
          cs_min = cs[1];
          cs_id_min1 = cs_id[1][0];
          cs_id_min2 = cs_id[1][1];
          bb_id = 1;
          }
          if ( cs[2] < cs_min ) {
          cs_min = cs[2];
          cs_id_min1 = cs_id[2][0];
          cs_id_min2 = cs_id[2][1];
          bb_id = 2;
          }       
   //   cout << "cs : " << cs_min <<  " : " << cs[0] << " : " << cs[1] << " : " << cs[2] << endl;      
        //     cout << "DDDD" << endl;                                          
   dihedrallist[kkkkk][0] = cs_id_min1;
   dihedrallist[kkkkk][1] = i2;
   dihedrallist[kkkkk][2] = i3;
   dihedrallist[kkkkk][3] = cs_id_min2;
   
   rasq = ax1[0]*ax1[0] + ay1[0]*ay1[0] + az1[0]*az1[0];
   rbsq = ax1[1]*ax1[1] + ay1[1]*ay1[1] + az1[1]*az1[1];
   rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
   rg = sqrt(rgsq);
   rginv = ra2inv = rb2inv = 0.0;
   rginv = 1.0/rg;
   ra2inv = 1.0/rasq;
   rb2inv = 1.0/rbsq;
   rabinv = sqrt(ra2inv*rb2inv); 
   c1 = (ax1[0]*ax1[1] + ay1[0]*ay1[1] + az1[0]*az1[1])*rabinv;
   s1 = rg*rabinv*(ax1[0]*vb1x1[1] + ay1[0]*vb1y1[1] + az1[0]*vb1z1[1]);
   
   theta1 = acos(c1)*180/3.141592;
   if ( s1 < 0 ) {
        theta1 = 360 - theta1;
        }
   
   rasq = ax1[0]*ax1[0] + ay1[0]*ay1[0] + az1[0]*az1[0];
   rbsq = ax1[2]*ax1[2] + ay1[2]*ay1[2] + az1[2]*az1[2];
   rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
   rg = sqrt(rgsq);
   rginv = ra2inv = rb2inv = 0.0;
   rginv = 1.0/rg;
   ra2inv = 1.0/rasq;
   rb2inv = 1.0/rbsq;
   rabinv = sqrt(ra2inv*rb2inv); 
   c2 = (ax1[0]*ax1[2] + ay1[0]*ay1[2] + az1[0]*az1[2])*rabinv;
   s2 = rg*rabinv*(ax1[0]*vb1x1[2] + ay1[0]*vb1y1[2] + az1[0]*vb1z1[2]);
   
   theta2 = acos(c2)*180/3.141592;
   if ( s2 < 0 ) {
        theta2 = 360 - theta2;
        }       
  
   //  cout << c2 << endl;  

    
  bbb_id = bb_id + 1;
  if ( bbb_id > 2 ) {
   bbb_id = bbb_id - 3;
  }

   rasq = bx1[bb_id]*bx1[bb_id] + by1[bb_id]*by1[bb_id] + bz1[bb_id]*bz1[bb_id];
   rbsq = bx1[bbb_id]*bx1[bbb_id] + by1[bbb_id]*by1[bbb_id] + bz1[bbb_id]*bz1[bbb_id];
   rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
   rg = sqrt(rgsq);
   rginv = ra2inv = rb2inv = 0.0;
   rginv = 1.0/rg;
   ra2inv = 1.0/rasq;
   rb2inv = 1.0/rbsq;
   rabinv = sqrt(ra2inv*rb2inv); 
   c3 = (bx1[bb_id]*bx1[bbb_id] + by1[bb_id]*by1[bbb_id] + bz1[bb_id]*bz1[bbb_id])*rabinv;
   s3 = rg*rabinv*(bx1[bb_id]*vb3x1[bbb_id] + by1[bb_id]*vb3y1[bbb_id] + bz1[bb_id]*vb3z1[bbb_id]);
   
   theta3 = acos(c3)*180/3.141592;
   if ( s3 < 0 ) {
        theta3 = 360 - theta3;
        }
   //  cout << theta3 << endl;

   bbb_id = bbb_id + 1;
   if ( bbb_id > 2 ) { 
   bbb_id = bbb_id - 3;
   }

   rasq = bx1[bb_id]*bx1[bb_id] + by1[bb_id]*by1[bb_id] + bz1[bb_id]*bz1[bb_id];
   rbsq = bx1[bbb_id]*bx1[bbb_id] + by1[bbb_id]*by1[bbb_id] + bz1[bbb_id]*bz1[bbb_id];
   rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
   rg = sqrt(rgsq);
   rginv = ra2inv = rb2inv = 0.0;
   rginv = 1.0/rg;
   ra2inv = 1.0/rasq;
   rb2inv = 1.0/rbsq;
   rabinv = sqrt(ra2inv*rb2inv); 
   c4 = (bx1[bb_id]*bx1[bbb_id] + by1[bb_id]*by1[bbb_id] + bz1[bb_id]*bz1[bbb_id])*rabinv;
   s4 = rg*rabinv*(bx1[bb_id]*vb3x1[bbb_id] + by1[bb_id]*vb3y1[bbb_id] + bz1[bb_id]*vb3z1[bbb_id]);
    
   theta4 = acos(c4)*180/3.141592;
   // cout << c1 << " :"  << c2 << " : " << c3 << " : " << c4 << endl;
   if ( s4 < 0 ) {
        theta4 = 360 - theta4;
        }    
   
   if ( theta2 > theta1 ) {
        if ( theta4 > theta3 ) {
         dihedrallist[kkkkk+1][0] = a_id[1];
         dihedrallist[kkkkk+1][1] = i2;
         dihedrallist[kkkkk+1][2] = i3;
         bb_id = bb_id + 1;
         if ( bb_id > 2 ) {
	   bb_id = bb_id - 3;
         }
         dihedrallist[kkkkk+1][3] = b_id[bb_id];
         
         dihedrallist[kkkkk+2][0] = a_id[2];
         dihedrallist[kkkkk+2][1] = i2;
         dihedrallist[kkkkk+2][2] = i3;
         bb_id = bb_id + 1; 
         if ( bb_id > 2 ) {
	   bb_id = bb_id - 3;
         }
         dihedrallist[kkkkk+2][3] = b_id[bb_id];            
         }
         if ( theta4 < theta3 ) {
         dihedrallist[kkkkk+1][0] = a_id[1];
         dihedrallist[kkkkk+1][1] = i2;
         dihedrallist[kkkkk+1][2] = i3;
         bb_id = bb_id - 1;
         if ( bb_id < 0 ) {
	   bb_id = bb_id + 3;
         }
         dihedrallist[kkkkk+1][3] = b_id[bb_id];
         
         dihedrallist[kkkkk+2][0] = a_id[2];
         dihedrallist[kkkkk+2][1] = i2;
         dihedrallist[kkkkk+2][2] = i3;
         bb_id = bb_id - 1;
         if ( bb_id < 0 ) {
	   bb_id = bb_id + 3;
         }
         dihedrallist[kkkkk+2][3] = b_id[bb_id];            
         }
   }       
     if ( theta1 > theta2 ) {
        if ( theta4 > theta3 ) {
         dihedrallist[kkkkk+1][0] = a_id[2];
         dihedrallist[kkkkk+1][1] = i2;
         dihedrallist[kkkkk+1][2] = i3;
         bb_id = bb_id + 1;
         if ( bb_id > 2 ) {
	   bb_id = bb_id - 3;
         }
         dihedrallist[kkkkk+1][3] = b_id[bb_id];
         
         dihedrallist[kkkkk+2][0] = a_id[1];
         dihedrallist[kkkkk+2][1] = i2;
         dihedrallist[kkkkk+2][2] = i3;
         bb_id = bb_id + 1;
         if ( bb_id > 2 ) {
	   bb_id = bb_id - 3;
         }
         dihedrallist[kkkkk+2][3] = b_id[bb_id];            
         }
         if ( theta4 < theta3 ) {
         dihedrallist[kkkkk+1][0] = a_id[2];
         dihedrallist[kkkkk+1][1] = i2;
         dihedrallist[kkkkk+1][2] = i3;
         bb_id = bb_id - 1;
         if ( bb_id < 0 ) {
	   bb_id = bb_id + 3;
         }
         dihedrallist[kkkkk+1][3] = b_id[bb_id];
         
         dihedrallist[kkkkk+2][0] = a_id[1];
         dihedrallist[kkkkk+2][1] = i2;
         dihedrallist[kkkkk+2][2] = i3;
         bb_id = bb_id - 1;
         if ( bb_id < 0 ) {
	   bb_id = bb_id + 3;
         }
         dihedrallist[kkkkk+2][3] = b_id[bb_id];            
         }
     }   
  
// cout << iii << " : " << theta1 << " : " << theta2 << " : " << theta3 << " : " << theta4 << endl; 
      kkkkk = kkkkk + 3;   
   //  cout << dihedrallist[0][0] << << endl;
     // cout << dihedrallist[0][0] << " : " << dihedrallist[0][1] << " : " << dihedrallist[0][2] << " : " << dihedrallist[0][3];
          
   }   // k;

//   cout << "DDDD" << " : " << n_dihedral << endl;  
//   cout << dihedrallist[0][0] << " : " << dihedrallist[0][1] << " : " << dihedrallist[0][2] << " : " << dihedrallist[0][3];
//cout << dihedrallist[1][0] << " : " << dihedrallist[1][1] << " : " << dihedrallist[1][2] << " : " << dihedrallist[1][3];
//cout << dihedrallist[2][0] << " : " << dihedrallist[2][1] << " : " << dihedrallist[2][2] << " : " << dihedrallist[2][3] << endl;
//   cout << dihedrallist[3][0] << " : " << dihedrallist[3][1] << " : " << dihedrallist[3][2] << " : " << dihedrallist[3][3];
//cout << dihedrallist[4][0] << " : " << dihedrallist[4][1] << " : " << dihedrallist[4][2] << " : " << dihedrallist[4][3];
//cout << dihedrallist[5][0] << " : " << dihedrallist[5][1] << " : " << dihedrallist[5][2] << " : " << dihedrallist[5][3] << endl;
//    cout << "AAAA" << endl;
  
// cout << "AA2" << endl;
 
  for (i = 0; i < kkkkk; i++) {
    i1 = dihedrallist[i][0];
    i2 = dihedrallist[i][1];
    i3 = dihedrallist[i][2];
    i4 = dihedrallist[i][3];
 
             
    vb1x = x[i1][0] - x[i2][0];
    vb1y = x[i1][1] - x[i2][1];
    vb1z = x[i1][2] - x[i2][2];
         if (vb1x >= -hlen[0] && vb1x < hlen[0])	
       ;
      else if (vb1x < -hlen[0]) {
	vb1x += len[0];
      }
      else {
	vb1x -= len[0];  
      }

      if (vb1y  >= -hlen[1] && vb1y  < hlen[1])
	;
      else if (vb1y  < -hlen[1]) {
	vb1y  += len[1];
      }
      else {
	vb1y  -= len[1];
      }
           if (vb1z >= -hlen[2] && vb1z < hlen[2])
	;
      else if (vb1z < -hlen[2]) {
	vb1z += len[2];
      }
      else {
	vb1z -= len[2];
      }

    // 2nd bond

    vb2x = x[i3][0] - x[i2][0];
    vb2y = x[i3][1] - x[i2][1];
    vb2z = x[i3][2] - x[i2][2];
    if (vb2x >= -hlen[0] && vb2x < hlen[0])	
       ;
      else if (vb2x < -hlen[0]) {
	vb2x += len[0];
      }
      else {
	vb2x -= len[0];  
      }

      if (vb2y  >= -hlen[1] && vb2y  < hlen[1])
	;
      else if (vb2y  < -hlen[1]) {
	vb2y  += len[1];
      }
      else {
	vb2y  -= len[1];
      }
           if (vb2z >= -hlen[2] && vb2z < hlen[2])
	;
      else if (vb2z < -hlen[2]) {
	vb2z += len[2];
      }
      else {
	vb2z -= len[2];
      }

    vb2xm = -vb2x;
    vb2ym = -vb2y;
    vb2zm = -vb2z;

    // 3rd bond

    vb3x = x[i4][0] - x[i3][0];
    vb3y = x[i4][1] - x[i3][1];
    vb3z = x[i4][2] - x[i3][2];
    if (vb3x >= -hlen[0] && vb3x < hlen[0])	
       ;
      else if (vb3x < -hlen[0]) {
	vb3x += len[0];
      }
      else {
	vb3x -= len[0];  
      }

      if (vb3y  >= -hlen[1] && vb3y  < hlen[1])
	;
      else if (vb3y  < -hlen[1]) {
	vb3y  += len[1];
      }
      else {
	vb3y  -= len[1];
      }
           if (vb3z >= -hlen[2] && vb3z < hlen[2])
	;
      else if (vb3z < -hlen[2]) {
	vb3z += len[2];
      }
      else {
	vb3z -= len[2];
      }
    
    // c,s calculation

    ax = vb1y*vb2zm - vb1z*vb2ym;
    ay = vb1z*vb2xm - vb1x*vb2zm;
    az = vb1x*vb2ym - vb1y*vb2xm;
    bx = vb3y*vb2zm - vb3z*vb2ym;
    by = vb3z*vb2xm - vb3x*vb2zm;
    bz = vb3x*vb2ym - vb3y*vb2xm;
    
    rasq = ax*ax + ay*ay + az*az;
    rbsq = bx*bx + by*by + bz*bz;
    rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
    rg = sqrt(rgsq);
    
    rginv = ra2inv = rb2inv = 0.0;
    rginv = 1.0/rg;
    ra2inv = 1.0/rasq;
    rb2inv = 1.0/rbsq;
    rabinv = sqrt(ra2inv*rb2inv);

    c = (ax*bx + ay*by + az*bz)*rabinv;
    s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);
    
  //  cout << "c : " << c << endl;
    
    E_tor = E_tor + 5 *(1+c);
    
}
//  cout << "AA3" << endl;                                                 
   if ( penalty == 1) {
        E_tor = 0;
  }                                   
    
    E_sys = E_sys+ E_tor ;
//cout << "TOR" << endl;

  E_sys = E_sys;

  // + (0.5*double(nSiO)*E_Si_O+0.5*double(nSiSi)*E_Si_Si) 
  //    + (0.5*double(nGeO)*E_Ge_O+0.5*double(nGeGe)*E_Ge_Ge) 
  //    + (0.5*double(nSiGe)*E_Si_Ge);      

  return E_sys;
}










































































































