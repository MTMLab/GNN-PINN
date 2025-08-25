/************************************************************
 * January 03, 2006 
 * Name: Utility (xmol --> POSCAR)
 * File name: MAIN.cpp 
 * Developers: Sangheon Lee and Gyeong S. Hwang
 * 
 * 
 *************************************************************/

#include "nr.h"
#include <iostream>
#include <cmath>
using std::cout;
using std::cin;
using std::endl;
using std::ios;
using std::cerr; 
#include <fstream>
using std::ofstream;
using std::ifstream;
#include <iomanip>
using std::setw;

#include <cstdlib>
#include <ctime>
#include <string>
using std::string;

using namespace std;

const int numC = 3;
const int numN = 4;
const DP pi = 3.141592654;
const int aSi = 0;
const int aC = 1;
const int aO = 2;
const int aH = 3;
const int aCu = 4;
const int aB = 5;
const int aP = 6;
const int aS = 7;
const int aCl = 8;
const int aN = 9;

//
// Global variables.
//
int namax; // Maximum number of atoms.
int nemax = 10; // Maximum number of elements.
int na; // Number of total atoms in the system.
int nnmax = 4; // Maximum number of neighbors.
DP cutoff1;
int penalty;
int npmax; 

int main()  
{
  //  cout << "Enter number of atoms:" << endl; 
  //  cin >> namax;
  namax = 300000;
  na = namax;
  Mat_IO_INT nn(namax, numN), b_id0(namax, 4), a_id0(namax, 5), b_id1(namax, 4), a_id1(namax, 5), a_dih0(namax,6), a_dih1(namax,6), a_im0(namax,6), a_im1(namax,6);
  Mat_IO_DP x(namax, numC), len(3,3), x0(namax, numC), x1(namax, numC);
  Mat_IO_DP b_coeff0(namax, 5), a_coeff0(namax, 5), d_coeff0(namax, 5), i_coeff0(namax, 5);
  Mat_IO_DP b_coeff1(namax, 5), a_coeff1(namax, 5), d_coeff1(namax, 5), i_coeff1(namax, 5);
  Vec_IO_INT elem(namax), cnt(namax), idb1(namax), idr1(namax), ide1(namax), id0(namax), id1(namax), atom_type0(namax), atom_type1(namax), mol_id0(namax), mol_id1(namax);
  Vec_IO_DP hlen(numC), cen(numC), q0(namax), q1(namax), mass1(namax), mass0(namax), l_lo(numC), l_hi(numC);
  Vec_IO_DP lj0_e(namax), lj0_s(namax), lj1_e(namax), lj1_s(namax), lj2_e(namax), lj2_s(namax);
  int i, j, k, l, m, m1, m2, m3, m4, m5, m6, kkk, kkkk, mol_id_0, mol_id_1, klkl;
  int idum;
  int logic1;
  int charge1;
  int nE1, nE2, nE3, nE4, nE5, nE6, nE7, nE8, nE9, nE10;
  DP ratio1, pp, sum_x, sum_y, sum_z, sum_x_m, sum_y_m, sum_z_m, avg_x1, avg_x2, avg_y1, avg_y2, avg_z1, avg_z2, avg_x_m, avg_y_m, avg_z_m;
  DP min_x1, max_x1, min_x_m, max_x_m;
  int n_type0, n_type1, n_bond0, n_bond1, n_bond_type0, n_bond_type1, n_angle0, n_angle1, n_angle_type0, n_angle_type1;
  int n_dihedral0, n_dihedral_type0, n_dihedral1, n_dihedral_type1, n_improper0, n_improper_type0, n_improper1, n_improper_type1;

  string Si("Si");
  string C("C");
  string Cu("Cu");
  string O("O");
  string H("H");
  string N("N");
  string B("B");
  string P("P");
  string S("S");
  string Cl("Cl");
  string string1;
  string system1;
  char char1,char2,char3,char4,char5,char6,char7,char8,char9,char10,char11,char12;
  DP freq, len_x, len_y, len_z;


  elem = -1;
  
  //  cout << "Box length: " << endl;
  //  cin >> len[0] >> len[1] >> len[2];
  

  n_improper0 = 0;
  n_improper_type0 = 0;
  // 
  // Read an initial configuration.
  //

  ifstream inputFile3333("system.data",ios::in);

  inputFile3333 >> string1 >> string1;
  inputFile3333 >> kkk >> string1;
    inputFile3333 >> n_bond0 >> string1;
  inputFile3333 >> n_angle0 >> string1;
    inputFile3333 >> n_dihedral0 >> string1;
    inputFile3333 >> n_improper0 >> string1;
    
  inputFile3333 >> n_type0 >> string1 >> string1;
  inputFile3333 >> n_bond_type0 >> string1 >> string1;
  inputFile3333 >> n_angle_type0 >> string1 >> string1;
  inputFile3333 >> n_dihedral_type0 >> string1 >> string1;
  inputFile3333 >> n_improper_type0 >> string1 >> string1;

  inputFile3333 >> l_lo[0] >> l_hi[0] >> string1 >> string1;
  inputFile3333 >> l_lo[1] >> l_hi[1]  >> string1 >> string1;
  inputFile3333 >> l_lo[2] >> l_hi[2]  >> string1 >> string1;

  inputFile3333 >> string1;
  
    for ( i = 0; i < n_type0; i++ ) {
  inputFile3333 >> string1 >> mass0[i] >> string1 >> string1;
  }

  inputFile3333 >> string1 >> string1 >> string1;
 
  for (i = 0; i < kkk; i++ ) {

    inputFile3333 >> id0[i] >> string1 >> atom_type0[i] >> q0[i] >> x0[i][0] >> x0[i][1] >> x0[i][2];

   }

  inputFile3333 >> string1;

  for ( i = 0; i < n_bond0; i++ ) {
  
   inputFile3333 >> b_id0[i][0] >> b_id0[i][1] >> b_id0[i][2] >> b_id0[i][3];
 
  }


  inputFile3333 >> string1;

  for ( i = 0; i < n_angle0; i++ ) {

   inputFile3333 >> a_id0[i][0] >> a_id0[i][1] >> a_id0[i][2] >> a_id0[i][3] >> a_id0[i][4];

  }
   
    inputFile3333 >> string1;

  for ( i = 0; i < n_dihedral0; i++ ) {

   inputFile3333 >> a_dih0[i][0] >> a_dih0[i][1] >> a_dih0[i][2] >> a_dih0[i][3] >> a_dih0[i][4] >> a_dih0[i][5];

  }
  
    inputFile3333 >> string1;

  for ( i = 0; i < n_improper0; i++ ) {

   inputFile3333 >> a_im0[i][0] >> a_im0[i][1] >> a_im0[i][2] >> a_im0[i][3] >> a_im0[i][4] >> a_im0[i][5];

  }

    
    for ( i = 0; i < kkk/5; i++) {
        mol_id0[i] = 1;
    }

    for ( i = kkk/5; i < 2*kkk/5; i++) {
        mol_id0[i] = 2;
    }
    
    for ( i = 2*kkk/5; i < 3*kkk/5; i++) {
        mol_id0[i] = 3;
    }
    
    for ( i = 3*kkk/5; i < 4*kkk/5; i++) {
        mol_id0[i] = 4;
    }
    
    for ( i = 4*kkk/5; i < kkk; i++) {
        mol_id0[i] = 5;
    }
    
    

  ofstream outputFile5("system2.data", ios::out);
  if (!outputFile5) {
    cerr << "POSCAR could not be opened" << endl;
    exit(1);
  }

outputFile5 << "LAMMPS Description" << endl;
outputFile5 << endl;
outputFile5 << kkk<< " atoms" << endl;
    outputFile5 << n_bond0 << " bonds" << endl;
    outputFile5 << n_angle0 << " angles" << endl;
    outputFile5 << n_dihedral0 << " dihedrals" << endl;
    outputFile5 << n_improper0 << " impropers" << endl;
    
outputFile5 << n_type0 <<" atom types" << endl;
outputFile5 << n_bond_type0 << " bond types" << endl;
outputFile5 << n_angle_type0 << " angle types" << endl;
outputFile5 << n_dihedral_type0 << " dihedral types" << endl;
outputFile5 << n_improper_type0 << " improper types" << endl;
    
outputFile5 << endl;
    
outputFile5 << l_lo[0] << " " << l_hi[0] << " xlo xhi" << endl;
outputFile5 << l_lo[1] << " " << l_hi[1] << " ylo yhi" << endl;
outputFile5 << l_lo[2] << " " << l_hi[2] << " zlo zhi" << endl;
outputFile5 << endl;
outputFile5 << "Masses" << endl;
outputFile5 << endl;
    
for ( i = 0; i < n_type0; i++ ) {
outputFile5 << i+1 << " " << mass0[i] << endl;
}

outputFile5 << endl;

outputFile5 << "Atoms # full" << endl;
outputFile5 << endl;

  for ( i = 0; i < kkk; i++ ) {

    outputFile5 << id0[i] << "   " << mol_id0[i] << "   " << atom_type0[i] << "   " << q0[i] << "   " << x0[i][0] << "   " << x0[i][1] << "   " << x0[i][2] << endl;

  }

outputFile5 << endl;
outputFile5 << "Bonds" << endl;
outputFile5 << endl;
  for ( i = 0; i < n_bond0; i++ ) {
  outputFile5 << b_id0[i][0] << "  " << b_id0[i][1] << "  " << b_id0[i][2] << "  " << b_id0[i][3] << endl;
  }
outputFile5 << endl;

outputFile5 << "Angles" << endl;
outputFile5 << endl;
  for ( i = 0; i < n_angle0; i++ ) {
  outputFile5 << a_id0[i][0] << "  " << a_id0[i][1] << "  " << a_id0[i][2] << "  " << a_id0[i][3] << "  " << a_id0[i][4] << endl;
  }
outputFile5 << endl;

outputFile5 << "Dihedrals" << endl;
outputFile5 << endl;
  for ( i = 0; i < n_dihedral0; i++ ) {
  outputFile5 << a_dih0[i][0] << "  " << a_dih0[i][1] << "  " << a_dih0[i][2] << "  " << a_dih0[i][3] << "  " << a_dih0[i][4] << "  " << a_dih0[i][5] << endl;
  }
outputFile5 << endl;


outputFile5 << "Impropers" << endl;
outputFile5 << endl;
  for ( i = 0; i < n_improper0; i++ ) {
  outputFile5 << a_im0[i][0] << "  " << a_im0[i][1] << "  " << a_im0[i][2] << "  " << a_im0[i][3] << "  " << a_im0[i][4] << "  " << a_im0[i][5] << endl;
  }
outputFile5 << endl;


  return 0;
}


