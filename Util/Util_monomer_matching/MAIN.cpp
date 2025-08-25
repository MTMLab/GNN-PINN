/************************************************************
 * Developers: Modified for -COOH group reordering
 * Original by: Sangheon Lee and Gyeong S. Hwang
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

int numC = 3;
int numN = 4;
DP pi = 3.141592654;
int aC = 0;
int aCL = 1;
int aO = 2;
int aH = 3;
int aN = 4;

//
// Global variables.
//
int namax; // Maximum number of atoms.
int nemax = 5; // Maximum number of elements.
int na; // Number of total atoms in the system.
int nnmax = 12; // Maximum number of neighbors.
DP cutoff1 = 3.0;
int penalty;
int npmax = 10; 


int main()  
{
  //  cout << "Enter number of atoms:" << endl; 
  //  cin >> namax;
  namax = 20000;
  na = namax;
  Mat_IO_INT nn(namax, nnmax);
  Mat_IO_DP x(namax, numC), r(namax, numN);
  Vec_IO_INT elem(namax), cnt(namax), idb1(namax), idr1(namax), ide1(namax), b_id0(namax), b_id1(namax), ad0(namax);
  Mat_IO_DP len(numC, numC), hlen(numC, numC);
  Vec_IO_DP v_d_ij(numC), q(namax);
  int i, j, k, l, ii, jj, kk, jjj, nb, kkk, iii, nnn;
  int idum;
  int logic1;
  int charge1;
  int nE1, nE2;
  int cnt1;
  int nCC, nCH;
  int nangles, nbonds, ndihedrals;
  int b_id0_up, b_id1_up;
  DP sum0, sum1, sum2, sum3, sum4, sum5, dd;
  DP R0, R1, R2, R3, R4, R5;
  DP n0, n1, n2, n3, n4, n5;
  DP dx, dy, dz, rij, rr;

  string C("C");
  string CL("CL");
  string O("O");
  string H("H");
  string N("N");
  string string1;
  string system1;
  string system2;

  elem = -1;
  
  // 
  // Read an initial configuration.
  //
  
    
  ifstream inputFile3("monomer.mol2",ios::in);
  if (!inputFile3) {
    cerr << "inputFile3 could not be opened" << endl;
    exit(1);
  }
  inputFile3 >> string1;
  inputFile3 >> string1;
  inputFile3 >> na >> nb >> string1 >> string1 >> string1;
  inputFile3 >> string1;
  inputFile3 >> string1;
  inputFile3 >> string1;

  for (ii = 0; ii < na; ii++) {
    inputFile3 >> string1 >> string1 >> x[ii][0] >> x[ii][1] >> x[ii][2] >> system1 >> system1 >> system1 >> q[ii];
    logic1 = string1.compare(C);
    if (logic1 == 0) {
      elem[ii] = aC;
    }
    else {
      logic1 = string1.compare(CL);
      if (logic1 == 0) {
        elem[ii] = aCL;
      }
      else {
        logic1 = string1.compare(O);
        if (logic1 == 0) {
          elem[ii] = aO;
        }
        else {
          logic1 = string1.compare(H);
          if (logic1 == 0) {
            elem[ii] = aH;
          }
          else {
            logic1 = string1.compare(N);
            if (logic1 == 0) {
              elem[ii] = aN;
            }
            else {
              elem[ii] = 10;
            }
          }
        }
      }
    }
  }
    
  inputFile3 >> string1;
    
  for (jj = 0; jj < nb; jj++) {
    inputFile3 >> rr >> b_id0[jj] >> b_id1[jj] >> string1;
  }
          
  // Find carbon atoms that are bonded to both O and OH (-COOH group)
  kkk = 0;
  for (iii = 0; iii < na; iii++) {
    if (elem[iii] == aC) { // If this is a carbon atom
      int hasO = 0; // Flag for oxygen bond (C=O)
      int hasOH = 0; // Flag for OH group
      
      // Check all bonds to see if this carbon is bonded to oxygen atoms
      for (jjj = 0; jjj < nb; jjj++) {
        int oxygenIndex = -1;
        
        // Check if carbon is bonded to oxygen
        if (b_id0[jjj] == iii + 1 && elem[b_id1[jjj]-1] == aO) {
          oxygenIndex = b_id1[jjj]-1;
        }
        if (b_id1[jjj] == iii + 1 && elem[b_id0[jjj]-1] == aO) {
          oxygenIndex = b_id0[jjj]-1;
        }
        
        if (oxygenIndex != -1) {
          // Check if this oxygen is part of OH group (bonded to hydrogen)
          int isOH = 0;
          for (int bondCheck = 0; bondCheck < nb; bondCheck++) {
            if ((b_id0[bondCheck] == oxygenIndex + 1 && elem[b_id1[bondCheck]-1] == aH) ||
                (b_id1[bondCheck] == oxygenIndex + 1 && elem[b_id0[bondCheck]-1] == aH)) {
              isOH = 1;
              break;
            }
          }
          
          if (isOH == 1) {
            hasOH = 1; // This carbon is bonded to OH group
          } else {
            hasO = 1; // This carbon is bonded to oxygen (likely C=O)
          }
        }
      }
      
      // If this carbon is bonded to both C=O and C-OH, it's part of -COOH group
      if (hasO == 1 && hasOH == 1) {
        ad0[kkk] = iii;
        kkk++;
      }
    }
  }
  
  // cout << "Found " << kkk << " -COOH groups" << endl;
  if (kkk == 0) {
    cerr << "No -COOH groups found in the molecule!" << endl;
    exit(1);
  }
 
  ofstream outputFile4("monomer_reorder.mol2",ios::out);
  if (!outputFile4) {
    cerr << "outputFile4 could not be opened" << endl;
    exit(1);
  }
  outputFile4 << "@<TRIPOS>MOLECULE" << endl;
  outputFile4 << "*****" << endl;
  outputFile4 << na << " " << nb << " " << "0 0 0" << endl;
  outputFile4 << "SMALL" << endl;
  outputFile4 << "GASTEIGER" << endl;
  outputFile4 << endl;
  outputFile4 << "@<TRIPOS>ATOM" << endl;
    
  ifstream inputFile33("monomer.mol2",ios::in);
  if (!inputFile33) {
    cerr << "inputFile33 could not be opened" << endl;
    exit(1);
  }
  inputFile33 >> string1;
  inputFile33 >> string1;
  inputFile33 >> na >> nb >> string1 >> string1 >> string1;
  inputFile33 >> string1;
  inputFile33 >> string1;
  inputFile33 >> string1;
  
  // Write the first -COOH carbon atom first
  for (i = 0; i < na; i++) {
    inputFile33 >> string1 >> string1 >> dd >> dd >> dd >> system1 >> system2 >> system2 >> dd;
    if(ad0[0] != i ) continue;
    outputFile4 << "1" << " " << string1 << "    " << x[ad0[0]][0] << "   " << x[ad0[0]][1] << "    " << x[ad0[0]][2] << " " << system1 << "  " << "1" << "  " << "UNL1" << "    " << q[ad0[0]] << endl;
  }
    
  ifstream inputFile331("monomer.mol2",ios::in);
  if (!inputFile331) {
    cerr << "inputFile331 could not be opened" << endl;
    exit(1);
  }
  inputFile331 >> string1;
  inputFile331 >> string1;
  inputFile331 >> na >> nb >> string1 >> string1 >> string1;
  inputFile331 >> string1;
  inputFile331 >> string1;
  inputFile331 >> string1;
    
  // Write atoms before the first -COOH carbon
  for (i = 0; i < ad0[0]; i++) {
    inputFile331 >> string1 >> string1 >> dd >> dd >> dd >> system1 >> system2 >> system2 >> dd;
    outputFile4 << i+2 << " " << string1 << "    " << x[i][0] << "   " << x[i][1] << "    " << x[i][2] << " " << system1 << "  " << "1" << "  " << "UNL1" << "    " << q[i] << endl;
  }
    
  ifstream inputFile332("monomer.mol2",ios::in);
  if (!inputFile332) {
    cerr << "inputFile332 could not be opened" << endl;
    exit(1);
  }
  inputFile332 >> string1;
  inputFile332 >> string1;
  inputFile332 >> na >> nb >> string1 >> string1 >> string1;
  inputFile332 >> string1;
  inputFile332 >> string1;
  inputFile332 >> string1;
    
  // Write atoms after the first -COOH carbon
  for (i = 0; i < na; i++) {
    inputFile332 >> string1 >> string1 >> dd >> dd >> dd >> system1 >> system2 >> system2 >> dd;
    if(i <= ad0[0]) continue;
    outputFile4 << i+1 << " " << string1 << "    " << x[i][0] << "   " << x[i][1] << "    " << x[i][2] << " " << system1 << "  " << "1" << "  " << "UNL1" << "    " << q[i] << endl;
  }

  outputFile4 << "@<TRIPOS>BOND" << endl;
  
  // Update bond indices according to the new atom ordering
  for (j = 0; j < nb; j++) {
    if(b_id0[j] < ad0[0]+1) {
      b_id0_up = b_id0[j] +1;
    }
    if(b_id1[j] < ad0[0]+1) {
      b_id1_up = b_id1[j] +1;
    }
    if(b_id0[j] == ad0[0]+1) {
      b_id0_up = 1;
    }
    if(b_id1[j] == ad0[0]+1) {
      b_id1_up = 1;
    }
    if(b_id0[j] > ad0[0]+1) {
      b_id0_up = b_id0[j];
    }
    if(b_id1[j] > ad0[0]+1) {
      b_id1_up = b_id1[j];
    }

    b_id0[j] = b_id0_up;
    b_id1[j] = b_id1_up;
  }

  ifstream inputFile333("monomer.mol2",ios::in);
  if (!inputFile333) {
    cerr << "inputFile333 could not be opened" << endl;
    exit(1);
  }
  inputFile333 >> string1;
  inputFile333 >> string1;
  inputFile333 >> na >> nb >> string1 >> string1 >> string1;
  inputFile333 >> string1;
  inputFile333 >> string1;
  inputFile333 >> string1;

  for (ii = 0; ii < na; ii++) {
    inputFile333 >> string1 >> string1 >> dd >> dd >> dd >> system1 >> system2 >> system2 >> dd;
  }
  inputFile333 >> string1;
    
  // Write bond information with updated indices
  for (jj = 0; jj < nb; jj++) {
    inputFile333 >> rr >> dd >> dd >> string1;
    outputFile4 << jj+1 << "    " << b_id0[jj] << "   " << b_id1[jj] << "   " << string1 << endl;
  }
    
  //cout << "Successfully reordered molecule with -COOH group at position 1" << endl;
  //cout << "Output written to monomer_reorder.mol2" << endl;

  return 0;
}
