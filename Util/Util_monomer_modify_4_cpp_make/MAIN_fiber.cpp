/************************************************************
 * Developers: Sangheon Lee and Gyeong S. Hwang
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
    Vec_IO_INT ad0_del(namax), ad0_cl(namax), b_del(namax), ad0_cl_del(namax), b_cl_del(namax);
    Vec_IO_INT b_id0_2(namax), b_id1_2(namax);
    Vec_IO_INT b_id0_3(namax), b_id1_3(namax), ad0_cl_del_2(namax), b_cl_del_2(namax), ad0_del_2(namax);
  Mat_IO_DP len(numC, numC), hlen(numC, numC);
  Vec_IO_DP v_d_ij(numC), q(namax);
  int i, j, k, l, ii, jj, kk, jjj, nb, kkk, iii, nnn, jjjj;
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
    DP deldel_x, deldel_y, deldel_z, x_back;
    DP sum_q, x_min, interval, d0, x_max;
    DP sum_xx, avg_xx;

  string C("C");
  string CL("CL");
  string O("O");
  string H("H");
  string N("N");
  string string1;
  string system1;
    string system2;

  elem = -1;
  
    nE1 = 0;
  // 
  // Read an initial configuration.
  //
  
    
    ifstream inputFile3("monomer_com1.mol2",ios::in);
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
          logic1 = string1.compare(N);
          if (logic1 == 0) {
              elem[ii] = aN;
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
              nE1++;
          }

             else {
           logic1 = string1.compare(H);
           if (logic1 == 0) {
             elem[ii] = aH;
               nE1++;
           }

          else {
        logic1 = string1.compare(C);
        if (logic1 == 0) {
          elem[ii] = aC;
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
          
 //   inputFile3 >> string1;
//  inputFile3 >> string1 >> string1 >> string1 >> string1 >> string1 >> string1 >> string1 >> string1 >> string1;
    
    kkk = 0;
          for (iii = 0; iii < na; iii++) {
            if (elem[iii] == aN) {
                nnn = 0;
                for ( jjj = 0; jjj < nb; jjj++ ) {
                    if ( b_id0[jjj] == iii + 1  && elem[b_id1[jjj]-1] == aH ) {
                        nnn++;
                    }
                }
                if ( nnn == 2) {
                    ad0[kkk] = iii;
                 //   ad0_del[kkk] = b_id1[jjj]-1;
                //    b_del[kkk] = jjj;
                    kkk++;
                }
            }
            }
    
    
    kkk = 0;
   

                for ( jjj = 0; jjj < nb; jjj++ ) {
                    if ( b_id0[jjj] == ad0[0] + 1  && elem[b_id1[jjj]-1] == aH ) {
                       ad0_del[kkk] = b_id1[jjj]-1;
                           b_del[kkk] = jjj;
                }
            }
    
    // cout << ad0[0] << endl;
    // cout << ad0_del[0] << endl;
    
    kkk = 0;
    for (iii = 0; iii < na; iii++) {
        if (elem[iii] != aC) continue;
     //     nnn = 0;
     //
     
        for ( jjjj = 0; jjjj < nb; jjjj++ ) {
            if ( (b_id0[jjjj] == iii + 1 && elem[b_id1[jjjj]-1] == aO  ) || (b_id1[jjjj] == iii + 1 && elem[b_id0[jjjj]-1] == aO )  ) {

          for ( jjj = 0; jjj < nb; jjj++ ) {
              if ( (b_id0[jjj] == iii + 1 && elem[b_id1[jjj]-1] == aCL  ) || (b_id1[jjj] == iii + 1 && elem[b_id0[jjj]-1] == aCL )  ) {
  //                nnn++;
                  ad0_cl[kkk] = iii;
                  ad0_cl_del[kkk] = b_id1[jjj]-1;
                  b_cl_del[kkk] = jjj;
                  kkk++;
              }
          }
      
          
      }
        }
    }
    
    // cout << ad0_cl[0] << endl;
    // cout << ad0_cl_del[0] << endl;
    // cout << b_cl_del[0] << endl;
    
    
    sum_q = 0;
    
    for (j = 0; j < na; j++) {
        if (j == ad0_del[0] || j == ad0_cl_del[0]) continue;
        sum_q += q[j];
    }
    
    // Create a new variable to count the number of hydrogen atoms.
    int hydrogen_count = 0;
    
    // Calculate the number of hydrogen atoms (excluding hydrogen that can be removed)
    for (j = 0; j < na; j++) {
    if (elem[j] == aH && j != ad0_del[0]) {
        hydrogen_count++;
    }
    }
    
    
    // cout << sum_q << endl;
    // cout << sum_q/hydrogen_count << endl;
    
    for (j = 0; j < na; j++) {
        if (elem[j] != aH) continue;
        q[j] = q[j] - sum_q/hydrogen_count;
    }
        
    
    deldel_x = abs(x[ad0[0]][0] - x[ad0_cl[0]][0]);
    deldel_y = abs(x[ad0[0]][1] - x[ad0_cl[0]][1]);
    deldel_z = abs(x[ad0[0]][2] - x[ad0_cl[0]][2]);
    
    
    if (deldel_y > deldel_x && deldel_y > deldel_z) {
        for (ii = 0; ii < na; ii++) {
            x_back = x[ii][0];
            x[ii][0] = x[ii][1];
            x[ii][1] = x_back;
        }
}
    
    if (deldel_z > deldel_x && deldel_z > deldel_y) {
        for (ii = 0; ii < na; ii++) {
            x_back = x[ii][0];
            x[ii][0] = x[ii][2];
            x[ii][2] = x_back;
        }
}
 
    
    sum_xx = 0;
    
    for ( ii = 0; ii < na; ii++) {
        sum_xx = sum_xx + x[ii][0];
    }
    
    avg_xx = sum_xx/na;
    
    if ((x[ad0[0]][0] - x[ad0_cl[0]][0]) > 0 ) {
        
        for ( ii = 0; ii < na; ii++) {
            x[ii][0] = 2*avg_xx - x[ii][0];
        }
        
    }
    
    
    
    x_min = 100000;
    for ( ii = 0; ii < na; ii++ ) {
        if ( x[ii][0] < x_min ) {
            x_min = x[ii][0];
        }
    }
    
    for ( ii = 0; ii < na; ii++ ) {
           x[ii][0] = abs(x_min) + x[ii][0];
    }
    
    
    x_max = -100000;
    for ( ii = 0; ii < na; ii++ ) {
        if ( x[ii][0] > x_max ) {
            x_max = x[ii][0];
        }
    }
    
    interval = x_max - x_min +1;
    
  ofstream outputFile4("monomer_reorder2.mol2",ios::out);
  if (!outputFile4) {
    cerr << "outputFile4 could not be opened" << endl;
    exit(1);
  }
  outputFile4 << "@<TRIPOS>MOLECULE" << endl;
  outputFile4 << "*****" << endl;
    outputFile4 << " " << na-2 << " " << nb-2 << " " << "0 0 0" << endl;
    outputFile4 << "SMALL" << endl;
    outputFile4 << "GASTEIGER" << endl;
    outputFile4 << endl;
    outputFile4 << "@<TRIPOS>ATOM" << endl;
    
    ifstream inputFile33("monomer_com2.mol2",ios::in);
    if (!inputFile33) {
      cerr << "inputFile33 could not be opened" << endl;
      exit(1);
    }
    inputFile33 >> string1;
        inputFile33 >> string1;
        inputFile33 >> na >> nb >> string1 >> string1 >> string1;
    inputFile33 >> string1;
    inputFile33 >> string1 >> string1 >> string1 >> string1 >> string1;
        
        inputFile33 >> string1;
  
    for (i = 0; i < na; i++) {
      inputFile33 >> string1 >> string1 >> dd >> dd >> dd >> system1 >> system2 >> system2 >> dd;
        if(ad0[0] != i ) continue;
        outputFile4 << "      " <<  "1" << " " << "N2000" << "           " << x[ad0[0]][0] << "   " << x[ad0[0]][1] << "    " << x[ad0[0]][2] << " " << "n" << "  " << "1" << "  " << "UNL1" << "    " << q[ad0[0]] << endl;
    }
    
    
    ifstream inputFile303("monomer_com2.mol2",ios::in);
    if (!inputFile303) {
      cerr << "inputFile303 could not be opened" << endl;
      exit(1);
    }
    inputFile303 >> string1;
        inputFile303 >> string1;
        inputFile303 >> na >> nb >> string1 >> string1 >> string1;
    inputFile303 >> string1;
    inputFile303 >> string1 >> string1 >> string1 >> string1 >> string1;
        
        inputFile303 >> string1;
  
    for (i = 0; i < na; i++) {
      inputFile303 >> string1 >> string1 >> dd >> dd >> dd >> system1 >> system2 >> system2 >> dd;
        if(ad0_cl[0] != i ) continue;
        outputFile4 << "      " <<  "2" << " " << "C2000" << "           " << x[ad0_cl[0]][0] << "   " << x[ad0_cl[0]][1] << "    " << x[ad0_cl[0]][2] << " " << system1 << "  " << "1" << "  " << "UNL1" << "    " << q[ad0_cl[0]] << endl;
    }
    
    
    ifstream inputFile331("monomer_com2.mol2",ios::in);
    if (!inputFile331) {
      cerr << "inputFile331 could not be opened" << endl;
      exit(1);
    }
    inputFile331 >> string1;
        inputFile331 >> string1;
        inputFile331 >> na >> nb >> string1 >> string1 >> string1;
    inputFile331 >> string1;
    inputFile331 >> string1 >> string1 >> string1 >> string1 >> string1;
        
        inputFile331 >> string1;
    
    kkk = 3;
    for (i = 0; i < na; i++) {
        inputFile331 >> string1 >> string1 >> dd >> dd >> dd >> system1 >> system2 >> system2 >> dd;
        if (ad0[0] == i || ad0_cl[0] == i || ad0_del[0] == i || ad0_cl_del[0] == i) continue;
        outputFile4 << "      " << kkk << " " << string1 << "           " << x[i][0] << "   " << x[i][1] << "    " << x[i][2] << " " << system1 << "  " << "1" << "  " << "UNL1" << "    " << q[i] << endl;
        kkk++;
   }

    outputFile4 << "@<TRIPOS>BOND" << endl;
  
    b_id0_2 = -1;
    b_id1_2 = -1;
    
    
    for (j = 0; j < nb; j++) {
        
        b_id0_up = b_id0[j];
        b_id1_up = b_id1[j];
        
        if(b_id0[j] < ad0[0]+1 && b_id0[j] < ad0_cl[0]+1 ) {
            b_id0_up = b_id0[j] +2;
        }
        if(b_id1[j] < ad0[0]+1 && b_id1[j] < ad0_cl[0]+1) {
            b_id1_up = b_id1[j] +2;
        }
        if( (b_id0[j] < ad0[0]+1 && b_id0[j] > ad0_cl[0]+1) || (b_id0[j] > ad0[0]+1 && b_id0[j] < ad0_cl[0]+1) ) {
            b_id0_up = b_id0[j] +1;
        }
        if( (b_id1[j] < ad0[0]+1 && b_id1[j] > ad0_cl[0]+1) || (b_id1[j] > ad0[0]+1 && b_id1[j] < ad0_cl[0]+1)) {
            b_id1_up = b_id1[j] +1;
        }
            
        if(b_id0[j] == ad0[0]+1) {
            b_id0_up = 1;
        }
        if(b_id1[j] == ad0[0]+1) {
            b_id1_up = 1;
        }
        if(b_id0[j] == ad0_cl[0]+1) {
            b_id0_up = 2;
        }
        if(b_id1[j] == ad0_cl[0]+1) {
            b_id1_up = 2;
        }

        b_id0_2[j] = b_id0_up;
        b_id1_2[j] = b_id1_up;
    }
    
    
      b_id0_3 = -1;
      b_id1_3 = -1;
    
    ad0_del_2 = -1;
    ad0_cl_del_2 = -1;
      
    if( ad0_del[0] < ad0[0] && ad0_del[0] < ad0_cl[0]) {
        ad0_del_2[0] = ad0_del[0] + 2;
    }
    
    if( (ad0_del[0] < ad0[0] && ad0_del[0] > ad0_cl[0]) || (ad0_del[0] > ad0[0] && ad0_del[0] < ad0_cl[0])) {
        ad0_del_2[0] = ad0_del[0] + 1;
    }
    
    if( ad0_del[0] > ad0[0] && ad0_del[0] > ad0_cl[0]) {
        ad0_del_2[0] = ad0_del[0];
    }
    
    
    if( ad0_cl_del[0] < ad0[0] && ad0_cl_del[0] < ad0_cl[0]) {
        ad0_cl_del_2[0] = ad0_cl_del[0] + 2;
    }
    
    if( (ad0_cl_del[0] < ad0[0] && ad0_cl_del[0] > ad0_cl[0]) || (ad0_cl_del[0] > ad0[0] && ad0_cl_del[0] < ad0_cl[0])) {
        ad0_cl_del_2[0] = ad0_cl_del[0] + 1;
    }
    
    if( ad0_cl_del[0] > ad0[0] && ad0_cl_del[0] > ad0_cl[0]) {
        ad0_cl_del_2[0] = ad0_cl_del[0];
    }
    
    
    for (j = 0; j < nb; j++) {
        
        b_id0_up = b_id0_2[j];
        b_id1_up = b_id1_2[j];
        
        if(b_id0_2[j] < ad0_del_2[0]+1 && b_id0_2[j] < ad0_cl_del_2[0]+1 ) {
            b_id0_up = b_id0_2[j];
        }
        if(b_id1_2[j] < ad0_del_2[0]+1 && b_id1_2[j] < ad0_cl_del_2[0]+1) {
            b_id1_up = b_id1_2[j];
        }
        if( (b_id0_2[j] < ad0_del_2[0]+1 && b_id0_2[j] > ad0_cl_del_2[0]+1) || (b_id0_2[j] > ad0_del_2[0]+1 && b_id0_2[j] < ad0_cl_del_2[0]+1) ) {
            b_id0_up = b_id0_2[j] -1;
        }
        if( (b_id1_2[j] < ad0_del_2[0]+1 && b_id1_2[j] > ad0_cl_del_2[0]+1) || (b_id1_2[j] > ad0_del_2[0]+1 && b_id1_2[j] < ad0_cl_del_2[0]+1)) {
            b_id1_up = b_id1_2[j] -1;
        }
    
        if( (b_id0_2[j] > ad0_del_2[0]+1 && b_id0_2[j] > ad0_cl_del_2[0]+1) ) {
            b_id0_up = b_id0_2[j] -2;
        }
        if( (b_id1_2[j] > ad0_del_2[0]+1 && b_id1_2[j] > ad0_cl_del_2[0]+1) ) {
            b_id1_up = b_id1_2[j] -2;
        }

        b_id0_3[j] = b_id0_up;
        b_id1_3[j] = b_id1_up;
    }
    


    ifstream inputFile333("monomer_com2.mol2",ios::in);
    if (!inputFile333) {
      cerr << "inputFile333 could not be opened" << endl;
      exit(1);
    }
    inputFile333 >> string1;
        inputFile333 >> string1;
        inputFile333 >> na >> nb >> string1 >> string1 >> string1;
    inputFile333 >> string1;
    inputFile333 >> string1 >> string1 >> string1 >> string1 >> string1;
        
        inputFile333 >> string1;

        for (ii = 0; ii < na; ii++) {
            inputFile333 >> string1 >> string1 >> dd >> dd >> dd >> system1 >> system2 >> system2 >> dd;
        }
    inputFile333 >> string1;
    
    jjjj = 0;
    for (jj = 0; jj < nb; jj++) {
        if (b_del[0] == jj || b_cl_del[0] == jj) continue;
        inputFile333 >> rr >> dd >> dd >> string1;
        outputFile4 << "     " << jjjj+1 << "    " << b_id0_3[jj] << "    " << b_id1_3[jj] << "    " << string1 << endl;
        jjjj++;
   }
    
    outputFile4 << "@<TRIPOS>SUBSTRUCTURE" << endl;
    outputFile4 << "      1 ***         1 TEMP              0 ****  ****    0 ROOT" << endl;
            

    ofstream outputFile40("polymer_new.lt",ios::out);
    if (!outputFile40) {
      cerr << "outputFile40 could not be opened" << endl;
      exit(1);
    }
    outputFile40 << "import monomer_add.lt" << endl;
outputFile40 << "import PPTA.lt" << endl;
    outputFile40 << endl;
    outputFile40 << "aramid inherits GAFF {" << endl;
    outputFile40 << endl;
    outputFile40 << "   create_var {$mol}" << endl;
    outputFile40 << endl;
    
// Polymer chain generation code start
outputFile40 << "monomers[0] = new cation.move(0,0,0)" << endl;
d0 = interval;
outputFile40 << "monomers[1] = new PPTA.move(" << d0 << ",0,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[2] = new cation.move(" << d0 << ",0,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[3] = new cation.move(" << d0 << ",0,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[4] = new PPTA.move(" << d0 << ",0,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[5] = new PPTA.move(" << d0 << ",0,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[6] = new PPTA.move(" << d0 << ",0,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[7] = new cation.move(" << d0 << ",0,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[8] = new cation.move(" << d0 << ",0,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[9] = new PPTA.move(" << d0 << ",0,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[10] = new cation.move(0,40,0)" << endl;
d0 = interval;
outputFile40 << "monomers[11] = new cation.move(" << d0 << ",40,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[12] = new cation.move(" << d0 << ",40,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[13] = new cation.move(" << d0 << ",40,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[14] = new PPTA.move(" << d0 << ",40,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[15] = new PPTA.move(" << d0 << ",40,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[16] = new PPTA.move(" << d0 << ",40,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[17] = new cation.move(" << d0 << ",40,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[18] = new PPTA.move(" << d0 << ",40,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[19] = new PPTA.move(" << d0 << ",40,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[20] = new PPTA.move(0,80,0)" << endl;
d0 = 12;
outputFile40 << "monomers[21] = new PPTA.move(" << d0 << ",80,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[22] = new cation.move(" << d0 << ",80,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[23] = new cation.move(" << d0 << ",80,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[24] = new PPTA.move(" << d0 << ",80,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[25] = new PPTA.move(" << d0 << ",80,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[26] = new cation.move(" << d0 << ",80,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[27] = new cation.move(" << d0 << ",80,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[28] = new cation.move(" << d0 << ",80,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[29] = new PPTA.move(" << d0 << ",80,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[30] = new PPTA.move(0,120,0)" << endl;
d0 = 12;
outputFile40 << "monomers[31] = new cation.move(" << d0 << ",120,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[32] = new PPTA.move(" << d0 << ",120,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[33] = new cation.move(" << d0 << ",120,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[34] = new PPTA.move(" << d0 << ",120,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[35] = new PPTA.move(" << d0 << ",120,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[36] = new cation.move(" << d0 << ",120,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[37] = new cation.move(" << d0 << ",120,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[38] = new cation.move(" << d0 << ",120,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[39] = new PPTA.move(" << d0 << ",120,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[40] = new cation.move(0,160,0)" << endl;
d0 = interval;
outputFile40 << "monomers[41] = new cation.move(" << d0 << ",160,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[42] = new PPTA.move(" << d0 << ",160,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[43] = new cation.move(" << d0 << ",160,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[44] = new PPTA.move(" << d0 << ",160,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[45] = new PPTA.move(" << d0 << ",160,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[46] = new cation.move(" << d0 << ",160,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[47] = new cation.move(" << d0 << ",160,0)" << endl;  // cation
d0 = d0 + interval;
outputFile40 << "monomers[48] = new PPTA.move(" << d0 << ",160,0)" << endl;  // PPTA
d0 = d0 + 12;
outputFile40 << "monomers[49] = new PPTA.move(" << d0 << ",160,0)" << endl;  // PPTA
d0 = d0 + 12;
// Polymer chain generation code end

// Bond structure definition start
outputFile40 << "   write('Data Bond List') {" << endl;
// First chain bonds (0-9)
outputFile40 << "    $bond:b1  $atom:monomers[0]/atom2 $atom:monomers[1]/atom1" << endl;
outputFile40 << "    $bond:b2  $atom:monomers[1]/atom2 $atom:monomers[2]/atom1" << endl;
outputFile40 << "    $bond:b3  $atom:monomers[2]/atom2 $atom:monomers[3]/atom1" << endl;
outputFile40 << "    $bond:b4  $atom:monomers[3]/atom2 $atom:monomers[4]/atom1" << endl;
outputFile40 << "    $bond:b5  $atom:monomers[4]/atom2 $atom:monomers[5]/atom1" << endl;
outputFile40 << "    $bond:b6  $atom:monomers[5]/atom2 $atom:monomers[6]/atom1" << endl;
outputFile40 << "    $bond:b7  $atom:monomers[6]/atom2 $atom:monomers[7]/atom1" << endl;
outputFile40 << "    $bond:b8  $atom:monomers[7]/atom2 $atom:monomers[8]/atom1" << endl;
outputFile40 << "    $bond:b9  $atom:monomers[8]/atom2 $atom:monomers[9]/atom1" << endl;
outputFile40 << "    $bond:b10  $atom:monomers[9]/atom2 $atom:monomers[0]/atom1" << endl;
// Second chain bonds (10-19)
outputFile40 << "    $bond:b11  $atom:monomers[10]/atom2 $atom:monomers[11]/atom1" << endl;
outputFile40 << "    $bond:b12  $atom:monomers[11]/atom2 $atom:monomers[12]/atom1" << endl;
outputFile40 << "    $bond:b13  $atom:monomers[12]/atom2 $atom:monomers[13]/atom1" << endl;
outputFile40 << "    $bond:b14  $atom:monomers[13]/atom2 $atom:monomers[14]/atom1" << endl;
outputFile40 << "    $bond:b15  $atom:monomers[14]/atom2 $atom:monomers[15]/atom1" << endl;
outputFile40 << "    $bond:b16  $atom:monomers[15]/atom2 $atom:monomers[16]/atom1" << endl;
outputFile40 << "    $bond:b17  $atom:monomers[16]/atom2 $atom:monomers[17]/atom1" << endl;
outputFile40 << "    $bond:b18  $atom:monomers[17]/atom2 $atom:monomers[18]/atom1" << endl;
outputFile40 << "    $bond:b19  $atom:monomers[18]/atom2 $atom:monomers[19]/atom1" << endl;
outputFile40 << "    $bond:b20  $atom:monomers[19]/atom2 $atom:monomers[10]/atom1" << endl;
// Third chain bonds (20-29)
outputFile40 << "    $bond:b21  $atom:monomers[20]/atom2 $atom:monomers[21]/atom1" << endl;
outputFile40 << "    $bond:b22  $atom:monomers[21]/atom2 $atom:monomers[22]/atom1" << endl;
outputFile40 << "    $bond:b23  $atom:monomers[22]/atom2 $atom:monomers[23]/atom1" << endl;
outputFile40 << "    $bond:b24  $atom:monomers[23]/atom2 $atom:monomers[24]/atom1" << endl;
outputFile40 << "    $bond:b25  $atom:monomers[24]/atom2 $atom:monomers[25]/atom1" << endl;
outputFile40 << "    $bond:b26  $atom:monomers[25]/atom2 $atom:monomers[26]/atom1" << endl;
outputFile40 << "    $bond:b27  $atom:monomers[26]/atom2 $atom:monomers[27]/atom1" << endl;
outputFile40 << "    $bond:b28  $atom:monomers[27]/atom2 $atom:monomers[28]/atom1" << endl;
outputFile40 << "    $bond:b29  $atom:monomers[28]/atom2 $atom:monomers[29]/atom1" << endl;
outputFile40 << "    $bond:b30  $atom:monomers[29]/atom2 $atom:monomers[20]/atom1" << endl;
// Fourth chain bonds (30-39)
outputFile40 << "    $bond:b31  $atom:monomers[30]/atom2 $atom:monomers[31]/atom1" << endl;
outputFile40 << "    $bond:b32  $atom:monomers[31]/atom2 $atom:monomers[32]/atom1" << endl;
outputFile40 << "    $bond:b33  $atom:monomers[32]/atom2 $atom:monomers[33]/atom1" << endl;
outputFile40 << "    $bond:b34  $atom:monomers[33]/atom2 $atom:monomers[34]/atom1" << endl;
outputFile40 << "    $bond:b35  $atom:monomers[34]/atom2 $atom:monomers[35]/atom1" << endl;
outputFile40 << "    $bond:b36  $atom:monomers[35]/atom2 $atom:monomers[36]/atom1" << endl;
outputFile40 << "    $bond:b37  $atom:monomers[36]/atom2 $atom:monomers[37]/atom1" << endl;
outputFile40 << "    $bond:b38  $atom:monomers[37]/atom2 $atom:monomers[38]/atom1" << endl;
outputFile40 << "    $bond:b39  $atom:monomers[38]/atom2 $atom:monomers[39]/atom1" << endl;
outputFile40 << "    $bond:b40  $atom:monomers[39]/atom2 $atom:monomers[30]/atom1" << endl;
// Fifth chain bonds (40-49)
outputFile40 << "    $bond:b41  $atom:monomers[40]/atom2 $atom:monomers[41]/atom1" << endl;
outputFile40 << "    $bond:b42  $atom:monomers[41]/atom2 $atom:monomers[42]/atom1" << endl;
outputFile40 << "    $bond:b43  $atom:monomers[42]/atom2 $atom:monomers[43]/atom1" << endl;
outputFile40 << "    $bond:b44  $atom:monomers[43]/atom2 $atom:monomers[44]/atom1" << endl;
outputFile40 << "    $bond:b45  $atom:monomers[44]/atom2 $atom:monomers[45]/atom1" << endl;
outputFile40 << "    $bond:b46  $atom:monomers[45]/atom2 $atom:monomers[46]/atom1" << endl;
outputFile40 << "    $bond:b47  $atom:monomers[46]/atom2 $atom:monomers[47]/atom1" << endl;
outputFile40 << "    $bond:b48  $atom:monomers[47]/atom2 $atom:monomers[48]/atom1" << endl;
outputFile40 << "    $bond:b49  $atom:monomers[48]/atom2 $atom:monomers[49]/atom1" << endl;
outputFile40 << "    $bond:b50  $atom:monomers[49]/atom2 $atom:monomers[40]/atom1" << endl;
// Bond structure definition end

    outputFile40 << "}" << endl;
    outputFile40 << "}" << endl;
    
    ofstream outputFile401("system_new.lt",ios::out);
    if (!outputFile401) {
      cerr << "outputFile401 could not be opened" << endl;
      exit(1);
    }
    outputFile401 << "import polymer.lt" << endl;
    outputFile401 << "polymers = new aramid" << endl;
    outputFile401 << "write_once(" << "\"Data Boundary\"" << ") {" << endl;
    outputFile401 << "0.0    " << d0 << "  xlo xhi" << endl;
    outputFile401 << "0.0    200"<< "  ylo yhi" << endl;
    outputFile401 << "0.0    40" << "  zlo zhi" << endl;
    outputFile401 << "}" << endl;
    

  return 0;
}
