#include "nr.h"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
using namespace std;


const int aC = 0;
const int aSi = 1;
const int aO = 2;
const int aH = 3;
const int aZn = 4;
const int aMg = 5;
const int aN = 6;
const int aTi = 7;
const int aI = 8;
const int aPb = 9;
const int aLa = 10;
const int aNa = 11;
const int aP = 12;
const int aF = 13;

const int numC = 3;

int namax;
int na;

int main() {

        std::sprintf;
	        std::ofstream outfile;
		std::ofstream outfile2;
        char command0[4500]; 
  char command1[4500]; 
  char command2[4500]; 
  char command3[4500]; 
  char command4[4500]; 
  char command5[4500];
    namax = 500000;
        string data0;
        string data1;
	string string1, string2, system1;
        DP N_node, rr;
        DP v, sin_a, sin_b, sin_c, cos_a, cos_b, cos_c;
    DP avg_x0, avg_x1, avg_x, avg_y0, avg_y1, avg_y, avg_z0, avg_z1, avg_z;
    DP sumx, sumy, sumz, del_x0, del_x1, del_y0, del_y1, del_z0, del_z1;
    DP len0, len1, dd, dd2;
    DP avg_n_x, avg_n_y, avg_n_z;
    DP del_xx, del_yy, del_zz;
    DP dis;
    Mat_IO_DP x(namax, numC);
    Vec_IO_INT elem(namax), b_id0(namax), b_id1(namax), ad0(namax), ad1(namax), h0(namax), h1(namax);
    
    string C("C");
    string Si("Si");
    string O("O");
    string H("H");
    string Zn("Zn");
    string Mg("Mg");
    string N("N");
    string Ti("Ti");
    string I("I");
    string Pb("Pb");
    string La("La");
    string Na("Na");
    string F("F");
    string P("P");
    string str1 = ".";
    string str2 = "+";
    string str3 = "-";
    string str4 = "C1";
    char* ptr;
    
     int numC = 3;
     int na, n_N, ii, jj, iii, jjj, nb, nnn, kkk, ik, naa, ikik;
    int logic1;
    DP penalty, penalty2, penalty3, ddddd;
    DP del_xxx, del_yyy, del_zzz, EE, HH, Pi, HE;
    
    elem = -1;
     DP pi;
      Vec_IO_DP len_c(numC), hlen(numC), cen(numC), angle(numC), comm(namax);
    pi = 3.141592;
    //    int system;
   //     int m,mm,mmm,mmm,mmmm;
   
    // Empty vector holding all names from file
         vector<string> names;
         
    //         // Read names from file LineUp.txt
                 ifstream in("name.txt");
                     if(!in.is_open())
                             cout << "Unable to open file\n";
    
    //                             // this is wrong, by the way: while(in.good()){
                                     string word;
                                        while(getline(in, word))
                                                     names.push_back(word);
    //
     //                                                    sort(names.begin(), names.end());
    //
    //                                                         // Loop to print names

    naa = 0;
   

    

    for (size_t i = 0; i < names.size(); i++)  {
        naa++;
    }
    
    
  //  cout << naa << endl;
    
  //  cout << com[5] << endl;
    
    ikik = 0;
    
    
    sprintf (command1, "mv E_h_bond.txt E_h_bond_back.txt");
    int systemRet1020 =  system(command1);
    
   //     sprintf (command1, "rm output_all.txt");
  //  int systemRet1025 =  system(command1);


    
    outfile.open("E_h_bond.txt", std::ios_base::app);
    
  for (size_t i = 0; i < names.size(); i++) {
      
      ofstream outputFile10077("run", ios::out);

      outputFile10077 << "#!/bin/bash" << endl;
      outputFile10077 << "rm -rf lammps/*" << endl;
      outputFile10077 << "mkdir lammps/" << names[i].c_str() << endl;
      outputFile10077 << "cp back_S/* lammps/" << names[i].c_str() << endl;
      
      outputFile10077 << "if ! command -v conda &> /dev/null; then" << endl;
      outputFile10077 << "    echo \"❌ Conda not found. Please install Anaconda or Miniconda.\"" << endl;
      outputFile10077 << "    exit 1" << endl;
      outputFile10077 << "fi" << endl;
      
      outputFile10077 << "source \"$(conda info --base)/etc/profile.d/conda.sh\"" << endl;
      outputFile10077 << "conda activate AmberTools23 || {" << endl;
      outputFile10077 << "    echo \"❌ AmberTools23 environment not found.\"" << endl;
      outputFile10077 << "    exit 1" << endl;
      outputFile10077 << "}" << endl;

      outputFile10077 << "cp structures/" << names[i].c_str() << ".smi test.smi" << endl;
      outputFile10077 << "taskset --cpu-list 0 obabel -ismi test.smi -omol2 --gen3d --partialcharge eem -O monomer_com1.mol2" << endl;
      outputFile10077 << "antechamber -i monomer_com1.mol2 -fi mol2 -fo mol2 -o monomer_com2.mol2 -pf y -at gaff" << endl;
      outputFile10077 << "../Util/Util_monomer_modify_4_new_make/exe" << endl;
      outputFile10077 << "cp monomer_reorder2.mol2 mol2tolt/test/monomer.mol2" << endl;

      outputFile10077 << "cp polymer_new.lt moltemplates/polymer.lt" << endl;
      outputFile10077 << "cp system_new.lt moltemplates/system.lt" << endl;

      outputFile10077 << "cd mol2tolt/" << endl;
      outputFile10077 << "./run.sh" << endl;

      outputFile10077 << "cp test/monomer.lt ../moltemplates/monomer_add.lt" << endl;
      outputFile10077 << "cd .." << endl;

      outputFile10077 << "cd moltemplates" << endl;
      outputFile10077 << "../../Util/moltemplate-master/moltemplate/scripts/moltemplate.sh system.lt > moltemplate.log 2>&1" << endl;

      outputFile10077 << "../../Util/Util_data_mol_modify/exe" << endl;


      outputFile10077 << "cp system2.data ../lammps/" << names[i].c_str() << "/system.data" << endl;
      outputFile10077 << "cd .." << endl;


      outputFile10077 << "cd lammps/" << names[i].c_str() << endl;
      outputFile10077 << "../../../Util/lammps-2Aug2023/src/lmp_serial -i run.in.npt2" << endl;
      outputFile10077 << "../../../Util/lammps-2Aug2023/src/lmp_serial -i run.in.npt2_pppm" << endl;

      outputFile10077 << "../../../Util/Util_Polymer_Output2/exe" << endl;

      sprintf (command1, "cp run run_exe");
      int systemRet1020 =  system(command1);


      sprintf (command1, "chmod 770 run_exe");
      int systemRet1030 =  system(command1);

      sprintf (command1, "./run_exe");
      int systemRet1010 =  system(command1);

      
                 sprintf (command1, "cp output_all_back.txt output_all.txt");
    int systemRet1090 =  system(command1);


      sprintf (command1, "cp lammps/%s/output_all.txt ./", names[i].c_str());
      int systemRet2010 =  system(command1);


      ifstream inputFile309("output_all.txt",ios::in);
      if (!inputFile309) {
        cerr << "inputFile309 could not be opened" << endl;
        exit(1);
      }
      
      inputFile309 >> EE >> HH >> Pi >> HE;
      
      outfile << names[i].c_str() <<  "    "  << EE << "    " << HH << "    " << Pi << "    " << HE << endl;


            sprintf (command1, "rm output_all.txt");
    int systemRet1050 =  system(command1);



      
}

 

  return 0;


                                                         
}
                                       
