/************************************ Libraries *****************************************/
#include <iostream>
using namespace std;
using std::cout;
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iomanip>

#include "./FunctDeclar.h"
#include "./Comparison.h"

/********************** Global variables for convenience ******************************/

string temp; int st, st2, st3; double t1, t2;

/***************************************************************************************/
/********************************** Main Program ***************************************/
/***************************************************************************************/

int main(int argc, char* argv[])
{
   // Checking the inline parameters given by the user
   if (argc != 2)   {
      cout << "\nPlease, write the command to execute Lammps. This program needs it together with the BH_Inputfile.txt that should be in this folder. So, as an example, to run BH you should type:\n mpirun -np 4 lammps \nAborting...\n";
      exit(1);
   }

   // Reading the inline input
   char *cpath = argv[1]; string const exec_path1 = cpath;

   // Reading the input file and writing some data to the output file
   #include "./Input.h"

   // Creating the required folders if they don't exist
   temp = "mkdir -p BH_Structures"; system(temp.c_str());

   // Checking the database of structures
   #include "Database.h"

   // Reading initial geometry and putting on a matrix - old style, to use old pieces of code
   Cluster Str("initial.xyz");
   int const nat = total_nat;
   string Species[nat];
   double Matrix[nat][3];
   double en_i=Str.get_energy();
   for(st=0;st<nat;st++)   {
      Species[st]=Str.get_species(st);
      for (int i = 0; i < 3; i++)
         Matrix[st][i]=Str.get_coord(st,i);                    
   }

   // Used variables
   double q2[nat][3]; double en_t;

   srand((unsigned)time(0)); // Seed for random number generation is read fom the time
   gout << "Starting from step 0...\n";
   for (int lo = 0; lo < maxdisp; lo++)   {
      gout << "Step " << lo << " -> "; gout.flush();
      bool accept = true;

      // Saving the coordinates
      for (st = 0; st < nat;st++)
         for (st2 = 0; st2 < 3; st2++)
            q2[st][st2]= Matrix[st][st2];

      // Changing the coordinates
      if (lo != 0)   {
         for (st2 = 0; st2 < nat; st2++)   {
            for (int i = 0; i < 3; i++)   {
               t1 = 1.0-0.001*(random() % 2000); 
               Matrix[st2][i] += step*t1;
            }
         }
      }

      // Evaluating the energy (relaxation)
      en_t = opt_lm(LABEL.c_str(),exec_path1,nat,Species,Matrix);
      string filename = "./BH_Structures/"+ LABEL + calc_zeros(lo) + ".xyz";
      write_to_file(filename.c_str(),nat,Species,Matrix,en_t); 

      // BH conditions - using KT value chosen by user
      double random_number = 0.001*(random() % 1000);
      gout << " new energy " << setprecision(12) << en_t;
      if ( random_number < exp((en_i-en_t)*KTi) )    {
         gout << ", structure accepted .\n";
         en_i=en_t;
      }
      else   {
         gout << ", structure not accepted .\n";
         accept=false;
         // Restore the coordinates
         for (st = 0; st < nat;st++)
            for (st2 = 0; st2 < 3; st2++)
               Matrix[st][st2]= q2[st][st2];
      }
      gout.flush();
  
   }
   // Cleaning
   gout << "Maximum number of steps!\n";
   gout << "*********************************************Program finished successfully! \n";

   temp = "rm -rf Database lammps.out log.lammps relaxed.xyz " + LABEL; system(temp.c_str());

   return 0;
}

/************************************ Functions ************************************/

#include "./Lammps.h"
#include "./Tools.h"
