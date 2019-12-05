/****************************************************************/
/* PW funtion to optimize the structure and evaluate the energy */
/****************************************************************/

double opt_lm(string const temporary_file, string const path, int const nat, string Species[], double Matrix[][3])   {
   // Building the Structure.dat file
   temp = "cp Model.dat Structure.dat"; system(temp.c_str());                
   ofstream str("Structure.dat",ios::app);
   for ( int i = 0; i < nat; i++)
      str << (i+1) << "   1   0   " << Matrix[i][0] << " " << Matrix[i][1] << " " << Matrix[i][2] << "\n";
   str.close();
   
   // Running the input
   temp = path+" -in lammps.in > lammps.out"; system(temp.c_str());    

   // Reading the output - energy
   double energy; 
   temp = "grep -A 1 'Energy initial, next-to-last' lammps.out | tail -n1 | awk '{print $3}' > " + temporary_file;
   system(temp.c_str());
   ifstream bar(temporary_file.c_str(),ios::in);  
   bar >> energy;
   bar.close();

   // Reading the output - geometry
   ifstream barr("relaxed.xyz",ios::in);  
   while(barr >> temp)   {
      getline(barr,temp); getline(barr,temp);
      for (int i = 0; i < nat; i++)   {
         barr >> temp;
         for (int j = 0; j < 3; j++)
            barr >> Matrix[i][j];
      }
   }

   return energy;
}
