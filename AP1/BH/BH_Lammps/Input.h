/*** Reading the input file ***/
   ifstream input("BH_Inputfile.txt",ios::in);
   if (!input)   {
      cout << "\nCan't find the BH_Inputfile.txt. Aborting.\n";
      exit(1);
   }

/*** Reading each parameter ***/

   // LABEL
   getline(input,temp);
   input >> temp;
   string const LABEL = temp;

   // Step-size and number of steps
   getline(input,temp);
   input >> t1;
   if (t1 < 0.1 || t1 > 100.0)   {
      cout << "**Warning** - Step size is out of range. Setting it as 1.0.\n"; t1 = 1.0;
   }
   double const step = t1;
   input >> st;
   if (st < 1 || st > 1000000)   {
      cout << "**Warning** - Number of steps is out of range. Setting it as 1000.\n"; st = 1000;
   }
   int const maxdisp = st;

   // KT                
   getline(input,temp);
   input >> t1;
   if (t1 < 0.0 || t1 > 100.0)   {
      cout << "**Warning** - KT is out of range. Setting it as 0.035Ry(0.5eV).\n"; t1 = 0.035;
   }
   double const KT = t1;  // KT
   double const KTi = 1.0/KT;  // KT-1

   // Total number of atoms and number of frozen atoms       
   getline(input,temp);
   input >> st;
   int const total_nat = st;

   input.close();

/*** Returning to user what was read from input ***/
   ofstream gout("BH_output.txt",ios::app);
   gout << "\nProgram has started. Reading the input.\n\n";
   gout << "The label used for temporary files is " << LABEL << ".\n";
   gout << "The step size in each x-y-z direction and number of steps are: " << step << " and " << maxdisp << ".\n";
   gout << "The KT value for the Metropolis criterion is " << KT << ".\n"; 
   gout << "Total number of atoms is: " << total_nat << ".\n"; 
   gout << "\nEnd of input reading.\n"; gout.flush();
