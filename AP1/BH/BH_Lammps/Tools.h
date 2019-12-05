// Function to round real numbers
double round(double x)   {
   if ( (x - floor(x)) < 0.5) return floor(x);
   else                       return (floor(x)+1.);
}

// Transformation of a number to the required format (like 1 -> _0001, 15 -> _0015)
string calc_zeros(unsigned int n)   {
   if (n == 0)   
      return "_0000";
   else   {
      string zeros = "_";
      int nt = 100*n;
      while (nt < 100000)   {
         zeros += "0";
         nt *= 10;
      }
      stringstream num;
      num << n;
      zeros += num.str();
      return zeros;
   }
}

// Function that writes a number of atoms, an energy, and a structure to a file named 'name'
void write_to_file(string name, int const nat, string Species[], double Matrix[][3], double en)   {
   ofstream output(name.c_str(),ios::out);   
   output << nat << "\n  " << setprecision(10) << en << "  \n";
   for (int i = 0; i < nat;i++)
      output << " " << Species[i] << " " << Matrix[i][0] << " " << Matrix[i][1] << " " << Matrix[i][2] <<"\n";
   output.close();
}
string convertInt(int number)   {
   stringstream ss;
   ss << number;
   return ss.str();
}
