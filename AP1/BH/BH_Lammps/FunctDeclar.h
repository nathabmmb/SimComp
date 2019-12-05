/**************************** Functions Declaration ************************************/
// Inline functions
inline double distance(double x1,double y1,double z1,double x2,double y2,double z2)
{return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));}
inline double distance_sq(double x1,double y1,double z1,double x2,double y2,double z2)
{return ((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));}

// Functions related to dft calculations - 'scf/relax', 'phonon', 'neb'
double opt_lm(string const,string const,int const,string [],double [][3]);

// Other functions
double round(double);
string calc_zeros(unsigned int);
void write_to_file(string, int const,string [],double [][3],double) ;
string convertInt(int number);
