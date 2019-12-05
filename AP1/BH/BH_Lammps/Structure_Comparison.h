   Cluster NEW(filename);
   NEW.center(); NEW.orient();
   bool is_new = true;
   string reason; Cluster Str3;
   for (int i = 0; i < count_db;i++)   {
      // First compare the energies, if they are similar, proceed
      if ( fabs(NEW.get_energy() - DB_Str[i].get_energy()) < energy_comparison )   {
         // Orient
         DB_Str[i].center(); DB_Str[i].orient();
         Cluster Sorted_NEW = sorting_comp(NEW,DB_Str[i]);

         // Evalute the smallest distance between the 2 structures
         double difference = 1000.0, temp_diff;
         // Apply rotation and mirror
         for (int axes = 0; axes < 3; axes++)   {
            for (int angle = 0; angle < 4; angle++)   {
               for (int mir = 0; mir < 8; mir++)   {
                  Str3 = DB_Str[i];
                  Str3.rotate(axes,angle*90); Str3.mirror(mir);
                  temp_diff = evaluate_diff(Str3,Sorted_NEW);
                  if (temp_diff < difference)
                     difference = temp_diff;
                }
             }
          }
          // If smallest difference is lower than the threshold (distance per atom), it is not new
          if (difference < distance_comparison)   {
             is_new = false;
             reason=DB_names[i];
             i=count_db;
          }
      }
   }
   if (is_new)   {
      temp = "cp " + filename + " " + database_path; system(temp.c_str());
      gout << "Structure " << filename << " is new and it was added to the database.\n";
      temp=database_path+"/"+LABEL+calc_zeros(lo);
      count_db++;
      DB_Str.resize(count_db,temp);
   }
   else
      gout << "Structure " << filename << " is similar to " << reason << ", so it will not be added to the database.\n";
