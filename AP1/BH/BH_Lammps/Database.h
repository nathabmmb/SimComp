   // Creating the required folder if they don't exist
   // temp = "mkdir -p " + database_path; system(temp.c_str());
/*
   // Checking what already exists - determining the step
   int count_db = 0;
   vector<string> DB_names;
   temp = "ls " + database_path + " 2> /dev/null > " + LABEL; system(temp.c_str());
   ifstream database(LABEL.c_str(),ios::in); 
   while (database >> temp)   {
      count_db++;
      temp=database_path+"/"+temp;
      DB_names.resize(count_db,temp);
   }
   database.close();
   gout << "\n" << count_db << " files read on the database folder given by the user. Storing them...\n";

   // Constructor: putting every structure from the database on the class structures
   Cluster TEMP;
   vector <Cluster> DB_Str;
   for (int i = 0; i < count_db;i++)
      DB_Str.resize(i+1,DB_names[i]);   
   gout << "Success. Proceeding.\n";*/
