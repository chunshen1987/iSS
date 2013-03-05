//===============================================================================
//  iSpectraSampler
//===============================================================================
//
//  Zhi Qiu & Chun Shen
//    qiu.24@asc.ohio-state.edu
//    shen.201@asc.ohio-state.edu
//
//        Date: 03/2012
// Version info written in the main function.



#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
#include <sys/time.h>
#include <stdlib.h>

#include "main.h"
#include "ParameterReader.h"
#include "Table.h"
#include "readindata.h"
#include "emissionfunction.h"
#include "arsenal.h"

using namespace std;


int main(int argc, char** argv)
{
   cout << endl
        << "                  iSpectraSampler                " << endl
        << endl
        << "  Ver 2.0.0.0  ---- Zhi Qiu & Chun Shen, 04/2012 " << endl;
   cout << endl << "**********************************************************" << endl;
   display_logo(3); // Hail to the king~
   cout << endl << "**********************************************************" << endl << endl;

   // Chun's input reading process
   string path="results";

   //load freeze out information
   int FO_length = 0;
   ostringstream decdatfile;
   cout << " -- Loading the decoupling data....";
   decdatfile << path << "/decdat2.dat";
   FO_length=get_filelength(decdatfile.str().c_str());
   cout <<" cells number: " <<  FO_length << endl;

   //read the data arrays for the decoupling information
   FO_surf* FOsurf_ptr = new FO_surf[FO_length];
   for(int i=0; i<FO_length; i++)
     for(int j=0; j<Maxparticle; j++)
         FOsurf_ptr[i].particle_mu[j] = 0.0e0;
   read_decdat(path, FO_length, FOsurf_ptr);

   //read the positions of the freeze out surface
   read_surfdat(path, FO_length, FOsurf_ptr);

   //read the chemical potential on the freeze out surface
   int N_stableparticle;
   ifstream particletable("EOS/EOS_particletable.dat");
   particletable >> N_stableparticle;
   double** particle_mu = new double* [N_stableparticle];
   for(int i=0; i<N_stableparticle; i++)
     particle_mu[i] = new double [FO_length];
   if(N_stableparticle >0)
   {
      read_decdat_mu(path, FO_length, N_stableparticle, particle_mu);
   }

   //read particle resonance decay table
   particle_info *particle = new particle_info [Maxparticle];
   for (int i=0; i<Maxparticle; i++) particle[i].decays=0; // to avoid infinite loop
   int Nparticle=read_resonance(particle);
   cout <<"particle number: " << Nparticle << endl;
   if(N_stableparticle >0)
   {
      cout << " -- EOS is partically chemical equilibrium " << endl;
      calculate_particle_mu(Nparticle, FOsurf_ptr, FO_length, particle, particle_mu);
   }
   else
   {
      cout << " -- EOS is chemical equilibrium. " << endl;
      for(int i=0; i<Nparticle; i++)
        for(int j=0; j<FO_length; j++)
           FOsurf_ptr[i].particle_mu[j] = 0.0e0;
   }
   cout << endl << " -- Read in data finished!" << endl << endl;


   // Next, Zhi's turn...

   // First other parameters; Zhi's style
   ParameterReader paraRdr;
   paraRdr.readFromFile("parameters.dat");
   paraRdr.readFromArguments(argc, argv);
   //paraRdr.echo();

   // init random seed from system time
   timeval a;
   gettimeofday(&a, 0);
   long randomSeed=paraRdr.getVal("randomSeed");
   if (randomSeed<0) randomSeed=a.tv_usec; // randomSeed<0 means to use CPU clock
   srand48(randomSeed);

   Table chosen_particles("EOS/chosen_particles.dat"); // skip others except for these particle
   Table pT_tab("tables/pT_gauss_table.dat"); // pt position and weight table
   Table phi_tab("tables/phi_gauss_table.dat"); // phi position and weight table
   Table eta_tab("tables/eta_gauss_table_20_full.dat"); // eta uniform dist table
   EmissionFunctionArray efa(&chosen_particles, &pT_tab, &phi_tab, &eta_tab, particle, Nparticle, FOsurf_ptr, FO_length, &paraRdr);

   efa.shell();
   //efa.combine_samples_to_OSCAR();
}
