
#include <sys/time.h>
#include "./Table.h"
#include "./arsenal.h"
#include "./iSS.h"

using namespace std;

iSS::iSS(ParameterReader *paraRdr_in) {
    path = "results";
    FO_length = 0;
    paraRdr_ptr = paraRdr_in;
    Nparticle = 0;
    flag_PCE = 0;

    FOsurf_ptr = NULL;
    particle = NULL;
}

iSS::~iSS() {
    if (FOsurf_ptr != NULL) {
        delete[] FOsurf_ptr;
    }

    if (particle != NULL) {
        delete[] particle;
    }
}

int iSS::shell() {
    read_in_FO_surface();
    generate_samples();
    return(0);
}

int iSS::read_in_FO_surface() {
    read_FOdata freeze_out_data(paraRdr_ptr, path);
    FO_length = freeze_out_data.get_number_of_freezeout_cells();
    cout << "total number of cells: " <<  FO_length << endl;
    FOsurf_ptr = new FO_surf[FO_length];
    freeze_out_data.read_in_freeze_out_data(FO_length, FOsurf_ptr);
    particle = new particle_info[Maxparticle];
    Nparticle = freeze_out_data.read_in_chemical_potentials(
                                    path, FO_length, FOsurf_ptr, particle);
    flag_PCE = freeze_out_data.get_flag_PCE();
    cout << endl << " -- Read in data finished!" << endl << endl;
    return(0);
}

int iSS::generate_samples() {
    timeval a;
    gettimeofday(&a, 0);
    long randomSeed = paraRdr_ptr->getVal("randomSeed");
    // randomSeed<0 means to use CPU clock
    if (randomSeed < 0) {
        randomSeed = a.tv_usec;
    }
    srand48(randomSeed);
    
    // skip others except for these particle
    Table chosen_particles("EOS/chosen_particles.dat");

    Table pT_tab("tables/pT_gauss_table.dat");    // pt tables
    Table phi_tab("tables/phi_gauss_table.dat");  // phi tables
    // eta uniform dist table
    Table eta_tab("tables/eta_uni_table.dat");
    // Table eta_tab("tables/eta_gauss_table_30_full.dat");

    EmissionFunctionArray efa(&chosen_particles, &pT_tab, &phi_tab, &eta_tab,
                              particle, Nparticle, FOsurf_ptr, FO_length,
                              flag_PCE, paraRdr_ptr);
    efa.shell();
    return(0);
}
