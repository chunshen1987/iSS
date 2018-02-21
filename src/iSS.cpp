
#include <sys/time.h>
#include "./Table.h"
#include "./arsenal.h"
#include "./iSS.h"

using namespace std;

iSS::iSS(string path_in) {
    path = path_in;
    FO_length = 0;
    Nparticle = 0;
    flag_PCE = 0;

    FOsurf_ptr = nullptr;
    particle = nullptr;
    efa = nullptr;

    paraRdr_ptr = new ParameterReader;
}

iSS::~iSS() {
    if (FOsurf_ptr != nullptr) {
        delete[] FOsurf_ptr;
    }

    if (particle != nullptr) {
        delete[] particle;
    }

    if (efa != nullptr) {
        delete efa;
    }
}

int iSS::shell() {
    int status = read_in_FO_surface();
    if (status != 0) {
        cout << "Some errors happened in reading in the hyper-surface" << endl;
        exit(-1);
    }
    set_random_seed();
    status = generate_samples();
    if (status != 0) {
        cout << "Some errors happened in generating particle samples" << endl;
        exit(-1);
    }
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

void iSS::set_random_seed(int randomSeed_in) {
    randomSeed = randomSeed_in;
    srand48(randomSeed);
}

void iSS::set_random_seed() {
    randomSeed = paraRdr_ptr->getVal("randomSeed");
    if (randomSeed < 0) {
        timeval a;
        gettimeofday(&a, 0);
        randomSeed = a.tv_usec;
    }
    srand48(randomSeed);
}

int iSS::generate_samples() {
    // skip others except for these particle
    Table chosen_particles(table_path + "/chosen_particles.dat");

    Table pT_tab(table_path + "/bin_tables/pT_gauss_table.dat");
    Table phi_tab(table_path + "/bin_tables/phi_gauss_table.dat");
    // eta uniform dist table
    Table eta_tab(table_path + "/bin_tables/eta_uni_table.dat");
    // Table eta_tab("tables/eta_gauss_table_30_full.dat");

    efa = new EmissionFunctionArray(
            &chosen_particles, &pT_tab, &phi_tab, &eta_tab,
            particle, Nparticle, FOsurf_ptr, FO_length,
            flag_PCE, paraRdr_ptr, path);
    efa->shell();

    return(0);
}

