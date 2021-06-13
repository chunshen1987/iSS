
#include <memory>

#include "Table.h"
#include "arsenal.h"
#include "iSS.h"
#include "Random.h"
#include "data_struct.h"
#include "Histogram.h"

//using namespace std;
using iSS_data::Vec4;

iSS::iSS(std::string path, std::string table_path,
         std::string particle_table_path, std::string inputfile,
         std::string surface_filename) :
        path_(path), table_path_(table_path),
        particle_table_path_(particle_table_path),
        surface_filename_(surface_filename) {

    flag_PCE_ = 0;
    paraRdr_ptr = new ParameterReader;
    paraRdr_ptr->readFromFile(inputfile);

}


iSS::~iSS() {
    clear();
    delete paraRdr_ptr;
}


void iSS::clear() {
    FOsurf_array_.clear();
    FOsurf_LRF_array_.clear();
    particle.clear();
}


int iSS::shell() {
    int status = read_in_FO_surface();
    if (status != 0) {
        messager << "Some errors happened in reading in the hyper-surface";
        messager.flush("error");
        exit(-1);
    }
    set_random_seed();
    status = generate_samples();
    if (status != 0) {
        messager << "Some errors happened in generating particle samples";
        messager.flush("error");
        exit(-1);
    }
    return(0);
}


void iSS::perform_checks() {
    messager.info("Performing checks for the samples ...");
    construct_Tmunu();
    Histogram hist_pion(0., 5., 100);
    Histogram hist_proton(0., 5., 100);
    int nev = get_number_of_sampled_events();
    for (int iev = 0; iev < nev; iev++) {
        int npart = get_number_of_particles(iev);
        for (int ipart = 0; ipart < npart; ipart++) {
            iSS_Hadron part_i = get_hadron(iev, ipart);
            //double p[4] = {part_i.E, part_i.px, part_i.py, part_i.pz};
            if (part_i.pid == 211) {
                double pT = sqrt(part_i.px*part_i.px + part_i.py*part_i.py);
                hist_pion.fill(pT, 1.);
            } else if (part_i.pid == 2212) {
                double pT = sqrt(part_i.px*part_i.px + part_i.py*part_i.py);
                hist_proton.fill(pT, 1.);
            }
        }
        hist_pion.add_an_event();
        hist_proton.add_an_event();
    }
    hist_pion.output_histogram("check_211_spectra.dat");
    hist_proton.output_histogram("check_2212_spectra.dat");
}


int iSS::read_in_FO_surface() {
    std::vector<FO_surf> FOsurf_temp;
    read_FOdata freeze_out_data(paraRdr_ptr, path_, table_path_,
                                particle_table_path_);
    //freeze_out_data.read_in_freeze_out_data(FOsurf_array_);
    freeze_out_data.read_in_freeze_out_data(FOsurf_temp, surface_filename_);
    messager << "total number of cells: " <<  FOsurf_temp.size();
    messager.flush("info");

    afterburner_type_ = freeze_out_data.get_afterburner_type();
    freeze_out_data.read_in_chemical_potentials(FOsurf_temp, particle);
    flag_PCE_ = freeze_out_data.get_flag_PCE();

    if (paraRdr_ptr->getVal("MC_sampling") == 4) {
        transform_to_local_rest_frame(FOsurf_temp, FOsurf_LRF_array_);
    } else {
        FOsurf_array_ = FOsurf_temp;
    }
    messager.info(" -- Read in data finished!");
    FOsurf_temp.clear();
    return(0);
}


void iSS::set_random_seed() {
    randomSeed_  = paraRdr_ptr->getVal("randomSeed");
    ran_gen_ptr_ = std::shared_ptr<RandomUtil::Random>(
                                        new RandomUtil::Random(randomSeed_));
}


void iSS::set_random_seed(int randomSeed_in) {
    randomSeed_ = randomSeed_in;
    ran_gen_ptr_ = std::shared_ptr<RandomUtil::Random>(
                                        new RandomUtil::Random(randomSeed_));
}


int iSS::generate_samples() {
    // skip others except for these particle
    Table chosen_particles;
    if (afterburner_type_ == AfterburnerType::SMASH) {
        chosen_particles.loadTableFromFile(
            particle_table_path_ + "/chosen_particles_SMASH.dat");
    } else if (afterburner_type_ == AfterburnerType::UrQMD) {
        chosen_particles.loadTableFromFile(
            particle_table_path_ + "/chosen_particles_urqmd_v3.3+.dat");
    } else {
        chosen_particles.loadTableFromFile(particle_table_path_
                                           + "/chosen_particles_s95p-v1.dat");
    }

    messager.info("Start computation and generating samples ...");
    if (paraRdr_ptr->getVal("MC_sampling") == 4) {
        spectra_sampler_ = std::unique_ptr<FSSW> (new FSSW(
                    ran_gen_ptr_, &chosen_particles,
                    particle, FOsurf_LRF_array_, flag_PCE_, paraRdr_ptr,
                    path_, table_path_, afterburner_type_));
        spectra_sampler_->shell();
    } else {
        Table pT_tab(table_path_ + "/bin_tables/pT_gauss_table.dat");
        Table phi_tab(table_path_ + "/bin_tables/phi_gauss_table.dat");
        // eta uniform dist table
        Table eta_tab(table_path_ + "/bin_tables/eta_uni_table.dat");
        //Table eta_tab("tables/eta_gauss_table_30_full.dat");
        efa_ = std::unique_ptr<EmissionFunctionArray> (
            new EmissionFunctionArray(
                ran_gen_ptr_, &chosen_particles, &pT_tab, &phi_tab, &eta_tab,
                particle, FOsurf_array_, flag_PCE_, paraRdr_ptr,
                path_, table_path_, afterburner_type_));
        efa_->shell();
    }
    return(0);
}


// this function transform all the variables to local rest frame of the fluid
// cell and trasform them to the t-z coordinate
void iSS::transform_to_local_rest_frame(
        std::vector<FO_surf> &FOsurf_ptr,
        std::vector<FO_surf_LRF> &FOsurf_LRF_ptr) {
    messager.info("Transforming fluid cells to their local rest frame ...");
    for (auto &surf_i: FOsurf_ptr) {
        FO_surf_LRF surf_LRF_i;
        surf_LRF_i.tau = surf_i.tau;
        surf_LRF_i.xpt = surf_i.xpt;
        surf_LRF_i.ypt = surf_i.ypt;
        surf_LRF_i.eta = surf_i.eta;
        surf_LRF_i.Edec = surf_i.Edec;
        surf_LRF_i.Tdec = surf_i.Tdec;
        surf_LRF_i.Pdec = surf_i.Pdec;
        surf_LRF_i.Bn = surf_i.Bn;
        surf_LRF_i.muB = surf_i.muB;
        surf_LRF_i.muS = surf_i.muS;
        surf_LRF_i.muQ = surf_i.muC;
        surf_LRF_i.bulkPi = surf_i.bulkPi;
        surf_LRF_i.particle_mu_PCE = surf_i.particle_mu_PCE;
        float cosh_eta = cosh(surf_i.eta);
        float sinh_eta = sinh(surf_i.eta);
        float ut = surf_i.u0*cosh_eta + surf_i.u3*sinh_eta;
        float uz = surf_i.u3*cosh_eta + surf_i.u0*sinh_eta;
        float ux = surf_i.u1;
        float uy = surf_i.u2;
        surf_LRF_i.u_tz[0] = ut;
        surf_LRF_i.u_tz[1] = ux;
        surf_LRF_i.u_tz[2] = uy;
        surf_LRF_i.u_tz[3] = uz;
        double LorentzBoost[4][4] = {
            {ut, -ux, -uy, -uz},
            {-ux, 1. + ux*ux/(ut + 1.), ux*uy/(ut + 1.), ux*uz/(ut + 1.)},
            {-uy, ux*uy/(ut + 1.), 1. + uy*uy/(ut + 1.), uy*uz/(ut + 1.)},
            {-uz, ux*uz/(ut + 1.), uy*uz/(ut + 1.), 1. + uz*uz/(ut + 1.)}
        };
        Vec4 da_upper = {
            surf_i.tau*surf_i.da0*cosh_eta - surf_i.da3*sinh_eta,
            -surf_i.tau*surf_i.da1,
            -surf_i.tau*surf_i.da2,
            -surf_i.da3*cosh_eta + surf_i.tau*surf_i.da0*sinh_eta};
        Vec4 da_LRF = {0., 0., 0., 0.};
        //double udotdsimga = surf_i.tau*(
        //    surf_i.u0*surf_i.da0 + surf_i.u1*surf_i.da1 + surf_i.u2*surf_i.da2
        //    + surf_i.u3*surf_i.da3/surf_i.tau);
        //double udotdsimga2 = (ut*da_upper[0] - ux*da_upper[1]
        //                      - uy*da_upper[2] - uz*da_upper[3]);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                da_LRF[i] += LorentzBoost[i][j]*da_upper[j];
            }
            if (i == 0) {
                surf_LRF_i.da_mu_LRF[i] = da_LRF[i];
            } else {
                surf_LRF_i.da_mu_LRF[i] = -da_LRF[i];
            }
        }
        if (surf_LRF_i.da_mu_LRF[0] < 0) continue;
        //cout << "check: " << udotdsimga << "  " << udotdsimga2 << "  "
        //     << surf_i.da_mu_LRF[0] << endl;

        // transform the diffusion current
        Vec4 qmu_tz = {
            surf_i.qmu0*cosh_eta + surf_i.qmu3*sinh_eta,
            surf_i.qmu1,
            surf_i.qmu2,
            surf_i.qmu3*cosh_eta + surf_i.qmu0*sinh_eta};
        Vec4 qmu_LRF = {0., 0., 0., 0.};
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                qmu_LRF[i] += LorentzBoost[i][j]*qmu_tz[j];
            }
        }
        //cout << qmu_LRF[0] << endl;
        surf_LRF_i.qmuLRF_x = qmu_LRF[1];
        surf_LRF_i.qmuLRF_y = qmu_LRF[2];
        surf_LRF_i.qmuLRF_z = qmu_LRF[3];

        // transform the shear stress tensor
        float pi_tz[4][4];
        pi_tz[0][0] = (  surf_i.pi00*cosh_eta*cosh_eta
                       + 2.*surf_i.pi03*cosh_eta*sinh_eta
                       + surf_i.pi33*sinh_eta*sinh_eta);
        pi_tz[0][1] = surf_i.pi01*cosh_eta + surf_i.pi13*sinh_eta;
        pi_tz[0][2] = surf_i.pi02*cosh_eta + surf_i.pi23*sinh_eta;
        pi_tz[0][3] = (  surf_i.pi00*cosh_eta*sinh_eta
                       + surf_i.pi03*(cosh_eta*cosh_eta + sinh_eta*sinh_eta)
                       + surf_i.pi33*sinh_eta*cosh_eta);
        pi_tz[1][0] = pi_tz[0][1];
        pi_tz[1][1] = surf_i.pi11;
        pi_tz[1][2] = surf_i.pi12;
        pi_tz[1][3] = surf_i.pi01*sinh_eta + surf_i.pi13*cosh_eta;
        pi_tz[2][0] = pi_tz[0][2];
        pi_tz[2][1] = surf_i.pi12;
        pi_tz[2][2] = surf_i.pi22;
        pi_tz[2][3] = surf_i.pi02*sinh_eta + surf_i.pi23*cosh_eta;
        pi_tz[3][0] = pi_tz[0][3];
        pi_tz[3][1] = pi_tz[1][3];
        pi_tz[3][2] = pi_tz[2][3];
        pi_tz[3][3] = (  surf_i.pi00*sinh_eta*sinh_eta
                       + 2.*surf_i.pi03*sinh_eta*cosh_eta
                       + surf_i.pi33*cosh_eta*cosh_eta);
        float pi_LRF[4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                pi_LRF[i][j] = 0.;
                for (int a = 0; a < 4; a++) {
                    for (int b = 0; b < 4; b++) {
                        pi_LRF[i][j] += (LorentzBoost[i][a]*pi_tz[a][b]
                                         *LorentzBoost[b][j]);
                    }
                }
            }
        }
        //std::cout << pi_LRF[0][0] << "  " << pi_LRF[0][1] << "  "
        //          << pi_LRF[0][2] << "  " << pi_LRF[0][3] << std::endl;
        surf_LRF_i.piLRF_xx = pi_LRF[1][1];
        surf_LRF_i.piLRF_xy = pi_LRF[1][2];
        surf_LRF_i.piLRF_xz = pi_LRF[1][3];
        surf_LRF_i.piLRF_yy = pi_LRF[2][2];
        surf_LRF_i.piLRF_yz = pi_LRF[2][3];

        FOsurf_LRF_ptr.push_back(surf_LRF_i);
    }
}


void iSS::construct_Tmunu() {
    messager.info("Constructing the fluid cell T^{mu nu} from samples ...");
    double volume = FOsurf_LRF_array_[0].da_mu_LRF[0];
    double T[4][4];
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            T[i][j] = 0.;
    int nev = get_number_of_sampled_events();
    for (int iev = 0; iev < nev; iev++) {
        int npart = get_number_of_particles(iev);
        for (int ipart = 0; ipart < npart; ipart++) {
            iSS_Hadron part_i = get_hadron(iev, ipart);
            double p[4] = {part_i.E, part_i.px, part_i.py, part_i.pz};
            for (int ii = 0; ii < 4; ii++)
                for (int jj = 0; jj < 4; jj++)
                    T[ii][jj] += p[ii]*p[jj]/p[0];
        }
    }
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            T[i][j] /= (nev*volume);
        }
    }
    messager << "check: e = " << T[0][0] << " GeV/fm^3, diff = "
             << T[0][0] - FOsurf_LRF_array_[0].Edec;
    messager.flush("info");
    double trace = T[1][1] + T[2][2] + T[3][3];
    double Pi = trace/3. - FOsurf_LRF_array_[0].Pdec;
    double pizz = T[3][3] - Pi - FOsurf_LRF_array_[0].Pdec;
    messager << "check: pi_zz = " << pizz << " GeV/fm^3, diff = "
              << pizz + (  FOsurf_LRF_array_[0].piLRF_xx
                         + FOsurf_LRF_array_[0].piLRF_yy);
    messager.flush("info");
    messager << "check: Bulk Pi = " << Pi << " GeV/fm^3, diff = "
             << Pi - FOsurf_LRF_array_[0].bulkPi;
    messager.flush("info");
}
