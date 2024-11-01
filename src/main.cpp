//=====================================================================
//  iSpectraSampler
//=====================================================================
//
//    Copyright @ 2012 Zhi Qiu & Chun Shen
//    qiu.24@asc.ohio-state.edu
//    shen.201@asc.ohio-state.edu
//
//        Date: 03/2012
// Version info written in the main function.
//=====================================================================

#include <iostream>
#include <string>

#include "arsenal.h"
#include "iSS.h"

using std::cout;
using std::endl;

int main(int argc, char** argv) {
    cout << endl
         << "                  iSpectraSampler                " << endl
         << endl
         << "  Ver 2.0.0.0  ---- Zhi Qiu & Chun Shen, 04/2012 " << endl;
    cout << endl
         << "**********************************************************"
         << endl;
    display_logo(3);  // Hail to the king~
    cout << endl
         << "**********************************************************" << endl
         << endl;

    const string table_path = "iSS_tables";
    const string particle_table_path = "iSS_tables";

    string path = "results";
    string input_file = "iSS_parameters.dat";
    string surface_filename = "surface.dat";
    if (argc > 1) {
        std::string s1 = argv[1];
        if (s1.find("=") == std::string::npos) input_file = s1;
    }
    cout << "input file : " << input_file << endl;
    if (argc > 2) {
        std::string s1 = argv[2];
        if (s1.find("=") == std::string::npos) path = s1;
    }
    cout << "work folder path : " << path << endl;
    if (argc > 3) {
        std::string s1 = argv[3];
        if (s1.find("=") == std::string::npos) surface_filename = s1;
    }
    cout << "surface filename : " << surface_filename << endl;

    iSS iSsampler(
        path, table_path, particle_table_path, input_file, surface_filename);
    // read in parameters
    iSsampler.paraRdr_ptr->readFromArguments(argc, argv);
    iSsampler.paraRdr_ptr->echo();

    int status = iSsampler.shell();

    int testFlag =
        static_cast<int>(iSsampler.paraRdr_ptr->getVal("perform_checks", 0));
    if (testFlag == 1) iSsampler.perform_checks();

    if (status == 0) cout << "Program executed normally." << endl;
    return (status);
}
