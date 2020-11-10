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

#include <stdlib.h>
#include <iostream>
#include <string>

#include "data_struct.h"
#include "iSS.h"
#include "arsenal.h"

using namespace std;

int main(int argc, char** argv) {
    cout << endl
         << "                  iSpectraSampler                " << endl
         << endl
         << "  Ver 2.0.0.0  ---- Zhi Qiu & Chun Shen, 04/2012 " << endl;
    cout << endl << "**********************************************************"
         << endl;
    display_logo(3);  // Hail to the king~
    cout << endl << "**********************************************************"
         << endl << endl;

    const string table_path = "iSS_tables";
    const string particle_table_path = "iSS_tables";

    string path = "results";
    string input_file = "iSS_parameters.dat";
    if (argc > 1)
        input_file = argv[1];
    cout << "input file : " << input_file << endl;
    if (argc > 2)
        path = argv[2];
    cout << "work folder path : " << path << endl;

    iSS iSsampler(path, table_path, particle_table_path, input_file);
    // read in parameters
    iSsampler.paraRdr_ptr->readFromArguments(argc, argv);
    iSsampler.paraRdr_ptr->echo();

    int status = iSsampler.shell();
    if (status == 0) {
        cout << "Program executed normally." << endl;
    }
    return(0);
}
