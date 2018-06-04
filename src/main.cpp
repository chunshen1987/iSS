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

    string path = "results";
    iSS iSsampler(path);
    // read in parameters
    iSsampler.paraRdr_ptr->readFromFile("iSS_parameters.dat");
    iSsampler.paraRdr_ptr->readFromArguments(argc, argv);
    iSsampler.paraRdr_ptr->echo();

    int status = iSsampler.shell();
    if (status == 0) {
        cout << "Program executed normally." << endl;
    }
}
