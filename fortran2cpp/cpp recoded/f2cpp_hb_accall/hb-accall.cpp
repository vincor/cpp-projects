#include <iostream>
#include "acc-process.h"

static void show_title() {
    cout << "**************************************************************************************" << endl
         << "Program Name: HB-ACCALL" << endl
         << "Description:  Accessibility calculations is the tool for the algorithm for analysis" << endl
         << "              of Asn, Gln and His side-chainscan be fully exploited." << endl
         << "**************************************************************************************" << endl << endl;
}
static void show_usage() {
    cout << "Usage:"<< endl 
         << "   hb-accall [options] <pdb_file> <vdw_radii_file>" << endl << endl;
    cout << "Options:"<< endl
         << "   -h, --help                                 Show this help message and exit." << endl
         << "   -p <probe_size>, --probe=<probe_size>      Set the probe size ( 1.4 by default )." << endl
         << "   -z <zslice_width>, --zslice=<zslice_width> Set the z-slice width (0.05 by default)." << endl
         << "   -het, --hetas                   Include HETATOMs." << endl
         << "   -hyd, --hydro                   Include HYDROGENs." << endl
         << "   -wat, --waters                  Include WATERs." << endl
         << "   -all, --all                     Include all: HETATOMs, HYDROGENs and WATERs." << endl << endl
         << "Arguments:" << endl
         << "   <pdb_file>                  The PDB file to be processed." << endl
         << "   <vdw_radii_file>            The Van der Waal radii file to be processed." << endl << endl
         << "Examples:" << endl
         << "   hb-accall pdbxxxx.new vdw.radii" << endl
         << "   hb-accall -p 1.7 pdbxxxx.new vdw.radii" << endl
         << "   hb-accall -p 1.7 --zslice=0.07 -het -wat pdbxxxx.new vdw.radii" << endl
         << "   hb-accall -p 1.7 -z 0.07 -all pdbxxxx.new vdw.radii" << endl
         << "   hb-accall -all pdbxxxx.new vdw.radii" << endl
         << "***********************************************************************************" << endl << endl;
}


int main(int argc, char *argv[]) {

    show_title();
    if (argc==1 || argc==2 && (string(argv[1])!="-h" || string(argv[1])!="--help")) {
        cout << "Incorrect args list! At least specify two input files." << endl << endl;
        show_usage();    
        return 1;
    }
    if (string(argv[1])=="-h" || string(argv[1])=="--help") {
        show_usage();
        return 0;
    }
    vector<string> args;
    for(unsigned int i=1;i<argc;args.push_back(argv[i++]));
    try {
        ProcessData pData = ProcessData(&args);
        string errMessage = pData.errstate? "with ERRORS" : "successfully";
        std::cout << "Program HB-ACCALL completed " + errMessage << endl 
                  << "Please check details in the <pdb>.log file" << endl << endl
                  << "******************************************" << endl;
    } catch (string const e) {cout << e << endl; return -1;} 

    return 0;
}
