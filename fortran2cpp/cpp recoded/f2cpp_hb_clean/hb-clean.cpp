#include <iostream>
#include "brk-process.h"

static void show_title() {
    cout << "****************************************************************************" << endl
         << "Program Name: HB-CLEAN" << endl
         << "Description:  Clean-up program for mistakes in brookhaven files (PDB files)." << endl 
         << "              Written by David Keith Smith, 1989." << endl 
         << "              Amended by Roman Laskowski, 1992." << endl 
         << "              Translated to C++ by Demid Vinnikov, 2024." << endl
         << "****************************************************************************" << endl << endl;
}
static void show_usage() {
    cout << "Usage:"<< endl 
         << "   hb-clean [options] <pdb_file>" << endl << endl;
    cout << "Options:"<< endl
         << "   -h, --help          Show this help message and exit." << endl << endl
         << "Arguments:" << endl
         << "   <pdb_file>          The PDB file to be processed." << endl << endl
         << "Examples:" << endl
         << "   hb-clean pdbxxxx.ent" << endl
         << "   hb-clean ./data/pdbxxxx.ent" << endl
         << "****************************************************************************" << endl << endl;
}


int main(int argc, char *argv[]) {
    show_title();
    if (argc<2) {
        cout << "Incorrect usage! Please specify Brookhaven Protein Data Bank (PDB) file as an argument." << endl << endl;
        show_usage();    
        return 1;
    }
    if (string(argv[1])=="-h" || string(argv[1])=="--help") {
        show_usage();
        return 0;
    }
    try {
        ProcessData pData = ProcessData(argv[1]);
        cout << "* Program completed" << endl;
    } catch (string const e) {cout << e << endl; return -1;} 

    return 0;
}
