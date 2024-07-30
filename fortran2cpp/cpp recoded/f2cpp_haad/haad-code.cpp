#include <iostream>
#include "datastructures.h"



int main(int argc, char *argv[]) {
    if (argc==1) {
        std::cout << "Use PDB file name as a parameter in command line (haad-code.exe <pdbfile>)" << std::endl;
        return -1;
    }
    try {
        ServiceData sData = ServiceData();
        ProcessData pData = ProcessData(argv[1], &sData);
        std::cout << "Atoms loaded: " << pData.inAtoms << std::endl;
        std::cout << "Output Atoms: " << pData.outAtoms << std::endl; 
// Example usage:
/*         std::cout << "Amino Acid: " << sData.aminoAcids[0].name << std::endl;
        std::cout << "Heavy Atoms: " << sData.aminoAcids[0].heavyAtoms << std::endl;
        std::cout << "Total Atoms: " << sData.aminoAcids[0].totalAtoms << std::endl;
        std::cout << "Atom Names: ";
        for (const string& atomName : sData.aminoAcids[0].atomNames) {
            std::cout << atomName << " ";
        }
        std::cout << std::endl;
        std::cout << "Atom Connections: ";
        for (const int conn : sData.aminoAcids[0].atomConnections) {
                std::cout << conn << " ";
        }
        std::cout << std::endl;
        std::cout << "bdc: " << endl;
        for (auto it = sData.aminoAcids[0].bdc.begin(); it!=sData.aminoAcids[0].bdc.end(); it++) {
            for (auto it1 = it->begin(); it1!=it->end(); it1++) {
                cout << (*it1) << " ";  
            }
            cout << endl;  
        }

        for (const string& conn : sData.atomConType) {
                std::cout << conn << " ";
        }
        std::cout << std::endl; */
        
    } catch (string const e) {std::cout << e << std::endl; return -1;} 

    return 0;
}
