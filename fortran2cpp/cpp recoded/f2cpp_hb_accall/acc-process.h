#ifndef ACCPROCESS_H
#define ACCPROCESS_H
#include <vector>
#include <map>
#include <string>
#include <iostream>

using namespace std;

struct ArgsSet
{
    string pdbfile;
    string vdwfile;
    double probe = 1.4;
    double zslice = 0.05;
    bool hetas = false, hydro = false, waters = false;
};
struct Cartesians
{
    double x = 999.9;
    double y = 999.9;
    double z = 999.9;
};

class ProcessData {
public:
    int errstate;
    ProcessData(vector<string> *args) {
        argSet = ArgsParsing(args);
        ReadVanDerWaalRadiiFile(argSet.vdwfile);
        ReadStandardDataFile();
        errstate = DataProcessing();

        int dotpos = argSet.pdbfile.find_last_of(".");
        string outfilename = dotpos>3? argSet.pdbfile.substr(0, dotpos):argSet.pdbfile;  

        if (logStream.size()>0 ) DataOutput(outfilename + ".log",&logStream);
        if (asaStream.size()>0 ) DataOutput(outfilename + ".asa",&asaStream);
        if (rsaStream.size()>0 ) DataOutput(outfilename + ".rsa",&rsaStream);
    }
private:    
// Constant data
    const double hyrad = 1.0;

//	Atoms classified as Non-Polar in sidechains
    const vector<string> Phobs = {" CA "," CB "," CD "," CD1"," CD2"," CG ", 
                                  " CG1"," CG2"," CE "," CE1"," CE2"," CE3",
                                  " CH2"," CZ "," CZ2"," CZ3"," SD "," SG "};
//	Atoms classified as Polar in sidechains
    const vector<string> Phils = {" AD1"," AD2"," AE1"," AE2"," ND1"," ND2"," NE ",
                                  " NE1"," NE2"," NH1"," NH2"," NZ "," OD1"," OD2",
                                  " OE3"," OE2"," OE1"," OG "," OG1"," OH "};
//	Main chain Atoms
    const vector<string> Mchain = {" N  "," C  "," O  "," OXT"};     


// ===========================
    ArgsSet argSet;
    vector<string> logStream;
    vector<string> asaStream;
    vector<string> rsaStream;

    
    vector<string> aacids, stacids;
    vector<unsigned int> numats;
    vector<vector<string>> anames;
    vector<vector<double>> vradii, standarea;

    ArgsSet ArgsParsing(vector<string> *args);
    void ReadVanDerWaalRadiiFile(string vdwname);
    void ReadStandardDataFile();
    int DataProcessing();
    void DataOutput(string outFileName, vector<string> *out_stream);

    vector<double> SOLVA(const vector<Cartesians> *sCoords, const vector<double> *sRads,
                         const double sProbe, const double sZslice);
    int getVecIndex(const vector<string> *invec, const string s2find);
    vector<unsigned int> SortIndex(vector<double> *arrIn);
    string sformat(double num, const unsigned int slen, const unsigned int digpos);

};

#endif