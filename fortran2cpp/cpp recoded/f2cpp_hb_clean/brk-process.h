#ifndef BRKPROCESS_H
#define BRKPROCESS_H
#include <vector>
#include <map>
#include <string>
#include "brkcln.h"

using namespace std;

class ProcessData {
public:
    ProcessData(const char *fn) {
        DataProcessing(fn);
        int dotpos = ((string)fn).find_last_of(".");
        string outfilename = dotpos>3? ((string)fn).substr(0, dotpos):(string)fn;  
        if (outStream.size()>0 ) DataOutput(outfilename + ".new",&outStream);
        if (altStream.size()>0 ) DataOutput(outfilename + ".alt",&altStream);
    }
private:    
    vector<string> outStream;
    vector<string> altStream;
    ConstSet bConst;    
    void DataProcessing(const char *fn);
    void DataOutput(string outFileName, vector<string> *out_stream);
    int getVecIndex(const vector<string> *invec, const string s2find);
    double getDist(const vector<double> *avec, const vector<double> *bvec);
    double diHed(const int angnum, const vector<double> *va, const vector<double> *vb, const vector<double> *vc, const vector<double> *vd);

};

#endif