#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include "brk-process.h"

using namespace std;

//************************ Brk file processing and output data processing *****************************************
void ProcessData::DataProcessing(const char *fn) {
    ifstream brkfi;
    string brkrec;
    // Open the file 
    brkfi.open(fn);
    if (!brkfi) throw "No Brookhaven file or directory! (" + (string)fn + ")";
    
    unsigned int cnt = 0; // total number of lines from the input file
    bool caonly = true, altern = false, startSeqFlag = true, specialCaseFlag = false, altset, endKey = false;
    double ave = -9.99, sd = -9.99, oldocc, oldthm, extocc, occup, therm;
    unsigned int seqlen = 0, aacnt = 0, scatct=0, thmcnt = 0, atmnum = 0, ninchn = 1, chnct = 1, chnst = 0;
    unsigned int natalt = 0, nhealt = 0;  
    int attype, sctype,  aacode = -1, nrec;
    string setatm, resnam, setloc, seqnum, oldatm, extatm, extnum, atmnam, altloc, altnum, wkres;
    vector<double> acoord = {0.0, 0.0, 0.0}, mcangs;
    vector<string> seqbcd, thmprs, chstrt;
    vector<string> inseqr, label, altmcc, altscc;
    vector<int> seqcod, chnsz;
    vector<vector<double>> mcocc, scocc, mcthrm, scthrm, scangs;
    vector<vector<int>> mcser, scser;
    vector<vector<string>> mcalt, scalt;
    vector<vector<vector<double>>> mcaaat, scaaat;
    
    vector<vector<bool>> atomin, scatin;
    vector<vector<int>> altmci = {{},{},{}}, altsci = {{},{},{}};
    vector<vector<double>> altmcr = {{},{}}, altmca = {{},{},{}}, altscr = {{},{}}, altsca = {{},{},{}};

    struct ExtDat {
        string dat;
        unsigned int res;
    };
    vector<ExtDat> extData;
    vector<ExtDat> altData;
    vector<int> totswp(bConst.sideaa, 0), restot(bConst.sideaa, 0); 

    for(unsigned int i=0; i<bConst.mnchna; i++) {
        mcocc.push_back({});
        mcthrm.push_back({});
        atomin.push_back({});
        mcser.push_back({});
        mcalt.push_back({});
        mcaaat.push_back({});
    }
    for(unsigned int i=0; i<bConst.sdchna; i++) {
        scocc.push_back({});
        scthrm.push_back({}); 
        scatin.push_back({});
        scser.push_back({});
        scalt.push_back({});
        scaaat.push_back({});
    }
    for(unsigned int i=0; i<bConst.maxchn; i++) chnsz.push_back(0);
    for(unsigned int i=0; i<bConst.nsdang; i++) scangs.push_back({});

// Routine to read the brookhaven data file and extract the sequence
//    and the atom names and coordinates. Checking is also done for
//    residue compatibility between sequence and atom. C-alpha only
//    files are flagged. The sequence residues are expected together,
//    then the secondary structure records (sheets together), and then
//    the atom records which are expected in residue order.

    while (getline(brkfi, brkrec) && !endKey) {
        bool nextRec = false;
        while (!nextRec) {
            nextRec = true;
            if (brkrec.length()<=80) {
                string keywd = brkrec.substr(0,bConst.keylen);
                if (keywd == bConst.endkey) {endKey = true; break;}
                if (brkrec.length()<66) throw "Broken PDB file! Too short record line. Line number: " + to_string(cnt+1); 
                if (keywd != bConst.terkey) {
                    if (keywd != bConst.maskey) {
                        if (keywd == bConst.modkey) throw "*** ERROR: MODEL record encountered unexpectedly in PDB file! Line number: " + to_string(cnt+1);
                        if (aacnt==0 && keywd != bConst.atmkey && keywd != bConst.htakey) outStream.push_back(brkrec);
                        if (keywd == bConst.seqkey) {
                            if (brkrec.length()>18) {
                                string aseq = brkrec.substr(18);
                                for(unsigned int i=0; i*4<aseq.length();i++) {
                                    string atcode = aseq.substr(i*4+1,3);
                                    if (atcode.size()==3 && atcode!="   " && getVecIndex(&inseqr, atcode) < 0) inseqr.push_back(atcode);
                                } 
                            }
                        }
                        if (keywd == bConst.atmkey || keywd == bConst.htakey) {
                            if (startSeqFlag || specialCaseFlag) {
                                if (!specialCaseFlag) {
                                    aacode = -1;
                                    altset = false;
                                    setatm = "    ";
                                    try {
                                        oldocc = brkrec.substr(54,6) == "      "? 1.0: stod(brkrec.substr(54,6));
                                        oldthm = stod(brkrec.substr(60,6));
                                    } catch (invalid_argument& e) {throw "Broken PDB file! Invalid numerical data. Line number: " + to_string(cnt+1);}
                                    oldatm = brkrec.substr(12,4);
                                    setloc = brkrec.substr(16,1);
                                    resnam = brkrec.substr(17,3);
                                    seqnum = brkrec.substr(21,6);
                                    aacode = getVecIndex(&bConst.aacids,resnam);
                                    if (aacode>=0) restot.at(aacode)++;
                                }
                                if (!specialCaseFlag && (inseqr.size()>0 && getVecIndex(&inseqr,resnam)<0 ||
                                                                            getVecIndex(&bConst.dudcod,resnam)>=0 || 
                                                                            getVecIndex(&bConst.modcod,resnam) >=0 || aacode < 0)) {
                                    if (keywd != bConst.htakey && aacode < 0) cout << " '" << resnam <<"': unknown amino acid code" << endl;
                                    if (setloc != " " && extData.size() > 0) {
                                        extatm = extData.back().dat.substr(12,4);
                                        extnum = extData.back().dat.substr(21,6);
                                        try {
                                            extocc = extData.back().dat.substr(54,6) == "      "? 1.0: stod(extData.back().dat.substr(54,6));
                                        } catch (invalid_argument& e) {throw "Broken PDB file (Ext Data)!  Invalid numerical data. Line number: " + to_string(cnt+1);}
                                        if (extatm == oldatm && extnum == seqnum) {
                                            if (extocc >= oldocc) {
                                                altData.push_back({brkrec.substr(0,66), aacnt});
                                            } else {
                                                altData.push_back({extData.back().dat, aacnt});
                                                extData.back().dat = brkrec.substr(0,66);
                                            }
                                        } else extData.push_back({brkrec.substr(0,66), aacnt});
                                    } else extData.push_back({brkrec.substr(0,66), aacnt});
                                    startSeqFlag = false;
                                    continue;
                                } else {
                                    if (!specialCaseFlag) {
                                        aacnt++;
                                        seqcod.push_back(aacode);
                                        seqbcd.push_back(seqnum);
                                        chstrt.push_back("  ");
                                        label.push_back(keywd);
                                        for(unsigned i=0;i<bConst.mnchna;i++) {
                                            mcocc.at(i).push_back(0.0);
                                            mcthrm.at(i).push_back(0.0);
                                            atomin.at(i).push_back(false);
                                            mcser.at(i).push_back(0);
                                            mcalt.at(i).push_back(" ");
                                            mcaaat.at(i).push_back({0.0, 0.0, 0.0}); 
                                        }
                                        for(unsigned i=0;i<bConst.sdchna;i++) {
                                            scocc.at(i).push_back(0.0);
                                            scthrm.at(i).push_back(0.0); 
                                            scatin.at(i).push_back(false);
                                            scser.at(i).push_back(0);
                                            scalt.at(i).push_back(" ");
                                            scaaat.at(i).push_back({0.0, 0.0, 0.0});
                                        }
                                        mcangs.push_back(bConst.null);
                                        for(unsigned int i=0; i<bConst.nsdang; i++) scangs.at(i).push_back(bConst.null);
                                    }
                                    specialCaseFlag = false;
                                    atmnam = brkrec.substr(12,4);
                                    altloc = brkrec.substr(16,1);
                                    altnum = brkrec.substr(21,6);
                                    try {
                                        atmnum = stoi(brkrec.substr(6,5));
                                        acoord.at(0) = stod(brkrec.substr(30,8));
                                        acoord.at(1) = stod(brkrec.substr(38,8));
                                        acoord.at(2) = stod(brkrec.substr(46,8));
                                        occup = brkrec.substr(54,6) == "      "? 1.0: stod(brkrec.substr(54,6));
                                        therm = stod(brkrec.substr(60,6)); 
                                    } catch (invalid_argument& e) {throw "Broken PDB file! Invalid numerical data. Line number: " + to_string(cnt+1);}
                                    
                                    if (atmnam==" OXT" || atmnam==" NXT" || atmnam.substr(1,1)=="H" || atmnam.substr(1,1)=="D" || atmnam.substr(1,1)=="Q" ) {
                                        if (altloc != " " && extData.size() > 0) {
                                            extatm = extData.back().dat.substr(12,4);
                                            extnum = extData.back().dat.substr(21,6);
                                            try {
                                                extocc = extData.back().dat.substr(54,6) == "      "? 1.0: stod(extData.back().dat.substr(54,6));
                                            } catch (invalid_argument& e) {throw "Broken PDB file (Ext Data)! Invalid numerical data. Line number: " + to_string(cnt+1);}
                                            if (extatm == atmnam && extnum == altnum) {
                                                if (extocc >= occup) {
                                                    altData.push_back({brkrec.substr(0,66), aacnt});
                                                } else {
                                                    altData.push_back({extData.back().dat, aacnt});
                                                    extData.back().dat = brkrec.substr(0,66);
                                                }
                                            } else extData.push_back({brkrec.substr(0,66), aacnt});
                                        } else extData.push_back({brkrec.substr(0,66), aacnt});
                                        startSeqFlag = false;
                                        continue;                        
                                    }
                                    attype = getVecIndex(&bConst.atname,atmnam);
                                    if (attype<0) {
                                        if (aacode > 43 && aacode <= 49) attype = getVecIndex(&bConst.batnam,atmnam);
                                        else if (aacode == 41) attype = getVecIndex(&bConst.tphnam,atmnam);
                                        else if (aacode == 40) attype = getVecIndex(&bConst.olenam,atmnam);
                                        else if (aacode == 35) attype = getVecIndex(&bConst.pyrnam,atmnam);
                                        else if (aacode == 38) attype = getVecIndex(&bConst.pphnam,atmnam);
                                    }
                                    if (altloc != " ") {
                                        if (setatm == "    ") {
                                            setatm = atmnam;
                                            setloc = altloc;
                                            oldocc = occup;
                                            oldthm = therm;
                                        }
                                        if (!altset)
                                            if (setatm == atmnam) {
                                                if (occup > oldocc) {
                                                    setloc = altloc;
                                                    oldocc = occup;
                                                    oldthm = therm;
                                                }
                                            } else altset = true;
                                    }
                                    if (attype >= 0 && (altloc == setloc || altloc == " ")) {
                                        label.back() = keywd;
                                        if (abs(mcocc.at(attype).back()) > 0.0001 && mcocc.at(attype).back() < occup) {
                                            altmci.at(0).push_back(aacnt);
                                            altmci.at(1).push_back(attype);
                                            altmci.at(2).push_back(mcser.at(attype).back());
                                            altmcr.at(0).push_back(mcocc.at(attype).back());
                                            altmcr.at(1).push_back(mcthrm.at(attype).back());
                                            altmcc.push_back(mcalt.at(attype).back());
                                            for(unsigned int i=0; i<3; i++) altmca.at(i).push_back(mcaaat.at(attype).back().at(i));
                                        }
                                        mcser.at(attype).back() = atmnum;
                                        mcocc.at(attype).back() = occup;
                                        mcthrm.at(attype).back() = therm;
                                        mcalt.at(attype).back() = altloc;
                                        for(unsigned int i=0; i<3; i++) mcaaat.at(attype).back().at(i) = acoord.at(i);
                                        atomin.at(attype).back() = true;
                                        if (attype != bConst.calph-1) caonly = false;

                                    } else if (attype>=0) {
                                        altmci.at(0).push_back(aacnt);
                                        altmci.at(1).push_back(attype);
                                        altmci.at(2).push_back(atmnum);
                                        altmcr.at(0).push_back(occup);
                                        altmcr.at(1).push_back(therm);
                                        altmcc.push_back(altloc);
                                        for(unsigned int i=0; i<3; i++) altmca.at(i).push_back(acoord.at(i));
                                    }
                                    if (attype == -1 && (altloc == setloc || altloc == " ")) {
                                        if (aacode != 26) {
                                            if (bConst.aacids.at(aacode) == "ILE" && atmnam == " CD ") atmnam = " CD1 ";
                                            if ((bConst.aacids.at(aacode) == "GLN" || bConst.aacids.at(aacode) == "GLU") && atmnam == " OE ") atmnam = " OE1 ";
                                            if (atmnam == bConst.cbet) sctype = 0;
                                            else {
                                                sctype = getVecIndex(&bConst.scname.at(aacode), atmnam);
                                                unsigned int altposcod = aacode == 7? 1: (aacode == 13? 2: (aacode == 17? 3: 0));
                                                if (sctype == -1 && altposcod >0) sctype = getVecIndex(&bConst.altnam.at(altposcod-1), atmnam);
                                                if (sctype>=0) sctype++;
                                            }
                                            if (sctype != -1) {
                                                if (abs(scocc.at(sctype).back()) > 0.0001 && scocc.at(sctype).back() < occup) {
                                                    altsci.at(0).push_back(aacnt);
                                                    altsci.at(1).push_back(sctype);
                                                    altsci.at(2).push_back(scser.at(sctype).back());
                                                    altscr.at(0).push_back(scocc.at(sctype).back());
                                                    altscr.at(1).push_back(scthrm.at(sctype).back());
                                                    altscc.push_back(scalt.at(sctype).back());
                                                    for(unsigned int i=0; i<3; i++) altsca.at(i).push_back(scaaat.at(sctype).back().at(i));
                                                }
                                                scatin.at(sctype).back() = true;
                                                scocc.at(sctype).back() = occup;
                                                scthrm.at(sctype).back() = therm;
                                                scser.at(sctype).back() = atmnum;
                                                scalt.at(sctype).back() = altloc;
                                                for(unsigned int i=0; i<3; i++) scaaat.at(sctype).back().at(i) = acoord.at(i);
                                            }
                                        } else extData.push_back({brkrec.substr(0,66), aacnt});
                                    } else {
                                        if (atmnam == bConst.cbet) sctype = 0;
                                        else {
                                            if (bConst.aacids.at(aacode) == "ILE" && atmnam == " CD ") atmnam = " CD1 ";
                                            if ((bConst.aacids.at(aacode) == "GLN" || bConst.aacids.at(aacode) == "GLU") && atmnam == " OE ") atmnam = " OE1 ";
                                            sctype = getVecIndex(&bConst.scname.at(aacode), atmnam);
                                            unsigned int altposcod = aacode == 7? 1: (aacode == 13? 2: (aacode == 17? 3: 0));
                                            if (sctype == -1 && altposcod >0) sctype = getVecIndex(&bConst.altnam.at(altposcod-1), atmnam);
                                            if (sctype>=0) sctype++;
                                        }
                                        if (sctype>=0) {
                                            altsci.at(0).push_back(aacnt);
                                            altsci.at(1).push_back(sctype);
                                            altsci.at(2).push_back(atmnum);
                                            altscr.at(0).push_back(occup);
                                            altscr.at(1).push_back(therm);
                                            altscc.push_back(altloc);
                                            for(unsigned int i=0; i<3; i++) altsca.at(i).push_back(acoord.at(i));
                                        }
                                    }

                                }
                                startSeqFlag = false;
                            } else {
                                atmnam = brkrec.substr(12,4);
                                altloc = brkrec.substr(16,1);
                                altnum = brkrec.substr(21,6);
                                try {
                                    occup = brkrec.substr(54,6) == "      "? 1.0: stod(brkrec.substr(54,6));
                                } catch (invalid_argument& e) {throw "Broken PDB file! Invalid numerical data. Line number: " + to_string(cnt+1);}
                                if (altnum == seqnum) {                         
                                    if (aacode == -1) {  
                                    
                                        if (altloc != " " && extData.size() > 0) {
                                            extatm = extData.back().dat.substr(12,4);
                                            extnum = extData.back().dat.substr(21,6);
                                            try {
                                                extocc = extData.back().dat.substr(54,6) == "      "? 1.0: stod(extData.back().dat.substr(54,6));
                                            } catch (invalid_argument& e) {throw "Broken PDB file (Ext Data)! Invalid numerical data. Line number: " + to_string(cnt+1);}
                                            if (extatm == atmnam && extnum == altnum) {
                                                if (extocc >= occup) {
                                                    altData.push_back({brkrec.substr(0,66), aacnt});
                                                } else {
                                                    altData.push_back({extData.back().dat, aacnt});
                                                    extData.back().dat = brkrec.substr(0,66);
                                                }
                                            } else extData.push_back({brkrec.substr(0,66), aacnt});
                                        } else extData.push_back({brkrec.substr(0,66), aacnt});
                                    } else {specialCaseFlag = true;nextRec = false;}
                                } else {
                                    startSeqFlag = true;
                                    nextRec = false;
                                }
                            }
                        }
                    } else {
                        try {
                            nrec = stoi(brkrec.substr(50,5));
                        } catch (invalid_argument& e) {throw "Broken PDB file! Invalid numerical data. Line number: " + to_string(cnt+1);}
                    }
                } else {
                    extData.push_back({brkrec.substr(0,66), aacnt});
                  }
            } else throw "Broken PDB file! Too long record line. Line number: " + to_string(cnt+1);
        }
        cnt++;
    }
    brkfi.close();

// Debugging ============================

/*     ofstream fout;
    fout.open((string ("debug.out")).c_str(),ifstream::out | ifstream::trunc);

    try {
        for(unsigned int j=0; j<aacnt; j++)
            for(unsigned int i=0; i<bConst.mnchna; i++) {
                char buffer[82];
                sprintf(buffer, "%5d%2d%5d %s %8.3f%8.3f%8.3f", j+1,
                                        atomin.at(i).at(j)? 1:0, mcser.at(i).at(j), mcalt.at(i).at(j).c_str(),
                                        mcaaat.at(i).at(j).at(0), mcaaat.at(i).at(j).at(1), mcaaat.at(i).at(j).at(2));                         
                string sformatted(buffer);
                fout << sformatted << endl;
            }
        fout << endl << endl;
        for(unsigned int j=0; j<aacnt; j++)
            for(unsigned int i=0; i<bConst.sdchna; i++) {
                char buffer[82];
                sprintf(buffer, "%5d%2d%5d %s %8.3f%8.3f%8.3f", j+1,
                                        scatin.at(i).at(j)? 1:0, scser.at(i).at(j), scalt.at(i).at(j).c_str(),
                                        scaaat.at(i).at(j).at(0), scaaat.at(i).at(j).at(1), scaaat.at(i).at(j).at(2));                         
                string sformatted(buffer);
                fout << sformatted << endl;
            }

    } catch (ifstream::failure e) {throw "File output Error: (debug.out)";}

    fout.close();    */ 

    

// ======================================

// Next part of data processing
    if (aacnt>0) {

// Processing routine 'chnbrk'
// Search through looking for distant peptide bonds. if any over the limit are found 
// or if a c or n atom doesn't exist then a chain break is deemed to occur   
        unsigned int i;
// chnset() -> lambda expression for repeating code ===|        
        auto chnset = [&]() {
            chnsz.at(chnct-1) = ninchn;
            ninchn = 1;
            if (i<aacnt) {
                cout << "Chain break between " << i << " (" << seqbcd.at(i-1) <<") and "<< i+1 << " (" << seqbcd.at(i) << ")" << endl;
                chnct++;
                chstrt.at(i-1)[0] = '>';
                chstrt.at(i)[1]   = '<';
            }
        };
//=====================================================|        
        for(i=1; i < aacnt;i++) { 
            if (caonly) {
                if (getDist(&mcaaat.at(bConst.calph-1).at(i-1), &mcaaat.at(bConst.calph-1).at(i)) > bConst.cadbnd)
                    chnset(); 
                else {
                    ninchn++;
                    if (seqbcd.at(i-1).at(0) != seqbcd.at(i).at(0))
                        cout << "Warning! Chain letter change, but no chain break between (" << i-1 << ": " << seqbcd.at(i-1) <<") and ("<< i << ": " << seqbcd.at(i) << ")" << endl;
                }
            } else 
                if (atomin.at(bConst.carb-1).at(i-1) && atomin.at(bConst.nmain-1).at(i)) {
                    if (getDist(&mcaaat.at(bConst.carb-1).at(i-1), &mcaaat.at(bConst.nmain-1).at(i)) > bConst.pepbnd)
                        chnset(); 
                    else {
                        ninchn++;
                        if (seqbcd.at(i-1).at(0) != seqbcd.at(i).at(0)) {
                            cout << "Warning! Chain letter change, but no chain break between (" << i-1 << ": " << seqbcd.at(i-1) <<") and ("<< i << ": " << seqbcd.at(i) << ")" << endl;
                            chstrt.at(i-1)[0] = '>';
                            chstrt.at(i)[1]   = '<';
                        }    
                    }
                } else chnset();
        }
        chnset();

        if (caonly) {
// Processing routine 'calout'
            unsigned int nathet = 0, istart = 0;
            for(unsigned int ii=0; ii < extData.size(); ii++)
                if (extData.at(ii).res == 1) {
                    outStream.push_back(extData.at(ii).dat);
                    if (extData.at(ii).dat.substr(0,bConst.keylen) == bConst.atmkey || extData.at(ii).dat.substr(0,bConst.keylen) == bConst.htakey)
                        nathet++;
                } else if (extData.at(ii).res > 1) {
                        istart = ii; break;
                    }
            for(unsigned int nxres=0; nxres < aacnt; nxres++) {
                wkres = bConst.aacids.at(seqcod.at(nxres));
                if (atomin.at(1).at(nxres)) {
                    char buffer[82];
                    sprintf(buffer, "ATOM  %5d  CA  %s %s   %8.3f%8.3f%8.3f%6.2f%6.2f %s",
                        mcser.at(bConst.calph-1).at(nxres), wkres.c_str(), seqbcd.at(nxres).c_str(),
                        mcaaat.at(bConst.calph-1).at(nxres).at(0), mcaaat.at(bConst.calph-1).at(nxres).at(1), mcaaat.at(bConst.calph-1).at(nxres).at(2),
                        mcocc.at(bConst.calph-1).at(nxres), mcthrm.at(bConst.calph-1).at(nxres), chstrt.at(nxres).c_str());                         
                    string sformatted(buffer);
                    outStream.push_back(sformatted);
                    nathet++;
                }
                for(unsigned int ii=istart; ii < extData.size(); ii++)
                    if (extData.at(ii).res == nxres + 1) {
                        outStream.push_back(extData.at(ii).dat);
                        if (extData.at(ii).dat.substr(0,bConst.keylen) == bConst.atmkey || extData.at(ii).dat.substr(0,bConst.keylen) == bConst.htakey)
                            nathet++;
                    } else if (extData.at(ii).res > nxres + 1) {
                        istart = ii; break;
                    }
            }
            string sformatted = to_string(nathet);
            if (sformatted.length()<6) sformatted = (string ("MASTER").append(6-sformatted.length(), ' ')) + sformatted; 
            outStream.push_back(sformatted);
            if (nathet != nrec) 
                cout <<"WARNING: Number of atoms/hetatms in new file differs from old file :" <<endl
                     << "           new number = "<< nathet << endl << "           old number = "<< nrec << endl;                    
            outStream.push_back("END");
            cout << "**** C-alpha only file" << endl;
// End of routine 'calout'

        } else {
// Processing routine 'mkangl'
            for(unsigned int nxchn = 0;nxchn<chnct; chnst = chnst + chnsz.at(nxchn), nxchn++) 
                if (chnsz.at(nxchn) > 0)
                    for(unsigned int nxres = chnst; nxres < chnst + chnsz.at(nxchn);nxres++) 
                        if (atomin.at(bConst.calph-1).at(nxres) &&
                            atomin.at(bConst.nmain-1).at(nxres) && 
                            atomin.at(bConst.carb-1).at(nxres) && 
                            scatin.at(bConst.cbeta-1).at(nxres)) 
                                mcangs.at(nxres) = diHed(1,&mcaaat.at(bConst.calph-1).at(nxres),
                                                           &mcaaat.at(bConst.nmain-1).at(nxres),
                                                           &mcaaat.at(bConst.carb-1).at(nxres),
                                                           &scaaat.at(bConst.cbeta-1).at(nxres)); 
// Call Avevar
            double aveangs = 0.0, sdevangs = 0.0;
            unsigned int jj = 0;
            for(unsigned int j=0;j<chnst;j++)
                if (abs(mcangs.at(j)) < bConst.null) {aveangs+= mcangs.at(j);jj++;}
            if (jj>0) aveangs/=jj;
            for(unsigned int j=0;j<chnst;j++)
                if (abs(mcangs.at(j)) < bConst.null) sdevangs+= (mcangs.at(j) - aveangs) * (mcangs.at(j) - aveangs);
            sdevangs = jj > 1? sqrt(sdevangs/(jj-1)) : 0.0;   
            cout << "Average value of CA-N-C-CB angle is " << aveangs << endl << "Standard deviation is               " << sdevangs << endl;

// Processing routine 'mksang'
            unsigned int astart = 0, istart = 0, nathet = 0;
            int nswap = 0;
            string outerr = " ";
            vector<vector<int>> nold;

            for(unsigned int i=0; i<bConst.sdchna; i++) nold.push_back({});

            for(istart = 0;istart<extData.size() && extData.at(istart).res == 0;istart++) {
                outStream.push_back(extData.at(istart).dat);
                if(extData.at(istart).dat.substr(0,bConst.keylen) == bConst.atmkey || extData.at(istart).dat.substr(0,bConst.keylen) == bConst.htakey) nathet++;
            }
            if(istart>=extData.size()) istart = 0;

            for(astart = 0;astart<altData.size() && altData.at(astart).res == 0;astart++) {
                altStream.push_back(altData.at(astart).dat);
                if(altData.at(astart).dat.substr(0,bConst.keylen) == bConst.atmkey || altData.at(astart).dat.substr(0,bConst.keylen) == bConst.htakey) nathet++;
            }
            if(astart>=altData.size()) astart = 0;
// 5000 start
            for(unsigned int nxres=0; nxres<aacnt; nxres++) {
                unsigned int aacod = seqcod.at(nxres);
                wkres = bConst.aacids.at(aacod);
                string ld, mcatom = "    ", scatom = "    ";
                unsigned int nchi = bConst.natsdc.at(aacod);
                unsigned int extatm = bConst.extprc.at(aacod) > 1 ? 1: 0;
                for(unsigned int i=0; i<bConst.sdchna; i++) nold.at(i).push_back(-1);

// chainswap() -> lambda expression for repeating code ===|        
                auto chainswap = [&]() {
                    map<unsigned int, vector<unsigned int>> sswap  = 
                    {{3, {3,0}}, {4, {4,0}}, {5, {3,5}}, {8, {2,0}}, {11, {3,0}}, {17, {6,0}}, {19, {2,0}}, {21, {2,0}}, {24, {3,5}},
                    {30, {3,5}}, {34, {3,5}}, {38, {3,5}}, {40, {3,0}}, {41, {3,5}}, {44, {2,0}}, {45, {2,0}}, {46, {3,5}}, {49, {3,5}}};
                    if (sswap.count(aacod) > 0) {
                        for(unsigned int k=0;k<2;k++) {
                            int nside = sswap.at(aacod).at(k)-1;
                            if (nside >= 0) {
                                double tdum;
                                for(unsigned int i=0; i<3; i++) {
                                    tdum = scaaat.at(nside).at(nxres).at(i);
                                    scaaat.at(nside).at(nxres).at(i) = scaaat.at(nside+1).at(nxres).at(i);
                                    scaaat.at(nside+1).at(nxres).at(i) = tdum; 
                                }
                                tdum = scocc.at(nside).at(nxres);
                                scocc.at(nside).at(nxres) = scocc.at(nside+1).at(nxres);
                                scocc.at(nside+1).at(nxres) = tdum;
                                tdum = scthrm.at(nside).at(nxres);
                                scthrm.at(nside).at(nxres) = scthrm.at(nside+1).at(nxres);
                                scthrm.at(nside+1).at(nxres) = tdum;
                                nold.at(nside).at(nxres) = nside + 1;
                                nold.at(nside+1).at(nxres) = nside;
                            }
                        }
                    }
                    totswp.at(aacod)++;
                };
//========================================================|  

                if(nchi > 0) {
                    bool swatm = false;
                    vector<bool> athere;
                    vector<vector<double>> relats;
                    for(unsigned int i=0; i<bConst.maxchi+3; i++) {
                        athere.push_back(false);
                        relats.push_back({0.0, 0.0, 0.0});
                    }
                    athere.at(0) = atomin.at(bConst.nmain-1).at(nxres);
                    athere.at(1) = atomin.at(bConst.calph-1).at(nxres);
                    for(unsigned int i=0; i<3; i++) {
                        relats.at(0).at(i) = mcaaat.at(bConst.nmain-1).at(nxres).at(i);
                        relats.at(1).at(i) = mcaaat.at(bConst.calph-1).at(nxres).at(i);
                    }
                    for(unsigned int i=0; i < nchi + 1 + extatm; i++) {
                        athere.at(2+i) = scatin.at(i).at(nxres);
                        for(unsigned int j=0; j<3; j++)
                            relats.at(2+i).at(j) = scaaat.at(i).at(nxres).at(j);
                    }
                    if (bConst.extprc.at(aacod) == 5) {
                        if (scatin.at(nchi-1).at(nxres) && scatin.at(nchi+1).at(nxres)) {
                            double distg1 = getDist(&scaaat.at(nchi-1).at(nxres), &scaaat.at(nchi+1).at(nxres));
                            if (distg1 > bConst.ccdst && scatin.at(nchi).at(nxres)) {
                                double distg2 = getDist(&scaaat.at(nchi).at(nxres), &scaaat.at(nchi+1).at(nxres));
                                if (distg2 < bConst.ccdst) {
                                    swatm = true;
                                    chainswap();
                                }  
                            }
                        }
                        double dihang = diHed(1, &scaaat.at(0).at(nxres),
                                                 &scaaat.at(1).at(nxres),
                                                 &scaaat.at(2).at(nxres),
                                                 &mcaaat.at(bConst.calph-1).at(nxres));
                        if (dihang < 0.0) cout << "**** ILE residue has wrong chirality at CB " << nxres+1 << endl;                        
                    }
                    for(unsigned int nxchi=0; nxchi < nchi; nxchi++)
                        if (athere.at(nxchi) && athere.at(nxchi+1) && athere.at(nxchi+2) && athere.at(nxchi+3))
                            scangs.at(nxchi).at(nxres) = diHed(bConst.nmnang + nxchi + 1, &relats.at(nxchi),
                                                                                          &relats.at(nxchi+1),
                                                                                          &relats.at(nxchi+2),
                                                                                          &relats.at(nxchi+3));
                        else scangs.at(nxchi).at(nxres) = bConst.null;                                                                                         
                    if (bConst.extprc.at(aacod) > 1 && abs(scangs.at(nchi-1).at(nxres)) < bConst.null) {
                        if (athere.at(nchi-1) && athere.at(nchi) && athere.at(nchi+1) && athere.at(nchi+2+extatm))
                            scangs.at(nchi).at(nxres) = diHed(bConst.nmnang + nchi + extatm, &relats.at(nchi-1),
                                                                                             &relats.at(nchi),
                                                                                             &relats.at(nchi+1),
                                                                                             &relats.at(nchi+2+extatm));
                        else scangs.at(nchi).at(nxres) = bConst.null; 
                        if (abs(scangs.at(nchi).at(nxres)) < bConst.null) {
                            if (bConst.extprc.at(aacod) == 2 || bConst.extprc.at(aacod) == 3) {
                                double ang1 = scangs.at(nchi-1).at(nxres),
                                       ang2 = scangs.at(nchi).at(nxres);
                                if (ang1 < 0.0) ang1 += 360.0;
                                if (ang2 < 0.0) ang2 += 360.0;
                                if(abs(ang1-ang2) > 180.0 && (bConst.extprc.at(aacod) == 2 && ang2 > ang1 || bConst.extprc.at(aacod) == 3 && ang2 < ang1) ||
                                   abs(ang1-ang2) <=180.0 && (bConst.extprc.at(aacod) == 2 && ang2 <=ang1 || bConst.extprc.at(aacod) == 3 && ang2 >=ang1)) {
                                    swatm = true;
                                    chainswap();
                                    scangs.at(nchi-1).at(nxres) = scangs.at(nchi).at(nxres);
                                }
                            }
                            if (bConst.extprc.at(aacod) == 4 && abs(scangs.at(nchi).at(nxres)) < abs(scangs.at(nchi-1).at(nxres))) {
                                swatm = true;
                                chainswap();
                                scangs.at(nchi-1).at(nxres) = scangs.at(nchi).at(nxres);    
                            }
                            scangs.at(nchi).at(nxres) = bConst.null;
                        }
                    }
                    if (swatm) {
                        if (nswap == 0) cout << "* Side chain atoms swapped for residues:" << endl;
                        int wkswap = nswap % 7 + 1;
                        outerr += wkres + " " + seqbcd.at(nxres) + " ";
                        if (wkswap == 7) {
                            cout << "* " << outerr << endl;
                            outerr = " ";
                        }
                        nswap++;
                    }
                }
// 3000 start       
                for(unsigned int j=0; j < 5; j++) {
                    if (aacod > 43 && aacod <= 49) mcatom = bConst.batnam.at(j);
                    else if (aacod == 41) mcatom = bConst.tphnam.at(j);
                    else if (aacod == 40) mcatom = bConst.olenam.at(j);
                    else if (aacod == 38) mcatom = bConst.pphnam.at(j);
                    else if (aacod == 35) mcatom = bConst.pyrnam.at(j);
                    else mcatom = bConst.atname.at(j);
                    ld = abs(mcangs.at(nxres)) < bConst.null && (mcangs.at(nxres) < 23.0 ||  mcangs.at(nxres) > 45.0)? "*" : " ";
                    if (atomin.at(j).at(nxres)) {
                        char buffer[82];
                        sprintf(buffer, "%s%5d %s%s%s %s   %8.3f%8.3f%8.3f%6.2f%6.2f %s    %s",
                            label.at(nxres).c_str(), mcser.at(j).at(nxres), mcatom.c_str(), mcalt.at(j).at(nxres).c_str(), wkres.c_str(), seqbcd.at(nxres).c_str(),
                            mcaaat.at(j).at(nxres).at(0), mcaaat.at(j).at(nxres).at(1), mcaaat.at(j).at(nxres).at(2),
                            mcocc.at(j).at(nxres), mcthrm.at(j).at(nxres), chstrt.at(nxres).c_str(), ld.c_str());                         
                        string sformatted(buffer);
                        outStream.push_back(sformatted);
                        nathet++;
                    }
                    if (altmcc.size() >0 && mcalt.at(j).at(nxres) != " ")
                        for(unsigned int i=0; i < altmcc.size(); i++)
                            if (altmci.at(1).at(i) == j && altmci.at(0).at(i) == nxres+1) {
                                char buffer[82];
                                sprintf(buffer, "%s%5d %s%s%s %s   %8.3f%8.3f%8.3f%6.2f%6.2f",
                                        label.at(nxres).c_str(), altmci.at(2).at(i), mcatom.c_str(), altmcc.at(i).c_str(), wkres.c_str(), seqbcd.at(nxres).c_str(),
                                        altmca.at(0).at(i), altmca.at(1).at(i), altmca.at(2).at(i), 
                                        altmcr.at(0).at(i), altmcr.at(1).at(i));                         
                                string sformatted(buffer);
                                altStream.push_back(sformatted);
// Most possible BUG in Fortran code
//                                char buffer1[82];
//                                string altlab = label.at(nxres).substr(0,2) + "ALT "; 
//                                sprintf(buffer1, "%s%5d %s%s%s %s   %8.3f%8.3f%8.3f%6.2f%6.2f",
//                                        altlab.c_str(), altsci.at(2).at(j), scatom.c_str(), altscc.at(j).c_str(), wkres.c_str(), seqbcd.at(nxres).c_str(),
//                                        altsca.at(0).at(j), altsca.at(1).at(j), altsca.at(2).at(j), 
//                                        altscr.at(0).at(j), altscr.at(1).at(j));                         
//                                sprintf(buffer1, "%s%5d %s%s%s %s   %8.3f%8.3f%8.3f%6.2f%6.2f",
//                                       altlab.c_str(), altmci.at(2).at(i), mcatom.c_str(), altmcc.at(i).c_str(), wkres.c_str(), seqbcd.at(nxres).c_str(),
//                                        altmca.at(0).at(i), altmca.at(1).at(i), altmca.at(2).at(i), 
//                                        altmcr.at(0).at(i), altmcr.at(1).at(i));                         
//                                string sformatted1(buffer1);
                                sformatted = label.at(nxres).substr(0,2) + "ALT " + sformatted.substr(6);    
                                outStream.push_back(sformatted);                                
                                if (label.at(nxres).substr(0,2) == "AT") natalt++;
                                else nhealt++; 
                                nathet++;
                            }
                }
// 3000 end                                        
// 4000 start
                for(unsigned int i=0; i < bConst.nscats.at(aacod); i++)
                    if(scatin.at(i).at(nxres)) {
                        scatom = i==0? bConst.cbet: bConst.scname.at(aacod).at(i-1);
                        string scold = nold.at(i).at(nxres) < 1 ? "    " : bConst.scname.at(aacod).at(nold.at(i).at(nxres)-1);
                        char buffer[82];
                        sprintf(buffer, "%s%5d %s%s%s %s   %8.3f%8.3f%8.3f%6.2f%6.2f %s%s%s",
                            label.at(nxres).c_str(), scser.at(i).at(nxres), scatom.c_str(), scalt.at(i).at(nxres).c_str(), wkres.c_str(), seqbcd.at(nxres).c_str(),
                            scaaat.at(i).at(nxres).at(0), scaaat.at(i).at(nxres).at(1), scaaat.at(i).at(nxres).at(2),
                            scocc.at(i).at(nxres), scthrm.at(i).at(nxres), chstrt.at(nxres).c_str(), scold.c_str(), ld.c_str());                         
                        string sformatted(buffer);
                        outStream.push_back(sformatted);
                        nathet++;
                        if (altscc.size() >0 && scalt.at(i).at(nxres) != " ")
                            for(unsigned int j=0; j < altscc.size(); j++)
                                if (altsci.at(1).at(j) == i && altsci.at(0).at(j) == nxres+1) {
                                    char buffer[82];
                                    sprintf(buffer, "%s%5d %s%s%s %s   %8.3f%8.3f%8.3f%6.2f%6.2f",
                                            label.at(nxres).c_str(), altsci.at(2).at(j), scatom.c_str(), altscc.at(j).c_str(), wkres.c_str(), seqbcd.at(nxres).c_str(),
                                            altsca.at(0).at(j), altsca.at(1).at(j), altsca.at(2).at(j), 
                                            altscr.at(0).at(j), altscr.at(1).at(j));                         
                                    string sformatted(buffer);
                                    altStream.push_back(sformatted);
                                    sformatted = label.at(nxres).substr(0,2) + "ALT " + sformatted.substr(6);    
                                    outStream.push_back(sformatted);                                
                                    if (label.at(nxres).substr(0,2) == "AT") natalt++;
                                    else nhealt++; 
                                    nathet++;
                                }                        
                    }
// 4000 end
                for(unsigned int ii=istart; ii < extData.size(); ii++)
                    if (extData.at(ii).res == nxres + 1) {
                        outStream.push_back(extData.at(ii).dat);
                        if (extData.at(ii).dat.substr(0,bConst.keylen) == bConst.atmkey || extData.at(ii).dat.substr(0,bConst.keylen) == bConst.htakey)
                            nathet++;
                    } else if (extData.at(ii).res > nxres + 1) {
                        istart = ii; break;
                    }
                for(unsigned int ii=astart; ii < altData.size(); ii++)
                    if (altData.at(ii).res == nxres + 1) {
                        altStream.push_back(altData.at(ii).dat);
                        if (altData.at(ii).dat.substr(0,bConst.keylen) == bConst.atmkey || altData.at(ii).dat.substr(0,bConst.keylen) == bConst.htakey)
                            nathet++;
                    } else if (altData.at(ii).res > nxres + 1) {
                        astart = ii; break;
                    }
            } 
// 5000 end
            string sformatted = to_string(nathet);
            if (sformatted.length()<6) sformatted = (string ("MASTER").append(6-sformatted.length(), ' ')) + sformatted; 
            outStream.push_back(sformatted);
            if (nathet != nrec) 
                cout <<"WARNING: Number of atoms/hetatms in new file differs from old file :" <<endl
                     << "           new number = "<< nathet << endl << "           old number = "<< nrec << endl;
            outStream.push_back("END");
            if (nswap%7 != 0) cout << "* " << outerr << endl;
            if (altmcc.size()>0) cout << "Number of alternate mainchain atoms is: "<< altmcc.size() << endl;
            if (altscc.size()>0) cout << "Number of alternate sidechain atoms is: "<< altscc.size() << endl;
// End of routine 'mksang'        
        }
        cout <<"Total number of residues changed and total number of residues:" << endl;
        for(unsigned k=0; k<bConst.sideaa;k++) {
            wkres = bConst.aacids.at(k);
            if (restot.at(k) > 0)
                cout << "     " << wkres << " " << totswp.at(k) << " -> " << restot.at(k) << endl;
        }
        if (natalt>0) cout << "* Number of alternate atoms (labelled ATALT)    " << natalt << endl; 
        if (nhealt>0) cout << "* Number of alternate hetatoms (labelled HEALT) " << nhealt << endl;

        
        cout << endl << "Number of atom processed: "<< aacnt << endl;

    }
    else throw string ("**** Error while reading file. No atoms encountered."); 

}


//************************ Processed data output *****************************************
void ProcessData::DataOutput(string outFileName, vector<string> *out_stream) {
    
    ofstream fout;
    fout.open(outFileName.c_str(),ifstream::out | ifstream::trunc);

    try {
        for(unsigned int i=0; i<out_stream->size();i++) {
            string paddedStr = out_stream->at(i);
            if (paddedStr.length()<80) paddedStr.append(80-paddedStr.length(), ' ');
            fout << paddedStr << endl;
        }    
    } catch (ifstream::failure e) {throw "File output Error: (" + outFileName + ") -> " + e.what();}

    fout.close();
} 


//////////// Analytical and service functions ///////////////////////////////
// Euclidian distance between pair of vectors
double ProcessData::getDist(const vector<double> *avec, const vector<double> *bvec) {
    double output = 0.0;
    unsigned int vsize = min(avec->size(),bvec->size());
    for(unsigned int ii=0;ii<vsize;ii++) 
        output += (avec->at(ii)-bvec->at(ii)) * (avec->at(ii)-bvec->at(ii));
    return sqrt(output);    
};

// Dihedral angle between the 4 vectors given
double ProcessData::diHed(const int angnum, const vector<double> *va, const vector<double> *vb, const vector<double> *vc, const vector<double> *vd) {
    const unsigned int vsize = min(min(va->size(),vb->size()),min(vc->size(),vd->size()));
    if (vsize<2) return 0.0;
    double cosang, output = 0.0;
    const vector<vector<double>> dihvec = {*va, *vb, *vc, *vd};
    vector<double> vecdst = {0.0, 0.0, 0.0};
    vector<vector<double>> codist = {{}, {}, {}};
    vector<vector<double>> dotprd = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    for(unsigned int i=0;i<3;i++) {
        vecdst.at(i) = getDist(&dihvec.at(i), &dihvec.at(i+1));
        for(unsigned int j=0;j<vsize;j++) 
            codist.at(i).push_back(dihvec.at(i+1).at(j) - dihvec.at(i).at(j));
    }
    for(unsigned int i=0;i<2;i++)
        for(unsigned int j=i+1;j<3;j++) {
            for(unsigned int k=0;k<vsize;k++)
                dotprd.at(i).at(j) += (codist.at(i).at(k) * codist.at(j).at(k));
            dotprd.at(i).at(j) = (abs(vecdst.at(i))<0.0001 || abs(vecdst.at(j))<0.0001) ? 1.0: dotprd.at(i).at(j) /= (vecdst.at(i) * vecdst.at(j));
        }
    cosang = (angnum > bConst.ndihed && angnum <= bConst.nmnang)? dotprd.at(0).at(2):
                 ((abs(dotprd.at(0).at(1)-1.0)<0.0001 || abs(dotprd.at(1).at(2)-1.0)<0.0001)? 1.0: 
                    (dotprd.at(0).at(1)*dotprd.at(1).at(2) - dotprd.at(0).at(2)) / (sqrt(1.0 - dotprd.at(0).at(1)*dotprd.at(0).at(1)) * sqrt(1.0 - dotprd.at(1).at(2)*dotprd.at(1).at(2))));
    if (abs(cosang) > 1.0) cosang = cosang < 0.0 ? -1.0: 1.0;
    
    output = acos(cosang) * bConst.radian;

    if (angnum <= bConst.ndihed || angnum > bConst.nmnang) {
        double detant = ((((((codist.at(0).at(0) * codist.at(1).at(1)) * codist.at(2).at(2)) - 
                            ((codist.at(0).at(0) * codist.at(1).at(2)) * codist.at(2).at(1))) + 
                            ((codist.at(0).at(1) * codist.at(1).at(2)) * codist.at(2).at(0))) - 
                            ((codist.at(0).at(1) * codist.at(1).at(0)) * codist.at(2).at(2))) + 
                            ((codist.at(0).at(2) * codist.at(1).at(0)) * codist.at(2).at(1))) -
                            ((codist.at(0).at(2) * codist.at(1).at(1)) * codist.at(2).at(0));
        if (detant < 0.0) output = -output;                    
    }
    return output;
};

int ProcessData::getVecIndex(const vector<string> *invec, const string s2find) {
    int idx=-1;
    if (invec->size()>0) {
        auto it = find(invec->begin(), invec->end(), s2find);
        idx = (it == invec->end()? -1: distance(invec->begin(), it));
    }
    return idx;
};
/////////////////////////////////////////////////////////////////////////    
