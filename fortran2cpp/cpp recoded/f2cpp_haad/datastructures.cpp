#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>
#include "datastructures.h"

using namespace std;


//************************ PDB data loading *****************************************
int ProcessData::DataLoad(char *fn) {
    // Creation of ifstream class object to read the file
    ifstream fin;
    string sline;
    // Open the file 
    fin.open(fn);
    if (!fin) throw "No such file or directory! (" + (string)fn + ")";
    unsigned int cnt = 0;
    LoadAtomInfo data2load;
    // Read the file line by line
    while (getline(fin, sline) && (sline.substr(0,3)!="END")  && (sline.substr(0,3)!="TER")) {
        if ((sline.substr(0,4) == "ATOM") && ((sline.substr(16,1) == "A") || (sline.substr(16,1) == " "))) {
            data2load.atomName = sline.substr(12,4);
            data2load.residue = sline.substr(17,3);
            try {
                data2load.seqNumber = stoi(sline.substr(22,4));
                data2load.x = stod(sline.substr(30,8));
                data2load.y = stod(sline.substr(38,8));
                data2load.z = stod(sline.substr(46,8));
            } catch (invalid_argument& e) {throw "Broken PDB file! Atom number: " + to_string(cnt+1);};
            
            data2load.kfp = sline.substr(0,30);
            if(resSeqIn.size()>0 && data2load.seqNumber!=resSeqIn.back()) resd2.push_back(loadedData.back().residue);
            if(resSeqIn.size()==0 || resSeqIn.back()!=data2load.seqNumber) resSeqIn.push_back(data2load.seqNumber);
            loadedData.push_back(data2load);
            cnt++;
        }    
    }
    if (cnt>0) resd2.push_back(loadedData.back().residue);
     // Close the file
    fin.close();
    if (cnt<10) throw  "Too short atom data file! Number of lines: " + to_string(cnt);

    return cnt;
} 

//************************ PDB processed data output *****************************************
void ProcessData::DataOutput(const char *fn, unsigned int outAtomsNumber) {
    // Creation of ofstream class object to write data into the file
    ofstream fin;
    // Open the file 
    fin.open(fn,ifstream::out | ifstream::trunc);
    try {
        for(unsigned int i=0; i<outAtomsNumber;i++) {
            try {
                char ch_int[70];
                if(assrd.at(i)>0) {
                    sprintf(ch_int, "ATOM      0 %s %s  %4d    %8.3f%8.3f%8.3f%3d", atomty.at(i).c_str(), resdOut.at(i).c_str(), resSeqOut.at(i), xx.at(i), yy.at(i), zz.at(i),assrd.at(i)); 
                } else {
                    sprintf(ch_int, "%30s%8.3f%8.3f%8.3f", loadedData.at(tkfp.at(i)).kfp.c_str(), xx.at(i), yy.at(i), zz.at(i)); 
                }
                fin << ch_int << endl;
            } catch (out_of_range& e) {throw "Index is Out of Range while writing an output file (assrd)#:" + to_string(i+1);}
        }
        fin << "TER" << endl;
    } catch (ifstream::failure e) {throw "File output Error: (" + (string)fn + ") -> " + e.what();}
    // Close the file
    fin.close();
} 

//************************ Get Residue Series and prepare data for an Addition ***************************************** 
int ProcessData::GetResidueSeries(unsigned short numRes) {
    unsigned int iAtom = 0; 
    unsigned int iResin = 0;
    map<string, unsigned short> otherResidues  = {{"HSD", 16}, {"HSE", 16}, {"HSP", 16}, {"HIE", 16}, {"MSE", 8}};
    vector<AcidInfo> amAcids = servData->aminoAcids;

    for (unsigned short ir=0; ir < numRes; ir++) {
        const string resN = resd2[ir];
        auto inList = [resN](AcidInfo el) { return el.name == resN; };
        if(auto iter = find_if(amAcids.begin(), amAcids.end(), inList); iter !=amAcids.end()) {
            iResin = distance(amAcids.begin(), iter);
        } else {
            if (auto otherIter = otherResidues.find(resN); otherIter !=otherResidues.end()) {
                iResin = (*otherIter).second;
            } else {
                iResin = 0; //GLY
                cout << "WARNING on reside abnormal: (" << ir+1 << ", " << resN << ")" << endl;
            }
        }
        ca.push_back({ iAtom+1, iResin});
        for(unsigned short j=0; j < amAcids[iResin].heavyAtoms;j++) {
            iAtom++;
            atomty.push_back(amAcids[iResin].atomNames[j]);
            atc.push_back(amAcids[iResin].atomConnections[j]);
            resdOut.push_back(amAcids[iResin].name);
            ress.push_back(ir);
            resSeqOut.push_back(resSeqIn.at(ir));
        }

        if (ir == 0) {
            atomty.push_back("1H  ");
            atomty.push_back("2H  ");
            atomty.push_back((iResin==7)? " HA " : "3H  "); // PRO -> HA : others -> 3H
            atc.push_back(20);
            atc.push_back(20);
            atc.push_back((iResin==7)? 19 : 20); // PRO -> 19 : others -> 20
            atc[0] = ((iResin==7)? 33 : 12); // PRO -> 33 : others -> 12
            resdOut.push_back(amAcids[iResin].name);
            resdOut.push_back(amAcids[iResin].name);
            resdOut.push_back(amAcids[iResin].name);
            ress.push_back(ir);
            ress.push_back(ir);
            ress.push_back(ir);
            resSeqOut.push_back(resSeqIn.at(ir));
            resSeqOut.push_back(resSeqIn.at(ir));
            resSeqOut.push_back(resSeqIn.at(ir));
            iAtom+=3;
        } else if (ir == nRes-1) {
            atc[ca[ir][0]+1] = 5; // change C to CC
            atomty.push_back(" OXT");
            atc.push_back(27);
            resdOut.push_back(amAcids[iResin].name);
            ress.push_back(ir);
            resSeqOut.push_back(resSeqIn.at(ir));
            iAtom++;
        }
        for(unsigned short j = amAcids[iResin].heavyAtoms + (( ir == 0 )? 1 : 0); j < amAcids[iResin].totalAtoms;j++) {
            iAtom++;
            atomty.push_back(amAcids[iResin].atomNames[j]);
            atc.push_back(amAcids[iResin].atomConnections[j]);
            resdOut.push_back(amAcids[iResin].name);
            ress.push_back(ir);
            resSeqOut.push_back(resSeqIn.at(ir));
        } 
    }
  
    // Assign the bond network
    // 0 - term
    iResin = ca[0][1];
    unsigned int itom;
    for (itom=0; itom < amAcids[iResin].heavyAtoms;itom++) {
        bonda.push_back({-1,-1,-1,-1,0});
        for (unsigned short j=0;j<4;j++) {
            if (amAcids[iResin].bdc[itom][j] != 0) {
                bonda[itom][j] = itom + amAcids[iResin].bdc[itom][j] + ((amAcids[iResin].bdc[itom][j] > 3) ? 2 : 0); 
            }
        }
    }
    bonda[0][2] = amAcids[iResin].heavyAtoms + 1 - 1;
    bonda[0][3] = amAcids[iResin].heavyAtoms + 2 - 1 ;
    bonda.push_back({-1,-1,-1,-1,0});
    bonda[itom][0] = 0;
    itom++;
    bonda.push_back({-1,-1,-1,-1,0});
    bonda[itom][0] = 0;
    for (unsigned int i=amAcids[iResin].heavyAtoms; i < amAcids[iResin].totalAtoms;i++) {
        itom++;
        bonda.push_back({-1,-1,-1,-1,0});
        bonda[itom][0] = itom + amAcids[iResin].bdc[i][0] - 2;
    }
    bonda[2][1] = ca[1][0] - 1;
    if (iResin == 7) bonda[0][1] = amAcids[iResin].heavyAtoms - 1; // PRO

    // {1..numRes-2} - term
    for (unsigned short k=1; k < numRes - 1; k++) {
        iResin = ca[k][1];
        for (unsigned int i=0; i < amAcids[iResin].totalAtoms;i++) {
            itom++;
            bonda.push_back({-1,-1,-1,-1,0});
            for (unsigned short j=0;j<4;j++) {
                if (amAcids[iResin].bdc[i][j] != 0) {
                    bonda[itom][j] = itom + amAcids[iResin].bdc[i][j]; 
                }
            }
        }
        bonda[ca[k][0]-1][2] = ca[k-1][0] + 1;
        bonda[ca[k][0]+1][1] = ca[k+1][0] - 1;
    }    

    // numRes-1 - term
    iResin = ca[numRes-1][1];
    for (unsigned int i=0; i < amAcids[iResin].heavyAtoms;i++) {
        itom++;
        bonda.push_back({-1,-1,-1,-1,0});
        for (unsigned short j=0;j<4;j++) {
            if (amAcids[iResin].bdc[i][j] != 0) {
                bonda[itom][j] = itom + amAcids[iResin].bdc[i][j] + ((amAcids[iResin].bdc[i][j] > 3) ? 1 : 0); 
            }
        }
    }
    itom++;
    bonda.push_back({-1,-1,-1,-1,0});
    bonda[itom][0] = ca[numRes-1][0] + 1;
    bonda[ca[numRes-1][0] + 1][1] = itom;

    for (unsigned int i=amAcids[iResin].heavyAtoms; i < amAcids[iResin].totalAtoms;i++) {
        itom++;
        bonda.push_back({-1,-1,-1,-1,0});
        bonda[itom][0] = itom + amAcids[iResin].bdc[i][0] - 1;
    }
    bonda[ca[numRes-1][0] - 1][2] = ca[numRes-2][0] + 1;
    if (iResin == 7) bonda[ca[numRes-1][0] - 1][1] -= 1;

    // Bond vector assign finish
    
    // Assign the bondlength 
    for (unsigned int i = 0; i < iAtom; i++) {
        assrd.push_back(0); // additional record flag initialisation
        bond.push_back({0.0,0.0,0.0,0.0});
        for (unsigned short j=0;j<4;j++) {
            if (bonda[i][j] >= 0) {
                bond[i][j] = (atc[i]<=atc[bonda[i][j]]) ? servData->bd0[atc[i]-1][atc[bonda[i][j]]-1]: servData->bd0[atc[bonda[i][j]]-1][atc[i]-1];
                if (bond[i][j] < 0.5) cout << "Low bond length: [" << i << "] -> " << bond[i][j] << endl;
            }
        }
    }
    // Final check
    for(unsigned int i=0;i<iAtom;i++) {
        int ip2 = -1;
        int icc = bonda[i][0];
        if (icc>=0) {
            for(unsigned int j=0;(j<4) && (bonda[icc][j] >=0);j++) {
                if ((i!=bonda[icc][j]) && (atomty[bonda[icc][j]].substr(1,1)  != "H")) {
                    ip2 = bonda[icc][j];
                }
            }
            if (ip2 >=0) {
                unsigned int i1 = (atc[i] > atc[ip2]) ? atc[ip2]-1 : atc[i]-1;
                unsigned int i3 = (atc[i] > atc[ip2]) ? atc[i]-1 : atc[ip2]-1;
                if (servData->ag0[i1][atc[icc]-1][i3] < 1.0) cout << "Warning! Abnormal angle on the output atom number : " << i+1 << endl;
            }
        }    
//        else cout << "Warning! Abnormal bond network on the output atom number: " << i+1 << endl;    

    }

    // Prefill x,y,z and number of heavy atoms bonded to each atom
    for(unsigned int i=0;i<iAtom;i++) {
        xx.push_back(999.0);
        yy.push_back(999.0);
        zz.push_back(999.0);
        tkfp.push_back(0);
        unsigned int j;
        for(j=0; j<inAtoms && ((atomty[i] != loadedData[j].atomName) || (resSeqOut[i] != loadedData[j].seqNumber)); j++);
        if (j<inAtoms) {
            xx[i] = loadedData[j].x;
            yy[i] = loadedData[j].y;
            zz[i] = loadedData[j].z;
            tkfp.at(i) = j;
        }

        for(unsigned int k=0;k<4;k++) 
            if (bonda[i][k]>=0 && atomty[bonda[i][k]].substr(1,1) != "H") bonda[i][4]++;

        assrd.at(i) += xx.at(i)>=990.0;    
    }
    return iAtom;
}

//======================================================================================================================    
//Function to add those misssing heavy atoms except CB, CG and CD
void ProcessData::Heavy_Add(unsigned int outAtomsNumber) {

    for(unsigned int i=0; i < outAtomsNumber; i++) {
        try { 
            if (atomty.at(i).substr(1,1) != "H" && xx.at(i)>=990.0) {
                double ax, ay, az, bd;
                int icc=-1, ip2=-1, ip3=-1; 
                
                for (unsigned int k=0; k<4 && bonda.at(i)[k]>=0 && icc==-1; k++) {
                    if (xx[bonda.at(i)[k]] <=990) {
                        icc = bonda.at(i)[k];
                        bd = bond.at(i)[k];
                    }    
                }
                if(icc>=0) { // This is redundant additional check, because icc has to be definitely found following previous logic
                    // Origin point
                    ax = xx.at(icc);
                    ay = yy.at(icc);
                    az = zz.at(icc);
                    for (unsigned int j=0; j<4 && bonda.at(icc)[j]>=0; j++) 
                        if(xx[bonda.at(icc)[j]] <= 990.0 && i!=bonda.at(icc)[j] && atomty[bonda.at(icc)[j]].substr(1,1) != "H") 
                            ip2=bonda.at(icc)[j];
    
                    for (unsigned int j=0; j<4 && ip3==-1; j++) 
                        if(bonda.at(ip2)[j]>=0 &&  icc!=bonda.at(ip2)[j] && atomty[bonda.at(ip2)[j]].substr(1,1) != "H" && xx[bonda.at(ip2)[j]] <= 990.0)
                            ip3=bonda.at(ip2)[j];
                    
                    // Now everything is ready for the addition proccess

                    //treat heavy atoms missing !ip3-ip2-icc-i (heavy atom)
                    //C-term ALA-Cb; ILE Cd; LEU Cd1,Cd2;MET Ce; Val CG1; THR CG2  sp3
                    //O-term O for all, ASN OD1; ASP OD1, OD2; GLN OE1; GLU OE1,OE2; TYR OH;sp2 || SER OG; THR OG1, OXT (sp3) 
                    //S-term CYS SG sp3
                    //N-term N-term LYS NZ sp3, ARG NH1,NH2; ASN ND2; GLN NE2 sp2
                    double x,y,z,r1,xi,yi,zi,fbd;
                    int ipp1 = ip3; // Just in case of abnormal record (assrd[i]=2) we'll need this
                    if (bonda.at(i)[4]==1) {
                        if((std::set<string> {" OD2"," OE2"," ND2"," NE2"," NH2"," OH "," OXT"}).count(atomty.at(i)) >0){
                            ip2 = icc-1; 
                            ip3 = (atomty.at(i) == " OH ")?  icc-2 : icc+1; 
                            x=0.5*(xx.at(ip2)+xx.at(ip3));
                            y=0.5*(yy.at(ip2)+yy.at(ip3));
                            z=0.5*(zz.at(ip2)+zz.at(ip3));
                            r1 = sqrt((x-ax)*(x-ax)+(y-ay)*(y-ay)+(z-az)*(z-az));
                            if (r1 > servData->minDist) {
                                xx.at(i) = ax + bd*(ax-x)/r1;
                                yy.at(i) = ay + bd*(ay-y)/r1;
                                zz.at(i) = az + bd*(az-z)/r1;
                            } else {
                                assrd.at(i) = 2; // Warning! Abnormal record
                                xx.at(i) = ax+xx.at(ip2)-xx.at(ipp1);
                                yy.at(i) = ay+yy.at(ip2)-yy.at(ipp1);
                                zz.at(i) = az+zz.at(ip2)-zz.at(ipp1);
                            } 
                        } else if ((atomty.at(i) == " CD2") && (resdOut.at(i) == "LEU")) {
                            xi = ax - xx.at(icc-1);
                            yi = ay - yy.at(icc-1);
                            zi = az - zz.at(icc-1);
                            x = xx.at(i-1);
                            y = yy.at(i-1);
                            z = zz.at(i-1);
                            Rotate_Matrix(x, y, z, ax, ay, az, xi, yi, zi, 2.18166);
                            xx.at(i) = x;
                            yy.at(i) = y;
                            zz.at(i) = z;
                        } else if (atomty.at(i) == " CG2") {
                            xi = ax - xx.at(icc-3);
                            yi = ay - yy.at(icc-3);
                            zi = az - zz.at(icc-3);
                            x = xx.at(i-1);
                            y = yy.at(i-1);
                            z = zz.at(i-1);
                            Rotate_Matrix(x, y, z, ax, ay, az, xi, yi, zi, 2.18166);
                            xx.at(i) = x;
                            yy.at(i) = y;
                            zz.at(i) = z;
                        } else if ((atomty.at(i) == " CD1") && (resdOut.at(i) == "ILE") || (atomty.at(i) == " CE ") && (resdOut.at(i) == "MET")) {
                            xi = xx.at(icc-1) - 0.5*(xx.at(i-1) + xx.at(icc-4));
                            yi = yy.at(icc-1) - 0.5*(yy.at(i-1) + yy.at(icc-4));
                            zi = zz.at(icc-1) - 0.5*(zz.at(i-1) + zz.at(icc-4));
                            x = xi - ax;
                            y = yi - ay;
                            z = zi - az;
                            r1 = sqrt(x*x + y*y + z*z);
                            if (r1<0.1) r1 = 0.1;
                            xx.at(i) = ax + bd*x/r1;
                            yy.at(i) = ay + bd*y/r1;
                            zz.at(i) = az + bd*z/r1;
                        } else { //!!sp3 and sp2 first atom  ?????????
                            if((std::set<string> {" CG1"," SG "," OG1"," OG "}).count(atomty.at(i)) >0) ip3 = ip2 +1;
                            x = xx.at(ip2) - xx.at(ip3);
                            y = yy.at(ip2) - yy.at(ip3);
                            z = zz.at(ip2) - zz.at(ip3);
                            r1 = sqrt(x*x + y*y + z*z);
                            if (r1 > servData->minDist) {
                                xx.at(i) = ax + bd*x/r1;
                                yy.at(i) = ay + bd*y/r1;
                                zz.at(i) = az + bd*z/r1;
                            } else {
                                assrd.at(i) = 2; // Warning! Abnormal record
                                xx.at(i) = 2.0*ax - xx.at(ip2);
                                yy.at(i) = 2.0*ay - yy.at(ip2);
                                zz.at(i) = 2.0*az - zz.at(ip2);
                            }
                        }
                    } else if(bonda.at(i)[4] >= 2) { //!!treats more atoms missing
                        if (atomty.at(i) == " CG ") ip3 = ip2 + 1;
                        x = xx.at(ip2) - xx.at(ip3);
                        y = yy.at(ip2) - yy.at(ip3);
                        z = zz.at(ip2) - zz.at(ip3);
                        r1 = sqrt(x*x + y*y + z*z);
                        if (r1 > servData->minDist) {
                            xx.at(i) = ax + bd*x/r1;
                            yy.at(i) = ay + bd*y/r1;
                            zz.at(i) = az + bd*z/r1;
                        } else {
                            assrd.at(i) = 2; // Warning! Abnormal record
                            xx.at(i) = ax + x;
                            yy.at(i) = ay + y;
                            zz.at(i) = az + z;
                        }

                        if (ca[ress.at(i)][1] == 7 || ca[ress.at(i)][1]>=16) {
                            if (atomty.at(i) == " CG ") {
                                if (resdOut.at(i) == "PRO") { //!!this part treats the ring PRO CG
                                    x = xx.at(i) - xx.at(ip2-1);
                                    y = yy.at(i) - yy.at(ip2-1);
                                    z = zz.at(i) - zz.at(ip2-1);
                                    r1 = sqrt(x*x + y*y + z*z);
                                    if (r1 < servData->minDist) r1 = servData->minDist;
                                    xx.at(i) = xx.at(ip2-1) + 2.22*x/r1;
                                    yy.at(i) = yy.at(ip2-1) + 2.22*y/r1;
                                    zz.at(i) = zz.at(ip2-1) + 2.22*z/r1;

                                    x = 0.5*(xx.at(ip2-1) + xx.at(i) - xx.at(ip2) - xx.at(i-1));
                                    y = 0.5*(yy.at(ip2-1) + yy.at(i) - yy.at(ip2) - yy.at(i-1));
                                    z = 0.5*(zz.at(ip2-1) + zz.at(i) - zz.at(ip2) - zz.at(i-1));
                                    
                                    assrd.at(i+1) = 1;
                                    xx.at(i+1) = 0.5*(xx.at(ip2) + xx.at(i-1)) + 1.618*x;
                                    yy.at(i+1) = 0.5*(yy.at(ip2) + yy.at(i-1)) + 1.618*y;
                                    zz.at(i+1) = 0.5*(zz.at(ip2) + zz.at(i-1)) + 1.618*z;
                                    
                                } else {
                                    double bx = xx.at(ip2+1) - xx.at(ip2-1);
                                    double by = yy.at(ip2+1) - yy.at(ip2-1);
                                    double bz = zz.at(ip2+1) - zz.at(ip2-1);
                                    double cx = xx.at(i) - xx.at(i-1);
                                    double cy = yy.at(i) - yy.at(i-1);
                                    double cz = zz.at(i) - zz.at(i-1);
                                    double ctbcx = by*cz - bz*cy;
                                    double ctbcy = bz*cx - bx*cz;
                                    double ctbcz = bx*cy - by*cx;
                                    r1 = sqrt(ctbcx*ctbcx + ctbcy*ctbcy + ctbcz*ctbcz);
                                    if (r1 < servData->minDist) r1 = servData->minDist;

                                    if (resdOut.at(i) == "PHE" || resdOut.at(i) == "TYR") {
                                        xx.at(i+5) = xx.at(i) + 1.788*cx;
                                        yy.at(i+5) = yy.at(i) + 1.788*cy;
                                        zz.at(i+5) = zz.at(i) + 1.788*cz;
                                        assrd.at(i+5) = 1;
                                        xi= xx.at(i+5) - xx.at(i);
                                        yi= yy.at(i+5) - yy.at(i);
                                        zi= zz.at(i+5) - zz.at(i);
                                        fbd = 1.191/r1;
                                        for (int k=1;k<=4;k++) {
                                            xx.at(i+k) = xx.at(i) + (0.25 + 0.5 * (k>2))*xi + (2*(k%2) - 1)*ctbcx*fbd;
                                            yy.at(i+k) = yy.at(i) + (0.25 + 0.5 * (k>2))*yi + (2*(k%2) - 1)*ctbcy*fbd;
                                            zz.at(i+k) = zz.at(i) + (0.25 + 0.5 * (k>2))*zi + (2*(k%2) - 1)*ctbcz*fbd;
                                            assrd.at(i+k) = 1;
                                        }
                                        if(resdOut.at(i) == "TYR") {
                                            xx.at(i+6) = xx.at(i+5) + 0.937*cx;
                                            yy.at(i+6) = yy.at(i+5) + 0.937*cy;
                                            zz.at(i+6) = zz.at(i+5) + 0.937*cz;
                                            assrd.at(i+6) = 1;
                                        }
                                    } else if (resdOut.at(i) == "TRP") {
                                        fbd = 1.19/r1;
                                        xx.at(i+1) = xx.at(i) + 0.63*cx + fbd*ctbcx;    //!CD1
                                        yy.at(i+1) = yy.at(i) + 0.63*cy + fbd*ctbcy;
                                        zz.at(i+1) = zz.at(i) + 0.63*cz + fbd*ctbcz;
                                        xx.at(i+2) = xx.at(i) + 0.448*cx - fbd*ctbcx;   //!!CD2
                                        yy.at(i+2) = yy.at(i) + 0.448*cy - fbd*ctbcy;
                                        zz.at(i+2) = zz.at(i) + 0.448*cz - fbd*ctbcz;
                                        fbd = 0.688/r1;
                                        xx.at(i+3) = xx.at(i+1) + 0.819*cx - fbd*ctbcx; //!!NE1
                                        yy.at(i+3) = yy.at(i+1) + 0.819*cy - fbd*ctbcy;
                                        zz.at(i+3) = zz.at(i+1) + 0.819*cz - fbd*ctbcz;
                                        xx.at(i+4) = xx.at(i+2) + 0.774*cx + fbd*ctbcx; //!!CE2
                                        yy.at(i+4) = yy.at(i+2) + 0.774*cy + fbd*ctbcy;
                                        zz.at(i+4) = zz.at(i+2) + 0.774*cz + fbd*ctbcz;
                                        fbd = 1.19/r1;
                                        xx.at(i+5) = xx.at(i+2) - fbd*ctbcx - 0.448*cx; //!CE3
                                        yy.at(i+5) = yy.at(i+2) - fbd*ctbcy - 0.448*cy;
                                        zz.at(i+5) = zz.at(i+2) - fbd*ctbcz - 0.448*cz; 
                                        xx.at(i+6) = xx.at(i+4) - fbd*ctbcx + 0.448*cx; //!CZ2
                                        yy.at(i+6) = yy.at(i+4) - fbd*ctbcy + 0.448*cy;
                                        zz.at(i+6) = zz.at(i+4) - fbd*ctbcz + 0.448*cz;      
                                        xx.at(i+7) = xx.at(i+5) - fbd*ctbcx + 0.448*cx; //!CZ3
                                        yy.at(i+7) = yy.at(i+5) - fbd*ctbcy + 0.448*cy;
                                        zz.at(i+7) = zz.at(i+5) - fbd*ctbcz + 0.448*cz;      
                                        xx.at(i+8) = xx.at(i+7) + xx.at(i+4) - xx.at(i+2);    //!CH2
                                        yy.at(i+8) = yy.at(i+7) + yy.at(i+4) - yy.at(i+2);
                                        zz.at(i+8) = zz.at(i+7) + zz.at(i+4) - zz.at(i+2); 
                                        for (unsigned int k=1;k<=8;k++) assrd.at(i+k) = 1;
                                    } else if (resdOut.at(i) == "HIS") {
                                        fbd = 1.19/r1;
                                        xx.at(i+1) = xx.at(i) + 0.6*cx + fbd*ctbcx;  //!!ND1
                                        yy.at(i+1) = yy.at(i) + 0.6*cy + fbd*ctbcy;
                                        zz.at(i+1) = zz.at(i) + 0.6*cz + fbd*ctbcz;
                                        xx.at(i+2) = xx.at(i) + 0.63*cx - fbd*ctbcx; //!!CD2
                                        yy.at(i+2) = yy.at(i) + 0.63*cy - fbd*ctbcy;
                                        zz.at(i+2) = zz.at(i) + 0.63*cz - fbd*ctbcz;
                                        fbd = 0.688/r1;
                                        xx.at(i+3) = xx.at(i) + 1.35*cx + fbd*ctbcx; //!!CE1
                                        yy.at(i+3) = yy.at(i) + 1.35*cy + fbd*ctbcy;
                                        zz.at(i+3) = zz.at(i) + 1.35*cz + fbd*ctbcz;
                                        xx.at(i+4) = xx.at(i) + 1.41*cx - fbd*ctbcx; //!!NE2
                                        yy.at(i+4) = yy.at(i) + 1.41*cy - fbd*ctbcy;
                                        zz.at(i+4) = zz.at(i) + 1.41*cz - fbd*ctbcz;
                                        for (unsigned int k=1;k<=4;k++) assrd.at(i+k) = 1;
                                    }
                                }       
                            } else if (atomty.at(i) == " CD " || atomty.at(i) == " CD1") {
                                if (resdOut.at(i) == "PRO") { //!!this part treats the ring PRO CG-CD
                                    ip3 = ca[ress.at(i)][0]; // !!Ca                               
                                    x = 0.5*(xx.at(ip3-1) + xx.at(i-1) - xx.at(ip3) - xx.at(i-2));
                                    y = 0.5*(yy.at(ip3-1) + yy.at(i-1) - yy.at(ip3) - yy.at(i-2));
                                    z = 0.5*(zz.at(ip3-1) + zz.at(i-1) - zz.at(ip3) - zz.at(i-2));
                                    xx.at(i+1) = 0.5*(xx.at(ip3) + xx.at(i-2)) + 1.618*x;
                                    yy.at(i+1) = 0.5*(yy.at(ip3) + yy.at(i-2)) + 1.618*y;
                                    zz.at(i+1) = 0.5*(zz.at(ip3) + zz.at(i-2)) + 1.618*z;
                                    assrd.at(i+1) = 1;
                                } else {
                                    double bx = xx.at(ip3+1) - xx.at(ip3-1);
                                    double by = yy.at(ip3+1) - yy.at(ip3-1);
                                    double bz = zz.at(ip3+1) - zz.at(ip3-1);
                                    double cx = xx.at(i-1) - xx.at(i-2);
                                    double cy = yy.at(i-1) - yy.at(i-2);
                                    double cz = zz.at(i-1) - zz.at(i-2);
                                    double ctbcx = by*cz - bz*cy;
                                    double ctbcy = bz*cx - bx*cz;
                                    double ctbcz = bx*cy - by*cx;
                                    r1 = sqrt(ctbcx*ctbcx + ctbcy*ctbcy + ctbcz*ctbcz);
                                    if (r1 < servData->minDist) r1 = servData->minDist;

                                    if (resdOut.at(i) == "PHE" || resdOut.at(i) == "TYR") {
                                        assrd.at(i+4) = 1;
                                        xx.at(i+4) = xx.at(i-1) + 1.788*cx;
                                        yy.at(i+4) = yy.at(i-1) + 1.788*cy;
                                        zz.at(i+4) = zz.at(i-1) + 1.788*cz;
                                        xi= xx.at(i+4) - xx.at(i-1);
                                        yi= yy.at(i+4) - yy.at(i-1);
                                        zi= zz.at(i+4) - zz.at(i-1);
                                        fbd = 1.191/r1;
                                        for (int k=0;k<=3;k++) {
                                            xx.at(i+k) = xx.at(i-1) + (0.25 + 0.5 * (k>1))*xi + (1 - 2*(k%2))*ctbcx*fbd;
                                            yy.at(i+k) = yy.at(i-1) + (0.25 + 0.5 * (k>1))*yi + (1 - 2*(k%2))*ctbcy*fbd;
                                            zz.at(i+k) = zz.at(i-1) + (0.25 + 0.5 * (k>1))*zi + (1 - 2*(k%2))*ctbcz*fbd;
                                            assrd.at(i+k) = 1;
                                        }
                                        if(resdOut.at(i) == "TYR") {
                                            xx.at(i+5) = xx.at(i+4) + 0.937*cx;
                                            yy.at(i+5) = yy.at(i+4) + 0.937*cy;
                                            zz.at(i+5) = zz.at(i+4) + 0.937*cz;
                                            assrd.at(i+5) = 1;
                                        }
                                    } else if (resdOut.at(i) == "TRP") {
                                        fbd = 1.19/r1;
                                        xx.at(i) = xx.at(i-1) + 0.63*cx + fbd*ctbcx;    //!CD1
                                        yy.at(i) = yy.at(i-1) + 0.63*cy + fbd*ctbcy;
                                        zz.at(i) = zz.at(i-1) + 0.63*cz + fbd*ctbcz;
                                        xx.at(i+1) = xx.at(i-1) + 0.448*cx - fbd*ctbcx;   //!!CD2
                                        yy.at(i+1) = yy.at(i-1) + 0.448*cy - fbd*ctbcy;
                                        zz.at(i+1) = zz.at(i-1) + 0.448*cz - fbd*ctbcz;
                                        fbd = 0.688/r1;
                                        xx.at(i+2) = xx.at(i) + 0.819*cx - fbd*ctbcx; //!!NE1
                                        yy.at(i+2) = yy.at(i) + 0.819*cy - fbd*ctbcy;
                                        zz.at(i+2) = zz.at(i) + 0.819*cz - fbd*ctbcz;
                                        xx.at(i+3) = xx.at(i+1) + 0.774*cx + fbd*ctbcx; //!!CE2
                                        yy.at(i+3) = yy.at(i+1) + 0.774*cy + fbd*ctbcy;
                                        zz.at(i+3) = zz.at(i+1) + 0.774*cz + fbd*ctbcz;
                                        fbd = 1.19/r1;
                                        xx.at(i+4) = xx.at(i+1) - fbd*ctbcx - 0.448*cx; //!CE3
                                        yy.at(i+4) = yy.at(i+1) - fbd*ctbcy - 0.448*cy;
                                        zz.at(i+4) = zz.at(i+1) - fbd*ctbcz - 0.448*cz; 
                                        xx.at(i+5) = xx.at(i+3) - fbd*ctbcx + 0.448*cx; //!CZ2
                                        yy.at(i+5) = yy.at(i+3) - fbd*ctbcy + 0.448*cy;
                                        zz.at(i+5) = zz.at(i+3) - fbd*ctbcz + 0.448*cz;      
                                        xx.at(i+6) = xx.at(i+4) - fbd*ctbcx + 0.448*cx; //!CZ3
                                        yy.at(i+6) = yy.at(i+4) - fbd*ctbcy + 0.448*cy;
                                        zz.at(i+6) = zz.at(i+4) - fbd*ctbcz + 0.448*cz;      
                                        xx.at(i+7) = xx.at(i+6) + xx.at(i+3) - xx.at(i+1);    //!CH2
                                        yy.at(i+7) = yy.at(i+6) + yy.at(i+3) - yy.at(i+1);
                                        zz.at(i+7) = zz.at(i+6) + zz.at(i+3) - zz.at(i+1); 
                                        for (unsigned int k=1;k<=7;k++) assrd.at(i+k) = 1;
                                    } else if (resdOut.at(i) == "HIS") {
                                        fbd = 1.19/r1;
                                        xx.at(i) = xx.at(i-1) + 0.6*cx + fbd*ctbcx;  //!!ND1
                                        yy.at(i) = yy.at(i-1) + 0.6*cy + fbd*ctbcy;
                                        zz.at(i) = zz.at(i-1) + 0.6*cz + fbd*ctbcz;
                                        xx.at(i+1) = xx.at(i-1) + 0.63*cx - fbd*ctbcx; //!!CD2
                                        yy.at(i+1) = yy.at(i-1) + 0.63*cy - fbd*ctbcy;
                                        zz.at(i+1) = zz.at(i-1) + 0.63*cz - fbd*ctbcz;
                                        fbd = 0.688/r1;
                                        xx.at(i+2) = xx.at(i-1) + 1.35*cx + fbd*ctbcx; //!!CE1
                                        yy.at(i+2) = yy.at(i-1) + 1.35*cy + fbd*ctbcy;
                                        zz.at(i+2) = zz.at(i-1) + 1.35*cz + fbd*ctbcz;
                                        xx.at(i+3) = xx.at(i-1) + 1.41*cx - fbd*ctbcx; //!!NE2
                                        yy.at(i+3) = yy.at(i-1) + 1.41*cy - fbd*ctbcy;
                                        zz.at(i+3) = zz.at(i-1) + 1.41*cz - fbd*ctbcz;
                                        for (unsigned int k=1;k<=3;k++) assrd.at(i+k) = 1;
                                    } else {
                                        x = 0.5*(xx.at(ip2) + xx.at(ip3));
                                        y = 0.5*(yy.at(ip2) + yy.at(ip3));
                                        z = 0.5*(zz.at(ip2) + zz.at(ip3));
                                        r1 = sqrt((x-ax)*(x-ax) + (y-ay)*(y-ay) + (z-az)*(z-az));
                                        if(r1>servData->minDist){
                                            fbd=bd/r1;
                                            xx.at(i) = ax + fbd*(ax-x);
                                            yy.at(i) = ay + fbd*(ay-y);
                                            zz.at(i) = az + fbd*(az-z);
                                        } else {
                                            assrd.at(i) = 2;           
                                            ip3 = ipp1;
                                            xx.at(i) = ax + xx.at(ip2) - xx.at(ip3);
                                            yy.at(i) = ay + yy.at(ip2) - yy.at(ip3);
                                            zz.at(i) = az + zz.at(ip2) - zz.at(ip3);
                                        }
                                    }
                                }       
                            }
                        }
                    }
                } else cout << "Warning! Abnormal index in addition of heavy atoms on the output atom number (icc): " << i+1 << endl;
            }
        } catch (out_of_range& e) {throw "Index is Out of Range on the output record #:" + to_string(i+1);}
    }
}

//======================================================================================================================    
//Function to add add hydrogen atoms
void ProcessData::HADD(unsigned int outAtomsNumber) {
    for(unsigned int i=outAtomsNumber-1;i>0;i--) {
        try {
            if (atomty.at(i).substr(1,1) == "H" && xx.at(i)>=990.0) {
                double ax, ay, az, bd, angg;
                int icc=-1, ip2=-1, ip3=-1; 
                unsigned short i1,i2,i3;
                vector<int> ipp0, ipp;
                icc = bonda.at(i)[0];  //!!bonded heavy atom
                bd = bond.at(i)[0];
                if(icc>=0) { // This is redundant additional check, because icc has to be definitely found following previous logic
                    assrd.at(i) = 1;
                    ax = xx.at(icc);  //!!origin point
                    ay = yy.at(icc);
                    az = zz.at(icc);
                    for (unsigned int j=0; j<4 && bonda.at(icc)[j]>=0; j++) 
                        if(xx[bonda.at(icc)[j]] <= 990.0 && i!=bonda.at(icc)[j] && atomty[bonda.at(icc)[j]].substr(1,1) != "H")
                            ipp0.push_back(bonda.at(icc)[j]);
                    if(ipp0.size()>0) ip2 = ipp0.at(0);

                    for (unsigned int j=0; j<4; j++) 
                        if(bonda.at(ip2)[j]>=0 && icc!=bonda.at(ip2)[j] && xx[bonda.at(ip2)[j]] <= 990.0 && atomty[bonda.at(ip2)[j]].substr(1,1) != "H")
                            ipp.push_back(bonda.at(ip2)[j]);
                    if(ipp.size()>0) ip3 = ipp.back();

                    i1 = atc.at(i)>atc.at(ip2) ? atc.at(ip2)-1:atc.at(i)-1;
                    i2 = atc.at(icc)-1;
                    i3 = atc.at(i)>atc.at(ip2) ? atc.at(i)-1:atc.at(ip2)-1;
                    angg = servData->HPI * servData->ag0[i1][i2][i3];
                    vector<double> vr0 = {xx.at(ip2)-ax, yy.at(ip2)-ay, zz.at(ip2)-az};
                    vector<vector<double>> vr1;
                    // Now everything is ready for the addition proccess
                    // !!!rotation axis , cross of vr0 and vr1, first H is transto ip2
                    double x,y,z,r1,xi,yi,zi,fbd;

                    if(atomty.at(i).substr(0,2) == "3H") {  // !!term     sp3-sp3  trans-, Gauche-, Gauche+
                        if(atomty.at(icc).substr(1,1) == "C") bd*=0.933;
                        assrd.at(i-1) = assrd.at(i-2) = 1;
                        for(unsigned int k=0;k<ipp.size();k++) //!!connect heavy atoms vectors
                            vr1.push_back({xx.at(ipp.at(k))-xx.at(ip2), yy.at(ipp.at(k))-yy.at(ip2), zz.at(ipp.at(k))-zz.at(ip2)});
                           
                        x = xx.at(ip2);
                        y = yy.at(ip2);
                        z = zz.at(ip2);
                        xi = vr1.at(0)[1]*vr0.at(2)-vr1.at(0)[2]*vr0.at(1);
                        yi = vr1.at(0)[2]*vr0.at(0)-vr1.at(0)[0]*vr0.at(2);
                        zi = vr1.at(0)[0]*vr0.at(1)-vr1.at(0)[1]*vr0.at(0);
                        Rotate_Matrix(x,y,z,ax,ay,az,xi,yi,zi,angg); //  !H on N, each H is trans    
                        r1 = sqrt((x-ax)*(x-ax) + (y-ay)*(y-ay) + (z-az)*(z-az));
                        if(r1<0.01 || r1>5.0) {
                            r1 = sqrt(pow(xx.at(ip3)-xx.at(ip2),2) + pow(yy.at(ip3)-yy.at(ip2),2) + pow(zz.at(ip3)-zz.at(ip2),2));
                            if(r1>10.0 || r1<0.001) {
                                xx.at(i) = ax + (xx.at(ip2) - xx.at(ip3))/bond.at(ip2)[1];
                                yy.at(i) = ay + (yy.at(ip2) - yy.at(ip3))/bond.at(ip2)[1];
                                zz.at(i) = az + (zz.at(ip2) - zz.at(ip3))/bond.at(ip2)[1];
                                assrd.at(i) = 2; // neighbour position severa wrong CH3-NH3
                            } else {
                                fbd = bd/r1;
                                xx.at(i) = ax + (xx.at(ip2) - xx.at(ip3))*fbd;
                                yy.at(i) = ay + (yy.at(ip2) - yy.at(ip3))*fbd;
                                zz.at(i) = az + (zz.at(ip2) - zz.at(ip3))*fbd;
                            }
                        } else {
                            fbd = bd/r1;  //!!adjust bond length
                            xx.at(i) = ax + fbd*(x-ax);
                            yy.at(i) = ay + fbd*(y-ay);
                            zz.at(i) = az + fbd*(z-az);
                        }
                        
                        angg = 2.0*servData->PI/3.0;
                        x = xx.at(i);
                        y = yy.at(i);
                        z = zz.at(i);
                        xi = vr0.at(0);
                        yi = vr0.at(1);
                        zi = vr0.at(2);
                        Rotate_Matrix(x,y,z,ax,ay,az,xi,yi,zi,angg); //  !H on N, each H is trans    
                        r1 = sqrt((x-ax)*(x-ax) + (y-ay)*(y-ay) + (z-az)*(z-az));
                        if(r1<0.01 || r1>100.0) {
                            xx.at(i-1) = ax + cos(servData->PI * 109.5/180)*(xx.at(i)-ax);
                            yy.at(i-1) = ay + sin(servData->PI * 109.5/180)*(yy.at(i)-ay);
                            zz.at(i-1) = az;   
                            assrd.at(i-1) = 2; //warning CH3,NH3-2
                        } else {
                            fbd = bd/r1;  // !!adjust bond length
                            xx.at(i-1) = ax + fbd*(x-ax);
                            yy.at(i-1) = ay + fbd*(y-ay);
                            zz.at(i-1) = az + fbd*(z-az);  
                        }
// Here is a bug in fortran code i-2 could be out of range, so I will implement adddititonal check
                        if (i>=2) {
                            x = xx.at(i);
                            y = yy.at(i);
                            z = zz.at(i);
                            Rotate_Matrix(x,y,z,ax,ay,az,xi,yi,zi,-angg); //  !H on N, each H is trans    
                            r1 = sqrt((x-ax)*(x-ax) + (y-ay)*(y-ay) + (z-az)*(z-az));
                            if(r1<0.01 || r1>10.0) {
                                xx.at(i-2) = ax + cos(servData->PI * 109.5/180)*(xx.at(i)-ax);
                                yy.at(i-2) = ay + sin(-servData->PI * 109.5/180)*(yy.at(i)-ay);
                                zz.at(i-2) = az;   
                                assrd.at(i-2) = 2; //warning CH3,NH3-2
                            } else {
                                fbd = bd/r1;  // !!adjust bond length
                                xx.at(i-2) = ax + fbd*(x-ax);
                                yy.at(i-2) = ay + fbd*(y-ay);
                                zz.at(i-2) = az + fbd*(z-az);  
                            }
                            i-=2;
                        } else i--;
                    // !!!!sp3 -term finish **********************************************************        

                    } else if(atomty.at(i).substr(0,2) == "2H") { //!!in term or middle
                        if(atomty.at(icc).substr(1,1) == "C") bd *= 0.933;
                        assrd.at(i-1) = 1;
                        if(atomty.at(icc).substr(1,1) == "N") {   // !!in term, sp2
                            ip3 = (atomty.at(ipp.at(0)).substr(1,1) == "C")? ipp.at(1) : ipp.at(0); // !!N=O -C
                            vr1.push_back({xx.at(ip3)-xx.at(ip2), yy.at(ip3)-yy.at(ip2), zz.at(ip3)-zz.at(ip2)});
                            x = xx.at(ip2);
                            y = yy.at(ip2);
                            z = zz.at(ip2);
                            xi = vr1.at(0)[1]*vr0.at(2)-vr1.at(0)[2]*vr0.at(1);
                            yi = vr1.at(0)[2]*vr0.at(0)-vr1.at(0)[0]*vr0.at(2);
                            zi = vr1.at(0)[0]*vr0.at(1)-vr1.at(0)[1]*vr0.at(0);
                            Rotate_Matrix(x,y,z,ax,ay,az,xi,yi,zi,angg); //  !H on N, each H is trans    
                            r1 = sqrt((x-ax)*(x-ax) + (y-ay)*(y-ay) + (z-az)*(z-az));
                            if(r1<0.01 || r1>5.0) {
                                assrd.at(i) = 2; // neighbour position wrong
                                r1 = sqrt(pow(xx.at(ip3)-xx.at(ip2),2) + pow(yy.at(ip3)-yy.at(ip2),2) + pow(zz.at(ip3)-zz.at(ip2),2));
                                if(r1>10.0 || r1<0.001) {
                                    xx.at(i) = ax + (xx.at(ip2) - xx.at(ip3))*0.6;
                                    yy.at(i) = ay + (yy.at(ip2) - yy.at(ip3))*0.6;
                                    zz.at(i) = az + (zz.at(ip2) - zz.at(ip3))*0.6;
                                } else {
                                    fbd = bd/r1;
                                    xx.at(i) = ax + (xx.at(ip2) - xx.at(ip3))*fbd;
                                    yy.at(i) = ay + (yy.at(ip2) - yy.at(ip3))*fbd;
                                    zz.at(i) = az + (zz.at(ip2) - zz.at(ip3))*fbd;
                                }
                            } else {
                                fbd = bd/r1;  //!!adjust bond length
                                xx.at(i) = ax + fbd*(x-ax);
                                yy.at(i) = ay + fbd*(y-ay);
                                zz.at(i) = az + fbd*(z-az);
                            }

                            x = xx.at(ip2);
                            y = yy.at(ip2);
                            z = zz.at(ip2);
                            Rotate_Matrix(x,y,z,ax,ay,az,xi,yi,zi,-angg); //  !H on N, each H is trans    
                            r1 = sqrt((x-ax)*(x-ax) + (y-ay)*(y-ay) + (z-az)*(z-az));
                            if(r1<0.01 || r1>5.0) {
                                r1 = sqrt(pow(xx.at(ip3)-xx.at(ip2),2) + pow(yy.at(ip3)-yy.at(ip2),2) + pow(zz.at(ip3)-zz.at(ip2),2));
                                if(r1>5.0 || r1<0.001) {
                                    xx.at(i-1) = ax + (xx.at(ip2) - xx.at(ip3))*0.6;
                                    yy.at(i-1) = ay + (yy.at(ip2) - yy.at(ip3))*0.6;
                                    zz.at(i-1) = az + (zz.at(ip2) - zz.at(ip3))*0.6;
                                    assrd.at(i-1) = 2; // neighbour position wrong NH2-2
                                } else {
                                    fbd = bd/r1;
                                    xx.at(i-1) = ax + (xx.at(ip3) - xx.at(ip2))*fbd;
                                    yy.at(i-1) = ay + (yy.at(ip3) - yy.at(ip2))*fbd;
                                    zz.at(i-1) = az + (zz.at(ip3) - zz.at(ip2))*fbd;
                                }
                            } else {
                                fbd = bd/r1;  //!!adjust bond length
                                xx.at(i-1) = ax + fbd*(x-ax);
                                yy.at(i-1) = ay + fbd*(y-ay);
                                zz.at(i-1) = az + fbd*(z-az);
                            }
                            if(resdOut.at(i) == "PRO" && ress.at(i)==0) {   //!!PRO N-term
                                xi = vr0.at(0);
                                yi = vr0.at(1);
                                zi = vr0.at(2);
                                for(int k=-1;k<=0;k++) {
                                    x = xx.at(i+k);
                                    y = yy.at(i+k);
                                    z = zz.at(i+k);
                                    Rotate_Matrix(x,y,z,ax,ay,az,xi,yi,zi, 0.5*servData->PI); // !H on N, each H is trans  
                                    xx.at(i+k) = x;
                                    yy.at(i+k) = y;
                                    zz.at(i+k) = z;
                                }
                            }
                        } else {  //  !!!C inter connected, not in term, rotate first H on of ip2 around icc-ip4
                            double ip4;

                            if(ipp0.at(1) == ip2) {
                                ip4 = ipp0.at(0);
                                if(abs(icc-ip4)>20) cout << "ip4 wrong!!! Atom number:"<< i+1 << endl;
                            } else ip4 = ipp0.at(1);
                            angg *=1.1;  // !!!120/110 degree

                            for(unsigned int ii=1;ii<=2;ii++) {
                                if(ii==2) {
                                    i2 = ip2;
                                    ip2 = ip4;
                                    ip4 = i2;
                                    angg=-angg;
                                }
                                x = xx.at(ip2);
                                y = yy.at(ip2);
                                z = zz.at(ip2);
                                xi = xx.at(ip4)-ax;
                                yi = yy.at(ip4)-ay;
                                zi = zz.at(ip4)-az;
                                Rotate_Matrix(x,y,z,ax,ay,az,xi,yi,zi,angg); //  !H on C, each H is trans   
                                r1 = sqrt((x-ax)*(x-ax) + (y-ay)*(y-ay) + (z-az)*(z-az));
                                if(r1<0.01 || r1>5.0) {
                                    r1 = sqrt(pow(xx.at(ip3)-xx.at(ip2),2) + pow(yy.at(ip3)-yy.at(ip2),2) + pow(zz.at(ip3)-zz.at(ip2),2));
                                    if(r1>10.0 || r1<0.001) {
                                        xx.at(i) = ax + (xx.at(ip2) - xx.at(ip3))/bond.at(ip2)[1];
                                        yy.at(i) = ay + (yy.at(ip2) - yy.at(ip3))/bond.at(ip2)[1];
                                        zz.at(i) = az + (zz.at(ip2) - zz.at(ip3))/bond.at(ip2)[1]; 
                                        assrd.at(i) = 2; //neighbour bond position severa wrong CH2-1
                                    } else {
                                        fbd = bd/r1;
                                        xx.at(i) = ax + (xx.at(ip3) - xx.at(ip2))*fbd;
                                        yy.at(i) = ay + (yy.at(ip3) - yy.at(ip2))*fbd;
                                        zz.at(i) = az + (zz.at(ip3) - zz.at(ip2))*fbd;  
                                    }
                                } else {
                                    fbd = bd/r1; //  !!adjust bond length
                                    if(ii==1) {
                                        xx.at(i) = ax + fbd*(x-ax);
                                        yy.at(i) = ay + fbd*(y-ay);
                                        zz.at(i) = az + fbd*(z-az);              
                                    } else {
                                        xx.at(i) = 0.5*(xx.at(i) + ax + fbd*(x-ax));
                                        yy.at(i) = 0.5*(yy.at(i) + ay + fbd*(y-ay));
                                        zz.at(i) = 0.5*(zz.at(i) + az + fbd*(z-az));
                                    }
                                }
                                x =xx.at(ip2);
                                y =yy.at(ip2);
                                z =zz.at(ip2);
                                Rotate_Matrix(x,y,z,ax,ay,az,xi,yi,zi,-angg); //  !H on C, each H is trans   
                                r1 = sqrt((x-ax)*(x-ax) + (y-ay)*(y-ay) + (z-az)*(z-az));
                                if(r1<0.01 || r1>5.0) {
                                    r1 = sqrt(pow(xx.at(ip3)-xx.at(ip2),2) + pow(yy.at(ip3)-yy.at(ip2),2) + pow(zz.at(ip3)-zz.at(ip2),2));
                                    if(r1>10.0 || r1<0.001) {
                                        xx.at(i) = ax + (xx.at(ip2) - xx.at(ip3))/bond.at(ip2)[1];
                                        yy.at(i) = ay + (yy.at(ip2) - yy.at(ip3))/bond.at(ip2)[1];
                                        zz.at(i) = az + (zz.at(ip2) - zz.at(ip3))/bond.at(ip2)[1]; 
                                        assrd.at(i) = 2; //neighbour bond position severa wrong CH2-2
                                    } else {
                                        fbd = bd/r1;
                                        xx.at(i-1) = ax + (xx.at(ip2) - xx.at(ip3))*fbd;
                                        yy.at(i-1) = ay + (yy.at(ip2) - yy.at(ip3))*fbd;
                                        zz.at(i-1) = az + (zz.at(ip2) - zz.at(ip3))*fbd; 
                                    }
                                } else {
                                    fbd = bd/r1; // !!adjust bond length
                                    if(ii<2) {
                                        xx.at(i-1) = ax + fbd*(x-ax);
                                        yy.at(i-1) = ay + fbd*(y-ay);
                                        zz.at(i-1) = az + fbd*(z-az);   
                                    } else {
                                        xx.at(i-1) = 0.5*(xx.at(i-1) + ax + fbd*(x-ax));
                                        yy.at(i-1) = 0.5*(yy.at(i-1) + ay + fbd*(y-ay));
                                        zz.at(i-1) = 0.5*(zz.at(i-1) + az + fbd*(z-az));
                                    }
                                }    
                            }

                        }
                        i--;
                    // !!  2H finished

                    } else { //!!1H case, sp2 has been assign, O-H and S-H also done trans-; only for sp3 C case
                        if(ipp0.size()==3) {  //   !!sp3 case
                            bd = bd*0.945;
                            x = (xx.at(ipp0.at(0)) + xx.at(ipp0.at(1)) + xx.at(ipp0.at(2)))/3.0;
                            y = (yy.at(ipp0.at(0)) + yy.at(ipp0.at(1)) + yy.at(ipp0.at(2)))/3.0;
                            z = (zz.at(ipp0.at(0)) + zz.at(ipp0.at(1)) + zz.at(ipp0.at(2)))/3.0;
                            r1 = sqrt((x-ax)*(x-ax) + (y-ay)*(y-ay) + (z-az)*(z-az));
                            if(r1<0.01 || r1>10.0) {
                                assrd.at(i) = 2; //warning H1-sp3
                                xx.at(i) = ax + 0.3;
                                yy.at(i) = ay + 0.6;
                                zz.at(i) = az + 0.7;
                            } else {
                                fbd = bd/r1; //  !!adjust bond length
                                xx.at(i) = ax + fbd*(ax-x);
                                yy.at(i) = ay + fbd*(ay-y);
                                zz.at(i) = az + fbd*(az-z);
                            }
                        } else if(ipp0.size()==2) { // !!sp2, aromatic ring H, H-N
                            bd = bd*0.95;
                            x = (xx.at(ipp0.at(0)) + xx.at(ipp0.at(1)))*0.5;
                            y = (yy.at(ipp0.at(0)) + yy.at(ipp0.at(1)))*0.5;
                            z = (zz.at(ipp0.at(0)) + zz.at(ipp0.at(1)))*0.5;
                            r1 = sqrt((x-ax)*(x-ax) + (y-ay)*(y-ay) + (z-az)*(z-az));
                            if(r1<0.01 || r1>5.0) {
                                assrd.at(i) = 2; //warning H1-sp3
                                xx.at(i) = ax + 0.3;
                                yy.at(i) = ay + 0.6;
                                zz.at(i) = az + 0.7;
                            } else {
                                fbd = bd/r1; //  !!adjust bond length
                                xx.at(i) = ax + fbd*(ax-x);
                                yy.at(i) = ay + fbd*(ay-y);
                                zz.at(i) = az + fbd*(az-z);
                            }
                        } else if(ipp0.size()==1) { // !!H-O,H-S
                            bd = bd*0.94;
                            vr1.push_back({xx.at(ipp.at(0)) - xx.at(ip2), yy.at(ipp.at(0)) - yy.at(ip2), zz.at(ipp.at(0)) - zz.at(ip2)});
                            x = xx.at(ip2);
                            y = yy.at(ip2);
                            z = zz.at(ip2);
                            xi = vr1.at(0)[1]*vr0.at(2) - vr1.at(0)[2]*vr0.at(1);
                            yi = vr1.at(0)[2]*vr0.at(0) - vr1.at(0)[0]*vr0.at(2);
                            zi = vr1.at(0)[0]*vr0.at(1) - vr1.at(0)[1]*vr0.at(0);
                            Rotate_Matrix(x,y,z,ax,ay,az,xi,yi,zi,angg); //  !H on N, each H is trans            
                            r1 = sqrt((x-ax)*(x-ax) + (y-ay)*(y-ay) + (z-az)*(z-az));
                            if(r1<0.01 || r1>5.0) {
                                assrd.at(i) = 2; // warning O-H
                                xx.at(i) = ax + 0.3;
                                yy.at(i) = ay + 0.6;
                                zz.at(i) = az + 0.7;  
                            } else {
                                fbd = bd/r1;  // !!adjust bond length
                                xx.at(i) = ax + fbd*(x-ax);
                                yy.at(i) = ay + fbd*(y-ay);
                                zz.at(i) = az + fbd*(z-az);
                            }
                        }       
                    }
                } else cout << "Warning! Abnormal index in addition of heavy atoms on the output atom number (icc): " << i+1 << endl;
            }
        } catch (out_of_range& e) {throw "Index is Out of Range on the output record #:" + to_string(i+1);}
    }
}

//======================================================================================================================    
//Function to adjust 1H case 
void ProcessData::AD_1H(unsigned int outAtomsNumber) {
    for(unsigned int i=0;i<outAtomsNumber;i++) {
        try {
            int icc=-1, ip2=-1;
            double ax, ay, az, bd, angg, x, y, z, r, xi, yi, zi;
            vector<int> HH1;
            vector<double> rHH1;
            if(atomty.at(i).substr(1,1) == "H") {
                icc = bonda.at(i)[0]; //  !!bonded heavy atom
                if(atomty.at(icc).substr(1,1)=="O" || atomty.at(icc).substr(1,1)=="S") {
                    bd = bond.at(i)[0];
                    ax = xx.at(icc); //  !!origin point
                    ay = yy.at(icc);
                    az = zz.at(icc);
                    ip2 = bonda.at(icc)[0]; // !ip2-icc-iH
//                    np0=0 !1record the number of atoms in contact
                    for(unsigned int j=0;j<outAtomsNumber;j++) {
                        if(ress.at(i) != ress.at(j)) {
                            r = pow(xx.at(i)-xx.at(j),2) + pow(yy.at(i)-yy.at(j),2) + pow(zz.at(i)-zz.at(j),2);

                            if(r>0.01 && r<25.0) {
                                HH1.push_back(j);
                                rHH1.push_back(r);
                            } else if(r>144.0) {
                                 j = (ress.at(j) < nRes-1)  ? ca.at(ress.at(j)+1)[0] - 2 : outAtomsNumber;
                            }
                        }
                    }

                    if(HH1.size()>0) {
                        vector<unsigned int> indx = SortIndex(&rHH1);
                        for(unsigned int j : indx) {
                            if(rHH1.at(j)>=0.01) {
                                const map<string, double> rSel  = {{"C", 9.0}, {"S", 9.0}, {"N", 8.12}, {"H", 4.0}, {"O", 7.67}};
                                if (auto rIt = rSel.find(atomty.at(HH1.at(j)).substr(1,1)); rIt !=rSel.end()) r = (*rIt).second;
                            }
                        }
                        x=y=z=0.0;
                        for(unsigned int j=0;j<HH1.size();j++) {
                            unsigned int ind = indx.at(j);
                            int ipin = HH1.at(ind);
                            if (rHH1.at(ind)>=0.01) {
                                int MH1=0;
                                if(atomty.at(ipin).substr(1,1)=="C" || atomty.at(ipin).substr(1,1)=="S") r = 9.0;
                                else if(atomty.at(ipin).substr(1,1)=="N") r = 8.12; // !2.85^2
                                else if(atomty.at(ipin).substr(1,1)=="H") r = 4.0;
                                else if(atomty.at(ipin).substr(1,1)=="O") r = 7.67; // !2.77^2
                                if(rHH1.at(ind)<r && atomty.at(ipin).substr(1,1)!="O") MH1 = -1; // !!repulsion and Hbond
                                else if(atomty.at(ipin).substr(1,1)=="O") MH1 = 7;
                                else if(atomty.at(ipin).substr(1,1)=="N") MH1 = 1;
                                x += MH1*(xx.at(ipin)-xx.at(i))/rHH1.at(ind);
                                y += MH1*(yy.at(ipin)-yy.at(i))/rHH1.at(ind);
                                z += MH1*(zz.at(ipin)-zz.at(i))/rHH1.at(ind);
                            }
                        }
                        
                        // !!the vector from current position
                        r = sqrt(x*x+y*y+z*z);
                        xi = xx.at(ip2) - ax;
                        yi = yy.at(ip2) - ay;
                        zi = zz.at(ip2) - az;
                        if(r<0.01) r = 0.01;
                        x/=r; y/=r; z/=r;
                        r = (x*(xx.at(i) - xx.at(icc)) + y*(yy.at(i) - yy.at(icc)) + z*(zz.at(i) - zz.at(icc)))/bond.at(i)[0];
                        r /= (atomty.at(icc).substr(1,1)=="O")? 0.956: 0.996;
                        if(r<1.0) {
                            angg = acos(r);
                            x = xx.at(i);
                            y = yy.at(i);
                            z = zz.at(i);
                            Rotate_Matrix(x,y,z,ax,ay,az,xi,yi,zi,angg);  // !H on N, each H is trans
                            xx.at(i) = x;
                            yy.at(i) = y;
                            zz.at(i) = z;
                        }
                    }
                }
            }
        } catch (out_of_range& e) {throw "Index is Out of Range on the output record #:" + to_string(i+1);}    
    }
}


//********** Service function Rotate_Matrix
//!rotate (x,y,z) with origin (ax,ay,az) around axis (xi,yi,zi) an
void ProcessData::Rotate_Matrix(double &x, double &y, double &z, double ax, double ay, double az, double &xi, double &yi, double &zi, double an) {
    double rf = sqrt(xi*xi + yi*yi + zi*zi);
    if(rf>=0.1) {
        xi = xi/rf;
        yi = yi/rf;
        zi = zi/rf;
        double r12=xi*xi+yi*yi;
        double r23=yi*yi+zi*zi;
        double r13=xi*xi+zi*zi;
        double rcor=xi*x+yi*y+zi*z;
        double xt = ax*r23+xi*(-ay*yi-az*zi+rcor)+((x-ax)*r23+xi*(ay*yi+az*zi-yi*y-zi*z))*cos(an)+(ay*zi-az*yi-zi*y+yi*z)*sin(an);
        double yt = ay*r13 + yi*(-ax*xi-az*zi+rcor) + ((y-ay)*r13+yi*(ax*xi+az*zi-xi*x-zi*z))*cos(an) + (-ax*zi+az*xi+zi*x-xi*z)*sin(an);
        z = az*r12+zi*(-ax*xi-ay*yi+rcor)+((z-az)*r12+zi*(ax*xi+ay*yi-xi*x-yi*y))*cos(an)+(ax*yi-ay*xi-yi*x+xi*y)*sin(an);
        x = xt;
        y = yt;
    }
}

//********** Sort index *****************
vector<unsigned int> ProcessData::SortIndex(vector<double> *arrIn) {
    map<double, unsigned int> map2sort;
    vector<unsigned int> result;
    for(unsigned int jj=0; jj<arrIn->size();jj++ ) {
        map2sort.insert (pair<double,unsigned int>(arrIn->at(jj),jj) );
    }
    for (auto it=map2sort.begin(); it!=map2sort.end(); it++) result.push_back((*it).second);
    return result;
}
