#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include "acc-process.h"

using namespace std;

//************************ Arguments Parsing ***************************************
ArgsSet ProcessData::ArgsParsing(vector<string> *args) {
    ArgsSet output;
    for(unsigned int i=0; i<args->size();i++) {
        string el = args->at(i), snum;
        if (el[0]=='-') {
            if (el.substr(0,2)=="-p"){
                if (i+1 == args->size()) throw string("Invalid value for the argument of the option '-p'");
                try { output.probe = stod(args->at(++i)); continue;
                } catch(invalid_argument& e) {throw string("Invalid value for the argument of the option '-p'");}
            } 
            if (el.substr(0,8)=="--probe=")            
                try { output.probe = stod(el.substr(el.find('=') + 1));
                } catch(invalid_argument& e) {throw string("Invalid value for the argument of the option '--probe='");}
            if (el.substr(0,2)=="-z"){
                if (i+1 == args->size()) throw string("Invalid value for the argument of the option '-z'");
                try { output.zslice = stod(args->at(++i)); continue;
                } catch(invalid_argument& e) {throw string("Invalid value for the argument of the option '-z'");}
            } 
            if (el.substr(0,9)=="--zslice=")            
                try { output.zslice = stod(el.substr(el.find('=') + 1));
                } catch(invalid_argument& e) {throw string("Invalid value for the argument of the option '--zslice='");}
            if(el=="-het" || el=="--hetas") output.hetas = true;
            if (el=="-hyd" || el=="--hydro") output.hydro = true;
            if (el=="-wat" || el=="--waters") output.waters = true;
            if (el=="-all" || el=="--all") output.hetas = output.hydro = output.waters = true;
        } else if (el.size()>2) if (output.pdbfile == "") output.pdbfile = el; else output.vdwfile = el; 
    }
    return output;
}

//************************ Read Van der Waal radii file data ***********************
void ProcessData::ReadVanDerWaalRadiiFile(string vdwname) {
    ifstream vdwfi;
    string vdwrec;
    // Open the file 
    vdwfi.open(vdwname.c_str());
    if (!vdwfi) throw "Fail to open Van der Waal radii file! -> " + vdwname ;
    
    unsigned int cnt = 0; // total number of lines from the input file
    while (getline(vdwfi, vdwrec)) {
        try {
            if (vdwrec.substr(0,7) == "RESIDUE") {
                aacids.push_back((vdwrec.substr(8,4) == "ATOM") ? vdwrec.substr(13,3):
                                 (vdwrec.substr(8,6) == "HETATM") ? vdwrec.substr(15,3): vdwrec.substr(14,3));
                numats.push_back(0);                                   
                anames.push_back({});
                vradii.push_back({});
            }
            if (vdwrec.substr(0,4) == "ATOM") {
                if (aacids.size()>0) {
                    numats.back()++;
                    anames.back().push_back(vdwrec.substr(5,4));
                    vradii.back().push_back(stod(vdwrec.substr(10)));
                } else throw "Fail to read Van der Waal radii file! Line number: " + cnt+1;
            } 
        } catch(logic_error& e) {throw "Fail to read Van der Waal radii file! Line number: " + to_string(cnt+1);}
        cnt++;
    }
    vdwfi.close();
}

//************************ Read standard.data file (Relative Accessibilites) if exists ***********************
void ProcessData::ReadStandardDataFile() {
    ifstream stdfi;
    string stdrec;
    // Open the file 
    stdfi.open("standard.data");
    if (!stdfi) return;
    
    unsigned int cnt = 0; // total number of lines from the input file
    while (getline(stdfi, stdrec)) {
        try {
            if (stdrec.substr(0,4) == "ATOM") {
                stacids.push_back(stdrec.substr(12,3));
                standarea.push_back({stod(stdrec.substr(16,7)), stod(stdrec.substr(29,7)), stod(stdrec.substr(42,7)),
                                     stod(stdrec.substr(55,7)), stod(stdrec.substr(68,7))});
            }
        } catch(logic_error& e) {
            cout << "WARNING! Fail to read standard.data file! Line number: " + to_string(cnt+1) << endl
                 << "Assuming from now there is no Relative Accessibilites data to use." << endl << endl;
            stacids.clear();
            standarea.clear();
            return;     
          }
        cnt++;
    }
    stdfi.close();
}

//************************ Data processing *****************************************
int ProcessData::DataProcessing() {
    ifstream pdbfi;
    string pdbrec;

    pdbfi.open(argSet.pdbfile);
    if (!pdbfi) throw "Fail to open PDB file! -> " + argSet.pdbfile ;

// Start Log records
    logStream.push_back(" ACCALL - Accessibility calculations");
    logStream.push_back(" PDB FILE INPUT " + argSet.pdbfile);
    logStream.push_back(" PROBE SIZE     " + sformat(argSet.probe, 6, 2));
    logStream.push_back(" Z-SLICE WIDTH  " + sformat(argSet.zslice, 6, 3));
    logStream.push_back(string ( argSet.hetas? " INCL":" EXCL") + " HETATOMS");
    logStream.push_back(string ( argSet.hydro? " INCL":" EXCL") + " HYDROGENS");
    logStream.push_back(string ( argSet.waters? " INCL":" EXCL") + " WATERS");
    logStream.push_back(" READVDW " + sformat(aacids.size(),3,0) + " residues input");


// Process PDB file
    unsigned int cnt = 0; // total number of lines from the input file
    
    unsigned int num_chains = 0;
    int resok = -1, aok = -1;
    char firstalt = '-', last_chain = '-';
    string atom;
    const vector<string> rLab = {"RES","HEM","HOH"};
    map<string, double> vguess = {{"C", 1.80}, {"N", 1.60}, {"S", 1.85}, {"O", 1.40}, {"P", 1.90},
                                        {"CA", 2.07}, {"FE", 1.47}, {"CU", 1.78}, {"ZN", 1.39}, {"MG", 1.73}};
    vector<unsigned int> rty, atomtype;
    vector<int> resindex;
    vector<string> resnam, label;
    vector<double> rads;
    vector<Cartesians> coords;

    while (getline(pdbfi, pdbrec)) {
        cnt++;
        if (pdbrec.size()<54) continue; // Ignore short non-standard records
        double vdw = 0.0;
        int atype = pdbrec.substr(17,3)=="HOH"? 2:
                    pdbrec.substr(0,6)=="HETATM"? 1:
                    pdbrec.substr(0,4)=="ATOM"? 0: -1;

        if (!(atype == 0 || atype == 1 && argSet.hetas || atype == 2 && argSet.waters)) continue; // Read only required records
        if (pdbrec[16] != ' ') {                        // Ignore Alternate positions, other than blanks or 1st encountered
            if (firstalt == '-') firstalt = pdbrec[16];
            else if (firstalt != pdbrec[16]) continue;
        }
        if (pdbrec[13] != 'H' && pdbrec[13] != 'D') {
// Residues processing          
            if (pdbrec[21] != last_chain) {last_chain = pdbrec[21]; num_chains++;}
            if (resnam.size() == 0 || pdbrec.substr(17,10) != resnam.back()) {
                resnam.push_back(pdbrec.substr(17,10));
                rty.push_back(atype);
                if ((resok = getVecIndex(&aacids, pdbrec.substr(17,3)))<0)
                    logStream.push_back(" UNKNOWN residue type.............> " + resnam.back());
            }
            atom = pdbrec.substr(12,4);
            atomtype.push_back(0);
            if (getVecIndex(&Mchain, atom)>=0) atomtype.back() = 3;
            else if (getVecIndex(&Phobs, atom)>=0) atomtype.back() = 1;
            else if (getVecIndex(&Phils, atom)>=0) atomtype.back() = 2;
            else if (atom[1] == 'O' || atom[1] == 'N') atomtype.back() = 2;
            else atomtype.back() = 3; 
            if (pdbrec.substr(12,4) != "OXT"){
                if (resok>=0) { // Known residue type
                    vdw = 0.0;
                    aok = getVecIndex(&anames.at(resok), atom);
                }
                unsigned int ires = resok;
                if (resok < 0 || aok < 0) {
                    // if Not found, try atoms in all residue types 
                    for(ires=0;ires < anames.size() && (aok = getVecIndex(&anames.at(ires), atom)) < 0;ires++);
                    logStream.push_back(" NON-STANDARD atom." + atom + " in residue> " + resnam.back());
                    if (aok>=0) {
                        logStream.push_back(" ASSUMED vdw of " + atom + " in " + resnam.back() + " = "
                                             + sformat(vradii.at(ires).at(aok), 5,2) + " (same as " + aacids.at(ires) + ")");
                    } 
                }
                if (aok>=0)
                    vdw = vradii.at(ires).at(aok); 
                else {
                    // Still Not found, make a guess                      
                    map<string,double>::iterator it = vguess.find(atom.substr(0,2));
                    vdw = it!=vguess.end()? it->second: (it = vguess.find(atom.substr(1,1)))!=vguess.end()? it->second : 1.8;   
                    logStream.push_back(" GUESSED vdw of " + atom + " in " + resnam.back() + " = " + sformat(vdw, 5,2));
                }
            } else vdw = 1.40; // Assign radius to atom. Special case: OXT
            // Finally store required data from the current PDB record
            rads.push_back(vdw);
            label.push_back(pdbrec.substr(0,30));
            resindex.push_back(resnam.size() - 1);
            coords.push_back({0.0, 0.0, 0.0}); // Will be updated a few lines later later to support hydro atoms case
        } else
            if (argSet.hydro) {
                vdw = hyrad;
                // Not adding a new record for hydro atoms. 
                rads.back() = vdw;
                label.back() = pdbrec.substr(0,30);
                resindex.back() = resnam.size() - 1;
            } else continue; // Get away if Hydro types are not required 
        try { coords.back() = {stod(pdbrec.substr(30,8)), stod(pdbrec.substr(38,8)), stod(pdbrec.substr(46,8))};
        } catch (invalid_argument& e) {throw "Broken PDB file! Invalid numerical data. Line number: " + to_string(cnt);}                  
    }
// End of PDB file      
    pdbfi.close();

    logStream.push_back(" ADDED VDW RADII");
    logStream.push_back(" CHAINS   " + sformat(num_chains, 5, 0));
    logStream.push_back(" RESIDUES " + sformat(resnam.size(), 5, 0));
    logStream.push_back(" ATOMS    " + sformat(atomtype.size(), 5, 0));

// Continue data processing

// Call SOLVA - LEE & RICHARDS TYPE ACCESSIBLITY CALCULATIONS
    vector<double> accs;
    try{
        accs = SOLVA(&coords, &rads, argSet.probe, argSet.zslice);
    } catch(char const* e) {
        logStream.push_back(string(e));
        return -1;
    }    

    logStream.push_back(" SOLVA: PROGRAM ENDS CORRECTLY");

    for(unsigned int i=0; i<atomtype.size(); i++) 
        asaStream.push_back(label.at(i) + sformat(coords.at(i).x, 8, 3)
                                        + sformat(coords.at(i).y, 8, 3)
                                        + sformat(coords.at(i).z, 8, 3)
                                        + "          " + sformat(accs.at(i),8 , 3));
    
    logStream.push_back(" CALCULATED ATOMIC ACCESSIBILITES");

// Sum atomic accessibilities by residue. Process ATOM and HETATM records.
// Use Relative Accessibilities for the 20 common aminos.
// Output is written to .rsa file

    struct AbsRel {
        double abs = 0.0;    
        double rel = 0.0;    
    };
    vector<double> tsums = {0.0, 0.0, 0.0, 0.0, 0.0};
    vector<vector<AbsRel>> ressums; 
    for (unsigned int i=0; i<atomtype.size(); i++) {
        unsigned int ir = resindex.at(i);
        int ires = getVecIndex(&stacids, resnam.at(ir).substr(0,3));
        if (ir >= ressums.size()) ressums.push_back({{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}});
        tsums.at(0) += accs.at(i);
        ressums.at(ir).at(0).abs += +accs.at(i);
        if (ires>=0) {
            if (atomtype.at(i) == 3) {
                ressums.at(ir).at(4).abs += accs.at(i);
                tsums.at(4) += accs.at(i);
            } else
                if (atomtype.at(i) != 0) {
                    ressums.at(ir).at(3).abs += accs.at(i);
                    tsums.at(3) += accs.at(i);
                    ressums.at(ir).at(atomtype.at(i)).abs += accs.at(i);
                    tsums.at(atomtype.at(i)) += accs.at(i);
                }
            // Calculate relative accessibilities
            if (i == atomtype.size()-1 || ir < resindex.at(i+1))
                for(unsigned int j=0;j<5;j++)
                    if (standarea.at(ires).at(j) > 0.0) 
                        ressums.at(ir).at(j).rel = 100.0*ressums.at(ir).at(j).abs/standarea.at(ires).at(j);

        } else for (unsigned int j=0;j<5;j++) ressums.at(ir).at(j).rel = -99.9;             
    }

// Final outputs    
    if (stacids.size()>0) {
        rsaStream.push_back("REM  Relative accessibilites read from external file \"standard.data\"");
        logStream.push_back(" RELATIVE (STANDARD) ACCESSIBILITIES READFOR " + sformat(stacids.size(), 3,0) + " AMINO ACIDS");
    } else
        logStream.push_back(" NO STANDARD VALUES INPUT");
    logStream.push_back(" SUMMED ACCESSIBILITIES OVER RESIDUES");

    rsaStream.push_back("REM  File of summed (Sum) and % (per.) accessibilities for ");
    rsaStream.push_back("REM RES _ NUM      All atoms   Non P side   Polar Side   Total Side   Main Chain");
    rsaStream.push_back("REM                ABS   REL    ABS   REL    ABS   REL    ABS   REL    ABS   REL");
    string output;
    for(unsigned int i=0; i<=resindex.back(); i++) { 
        output = rLab.at(rty.at(i)) + " " + resnam.at(i) + " ";
        for(unsigned int j=0; j<5; j++) 
            output += sformat(ressums.at(i).at(j).abs, 7, 2) + sformat(ressums.at(i).at(j).rel, 6, 1); 
        rsaStream.push_back(output);
    }
    rsaStream.push_back("END  Absolute sums over accessible surface ");
    output = "TOTAL   ";
    for(unsigned int j=0; j<5; j++) 
            output +=  "     " + sformat(tsums.at(j), 8, 1);    
    rsaStream.push_back(output);

    return 0;
// =================    
}

//************************ Processed data output *****************************************
void ProcessData::DataOutput(string outFileName, vector<string> *out_stream) {
    
    ofstream fout;
    fout.open(outFileName.c_str(),ifstream::out | ifstream::trunc);
    try {
        for(unsigned int i=0; i<out_stream->size();i++) {
            fout << out_stream->at(i) << endl;
        }    
    } catch (ifstream::failure e) {throw "File output Error: (" + outFileName + ") -> " + e.what();}

    fout.close();
} 


//////////// Analytical and service functions ///////////////////////////////

// SOLVA - LEE & RICHARDS TYPE ACCESSIBLITY CALCULATIONS
// 
//  Calculate accessible surface area for a group of atoms.
//  The accessible area for a given atom is calculated by the 
//  formula:
//      (arcsum) x (atom radius+probe radius) x (deltaz)
//  Numerical integration is carried out over z. in each z-
//  section, the arcsum for a given atom is the arclength of 
//  the circle (intersection of the atom sphere with the z-
//  section) that is not interior to any other atom circles
//  in the same z-section.
// 
// ================================================================
// 
//   error parameter (zslice) - this gives accuracy of calculation
//                            - suitable values are 0.01 (high 
//                              accuracy) to 0.1 (low accuracy)
//                            - in detail the z sections are spaced
//                              at about error*diameter of atom
// 
//   probe size       - radius of probe in angstroms
//                    - suitable value for water = 1.4
//
// ================================================================

vector<double> ProcessData::SOLVA(const vector<Cartesians> *sCoords, const vector<double> *sRads,
                                  const double sProbe, const double sZslice) {
    const unsigned int nInt = 2000;   // maximum number of sphere intersections
    const unsigned int nCube = 10000; // maximum number of cubes allowed for placing of atoms
    const unsigned int nAC = 150;     // maximum number of atoms per cube

    const double PI = acos(-1.0);
    const double PIx2 = 2.0*PI;  

    const unsigned int nats = sRads->size();
    const int nzp = floor(1.0/sZslice + 0.5);
    vector<double> accs = vector<double> (nats, 0.0); // <-Initialising an array for Output
    Cartesians coordMins = {9999.0, 9999.0, 9999.0};
    Cartesians coordMaxs = {-9999.0, -9999.0, -9999.0};
    vector<unsigned int> itab = vector<unsigned int> (nCube, 0);
    vector<int> cube;
    
    vector<vector<unsigned int>> natm; // natm[nAC][nCube]
        
    int karc = nInt, idim, jidim, kjidim, i, j, k, kji;
    double rmax = 0.0;

    for(i=0;i<nAC;i++) natm.push_back(vector<unsigned int>(nCube, 0)); // Cubes initialisation

//  Radius of an atom sphere = atom radius + probe radius
//  Find maxima and minima
    for(unsigned int ii = 0; ii < nats; ii++) {
        if (sRads->at(ii) + sProbe > rmax) rmax = sRads->at(ii) + sProbe;
        if (coordMins.x > sCoords->at(ii).x) coordMins.x = sCoords->at(ii).x;
        if (coordMins.y > sCoords->at(ii).y) coordMins.y = sCoords->at(ii).y;
        if (coordMins.z > sCoords->at(ii).z) coordMins.z = sCoords->at(ii).z;
        if (coordMaxs.x < sCoords->at(ii).x) coordMaxs.x = sCoords->at(ii).x;
        if (coordMaxs.y < sCoords->at(ii).y) coordMaxs.y = sCoords->at(ii).y;
        if (coordMaxs.z < sCoords->at(ii).z) coordMaxs.z = sCoords->at(ii).z;
    }
//  rmax = max diameter
    rmax *= 2.0;

// Cubicals containing the atoms are setup. 
// The dimension of an edge equals the largest atom sphere radius
// The cubes have a single index
// Minimum of 3 by 3 cubic grid
// EXIT if max cubes exceeded
    idim = (coordMaxs.x - coordMins.x)/rmax+1.0;
    if (idim < 3.0) idim = 3.0;
    jidim = (coordMaxs.y - coordMins.y)/rmax+1.0;
    if (jidim < 3.0) jidim = 3.0;
    jidim *= idim;
    kjidim = (coordMaxs.z - coordMins.z)/rmax+1.0;
    if (kjidim < 3.0) kjidim = 3.0;
    kjidim *= jidim;
    if(kjidim > nCube) throw " SOLVA ERROR: max cubes exceeded";    
// Prepare the grid upto <nCube> cubes each containing upto <nAC> atoms.
// The cube index is kji. The atom index for each cube is stored in itab[nCube]
    for(unsigned int n=0;n<nats;n++) {
        i = (sCoords->at(n).x-coordMins.x)/rmax + 1.0;
        j = (sCoords->at(n).y-coordMins.y)/rmax; 
        k = (sCoords->at(n).z-coordMins.z)/rmax;
        kji = k*jidim + j*idim + i;
        if (itab.at(kji-1) + 1 > nAC) throw " SOLVA ERROR: max atoms per cube exceeded";
        itab.at(kji-1)++;
        natm.at(itab.at(kji-1)-1).at(kji-1) = n;
        cube.push_back(kji);
    }

// Process each atom in turn
    for (unsigned int ir = 0; ir<nats; ir++) {
        vector<double> dx,dy;
        vector<unsigned int> inov;
        double xr = sCoords->at(ir).x, yr = sCoords->at(ir).y, zr = sCoords->at(ir).z;
        double area = 0.0, rr = sRads->at(ir) + sProbe, rrx2 = rr*2.0,
               rrsq = (sRads->at(ir) + sProbe)*(sRads->at(ir) + sProbe);
        
        kji = cube.at(ir);
        //Find the 'mkji' cubes neighboring the kji cube
        bool okflag = true;
        for(k=-1;k<=1 && okflag;k++)
            for(j=-1;j<=1 && okflag;j++)
                for(i=-1;i<=1 && okflag;i++) {
                    int mkji = kji + k*jidim + j*idim + i;
                    if (mkji >= 1) {
                        if (mkji <= kjidim) {
                            for (int m = 0; m < itab.at(mkji-1); m++) {
                                unsigned int in = natm.at(m).at(mkji-1);
                                if (in != ir) {
                                    if (inov.size() + 1 > nInt) throw " SOLVA ERROR: number of sphere intersections > Max"; 
                                    dx.push_back(xr - sCoords->at(in).x);
                                    dy.push_back(yr - sCoords->at(in).y);
                                    inov.push_back(in);
                                }
                            }
                        } else  okflag = false; 
                    }
                }
        if (inov.size()>0) {
            // z resolution determined
            double zres = rrx2/nzp,
                   zgrid = sCoords->at(ir).z - rr - zres/2.0;
            for (i = 0; i < nzp; i++) {
                bool skipFlag = false;
                zgrid += zres;
                // find the radius of the circle of intersection of the <ir> sphere on the current z-plane
                double rsec2r = rrsq - (zgrid - zr)*(zgrid - zr), rsecr=sqrt(rsec2r);
                vector<double> arci, arcf; 
                karc = 0;
                for (j = 0; j < inov.size(); j++) {
                    skipFlag = false;
                    unsigned int in = inov.at(j);
                    // find radius of circle locus
                    double rsec2n = (sRads->at(in) + sProbe)*(sRads->at(in) + sProbe)
                                  - (zgrid - sCoords->at(in).z)*(zgrid - sCoords->at(in).z);
                    if (rsec2n <= 0.0) continue;                                   
                    double rsecn = sqrt(rsec2n);
                    double dsqj = dx.at(j)*dx.at(j) + dy.at(j)*dy.at(j);
                    // find intersections of n.circles with <ir> circles in section
                    double dj = sqrt(dsqj);
                    if (dj >= rsecr + rsecn) continue;
                    // do the circles intersect, or is one circle completely inside the other?
                    if (dj > abs(rsecr-rsecn)) {
                        karc++;
                        arci.push_back(0.0);
                        arcf.push_back(0.0);
                        if (arci.size() >= nInt) throw " SOLVA ERROR: max intersections exceeded";
                        // Initial and final arc endpoints are found for the ir circle intersected
                        // by a neighboring circle contained in the same plane. The initial endpoint
                        // of the enclosed arc is stored in arci, and the final arc in arcf
                        // law of cosines
                        double trig_test = (dsqj + rsec2r - rsec2n)/(2.0*dj*rsecr); 
                        if (trig_test >= 1.0) trig_test = 0.99999;
                        if (trig_test <=-1.0) trig_test =-0.99999;
                        // <alpha> is the angle between a line containing a point of intersection and
                        // the reference circle center and the line containing both circle centers
                        double alpha = acos(trig_test);
                        // <beta> is the angle between the line containing both circle centers and the x-axis
                        double beta = atan2(dy.at(j), dx.at(j)) + PI;
                        double ti = beta - alpha, tf = beta + alpha;
                        if (ti < 0.0) ti +=PIx2;
                        if (tf > PIx2) tf -=PIx2;
                        arci.back() = ti;
                        if (tf < ti) {
                            // if the arc crosses zero, then it is broken into two segments.
                            // the first ends at PIx2 and the second begins at zero
                            arcf.back() = PIx2;
                            karc++;
                            arci.push_back(0.0);
                            arcf.push_back(0.0);
                        }
                        arcf.back() = tf;
                    } else if (rsecr <= rsecn) {skipFlag = true; break;}  
                }
                if (skipFlag) continue; // skip some calculations on this turn

                // find the accessible surface area for the sphere ir on this section
                double arcsum;
                if (arci.size()>0) {
                    // the arc endpoints are sorted on the value of the initial arc endpoint
                    vector<unsigned int> tag = SortIndex(&arci);    
                    // calculate the accessible area
                    arcsum = arci.at(0);
                    double t = arcf.at(tag.at(0));
                    if (arci.size()>1)
                        for (k=1; k < arci.size(); k++) {
                            if(t < arci.at(k)) arcsum += arci.at(k) - t;
                            double tt = arcf.at(tag.at(k));
                            if (tt > t) t = tt;    
                        }
                    arcsum += PIx2 - t;
                } else arcsum = PIx2;
                // The area/radius is equal to the accessible arc length x the section thickness.
                // Add the accessible area for this atom in this section to the area for this
                // atom for all the section encountered thus far
                area += arcsum * zres;

            }       
        } else area = PIx2*rrx2;
        // scale area to vdw shell
        accs.at(ir) = area * rr;
    }
// End of SOLVA calculations
    return accs;
}


// Find the string in the vector of strings and return an index or -1 if not found
int ProcessData::getVecIndex(const vector<string> *invec, const string s2find) {
    int idx=-1;
    if (invec->size()>0) {
        auto it = find(invec->begin(), invec->end(), s2find);
        idx = (it == invec->end()? -1: distance(invec->begin(), it));
    }
    return idx;
};

// Sort array of doubles and return array of indexes of the sorted array 
// index[] = SortIndex(&arrIn[]={2.0, 3.5, -7.3}) -> arrIn[]=={-7.3, 2.0, 3.5}; index[] == {2, 0, 1}  
vector<unsigned int> ProcessData::SortIndex(vector<double> *arrIn) {
    multimap<double, unsigned int> map2sort;
    vector<unsigned int> output;
    unsigned int jj;
    for(jj=0; jj<arrIn->size();jj++ ) {
        map2sort.insert (pair<double,unsigned int>(arrIn->at(jj),jj) );
    }
    jj = 0;
    for (auto it=map2sort.begin(); it!=map2sort.end(); it++) {
        arrIn->at(jj++) = it->first;
        output.push_back(it->second);
    }    
    return output;
}

// Format number to a string with the given length <slen> and number of digits <digpos>
string ProcessData::sformat(double num, const unsigned int slen, const unsigned int digpos) {
    string output = ""; 
    string tstr = to_string(round(num * pow(10, digpos))/pow(10, digpos));
    int dotpos = tstr.find(".");
    if (dotpos < 0) {
        if (digpos > 0) {tstr+=".";tstr.append(digpos, '0');}  
    } else
        if (digpos > 0) {
            unsigned int dsize = tstr.substr(dotpos).size();
            if (dsize<digpos+1) tstr.append(digpos+1-dsize, '0');
            else tstr = tstr.substr(0,dotpos+digpos+1);
        } else tstr = tstr.substr(0,dotpos);
    if (tstr.size()<slen) output.append(slen-tstr.size(), ' ');
    output += tstr;
    return output;
}
/////////////////////////////////////////////////////////////////////////    
