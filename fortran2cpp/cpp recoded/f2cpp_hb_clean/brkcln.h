//===============================================================================
//     brkcln header file (fortran legacy brkcln.par)
//===============================================================================


class ConstSet {
public:
    const short aaonln =    100;
    const short abetsz =     29;
    const short ahsz   =      4;
    const short alphah =      1;
    const short altfi  =     15;
    const short altmax =    300;
    const short authfi =     14;
    const short bdgsep =      2;
    const short bend   =     10;
    const short blgszl =      5;
    const short blgszs =      2;
    const short bridge =      4;
    const short bridg1 =      5;
    const short bridg2 =      6;
    const short brkfi  =     11;
    const short calph  =      2;
    const short carb   =      3;
    const short cbeta  =      1;
    const short chiral =      4;
    const short chisgn =     11;
    const short chnfi  =      9;
    const short cst    =      1;
    const short cyscod =      3;
    const short dihdat =      6;
    const short errfi  =     12;
    const short fnamln =     80;
    const short impld  =      1;
    const short kappa  =      5;
    const short keylen =      6;
    const short logfi  =     18;
    const short maxalt =    600;
    const short maxbnd =      4;
    const short maxchi =      6;
    const short maxchn =     50;
    const short maxlin =   8500;
    const short maxres =   5700; // no need
    const short mnchna =      6;
    const short mprcod =     34;
    const short namfi  =     10;
    const short ncoord =      3;
    const short ndihed =      5;
    const short nhelix =      3;
    const short nhyd   =      5;
    const short nmain  =      1;
    const short nmnang =      7;
    const short nmod   =     14;
    const short nsdang =      6;
    const short nst    =      3;
    const short nstchs =     26;
    const short nstruc =     11;
    const short nturns =      3;
    const short omega  =      3;
    const short ooi1   =      1;
    const short ooi2   =      2;
    const short outfi  =     13;
    const short oxyg   =      4;
    const short phi    =      1;
    const short pih    =      8;
    const short pihsz  =      5;
    const short procod =     16;
    const short protfi =      8;
    const short psi    =      2;
    const short sdchna =     12;
    const short seccod =     32;
    const short sgamma =      2;
    const short sheet  =      3;
    const short sheet1 =      2;
    const short sideaa =     50;
    const short stdaa  =     26;
    const short tco    =      6;
    const short thrtnh =      7;
    const short thrtsz =      3;
    const short turn   =      9;
    const short unkcod =     -1;
    
// Add-ons from subroutines
    const short aapln  =     13;  
//------------------------------    

    const int allatm = maxres*mnchna;
    const int allcrd = maxres*mnchna*ncoord;


    const double cadbnd =   5.0;
    const double ccdst  =   2.0;
    const double null   = 999.9;
    const double ooi1rd =   8.0;
    const double ooi2rd =  14.0;
    const double pepbnd =   2.5;
    const double radian = (180.0 / 3.14159265);
    const double sssep  =   3.0;

    const std::string atmkey = "ATOM  ";
    const std::string bendch = "S";
    const std::string bndbeg = ">";
    const std::string bndbth = "X";
    const std::string bndend = "<";
    const std::string brkch  = "!";
    const std::string bulgch = "*";
    const std::string calpha = " CA ";
    const std::string cmpkey = "COMPND";
    const std::string dsfkey = "SSBOND";
    const std::string endkey = "END   ";
    const std::string fmlkey = "FORMUL";
    const std::string hedkey = "HEADER";
    const std::string hetkey = "HET   ";
    const std::string hlxkey = "HELIX ";
    const std::string htakey = "HETATM";
    const std::string maskey = "MASTER";
    const std::string modkey = "MODEL ";
    const std::string nonatm = "ZZZZ";
    const std::string pbch   = "b";
    const std::string seqkey = "SEQRES";
    const std::string shtkey = "SHEET ";
    const std::string shtsma = "e";
    const std::string shtsym = "E";
    const std::string sitkey = "SITE  ";
    const std::string soukey = "SOURCE";
    const std::string space  = " ";
    const std::string terkey = "TER   ";
    const std::string trnch  = "T";
    const std::string trnkey = "TURN  ";
    const std::string allspc = "                                            ";
    
    
// Next is from const data of clean.f
    const std::string cbet = " CB ";
    // aacids is a list aminoacids -> aacode3 and hetcod combined in 1 list of 50 aacids.
    const std::vector<std::string> aacids = {"ALA","ASX","CYS","ASP","GLU","PHE","GLY",         // 29 - aacids from aacode3 
                                             "HIS","ILE","XXX","LYS","LEU","MET","ASN",
                                             "XXX","PRO","GLN","ARG","SER","THR","XXX",
                                             "VAL","TRP","XXX","TYR","GLX","UNK","PCA","INI",  
                                             "AIB","PHL","SEC","ALM","MPR","FRD","PYR",         // 21 - from hetcod 
                                             "LYM","GLM","PPH","PGL","OLE","TPH","ABA",
                                             "NLE","B2V","B2I","B1F","BNO","B2A","B2F"};
    const std::vector<std::string> dudcod = {"  A","  C","  G","  T","  U","1MA","5MC",
                                             "OMC","1MG","2MG","M2G","7MG","OMG"," YG",
                                             "  I"," +U","H2U","5MU","PSU"," +C"," +A",
                                             "GAL","AGL","GCU","ASG","MAN","GLC","CEG",
                                             "G4S","AGS","NAG","GLS","NGS","NAM","EXC"};
    const std::vector<std::string> modcod = {"ACE","ANI","BZO","CBZ","CLT","FOR","NH2",
                                             "PHO","RHA","TFA","TOS"," NH","MYR","BOC"};
    const std::vector<std::string> atname = {" N  "," CA "," C  "," O  ","    "};
    const std::vector<std::string> tphnam = {" N  "," CA "," C  "," O1 "," O2 "};
    const std::vector<std::string> pphnam = {" N  "," CA "," P  "," OP1"," OP2"};
    const std::vector<std::string> olenam = {" ON "," CA "," C  "," O  ","    "};
    const std::vector<std::string> pyrnam = {" ON "," CA "," C  "," O  ","    "};
    const std::vector<std::string> batnam = {" N  "," CA "," B  "," O1 "," O2 "};
    const std::vector<std::vector<std::string>> scname = {{"    ","    ","    ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," AD1"," AD2","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" SG ","    ","    ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," OD1"," OD2","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD "," OE1"," OE2","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD1"," CD2"," CE1"," CE2"," CZ ","    ","    ","    ","    ","    "},
                                                          {"    ","    ","    ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," ND1"," CD2"," CE1"," NE2"," AD1"," AD2"," AE1"," AE2","    ","    "},
                                                          {" CG1"," CG2"," CD1","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {"    ","    ","    ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD "," CE "," NZ ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD1"," CD2","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," SD "," CE ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," OD1"," ND2","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {"    ","    ","    ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD ","    ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD "," OE1"," NE2","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD "," NE "," CZ "," NH1"," NH2","    ","    ","    ","    ","    "},
                                                          {" OG ","    ","    ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" OG1"," CG2","    ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {"    ","    ","    ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG1"," CG2","    ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD1"," CD2"," NE1"," CE2"," CE3"," CZ2"," CZ3"," CH2","    ","    "},
                                                          {"    ","    ","    ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD1"," CD2"," CE1"," CE2"," CZ "," OH ","    ","    ","    ","    "},
                                                          {" CG "," CD "," AE1"," AE2","    ","    ","    ","    ","    ","    ","    "},
                                                          {"    ","    ","    ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD "," OE ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD "," CE "," NZ "," CI1"," CI2"," CI3"," CI4"," NI2"," CI5"," CI6"},
                                                          {" CB1"," CB2","    ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD1"," CD2"," CE1"," CE2"," CZ ","    ","    ","    ","    ","    "},
                                                          {"SEG "," OD1"," OD2","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CM ","    ","    ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" SG ","    ","    ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD1"," CD2"," CE1"," CE2"," CZ ","    ","    ","    ","    ","    "},
                                                          {"    ","    ","    ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD "," CE "," NZ "," CM ","    ","    ","    ","    ","    ","    "},
                                                          {" CM ","    ","    ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD1"," CD2"," CE1"," CE2"," CZ ","    ","    ","    ","    ","    "},
                                                          {" P  "," O1P"," O2P","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD1"," CD2","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD1"," CD2"," CE1"," CE2"," CZ ","    ","    ","    ","    ","    "},
                                                          {" CG ","    ","    ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD "," CE ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG1"," CG2","    ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG1"," CG2"," CD1","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD1"," CD2"," CE1"," CE2"," CZ ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD "," CE ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {"    ","    ","    ","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD1"," CD2"," CE1"," CE2"," CZ ","    ","    ","    ","    ","    "}};
    const std::vector<std::vector<std::string>> altnam = {{" CG "," AD1"," AD2","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," AD1"," AD2","    ","    ","    ","    ","    ","    ","    ","    "},
                                                          {" CG "," CD "," AE1"," AE2","    ","    ","    ","    ","    ","    ","    "}};

    const std::vector<unsigned int> natsdc = {0, 0, 1, 2, 3, 2, 0, 2, 2, 0, 4, 2, 3, 2, 0, 1, 3, 5, 1, 1, 0, 1, 2, 0, 2,
                                              0, 0, 3, 4, 0, 2, 1, 1, 1, 2, 0, 4, 0, 2, 0, 2, 2, 1, 3, 1, 2, 2, 3, 0, 2};  

    const std::vector<unsigned int> extprc = {0, 0, 0, 4, 4, 4, 0, 0, 5, 0, 0, 2, 0, 0, 0, 0, 0, 4, 0, 3, 0, 2, 0, 0, 4,
                                              0, 0, 0, 0, 0, 4, 0, 0, 0, 4, 0, 0, 0, 4, 0, 2, 4, 0, 0, 2, 5, 4, 0, 0, 4};

    const std::vector<int> nscats = {1, 4, 2, 4, 5, 7, 0, 10, 4,-1, 5, 4, 4, 4,-1, 3, 5, 7, 2, 3,-1, 3,10,-1, 8,
                                     5,-1, 4, 12, 3, 7, 4, 2, 2, 7, 8, 6, 2,10, 4, 5, 9, 2, 4, 3, 4, 7, 4, 1, 7};

    ConstSet() {};
};
//===============================================================================
