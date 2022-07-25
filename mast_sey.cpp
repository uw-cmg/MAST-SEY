#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <random>
#include <ctime>
#include <chrono>
#include <string.h>
#include <array>

#define PI 3.14159265359
#define HA2EV 27.21138602
#define EV2HA 1/27.21138602
#define ANG2BOHR 1/0.5291772109
#define BOHR2ANG 0.5291772109

using namespace std;
string code_version = "4.1 (04-05-2022)";
vector<string> elems = {"XX","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"};


random_device rd;
mt19937 mt(rd());

streambuf* orig_buf = cout.rdbuf();


vector<double> ie_arr, q_arr, de_arr;
vector<vector<array<double,2> > > inel_arr, phon_plus_arr, phon_minus_arr, elas_arr, ene_elf, jdos_arr;
vector<array<double,2> > elas_alloy_arr, dos_arr;

int dircos_algo = 0;
int spher_sec = 0;
char dirname[100];
vector <int> atnum;
vector <double> atcomp;
double mass, rho, ef, wf, u0, vol;
double ini_angle = 0.0;
double eb = 0.;
double erange = 1000.*EV2HA;
double ebeg = 0.0;
int egrid = 500, icsintgrid = 1000, qintgrid = 100, qdepgrid = 100, qdepomgrid = 1000;
bool lin_prep = false;
int mc_elec = 1000;
int qdep = 1;
string qdepname = "SPA";
string elfeps;
bool elsepa = false;
int es_nuc = 3, es_el = 1, es_ex = 1, es_muffin = 0, es_mcpol = 0;
bool emfp_only = false;
bool fermi_only = false;
bool prep = false;
bool preprange = false;
bool save_sumr = false;
bool save_qdep = false;
int sq_egrid, sq_qgrid;
double sq_qmax;
bool save_coords = false;
bool distrib = false;
bool noout = false;
bool classical_ang = false;
bool use_dos = false;
bool feg_dos = false;
bool notir = false;
bool test = false;
double cc = 137.;

bool ins = false;
double eg, evb, chi;
bool ph = false;
double eps_zero, eps_inf;
vector <double> ph_loss;

string getTime();
void getInput(int argc, char** argv);
void printVersion(char** argv);
void readMaterialFile(string filename="material.in");
void readPhononFile(string filename="ph.in");
void readEpsFile(string filename="eps.in");
vector<double> read1colFile(string filename="energies.in");
vector<array<double,2> > read2colFile(string filename="dos.in");
vector<vector<array<double,2> > > readScatteringFile(string filename);
void checkSumRules(bool save_sumr);
void printInput();
vector<double> range(double a, double b, int n);
vector<double> logrange(double a, double b, int n);
vector<double> getDeArr();

vector<array<double,2> > elas(double ie, int at=0, double comp=1.0);
vector<array<double,2> > inel(double ie);
vector<array<double,2> > phonon_plus(double ie);
vector<array<double,2> > phonon_minus(double ie);
vector<array<double,2> > int_elastic_ang(const vector<double> &xarr, const vector<double> &yarr, bool cumint=false);
vector<array<double,2> > int_inelastic_ene(double (*f)(double,double,int), double x0, double x1, int ndiv, double args, bool cumint=false);
vector<array<double,2> > int_inelastic_ang(double (*f)(double,double,int), double x0, double x1, int ndiv, double ie, double de, bool cumint=false);
vector<array<double,2> > cumintVect(const vector<array<double,2>> &xyarr);
array<double,3> f_rotdircos(array<double,3> uvw, double ang0, double ang1);
vector<array<double,2> > linterp2dline(double x0,double y0, const vector<double> &xarr, const vector<vector<array<double,2> > > &ytuparr);
double linterp(double x, const vector<array<double,2> > &xyarr, bool xfind=false);
double linterp2d(double x0,double y0, const vector<double> &xarr, const vector<vector<array<double,2> > > &ytuparr, bool total=false, bool xfind=false);
double linterp3d(double x0, double y0, double z0, const vector<double> &xarr, const vector<vector<vector<array<double,2> > > > &ytuparr, bool total, bool xfind);
double fzero(double (*f)(double,double,double), double x0, double x1, double ww, double qq, double tol=1e-6);
double jdos(double e, double de, double r);
double sspa_elf(double w, double q);
double spa_dispers(double w0, double w, double q);
double spa_elf(double w, double q);
double elfq(double q, double om, int dq=1);
double qIntFun(double om, double ie, int dummy);

void prepareJDOS(const vector<array<double,2> > &dos);
void saveVector(vector<double> arr, string filename);
void saveVector(vector<array<double,2> > arr, string filename);
void saveVector(vector<array<double,3> > arr, string filename);
void saveVector(vector<vector<double> > arr, string filename, int ncols=2);
void saveVector(vector<vector<array<double,2> > > arr, string filename);
void saveVector(vector<vector<vector<array<double,2> > > > arr, string filename);
void saveMFP(string filename);
void saveQdep(int n_ene, int n_q, double q_max);
void saveCoordVector(vector<vector<array<double,3> > > arr, vector<int> second, string filename);
void printVector(vector<double> &arr);
void printVector(vector<array<double,2> > &arr);
void printVector(vector<vector<array<double,2> > > &arr);
void marker(string num);
void print(string text);
void print(int num);
void print(double num);
void print(string text, int num);
void print(string text, double num);
void printProgress();
int printStars(int progr, int is, int size);
double rand01();
string checkName (string name);

class Electron
{
    public:
        double e,de,iimfp,iemfp,itmfp,iphmfp_plus,iphmfp_minus,s_ef,s_u0,s_vb,s_eg,pathlength;
        int secondary;
        double angles[2];
        double defl[2];
        bool sc_type_el;
        bool sc_type_ph;
        array<int,2> sc_type_elinel{0,0};
        array<int,1> sc_type_phon{0};
        array<double,3> xyz{0.,0.,0.};
        array<double,3> uvw{0.,0.,1.};
        vector<array<double,3> > coord;
        bool inside, dead, phonon;
    Electron(double ie, double x=0.0, double y=0.0, double z=0.0, double u=0.0, double v=0.0, double w=1.0, int sec=0, bool ph=false)
    {
        e = ie;
        de = 0.0;
        xyz[0] = x;
        xyz[1] = y;
        xyz[2] = z;
        uvw[0] = u;
        uvw[1] = v;
        uvw[2] = w;
        if (save_coords)
        {
            coord.push_back(xyz);
        }
        defl[0] = 0.0;
        defl[1] = 0.0;
        iimfp = IIMFP();
        iemfp = IEMFP();  
        phonon = ph;
        if (phonon)
        {
            iphmfp_plus = IPHMFP_plus();
            iphmfp_minus = IPHMFP_minus();
            itmfp = iemfp+iimfp+iphmfp_plus+iphmfp_minus;
        }
        else
        {
            itmfp = iemfp+iimfp;
        }  
        inside = true;
        dead = false;
        secondary = sec;

        pathlength = 0.0;
    }
    void travel_s()
    {
        double rn = random01();
        double s = -(1./itmfp)*log(rn);
        pathlength += s;
        xyz[0]=xyz[0]+uvw[0]*s;
        xyz[1]=xyz[1]+uvw[1]*s;
        xyz[2]=xyz[2]+uvw[2]*s;
        if (save_coords)
        {
            coord.push_back(xyz);
        }
    }

    void dircos2ang()
    {
        angles[0]=acos(uvw[2]);
        angles[1]=atan2(uvw[1],uvw[0]);
    }

    void determ_scatter()
    {
        double rn = random01();
        if (rn < iemfp/itmfp)
        {
            sc_type_el = true;
            sc_type_elinel[0]++;
        }
        else if (rn < (iemfp+iimfp)/itmfp)
        {
            sc_type_el = false;
            sc_type_elinel[1]++;
        }
        else
        {
            sc_type_el = false;
            sc_type_ph = true;
            sc_type_phon[0]++;
        }
    }

    bool scatter()
    {
        double rn = random01();
        double rn2 = random01();
        double rn3 = random01();
        double rn4 = random01();
        double rn5 = random01();
        double rn6 = random01();
        defl[1] = rn*2.*PI;
        if (sc_type_el)
        {
            double tot_elast_int = linterp2d(e,-1,ie_arr,elas_arr,true);
            defl[0] = linterp2d(e,rn2*tot_elast_int,ie_arr,elas_arr,false,true);
            return false;
        }
        else if (sc_type_ph)
        {
            double detot_ph_plus_int = linterp2d(e,-1,ie_arr,phon_plus_arr,true);
            double detot_ph_minus_int = linterp2d(e,-1,ie_arr,phon_minus_arr,true);
            if (rn5 < detot_ph_plus_int/(detot_ph_plus_int+detot_ph_minus_int))
            {
                // de = linterp2d(e,rn6*detot_ph_plus_int,ie_arr,phon_plus_arr,false,true);
                de = 0.0095*EV2HA;
                e = e-de;
                died();
                if (! dead)
                {
                    iimfp = IIMFP();
                    iemfp = IEMFP();
                    if (phonon)
                    {
                        iphmfp_plus = IPHMFP_plus();
                        iphmfp_minus = IPHMFP_minus();
                        itmfp = iemfp+iimfp+iphmfp_plus+iphmfp_minus;
                    }
                    else
                    {
                        itmfp = iemfp+iimfp;
                    }
                }
            }
            else
            {
                // de = linterp2d(e,rn6*detot_ph_minus_int,ie_arr,phon_minus_arr,false,true);
                de = 0.0095*EV2HA;
                e = e+de;
                died();
                if (! dead)
                {
                    iimfp = IIMFP();
                    iemfp = IEMFP();
                    if (phonon)
                    {
                        iphmfp_plus = IPHMFP_plus();
                        iphmfp_minus = IPHMFP_minus();
                        itmfp = iemfp+iimfp+iphmfp_plus+iphmfp_minus;
                    }
                    else
                    {
                        itmfp = iemfp+iimfp;
                    }
                }
            }
            return false;
        }
        else
        {
            double detot_inel_int = linterp2d(e,-1,ie_arr,inel_arr,true);
            de = linterp2d(e,rn3*detot_inel_int,ie_arr,inel_arr,false,true);
            if (phonon && de < eg)
            {
                while (de < eg)
                {
                    rn3 = random01();
                    de = linterp2d(e,rn3*detot_inel_int,ie_arr,inel_arr,false,true);
                }
            }
            if (classical_ang)
            {
                defl[0] = asin(sqrt(de/e));
            }
            else
            {
                vector<array<double,2> > int_inelas_ang = int_inelastic_ang(&elfq,0.0,PI/2.,100,e,de,true);
                defl[0] = linterp(rn*int_inelas_ang[int_inelas_ang.size()-1][1],int_inelas_ang,true);
            }
            e = e-de;
            died();
            if (! dead)
            {
                iimfp = IIMFP();
                iemfp = IEMFP();
                if (phonon)
                {
                    iphmfp_plus = IPHMFP_plus();
                    iphmfp_minus = IPHMFP_minus();
                    itmfp = iemfp+iimfp+iphmfp_plus+iphmfp_minus;
                }
                else
                {
                    itmfp = iemfp+iimfp;
                }
            }
            if (use_dos) 
            {
                if (feg_dos)
                {
                    if (ins) {
                        s_ef = fzero(&jdos,0.,evb,de,rn4);
                    } else {
                        s_ef = fzero(&jdos,0.,ef,de,rn4);
                    }
                } else {
                    double s_ef_int = linterp2d(de,-1,de_arr,jdos_arr,true);
                    s_ef = linterp2d(de,rn4*s_ef_int,de_arr,jdos_arr,false,true);
                }
            } else {
                if (ins)
                    s_ef = evb;
                else
                    s_ef = ef;
            }
            return true;
        }
    }

    double IEMFP()
    {
        double dcs = linterp2d(e,-1,ie_arr,elas_arr,true);
        return dcs/vol;
    }

    double IIMFP()
    {
        return linterp2d(e,-1,ie_arr,inel_arr,true);
    }

    double IPHMFP_minus()
    {
        return linterp2d(e,-1,ie_arr,phon_minus_arr,true);
    }

    double IPHMFP_plus()
    {
        return linterp2d(e,-1,ie_arr,phon_plus_arr,true);
    }

    void escaped()
    {
        double rn = random01();
        double t;
        dircos2ang();
        if (xyz[2]<0.0)
        {
            double beta = PI-angles[0];
            if (e*cos(beta)*cos(beta) > u0)
            {
                t = 4.*sqrt(1.-u0/(e*cos(beta)*cos(beta)))/pow((1.+sqrt(1.-u0/(e*cos(beta)*cos(beta)))),2);
            }
            else
            {
                t = 0.;
            }
            if (rn < t)
            {
                inside = false;
                xyz[0] = xyz[0]+sin(beta)*cos(angles[1])*xyz[2]/cos(beta);
                xyz[1] = xyz[1]+sin(beta)*sin(angles[1])*xyz[2]/cos(beta);
                xyz[2] = 0.0;
                angles[0] = PI-asin(sin(beta)*sqrt(e/(e-u0)));
                e = e-u0;
                if (save_coords)
                {
                    coord[coord.size()-1] = xyz;
                    coord.push_back({xyz[0]+1000.*sin(angles[0])*cos(angles[1]),xyz[1]+1000.*sin(angles[0])*sin(angles[1]),xyz[2]+1000.*cos(angles[0])});
                }
            }
            else
            {
                if (notir)
                {
                    dead = true;
                } else {
                    uvw[2] = -1.*uvw[2];
                    xyz[2] = -1.*xyz[2];
                    if (save_coords)
                    {
                        coord[coord.size()-1] = {xyz[0],xyz[1],-1.*xyz[2]};
                    }
                }
            }
        }
    }

    void died()
    {
        if (e<u0+eg) { dead = true; }
    }

    double random01()
    {
        uniform_real_distribution<double> dist(0.0, 1.0);
        return dist(mt);
    }

};

int main(int argc, char** argv)
{
    printVersion(argv);
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-ins") == 0) {
            ins = true;
        }
    }
    readMaterialFile();
    getInput(argc,argv);
    if (ph)
    {
        readPhononFile();
    }
    // ===================================
    readEpsFile(); // ??? we call it later for prep, delete?
    // ===================================
    if(noout) {
        cout.rdbuf(NULL);
    }
    printInput();
    if (prep)
    {
        readEpsFile();

        if (lin_prep)
        {
            ie_arr = range(ebeg,erange,egrid);
        }
        else
        {
            ie_arr = logrange(ebeg,erange,egrid);
        }
        if (save_qdep)
        {
            saveQdep(sq_egrid,sq_qgrid,sq_qmax);
            print("\n#\nSaved ELF(q,e) to file elf_qdep.out.");
            exit(1);
        }
        checkSumRules(save_sumr);
        print("# Starting mc_sey in \"prepare\" mode to get scattering properties");
        printProgress();
        int progress = 0;
        for (size_t i = 0; i < ie_arr.size(); i++)
        {
            progress = printStars(progress,i,ie_arr.size());

            elas_arr.push_back(elas(ie_arr[i],atnum[0],atcomp[0]));
            for (size_t ia = 1; ia < atnum.size(); ia++)
            {
                elas_alloy_arr = elas(ie_arr[i],atnum[ia],atcomp[ia]);
                for (int k = 0; k < 606; k++)
                {
                    elas_arr[i][k][1] = elas_arr[i][k][1]+elas_alloy_arr[k][1];
                }
            }
            elas_alloy_arr.clear();

            if(!emfp_only)
            {
                inel_arr.push_back(inel(ie_arr[i]));
                if (ph)
                {
                    phon_plus_arr.push_back(phonon_plus(ie_arr[i]));
                    phon_minus_arr.push_back(phonon_minus(ie_arr[i]));
                }
            }
        }
        if (preprange)
        {
            string er_str = to_string(ebeg*HA2EV)+"_"+to_string(erange*HA2EV);
            saveVector(ie_arr, "energies_"+er_str+".in");
            saveVector(elas_arr, "elastic_"+er_str+".in");
            if(!emfp_only)
            {
                saveVector(inel_arr, "inelastic_"+er_str+".in");
                if (ph)
                {
                    saveVector(phon_plus_arr, "phonon_plus_"+er_str+".in");
                    saveVector(phon_minus_arr, "phonon_minus_"+er_str+".in");
                }
            }
            saveMFP("mfp_"+er_str+".plot");
            cout << "\nFinished preparing input files, results in \"elastic_"+er_str+".in\", \"inelastic_"+er_str+".in\" and \"energies_"+er_str+".in\"" << endl;
        } else
        {
            saveVector(ie_arr, "energies.in");
            saveVector(elas_arr, "elastic.in");
            if(!emfp_only)
            {
                saveVector(inel_arr, "inelastic.in");
                if (ph)
                {
                    saveVector(phon_plus_arr, "phonon_plus.in");
                    saveVector(phon_minus_arr, "phonon_minus.in");
                }
            }
            saveMFP("mfp.plot");
            cout << "\nFinished preparing input files, results in \"elastic.in\", \"inelastic.in\" and \"energies.in\"" << endl;
        }

        system("rm *dat elsepa.in");
        return 0;
    }
    else
    {
        erange = erange+u0; 
        ie_arr = read1colFile();
        if (erange>=ie_arr[ie_arr.size()-1])
        {
            cerr << "Incident energy higher than range generated in preparation step." << endl;
            exit(1);
        } else if (erange<=ie_arr[0])
        {
            cerr << "Incident energy lower than range generated in preparation step." << endl;
            exit(1);
        }
        elas_arr = readScatteringFile("elastic.in");
        inel_arr = readScatteringFile("inelastic.in");
        if (ph)
        {
            phon_plus_arr = readScatteringFile("phonon_plus.in");
            phon_minus_arr = readScatteringFile("phonon_minus.in");
        }
        if (use_dos && !feg_dos)
        {
            ifstream dosfile("jdos.in");
            dos_arr = read2colFile("dos.in");
            de_arr = getDeArr();
            if (!dosfile)
            {
                print("# Could not find jdos.in - preparing...");
                prepareJDOS(dos_arr);
                print("# Prepared Joint DOS and saved to jdos.in\n#");
            }
            jdos_arr = readScatteringFile("jdos.in");
            print("# Secondaries generated from joint DOS provided by user in jdos.in");
        } else if (feg_dos){
            print("# Secondaries generated from joint DOS of a Free Electron Gas");
        }
        if (spher_sec==1)
        {
            print("# Secondaries generated with spherical symmetry");
        } else if (spher_sec==2) {
            print("# Secondaries generated randomly");
        }
        if (classical_ang)
        {
            print("# Inelastic scattering angle approximated: asin(sqrt(de/e))");
            print("#");
        } else {
            print("#");
        }
        print("# Input files read");
        print("# Starting mc_sey in \"Simulation\" mode to get SEY\n#");
        vector<Electron > elec_arr;
        ini_angle = asin(sin(ini_angle)*sqrt((erange-u0)/erange));
        int i, progress;
        double s_ene;
        array<double,3> s_uvw{0.0,0.0,0.0};
        array<double,3> s_xyz{0.0,0.0,0.0};
        i = -1;
        progress = 0;
        auto t1 = std::chrono::high_resolution_clock::now();

        for (int n_e = 0; n_e < mc_elec; n_e++)
        {
            if (n_e==9)
            {
                auto t2 = std::chrono::high_resolution_clock::now();
                auto ttt = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
                cout << "# Estimated time of execution: " << setprecision(2) << (ttt*mc_elec/10./1000000.)/60. << " min\n#" << setprecision(17) << endl;
                printProgress();
            }
            progress = printStars(progress,n_e,mc_elec);
            elec_arr.push_back(Electron(erange,0.0,0.0,0.0,sin(ini_angle),0.0,cos(ini_angle),0,ph));
            while (i < (int)elec_arr.size()-1)
            {
                i++;
                
                while (elec_arr[i].inside && ! elec_arr[i].dead)
                {
                    
                    elec_arr[i].travel_s();
                    
                    elec_arr[i].escaped();
                    
                    if (elec_arr[i].inside && ! elec_arr[i].dead)
                    {
                        elec_arr[i].determ_scatter();
                        if (elec_arr[i].scatter())
                        {
                            // ionisation
                            if (elec_arr[i].de-eb>u0 && eb>0.001)
                            {
                                s_ene = elec_arr[i].de-eb;
                                s_xyz[0] = elec_arr[i].xyz[0];
                                s_xyz[1] = elec_arr[i].xyz[1];
                                s_xyz[2] = elec_arr[i].xyz[2];
                                if (spher_sec==1) {
                                    s_uvw = { sin(acos(2.*rand01()-1.))*cos(2.*rand01()*PI),
                                              sin(acos(2.*rand01()-1.))*sin(2.*rand01()*PI),
                                              cos(acos(2.*rand01()-1.)) };
                                } else  if (spher_sec==2) {
                                    s_uvw = { sin(rand01()*PI)*cos(2.*rand01()*PI),
                                              sin(rand01()*PI)*sin(2.*rand01()*PI),
                                              cos(rand01()*PI) };
                                } else {
                                    s_uvw = f_rotdircos(elec_arr[i].uvw,asin(cos(elec_arr[i].defl[0])),elec_arr[i].defl[1]+PI);
                                }
                                elec_arr.push_back(Electron(s_ene,s_xyz[0],s_xyz[1],s_xyz[2],s_uvw[0],s_uvw[1],s_uvw[2],elec_arr[i].secondary+1,ph));
                            }
                            // valence band interaction
                            else if (elec_arr[i].de+elec_arr[i].s_ef>u0+eg)
                            {
                                s_ene = elec_arr[i].de+elec_arr[i].s_ef;
                                s_xyz[0] = elec_arr[i].xyz[0];
                                s_xyz[1] = elec_arr[i].xyz[1];
                                s_xyz[2] = elec_arr[i].xyz[2];
                                if (spher_sec==1) {
                                    s_uvw = { sin(acos(2.*rand01()-1.))*cos(2.*rand01()*PI),
                                              sin(acos(2.*rand01()-1.))*sin(2.*rand01()*PI),
                                              cos(acos(2.*rand01()-1.)) };
                                } else  if (spher_sec==2) {
                                    s_uvw = { sin(rand01()*PI)*cos(2.*rand01()*PI),
                                              sin(rand01()*PI)*sin(2.*rand01()*PI),
                                              cos(rand01()*PI) };
                                } else {
                                    s_uvw = f_rotdircos(elec_arr[i].uvw,asin(cos(elec_arr[i].defl[0])),elec_arr[i].defl[1]+PI);
                                }
                                elec_arr.push_back(Electron(s_ene,s_xyz[0],s_xyz[1],s_xyz[2],s_uvw[0],s_uvw[1],s_uvw[2],elec_arr[i].secondary+1,ph));
                            }
                        } 
                        else if (elec_arr[i].sc_type_el)
                        {
                            elec_arr[i].uvw = f_rotdircos(elec_arr[i].uvw,elec_arr[i].defl[0],elec_arr[i].defl[1]);
                        }
                    }
                }
            }
        }
        vector<vector<array<double,3> > > coord_vec;
        vector<int> secondary_ind;
        vector<vector<double> > ene_distrib;
        int em = 0, tem = 0, bsc = 0, nem = 0, d_prim = 0, e_bsc = 0;
        for (size_t ei = 0; ei < elec_arr.size()-1; ei++)
        {
            if (save_coords)
            {
                coord_vec.push_back(elec_arr[ei].coord);
                secondary_ind.push_back(elec_arr[ei].secondary);
            }
            if (! elec_arr[ei].inside)
            {
                em++; 
                if (elec_arr[ei].e < 50.*EV2HA)
                {
                    tem++; 
                    if (distrib)
                    {
                        ene_distrib.push_back({elec_arr[ei].e*HA2EV,elec_arr[ei].angles[0],elec_arr[ei].angles[1],elec_arr[ei].xyz[0],elec_arr[ei].xyz[1],(double)elec_arr[ei].secondary});
                    }
                } else {
                    bsc++; 
                }
                if (elec_arr[ei].secondary == 0 && elec_arr[ei].e > (erange-u0)-0.0001)
                {
                    e_bsc++; 
                }
                if (elec_arr[ei].secondary == 0 && elec_arr[ei].e <= (erange-u0)-0.0001)
                {
                    d_prim++; 
                }
            } else {
                nem++; 
            }
        }
        if (save_coords) { saveCoordVector(coord_vec,secondary_ind,checkName("mc_coords.plot")); }
        if (distrib)
        {
            saveVector(ene_distrib,checkName("mc_distrib.plot"),6);
        }
        print("\n#");
        cout << fixed << setprecision(4) << setfill(' ');
        cout << "# Energy[eV]     SEY TrueSEY   Bcksc DifPrim  eBcksc" << endl;
        if (noout) { cout.rdbuf(orig_buf); }
        cout << setw(12) << (erange-u0)*HA2EV;
        cout << setw(8) << (double)em/(double)mc_elec;
        cout << setw(8) << (double)tem/(double)mc_elec;
        cout << setw(8) << (double)bsc/(double)mc_elec;
        cout << setw(8) << (double)d_prim/(double)mc_elec;
        cout << setw(8) << (double)e_bsc/(double)mc_elec << endl;
        return 0;
    }
    return 0;
}

string getTime()
{
    time_t now = time(0);
    tm *ltm = localtime(&now);
    string dt = to_string(ltm->tm_mday)+"."+to_string(1+ltm->tm_mon)+"."+to_string(-100+ltm->tm_year)+" at "+to_string(ltm->tm_hour)+":"+to_string(ltm->tm_min)+":"+to_string(ltm->tm_sec);
    return dt;
}

void getInput(int argc, char** argv)
{
    double elow;
    if (argc == 1)
    {
        cerr << "No arguments specified, use \"-h\" flag for options.\nAt least the \"prepare\" keyword or \"-e\" flag is needed." << endl;
        exit(1);
    }
    if (strcmp(argv[1], "prepare") == 0) { prep = true; }

    if (ins) {
        elow = eg+eg+evb+1e-4;
    } else {
        elow = ef+1e-4;
    }
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-e") == 0 && prep)
        {
            if (argc > i+3)
            {
                if (argv[i+1][0] != '-' && argv[i+2][0] != '-' && argv[i+3][0] != '-')
                {
                    preprange = true;
                    ebeg = stod(argv[i+1])*EV2HA;
                    erange = stod(argv[i+2])*EV2HA;
                    egrid = stoi(argv[i+3]);
                    if (ebeg <= elow)
                    {
                        if (ins) {
                            cout << "# WARNING: Initial energy too low, setting to E_gap + E_vb = " << (eg+evb)*HA2EV << "\n#" << endl;
                        } else {
                            cout << "# WARNING: Initial energy too low, setting to E_fermi = " << ef*HA2EV << "\n#" << endl;
                        }
                        ebeg = elow;
                    }
                } else if (argv[i+1][0] != '-' && argv[i+2][0] != '-' && argv[i+3][0] == '-')
                {
                    ebeg = elow;
                    erange = stod(argv[i+1])*EV2HA;
                    egrid = stoi(argv[i+2]);
                    if (erange <= elow) {
                        if (ins)
                            cerr << "Energy range lower than E_gap + E_vb, stopping." << endl;
                        else
                            cerr << "Energy range lower than E_fermi, stopping." << endl;
                        exit(1);
                    }
                } else {
                    cerr << "Too few arguments for " << argv[i] << ", needs 2 or 3 arguments." << endl;
                    exit(1);
                }
            } else if (argc == i+3)
            {
                if (argv[i+1][0] != '-' && argv[i+2][0] != '-')
                {
                    ebeg = elow;
                    erange = stod(argv[i+1])*EV2HA;
                    egrid = stoi(argv[i+2]);
                } else {
                    cerr << "Too few arguments for " << argv[i] << ", needs 2 or 3 arguments." << endl;
                    exit(1);
                }
            } else {
                cerr << "Too few arguments for " << argv[i] << ", needs 2 or 3 arguments." << endl;
                exit(1);
            }
        } else if (strcmp(argv[i], "-e") == 0) {
            if (argc > i+1)
            {
                if (argv[i+1][0] != '-')
                {
                    erange = stod(argv[i+1])*EV2HA;
                } else {
                    cerr << "Too few arguments for " << argv[i] << ", needs 1 argument." << endl;
                    exit(1);
                }
            } else {
                cerr << "Too few arguments for " << argv[i] << ", needs 1 argument." << endl;
                exit(1);
            }
        }

        if (strcmp(argv[i], "-lin") == 0) { lin_prep = true; }
        if (strcmp(argv[i], "-coord") == 0) { save_coords = true; }
        if (strcmp(argv[i], "-distr") == 0) { distrib = true; }
        if (strcmp(argv[i], "-sumr") == 0) { save_sumr = true; }
        if (strcmp(argv[i], "-noout") == 0) { noout = true; }
        if (strcmp(argv[i], "-noang") == 0) { classical_ang = true; }
        if (strcmp(argv[i], "-sphsec") == 0) { spher_sec = 1; }
        if (strcmp(argv[i], "-rndsec") == 0) { spher_sec = 2; }
        if (strcmp(argv[i], "-notir") == 0) { notir = true; }
        if (strcmp(argv[i], "-test") == 0) { test = true; }
        if (strcmp(argv[i], "-emfp") == 0) { emfp_only = true; }
        if (strcmp(argv[i], "SOLID") == 0) { es_muffin = 1; }
        if (strcmp(argv[i], "LDA") == 0) { es_mcpol = 2; }

        if (strcmp(argv[i], "-ph") == 0)
        {
            ph = true;
        }

        if (strcmp(argv[i], "-dircos") == 0)
        {
            if (argc > i+1)
            {
                if (argv[i+1][0] != '-')
                {
                    dircos_algo = stoi(argv[i+1]);
                } else {
                    cerr << "Too few arguments for " << argv[i] << ", needs 1 argument." << endl;
                    exit(1);
                }
            } else {
                cerr << "Too few arguments for " << argv[i] << ", needs 1 argument." << endl;
                exit(1);
            }
        }
            

        if (strcmp(argv[i], "-dos") == 0)
        {
            use_dos = true;
            if (argc > i+1)
            {
                if (argv[i+1][0] != '-')
                {
                    if (strcmp(argv[i+1], "FEG") == 0)
                    {
                        feg_dos = true;
                    }
                }
            }
        }
        if (strcmp(argv[i], "-i") == 0)
        {
            if (argc > i+2)
            {
                if (argv[i+1][0] != '-' && argv[i+2][0] != '-')
                {
                    icsintgrid = stoi(argv[i+1]);
                    qintgrid = stoi(argv[i+2]);
                } else {
                    cerr << "Too few arguments for " << argv[i] << ", needs 2 arguments." << endl;
                    exit(1);
                }
            } else {
                cerr << "Too few arguments for " << argv[i] << ", needs 2 arguments." << endl;
                exit(1);
            }
        }
        if (strcmp(argv[i], "-q") == 0)
        {
            if (argc > i+2)
            {
                if (argv[i+1][0] != '-' && argv[i+2][0] != '-')
                {
                    qdepgrid = stoi(argv[i+1]);
                    qdepomgrid = stoi(argv[i+2]);
                } else {
                    cerr << "Too few arguments for " << argv[i] << ", needs 2 arguments." << endl;
                    exit(1);
                }
            } else {
                cerr << "Too few arguments for " << argv[i] << ", needs 2 arguments." << endl;
                exit(1);
            }
        }
        if (strcmp(argv[i], "-m") == 0)
        {
            if (argc > i+1)
            {
                if (argv[i+1][0] != '-')
                {
                    mc_elec = stoi(argv[i+1]);
                } else {
                    cerr << "Too few arguments for " << argv[i] << ", needs 1 argument." << endl;
                    exit(1);
                }
            } else {
                cerr << "Too few arguments for " << argv[i] << ", needs 1 argument." << endl;
                exit(1);
            }
            if (mc_elec<1000)
            {
                cerr << mc_elec << " MC e-'s is too few, results would make no sense." << endl;
                exit(1);
            }
        }
        if (strcmp(argv[i], "-core") == 0)
        {
            if (argc > i+1)
            {
                if (argv[i+1][0] != '-')
                {
                    eb = stod(argv[i+1])*EV2HA;
                } else {
                    cerr << "Too few arguments for " << argv[i] << ", needs 1 argument." << endl;
                    exit(1);
                }
            } else {
                cerr << "Too few arguments for " << argv[i] << ", needs 1 argument." << endl;
                exit(1);
            }
        }
        if (strcmp(argv[i], "-pa") == 0)
        {
            if (argc > i+1)
            {
                if (argv[i+1][0] != '-')
                {
                    ini_angle = stod(argv[i+1])*PI/180.;
                } else {
                    cerr << "Too few arguments for " << argv[i] << ", needs 1 argument." << endl;
                    exit(1);
                }
            } else {
                cerr << "Too few arguments for " << argv[i] << ", needs 1 argument." << endl;
                exit(1);
            }
        }
        if (strcmp(argv[i], "-saveq") == 0)
        {
            save_qdep = true;
            if (argc > i+3)
            {
                if (argv[i+1][0] != '-' && argv[i+2][0] != '-' && argv[i+3][0] != '-')
                {
                    sq_egrid = stoi(argv[i+1]);
                    sq_qgrid = stoi(argv[i+2]);
                    sq_qmax = stod(argv[i+3]);
                } else {
                    cerr << "Too few arguments for " << argv[i] << ", needs 3 arguments." << endl;
                    exit(1);
                }
            } else {
                cerr << "Too few arguments for " << argv[i] << ", needs 3 arguments." << endl;
                exit(1);
            }
        }
        if (strcmp(argv[i], "-qdep") == 0)
        {
            if (argc > i+1)
            {
                if (argv[i+1][0] != '-')
                {
                    qdepname=argv[i+1];
                } else {
                    cerr << "Too few arguments for " << argv[i] << ", needs 1 argument." << endl;
                    exit(1);
                }
            } else {
                cerr << "Too few arguments for " << argv[i] << ", needs 1 argument." << endl;
                exit(1);
            }
            if (strcmp(argv[i+1], "CUSTOM") == 0)
            {
                qdep = 0;
            } else if (strcmp(argv[i+1], "DFT") == 0)
            {
                qdep = 0;
            } else if (strcmp(argv[i+1], "SPA") == 0)
            {
                qdep = 1;
            } else if (strcmp(argv[i+1], "SSPA") == 0)
            {
                qdep = 2;
            } else if (strcmp(argv[i+1], "FPA") == 0)
            {
                qdep = 3;
            } else
            {
                cout << "# WARNING: option " << argv[i+1] << " for qdep not known, using SPA." << endl;
                qdep = 1;
            }
        }
        if (strcmp(argv[i], "-elastic") == 0)
        {
            elsepa = true;
            if (argc > i+3)
            {
                if (argv[i+1][0] != '-' && argv[i+2][0] != '-' && argv[i+3][0] != '-')
                {
                    if (strcmp(argv[i+1], "P") == 0)
                    {
                        es_nuc = 1;
                    } else if (strcmp(argv[i+1], "U") == 0)
                    {
                        es_nuc = 2;
                    } else if (strcmp(argv[i+1], "F") == 0)
                    {
                        es_nuc = 3;
                    } else if (strcmp(argv[i+1], "UU") == 0)
                    {
                        es_nuc = 4;
                    } else if (strcmp(argv[i+1], "radial") == 0)
                    {
                        elsepa = false;
                    } else {
                        cout << "# WARNING: option " << argv[i+1] << " for elsepa nucl. not known, using F." << endl;
                    }
                    if (strcmp(argv[i+2], "TFM") == 0)
                    {
                        es_el = 1;
                    } else if (strcmp(argv[i+2], "TFD") == 0)
                    {
                        es_el = 2;
                    } else if (strcmp(argv[i+2], "DHFS") == 0)
                    {
                        es_el = 3;
                    } else if (strcmp(argv[i+2], "DF") == 0)
                    {
                        es_el = 4;
                    } else {
                        cout << "# WARNING: option " << argv[i+2] << " for elsepa elec. not known, using TFM." << endl;
                    }
                    if (strcmp(argv[i+3], "NO") == 0)
                    {
                        es_ex = 0;
                    } else if (strcmp(argv[i+3], "FM") == 0)
                    {
                        es_ex = 1;
                    } else if (strcmp(argv[i+3], "TF") == 0)
                    {
                        es_ex = 2;
                    } else if (strcmp(argv[i+3], "RT") == 0)
                    {
                        es_ex = 3;
                    } else {
                        cout << "# WARNING: option " << argv[i+3] << " for exch. not known, using FM." << endl;
                    }
                } else {
                    cerr << "Too few arguments for " << argv[i] << ", needs 3 arguments." << endl;
                    exit(1);
                }
            } else {
                cerr << "Too few arguments for " << argv[i] << ", needs 3 arguments." << endl;
                exit(1);
            }
        }
        if (strcmp(argv[i], "-radial") == 0)
        {
            if (argc > i+1)
            {
                if (argv[i+1][0] != '-')
                {
                    if (strcmp(argv[i+1], "DHFS") == 0)
                    {
                        es_el = 1;
                    } else if (strcmp(argv[i+1], "TFM") == 0)
                    {
                        es_el = 2;
                    } else {
                        cout << "# WARNING: option " << argv[i+1] << " for radial not known, using DHFS." << endl;
                    }
                } else {
                    cerr << "Too few arguments for " << argv[i] << ", needs 1 argument." << endl;
                    exit(1);
                }
            } else {
                cerr << "Too few arguments for " << argv[i] << ", needs 1 argument." << endl;
                exit(1);
            }
        }
    }
}
void printVersion(char** argv)
{
    if (strcmp(argv[1], "-v") == 0)
    {
            cout << "MAterials Simulation Toolkit for Secondary Electron Emission (MAST-SEY)" << endl;
            cout << "Cite as: Comput. Mater. Sci. 193 (2021), 110281 (https://doi.org/10.1016/j.commatsci.2021.110281)" << endl;
            cerr << "(c) 2021 Maciej P. Polak (mppolak@wisc.edu) & Dane Morgan\n" << endl;
            cout << "Code version "<< code_version << endl;
            exit(0);
    }
    else if (strcmp(argv[1], "-h") == 0)
    {
        cout << "MAterials Simulation Toolkit for Secondary Electron Emission (MAST-SEY)" << endl;
        cout << "Cite as: Comput. Mater. Sci. 193 (2021), 110281 (https://doi.org/10.1016/j.commatsci.2021.110281)" << endl;
        cout << "(c) 2021 Maciej P. Polak (mppolak@wisc.edu) & Dane Morgan\n" << endl;
        cout << "\nInput options:\n" << endl;
        cout << "\"prepare\" as first argument will run input preparation from \"eps/elf.in\" and \"material.in\"" << endl;
        cout << "otherwise, the \"simulate\" version will be executed" << endl;
        cout << "\n\"prepare\" options:" << endl;
        cout << "-e       [iniE(eV,optional) range(eV) grid] energy range and grid (def: 1000 500)" << endl;
        cout << "-lin     generate the energy grid on a linear scale (default is logarithmic)" << endl;
        cout << "-i       [ICS q-int] grids for ICS and q integration (def: 1000 100)" << endl;
        cout << "-qdep    [SPA/SSPA/CUSTOM] specify type of q-dependence of ELF (def: SPA)" << endl;
        cout << "-sumr    output sum rules for plotting" << endl;
        cout << "-saveq   [E_grid q_grid q_max] save q-dependence for plotting" << endl;
        cout << "-elastic [nuclear electron exchange (SOLID LDA opt.)] models to use in elastic scattering (def: F TFM FM)" << endl;
        cout << "         nuclear: [P]oint/[U]niform/[F]ermi" << endl;
        cout << "         electron: [TFM]Thomas–Fermi–Moliere/[TFD]Thomas-Fermi-Dirac/[DHFS]Dirac–Hartree–Fock–Slater/[DF]Dirac-Fock" << endl;
        cout << "         exchange: [NO]/[FM]Furness–McCarthy/[TF]Thomas-Fermi/[RT]Riley–Truhlar" << endl;
        cout << "         (optional): [SOLID] muffin-tin model potential" << endl;
        cout << "         (optional): [LDA] LDA correlation–polarization potential model" << endl;
        cout << "\n\"simulate\" options:" << endl;
        cout << "-e       [incident_energy(eV)] energy of incident energy" << endl;
        cout << "-m       [number_of_e-] number of incident electrons (def: 1000)" << endl;
        cout << "-core    [energy(eV)] allow secondaries to come from bound states" << endl;
        cout << "-dos     [FEG (optional)] generate secondaries from joint DOS from prepared \"jdos.in\"" << endl;
        cout << "         or from parabolic free electron gas approximation (FEG)" << endl;
        cout << "-pa      [angle(deg)] angle of incident electrons with respect to surface normal" << endl;
        cout << "-coord   save travel paths of e-" << endl;
        cout << "-distr   save distribution of secondaries" << endl;
        cout << "-noang   use classical approach to inelastic angle scattering" << endl;
        cout << "-noout   supress all output" << endl;
        cout << "\n-v       display version of the code" << endl;
        cout << "-h       this message\n" << endl;
        cout << "\nPlease be careful when giving input arguments, there is no extensive input checks" << endl;
        cout << "Example executions:" << endl;
        cout << "./mast_sey prepare -e 700 1000 -i 200 100 -elastic F TFD TF" << endl;
        cout << "./mast_sey -e 350 -m 10000" << endl;
        exit(0);
    }
}

void readMaterialFile(string filename)
{
    ifstream infile(filename);
    if(!infile)
    {
        cerr << "MAterials Simulation Toolkit for Secondary Electron Emission (MAST-SEY)" << endl;
        cerr << "Cite as: Comput. Mater. Sci. 193 (2021), 110281 (https://doi.org/10.1016/j.commatsci.2021.110281)" << endl;
        cerr << "(c) 2021 Maciej P. Polak (mppolak@wisc.edu) & Dane Morgan\n" << endl;
        cerr << "Cannot open obligatory file " << filename << endl;
        cerr << "The file should have, line by line, the following information:" << endl;
        cerr << "Atomic Number, Unit Cell Volume [A^3], Fermi energy [eV], Work Function [eV]\n" << endl;

        exit(1);
    }
    int c_mat=1;
    double atnum_tmp;
    while ( infile.peek() != '\n' )
    {
        infile >> atnum_tmp;
        if (c_mat%2==0)
        {
            atcomp.push_back(atnum_tmp);
        } else {
            atnum.push_back(atnum_tmp);
        }
        c_mat=c_mat+1;
    }
    atnum_tmp=0.0;
    if (atcomp.size() > 0)
    {
        for (size_t i = 0; i < atcomp.size(); i++)
        {
            atnum_tmp=atnum_tmp+atcomp[i];
        }
    } else {
            atcomp.push_back(1.0);
    }
    if (abs(atnum_tmp-1.0)>0.0001 && c_mat>2) {
        cerr << "The sum of compositions does not equal 1. The first line of material.in should be:" << endl;
        cerr << "AtomicNumber_1 Composition_1 AtomicNumber_2 Composition_2 .... AtNum_N Comp_N" << endl;
        cerr << "Where Comp1+Comp2+...+CompN = 1.0" << endl;
        exit(1);
    }
    if (ins) {
        infile >> vol >> eg >> evb >> chi;
        eg = EV2HA*eg;
        evb = EV2HA*evb;
        chi = EV2HA*chi;
        u0 = eg+evb+chi;
        ebeg = eg+evb+1e-4;
    } else {
        infile >> vol >> ef >> wf;
        ef = EV2HA*ef;
        wf = EV2HA*wf;
        u0 = ef+wf;
        ebeg = ef+1e-4;
    }
    vol = vol*ANG2BOHR*ANG2BOHR*ANG2BOHR;
}

void readPhononFile(string filename)
{
    ifstream infile(filename);
    double ph_tmp;
    while ( infile.peek() != '\n' )
    {
        infile >> ph_tmp;
        ph_loss.push_back(ph_tmp);
    }
    infile >> eps_zero >> eps_inf;
}

void readEpsFile(string filename)
{
    ifstream epsfile("eps.in");
    ifstream elffile("elf.in");
    bool eps = true;
    if(epsfile && !elffile)
    {
        elfeps = "#  Using dielectric function from file \"eps.in\"";
        filename = "eps.in";
    } else if(!epsfile && elffile)
    {
        elfeps = "#  Using energy loss function function from file \"elf.in\"";
        filename = "elf.in";
        eps = false;
    } else if(epsfile && elffile)
    {
        elfeps = "# WARNING: Found both \"eps.in\" and \"elf.in\", using dielectric function";
        filename = "eps.in";
    } else
    {
        cerr << "Cannot find \"eps.in\" or \"elf.in\" files" << endl; exit(1);
    }   
    ifstream infile(filename);
    if(!infile) {cerr << "Cannot open file: " << filename << endl; exit(1);}
    double e, re, im;
    vector<array<double,2> > ene_elf_tmp;
    if (eps)
    {
        while (infile >> e >> re >> im)
        {
            ene_elf_tmp.push_back({e*EV2HA,im/(re*re+im*im)});
        }
        ene_elf.push_back(ene_elf_tmp);
    } else
    {
        if (qdep == 0)
        {
            while (infile >> e >> re)
            {
                ene_elf_tmp.push_back({e*EV2HA,re});
                if (e == 99999999 && re == 99999999)
                {
                    ene_elf_tmp.pop_back();
                    q_arr.push_back(ene_elf_tmp[ene_elf_tmp.size()-1][0]*HA2EV*BOHR2ANG);
                    ene_elf_tmp.pop_back();
                    ene_elf.push_back(ene_elf_tmp);
                    ene_elf_tmp.clear();
                }
            }
        } else {
            while (infile >> e >> re)
            {
                ene_elf_tmp.push_back({e*EV2HA,re});
                if (e == 99999999 && re == 99999999)
                {
                    ene_elf_tmp.pop_back();
                    ene_elf_tmp.pop_back();
                    break;
                }
            }
            ene_elf.push_back(ene_elf_tmp);
        }
    }
}

void checkSumRules(bool save_sumr)
{
    print("# Sum rules for q=0 Energy Loss Function:");
    vector<array<double,2> > cumsum_sumr2;
    vector<array<double,2> > cumsum_sumr3;
    double sum2 = 0.0;
    double sum3 = 0.0;
    double de,a2,a3,b2,b3;
    for (size_t i = 2; i < ene_elf[0].size(); i++)
    {
        a2 = ene_elf[0][i][1]*ene_elf[0][i][0];
        b2 = ene_elf[0][i-1][1]*ene_elf[0][i-1][0];
        a3 = ene_elf[0][i][1]*1./ene_elf[0][i][0];
        b3 = ene_elf[0][i-1][1]*1./ene_elf[0][i-1][0];
        de = ene_elf[0][i][0]-ene_elf[0][i-1][0];
        sum2 += (a2+b2)*de/2.;
        sum3 += (a3+b3)*de/2.;
        cumsum_sumr2.push_back({ene_elf[0][i][0],sum2});
        cumsum_sumr3.push_back({ene_elf[0][i][0],sum3});
    }
    if (save_sumr)
    {
        saveVector(cumsum_sumr2,"sumrule2.out");
        saveVector(cumsum_sumr3,"sumrule3.out");
    }
    print("# Neff = int(ELF*w*dw) = ",cumsum_sumr2[cumsum_sumr2.size()-1][1]);
    print("# int(ELF*1/w*dw)-pi/2 = ",cumsum_sumr3[cumsum_sumr3.size()-1][1]-PI/2.);
    print("#");
}

vector<vector<array<double,2> > > readScatteringFile(string filename)
{
    ifstream infile(filename);
    if(!infile) {cerr << "Cannot open file: " << filename << endl; exit(1);}
    double x,y;
    vector<array<double,2> > arr2d;
    vector<vector<array<double,2> > > arr3d;
    while (infile >> x >> y)
    {
        arr2d.push_back({x,y});
        if (x == 99999999 && y == 99999999)
        {
            arr2d.pop_back();
            arr3d.push_back(arr2d);
            arr2d.clear();
        }
    }
    return arr3d;
}

vector<double> read1colFile(string filename)
{
    ifstream infile(filename);
    if(!infile) {cerr << "Cannot open file: " << filename << endl; exit(1);}
    double x;
    vector<double> arr1d;
    while (infile >> x)
    {
        arr1d.push_back(x);
    }
    return arr1d;
}

vector<array<double,2> > read2colFile(string filename)
{
    ifstream infile(filename);
    if(!infile) {cerr << "Cannot open file: " << filename << endl; exit(1);}
    double x,y;
    vector<array<double,2> > arr2d;
    while (infile >> x >> y)
    {
        arr2d.push_back({x*EV2HA,y});
    }
    return arr2d;
}

void prepareJDOS(const vector<array<double,2> > &dos)
{
    vector<array<double,2> > arr2d;
    vector<array<double,2> > arr2dint;
    vector<vector<array<double,2> > > arr3d;
    arr2d.reserve(200);
    arr2dint.reserve(200);
    double d_ev = ef/200.;
    for (size_t di = 0; di < de_arr.size(); di++)
    {
        for (int n_ev = 0; n_ev < 200; n_ev++)
        {
            arr2d.push_back({n_ev*d_ev,linterp(de_arr[di]+n_ev*d_ev,dos)*linterp(n_ev*d_ev,dos)});
            arr2dint = cumintVect(arr2d);
        }
        arr3d.push_back(arr2dint);
        arr2d.clear();
        arr2dint.clear();
    }
    saveVector(arr3d,"jdos.in");                                                                                                                      
}

void printInput()
{
    cout << "# MAterials Simulation Toolkit for Secondary Electron Emission (MAST-SEY)" << endl;
    cout << "# Cite as: Comput. Mater. Sci. 193 (2021), 110281 (https://doi.org/10.1016/j.commatsci.2021.110281)" << endl;
    cout << "# (c) 2021 Maciej P. Polak (mppolak@wisc.edu) & Dane Morgan\n" << endl;
    cout << "# Code version " << code_version << "\n#" << endl;
    cout << "# Job started on " << getTime() << "\n#" << endl;
    cout << "#  Input feedback:\n#" << endl;
    if (prep)
    {
        cout << elfeps << "\n#" << endl;
        if (atnum.size()>1)
        {
            cout << "#" << "  Alloy: ";
            for (size_t i = 0; i < atnum.size(); i++)
            {
                cout << elems[atnum[i]] << atcomp[i];
            }
            cout << endl << "#" << endl;
            if (ins) {
                cout << "#" << setw(9) << "Volume" << setw(9) << "Eg" << setw(9) << "Evb" << endl;
                cout << "#" << setw(9) << "[A^3]" << setw(9) << "[eV]" << setw(9) << "[eV]" << endl;
                cout << "#" << setw(9) << vol*BOHR2ANG*BOHR2ANG*BOHR2ANG << setw(9) << eg*HA2EV << setw(9) << evb*HA2EV << setw(9) << "\n#" << endl;
            } else {
                cout << "#" << setw(9) << "Volume" << setw(9) << "EFermi" << endl;
                cout << "#" << setw(9) << "[A^3]" << setw(9) << "[eV]" << endl;
                cout << "#" << setw(9) << vol*BOHR2ANG*BOHR2ANG*BOHR2ANG << setw(9) << ef*HA2EV << setw(9) << "\n#" << endl;
            }
        } else {
            if (ins) {
                cout << "#" << setw(9) << "Elem" << setw(9) << "Volume" << setw(9) << "Eg" << setw(9) << "Evb" << endl;
                cout << "#" << setw(9) << " " << setw(9) << "[A^3]" << setw(9) << "[eV]" << setw(9) << "[eV]" << endl;
                cout << "#" << setw(9) << elems[atnum[0]] << setw(9) << vol*BOHR2ANG*BOHR2ANG*BOHR2ANG << setw(9) << eg*HA2EV << setw(9) << evb*HA2EV << setw(9) << "\n#" << endl;
            } else {
                cout << "#" << setw(9) << "Elem" << setw(9) << "Volume" << setw(9) << "EFermi" << endl;
                cout << "#" << setw(9) << " " << setw(9) << "[A^3]" << setw(9) << "[eV]" << endl;
                cout << "#" << setw(9) << elems[atnum[0]] << setw(9) << vol*BOHR2ANG*BOHR2ANG*BOHR2ANG << setw(9) << ef*HA2EV << setw(9) << "\n#" << endl;
            }
        }
    } else
    {
        if (eb>0.01)
        {
            if (ins) {
                cout << "#" << setw(9) << "Eg" << setw(9) << "Evb" << setw(9) << "Affinity" << setw(9) << "BindEne" << setw(9) << "IncAng "<< endl;
                cout << "#" << setw(9) << "[eV]" << setw(9) << "[eV]" << setw(9) << "[eV]" << setw(9) << "[eV]" << setw(9) << "[deg]" << endl;
                cout << "#" << setw(9) << eg*HA2EV << setw(9) << evb*HA2EV << setw(9) << chi*HA2EV << setw(9) << eb*HA2EV << setw(9) << ini_angle*180./PI << "\n#" << endl;
            } else {
                cout << "#" << setw(9) << "EFermi" << setw(9) << "WorkFun" << setw(9) << "BindEne" << setw(9) << "IncAng "<< endl;
                cout << "#" << setw(9) << "[eV]" << setw(9) << "[eV]" << setw(9) << "[eV]" << setw(9) << "[deg]" << endl;
                cout << "#" << setw(9) << ef*HA2EV << setw(9) << wf*HA2EV << setw(9) << eb*HA2EV << setw(9) << ini_angle*180./PI << "\n#" << endl;
            }
        } else
        {
            if (ins) {
                cout << "#" << setw(9) << "Eg" << setw(9) << "Evb" << setw(9) << "Affinity" << setw(9) << "IncAng" << endl;
                cout << "#" << setw(9) << "[eV]" << setw(9) << "[eV]" << setw(9) << "[eV]" << setw(9) << "[deg]" << endl;
                cout << "#" << setw(9) << eg*HA2EV << setw(9) << evb*HA2EV << setw(9) << chi*HA2EV << setw(9) << ini_angle*180./PI << "\n#" << endl;
            } else {
                cout << "#" << setw(9) << "EFermi" << setw(9) << "WorkFun" << setw(9) << "IncAng" << endl;
                cout << "#" << setw(9) << "[eV]" << setw(9) << "[eV]" << setw(9) << "[deg]" << endl;
                cout << "#" << setw(9) << ef*HA2EV << setw(9) << wf*HA2EV << setw(9) << ini_angle*180./PI << "\n#" << endl;
            }
        }
    }

    if (prep)
    {
        cout << "#" << setw(9) << "Eini[eV]" << setw(9) << "Eend[eV]" << setw(9) << "E-grid" << setw(9) << "ICS-int" << setw(9) << "q-int" << setw(9) << "q-dep" << endl;
        cout << "#" << setw(9) << ebeg*HA2EV << setw(9) << erange*HA2EV << setw(9) << egrid << setw(9) << icsintgrid << setw(9) << qintgrid << setw(9) << qdepname << "\n#" << endl;
        if (elsepa)
        {
            cout << "#" << setw(9) << "ElastSc" << setw(9) << "Nucl." << setw(9) << "Elect." << setw(9) << "Exch.";
            if (es_muffin > 0) { cout << setw(9) << "Solid"; }
            if (es_mcpol > 0) { cout << setw(9) << "CorPol"; }
            cout << endl;
            cout << "#" << setw(9) << "elsepa" << setw(9) << es_nuc << setw(9) << es_el << setw(9) << es_ex;
            if (es_muffin > 0) { cout << setw(9) << es_muffin; }
            if (es_mcpol > 0) { cout << setw(9) << es_mcpol; }
            cout << endl;
        } else {
            cout << "#" << setw(9) << "ElastSc" << setw(9) << "Elect." << endl;
            cout << "#" << setw(9) << "radial" << setw(9) << es_el << endl;
        }
    }
    else
    {
        cout << "#" << setw(9) << "Inc ene" << setw(9) << "No. e-" << endl;
        cout << "#" << setw(9) << erange*HA2EV << setw(9) << mc_elec << endl;
    }
    cout << "#" << endl;
}

vector<double> range(double a, double b, int n)
{
    vector<double> range;
    double d = (b-a)/n;
    for (int i = 0; i <= n; i++)
    {
        range.push_back(a+i*d);
    }
    return range;
}

vector<double> logrange(double a, double b, int n)
{
    vector<double> range;
    a = log10(a);
    b = log10(b);
    double d = (b-a)/(n-1);
    for (int i = 0; i < n; i++)
    {
        range.push_back(pow(10,a+i*d));
    }
    return range;
}

vector<double> getDeArr()
{
    vector<double> dearr;
    vector<double> dearr1 = range(0.00,0.75,100);
    vector<double> dearr2 = range(0.751,dos_arr[dos_arr.size()-1][0],100);
    dearr.reserve(dearr1.size()+dearr2.size());
    dearr.insert(dearr.end(), dearr1.begin(), dearr1.end());
    dearr.insert(dearr.end(), dearr2.begin(), dearr2.end());
    return dearr;
}

void saveVector(vector<double> arr, string filename)
{
    ofstream outfile(filename);
    for (size_t i = 0; i < arr.size(); i++)
    {
        outfile << setprecision(17) << arr[i] << endl;
    }
}
void saveVector(vector<array<double,2> > arr, string filename)
{
    ofstream outfile(filename);
    for (size_t i = 0; i < arr.size(); i++)
    {
        outfile << setprecision(17) << arr[i][0] << " " << setprecision(17) << arr[i][1] << endl;
    }
}
void saveVector(vector<array<double,3> > arr, string filename)
{
    ofstream outfile(filename);
    for (size_t i = 0; i < arr.size(); i++)
    {
        outfile << setprecision(17) << arr[i][0] << " " << setprecision(17) << arr[i][1] << " " << setprecision(17) << arr[i][2] << endl;
    }
}
void saveVector(vector<vector<double> > arr, string filename, int ncols)
{
    ofstream outfile(filename);
    if (ncols==2)
    {
        for (size_t i = 0; i < arr.size(); i++)
        {
            outfile << setprecision(17) << arr[i][0] << " " << setprecision(17) << arr[i][1] << endl;
        }
    } else if (ncols==3)
    {
        for (size_t i = 0; i < arr.size(); i++)
        {
            outfile << setprecision(17) << arr[i][0] << " " << setprecision(17) << arr[i][1] << " " << setprecision(17) << arr[i][2] << endl;
        }
    } else if (ncols==6)
    {
        for (size_t i = 0; i < arr.size(); i++)
        {
            outfile << setprecision(17) << arr[i][0] << " " << setprecision(17) << arr[i][1] << " " << setprecision(17) << arr[i][2] << " " << setprecision(17) << arr[i][3] << " " << setprecision(17) << arr[i][4] << " " << arr[i][5] << endl;
        }
    } else { cerr << "Something went wrong with saveVector" << endl; exit(1); }
}
void saveVector(vector<vector<array<double,2> > > arr, string filename)
{
    ofstream outfile(filename);
    for (size_t i = 0; i < arr.size(); i++)
    {
        for (size_t j = 0; j < arr[0].size(); j++)
        {
            outfile << setprecision(17) << arr[i][j][0] << " " << arr[i][j][1] << endl;
        }
        outfile << "99999999 99999999" << endl;
    }
}
void saveVector(vector<vector<vector<array<double,2> > > > arr, string filename)
{
    ofstream outfile(filename);
    for (size_t i = 0; i < arr.size(); i++)
    {
        for (size_t j = 0; j < arr[0].size(); j++)
        {
            for (size_t k = 0; k < arr[0].size(); k++)
            {
                outfile << setprecision(17) << arr[i][j][k][0] << " " << arr[i][j][k][1] << endl;
            }
            outfile << "11111111 11111111" << endl;
        }
        outfile << "99999999 99999999" << endl;
    }
}

void saveMFP(string filename)
{
    ofstream outfile(filename);
    if (emfp_only)
    {
        outfile << "#Energy[eV] EMFP[A]" << endl;
    } else if (ph) {
        outfile << "#Energy[eV] IMFP[A] EMFP[A] PHMFP_PLUS[A] PHMFP_MINUS[A]" << endl;
    } else {
        outfile << "#Energy[eV] IMFP[A] EMFP[A]" << endl;
    }
    for (size_t ee = 0; ee < ie_arr.size(); ee++)
    {
        if (emfp_only)
        {
            double iemfp_e = linterp2d(ie_arr[ee],-1,ie_arr,elas_arr,true)/vol;
            outfile << setprecision(17) << ie_arr[ee]*HA2EV << " " << setprecision(17) << 1./iemfp_e*BOHR2ANG << endl;
        } else {
            double iemfp_e = linterp2d(ie_arr[ee],-1,ie_arr,elas_arr,true)/vol;
            double iimfp_e = linterp2d(ie_arr[ee],-1,ie_arr,inel_arr,true);
            if (ph)
            {
                double iphmfp_plus_e = linterp2d(ie_arr[ee],-1,ie_arr,phon_plus_arr,true);
                double iphmfp_minus_e = linterp2d(ie_arr[ee],-1,ie_arr,phon_minus_arr,true);
                outfile << setprecision(17) << ie_arr[ee]*HA2EV << " " << setprecision(17) << 1./iimfp_e*BOHR2ANG << " " << setprecision(17) << 1./iemfp_e*BOHR2ANG << " " << setprecision(17) << 1./iphmfp_plus_e*BOHR2ANG << " " << setprecision(17) << 1./iphmfp_minus_e*BOHR2ANG << endl;
            }
            else {
                outfile << setprecision(17) << ie_arr[ee]*HA2EV << " " << setprecision(17) << 1./iimfp_e*BOHR2ANG << " " << setprecision(17) << 1./iemfp_e*BOHR2ANG << endl;
            }
        }
    }
}

void saveQdep(int n_ene, int n_q, double q_max)
{
    double de = ie_arr[ie_arr.size()-1]/n_ene;
    double dq = q_max/n_q;
    double ee;
    double qq = 0;
    int progress = 0;
    vector<array<double,3> > elfq_arr;
    printProgress();
    for (int i = 0; i <= n_q; i++)
    {
        progress = printStars(progress,i,n_q);
        ee = 0;
        for (int j = 0; j <= n_ene; j++)
        {
            elfq_arr.push_back({ee*HA2EV,qq,elfq(qq,ee,0)});
            ee += de;
        }
        qq += dq;
        elfq_arr.push_back({99999999,99999999,99999999});
    }
    saveVector(elfq_arr,"elf_qdep.out");
}

void saveCoordVector(vector<vector<array<double,3> > > arr, vector<int> second, string filename)
{
    ofstream outfile(filename);
    for (size_t i = 0; i < arr.size(); i++)
    {
        outfile << second[i] << endl;
        for (size_t j = 0; j < arr[i].size(); j++)
        {
            outfile << setprecision(7) << arr[i][j][0] << " " << arr[i][j][1] << " " << arr[i][j][2] << endl;
        }
        outfile << endl;
    }
}

vector<array<double,2> > cumintVect(const vector<array<double,2> > &xyarr)
{
    double intgrl = 0.;
    double a,b,dx;
    vector<array<double,2> > int_arr;
    int_arr.push_back({0.,0.});
    for (size_t i = 0; i < xyarr.size()-1; i++)
    {
        a = xyarr[i][1];
        b = xyarr[i+1][1];
        dx = xyarr[i+1][0]-xyarr[i][0];
        intgrl += (a+b)*dx*0.5;
        int_arr.push_back({xyarr[i+1][0],intgrl});
    }
    return int_arr;
}

array<double,3> f_rotdircos(array<double,3> uvw, double ang0, double ang1)
{
    array<double,3> base_uvw{0.,0.,0.};
    if (dircos_algo == 0) {
    double sdt = sin(ang0);
    double cdt = cos(ang0);
    double sdf = sin(ang1);
    double cdf = cos(ang1);
    double dxy = uvw[0]*uvw[0]+uvw[1]*uvw[1];
    double dxyz = dxy+uvw[2]*uvw[2];
    if (abs(dxyz-1.0) > 1e-9) {
        double fnorm = 1.0/sqrt(dxyz);
        uvw[0] = fnorm*uvw[0];
        uvw[1] = fnorm*uvw[1];
        uvw[2] = fnorm*uvw[2];
        dxy = uvw[0]*uvw[0]+uvw[1]*uvw[1];
    }
    if (dxy > 1e-9) {
        sdt = sqrt((1.0-cdt*cdt)/dxy);
        double up = uvw[0];
        base_uvw[0] = uvw[0]*cdt+sdt*(up*uvw[2]*cdf-uvw[1]*sdf);
        base_uvw[1] = uvw[1]*cdt+sdt*(uvw[1]*uvw[2]*cdf+up*sdf);
        base_uvw[2] = uvw[2]*cdt-dxy*sdt*cdf;
    } else {
        sdt = sqrt(1.0-cdt*cdt);
        base_uvw[1] = sdt*sdf;
        if (uvw[2]>0.0) {
            base_uvw[0] = sdt*cdf;
            base_uvw[2] = cdt;
        } else {
            base_uvw[0] = -1.0*sdt*cdf;
            base_uvw[2] = -1.0*cdt;
        }
    }
    }

    if (dircos_algo == 1) {

    if (abs(uvw[2])>0.99) {
        base_uvw[0]=sin(ang0)*cos(ang1);
        base_uvw[1]=sin(ang0)*sin(ang1);
        base_uvw[2]=cos(ang0);
    } else {
        double sq1=sqrt(1-pow(cos(ang0),2));
        double sqw=sqrt(1-pow(uvw[2],2));
        base_uvw[0]=uvw[0]*cos(ang0)+sq1/sqw*(uvw[0]*uvw[2]*cos(ang1)-uvw[1]*sin(ang1));
        base_uvw[1]=uvw[1]*cos(ang0)+sq1/sqw*(uvw[1]*uvw[2]*cos(ang1)+uvw[0]*sin(ang1));
        base_uvw[2]=uvw[2]*cos(ang0)-sq1*sqw*cos(ang1);
    }
    }

    if (dircos_algo == 2) {

    double th0 = acos(uvw[2]);
    double ph0 = atan2(uvw[1],uvw[0]);
    double theta = acos(cos(th0)*cos(ang0)-sin(th0)*sin(ang0)*cos(ang1));
    double phi = asin(sin(ang0)*sin(ang1)/sin(theta))+ph0;

    base_uvw[0]=sin(theta)*cos(phi);
    base_uvw[1]=sin(theta)*sin(phi);
    base_uvw[2]=cos(theta);
    }

    return base_uvw;
}

vector<array<double,2> > int_elastic_ang(const vector<double> &xarr, const vector<double> &yarr, bool cumint)
{
    double intgrl = 0.;
    double a,b,dx;
    vector<array<double,2> > int_arr;
    int_arr.push_back({0.,0.});
    for (size_t i = 0; i < xarr.size()-1; i++)
    {
        a = yarr[i]*sin(xarr[i]);
        b = yarr[i+1]*sin(xarr[i+1]);
        dx = xarr[i+1]-xarr[i];
        intgrl += (a+b)*dx*0.5;
        if (cumint) { int_arr.push_back({xarr[i+1],2.*PI*intgrl}); }
    }
    if (cumint) { return int_arr; }
    else
    {
        int_arr.clear();
        int_arr.push_back({xarr[xarr.size()],2.*PI*intgrl});
        return int_arr;
    }
}

vector<array<double,2> > int_inelastic_ene(double (*f)(double,double,int), double x0, double x1, int ndiv, double args, bool cumint)
{
    double x = x0;
    double dx = (x1-x0)/ndiv;
    double a,b;
    double intgrl = 0.;
    vector<array<double,2> > int_arr;
    int_arr.push_back({0.,0.});
    for (int i = 0; i < ndiv; i++)
    {
        a = (*f)(x,args,1);
        b = (*f)(x+dx,args,1);
        intgrl += (a+b)*dx*0.5;
        x = x+dx;
        if (cumint) { int_arr.push_back({x,intgrl}); }
    }
    if (cumint) { return int_arr; }
    else
    {
        int_arr.clear();
        int_arr.push_back({x1,intgrl});
        return int_arr;
    }
}

vector<array<double,2> > int_inelastic_ang(double (*f)(double,double,int), double x0, double x1, int ndiv, double ie, double de, bool cumint)
{
    double x = x0;
    double dx = (x1-x0)/ndiv;
    double q2,q,a,b;
    double intgrl = 0.;
    vector<array<double,2> > int_arr;
    int_arr.push_back({0.,0.});
    for (int i = 0; i < ndiv; i++)
    {
        if (de<1e-10)
        {
            x = x+dx;
            if (cumint) { int_arr.push_back({x,0.0}); }
        }
        q2 = 4.*ie-2.*de-4.*sqrt(ie*(ie-de))*cos(x);
        q = sqrt(q2);
        a = 1./(PI*PI*ie)*(*f)(q,de,2)*sqrt(ie*(ie-de))*sin(x);

        q2 = 4.*ie-2.*de-4.*sqrt(ie*(ie-de))*cos(x+dx);
        q = sqrt(q2);
        b = 1./(PI*PI*ie)*(*f)(q,de,2)*sqrt(ie*(ie-de))*sin(x+dx);

        intgrl += (a+b)*dx*0.5;
        x = x+dx;
        if (cumint) { int_arr.push_back({x,intgrl}); }
    }
    if (cumint) { return int_arr; }
    else
    {
        int_arr.clear();
        int_arr.push_back({x1,intgrl});
        return int_arr;
    }
}

vector<array<double,2> > inel(double ie)
{
    vector<array<double,2> > iimfpint;
    if (ie<1e-10)
    {
        for (int i = 0; i <= icsintgrid; i++)
        {
            iimfpint.push_back({0.0,0.0});
        }
    }
    else
    {
        if (ins)
        {
            iimfpint = int_inelastic_ene(&qIntFun,eg,ie-eg-evb,icsintgrid,ie,true);
        }
        else
        {
            iimfpint = int_inelastic_ene(&qIntFun,1e-5,ie-ef,icsintgrid,ie,true);
        }
    }
    return iimfpint;
}

vector<array<double,2> > phonon_plus(double ie)
{
    double t = 300.0;
    double k_b = 8.617e-5;
    double n_lo, ln_plus, term_plus;
    vector<array<double,2> > iphmfp_plus;
    if (ie<1e-10)
    {
        iphmfp_plus.push_back({0.0,0.0});
    }
    else
    {
        for (size_t i = 0; i < ph_loss.size(); i++)
        {
            n_lo = 1.0/(exp(ph_loss[i]/(k_b*t))-1.0);
            ln_plus = (1.0 + sqrt(abs(1.0 - ph_loss[i]/ie))) / (1.0 - sqrt(abs(1.0 - ph_loss[i]/ie)));
            term_plus = (n_lo + 1.0)*(1.0/eps_inf - 1.0/eps_zero)*ph_loss[i]/(2.0*ie);
            iphmfp_plus.push_back({ph_loss[i], term_plus*log(ln_plus)});
        }
    }
    return iphmfp_plus;
}

vector<array<double,2> > phonon_minus(double ie)
{
    double t = 300.0;
    double k_b = 8.617e-5;
    double n_lo, ln_minus, term_minus;
    vector<array<double,2> > iphmfp_minus;
    if (ie<1e-10)
    {
        iphmfp_minus.push_back({0.0,0.0});
    }
    else
    {
        for (size_t i = 0; i < ph_loss.size(); i++)
        {
            n_lo = 1.0/(exp(ph_loss[i]/(k_b*t))-1.0);
            ln_minus = (1.0 + sqrt(abs(1.0 + ph_loss[i]/ie))) / (-1.0 + sqrt(abs(1.0 + ph_loss[i]/ie)));
            term_minus = n_lo*(1.0/eps_inf - 1.0/eps_zero)*ph_loss[i]/(2.0*ie);
            iphmfp_minus.push_back({ph_loss[i], term_minus*log(ln_minus)});
        }
    }
    return iphmfp_minus;
}

vector<vector<array<double,2> > > inelang(double ie)
{
    vector<array<double,2> > inangint;
    vector<vector<array<double,2> > > inangint_de;
    int nde = 100;
    double de;
    for (int d = 0; d <= nde; d++)
    {
        if (ins) {
            de = eg + (double)d/(double)nde*(ie-eg-evb);
        } else {
            de = (double)d/(double)nde*(ie-ef);
        }
        inangint = int_inelastic_ang(&elfq,0.0,PI/2.,nde,ie,de,true);
        inangint_de.push_back(inangint);
        inangint.clear();
    }
    return inangint_de;
}

vector<array<double,2> > elas(double ie, int at, double comp)
{
    vector<array<double,2> > elasint;
    elasint.reserve(606);
    if (ie<1e-10)
    {
        for (int i = 0; i < 606; i++)
        {
            elasint.push_back({0.0,0.0});
        }
    }
    else
    {
        if (elsepa)
        {
            string command = "getDDCS "+to_string(es_nuc)+" "+to_string(es_el)+" "+to_string(es_ex)+" "+to_string(es_muffin)+" "+to_string(es_mcpol)+" "+to_string(at)+" "+to_string(ie*HA2EV)+" "+to_string(1);
            system(command.c_str());
        }
        else
        {
            if (es_el == 3) { es_el = 1; }
            else if (es_el == 1) { es_el = 2; }
            else { es_el = 1; }
            string command = "getDDCS "+to_string(es_nuc)+" "+to_string(es_el)+" "+to_string(es_ex)+" "+to_string(at)+" "+to_string(ie*HA2EV)+" "+to_string(0);
            system(command.c_str());
        }
        ifstream infile("elast_cs.dat");
        if(!infile) {cerr << "Cannot open file: elast_cs.dat" << endl; exit(1);}
        
        vector<double> ang, ddcs;
        double a, ecs;
        while (infile >> a >> ecs)
        {
            ang.push_back(a);
            ddcs.push_back(ecs*comp);
        }
        elasint = int_elastic_ang(ang,ddcs,true);
    }
    return elasint;
}

double linterp(double x, const vector<array<double,2> > &xyarr, bool xfind)
{
    for (size_t i = 0; i < xyarr.size()-1; i++)
    {
        if (xfind)
        {
            if(x<xyarr[i+1][1] && x>=xyarr[i][1])
            {
                return xyarr[i][0]+(xyarr[i+1][0]-xyarr[i][0])*(x-xyarr[i][1])/(xyarr[i+1][1]-xyarr[i][1]);
                break;
            }
        }
        else 
        {
            if(x<xyarr[i+1][0] && x>=xyarr[i][0])
            {
                return xyarr[i][1]+(xyarr[i+1][1]-xyarr[i][1])*(x-xyarr[i][0])/(xyarr[i+1][0]-xyarr[i][0]);
                break;
            }
        }
    }
    return 0.;
}

double linterp2d(double x0,double y0, const vector<double> &xarr, const vector<vector<array<double,2> > > &ytuparr, bool total, bool xfind)
{
    vector<array<double,2> > ytup_intrp;
    if (!total)
    {
        ytup_intrp.reserve(ytuparr[0].size());
    }
    double dx;
    for (size_t i = 0; i < xarr.size()-1; i++)
    {
        if (x0<=xarr[i+1] && x0>=xarr[i])
        {
            dx = (x0-xarr[i])/(xarr[i+1]-xarr[i]);
            if (total)
            {
                return ytuparr[i][ytuparr[i].size()-1][1]+(ytuparr[i+1][ytuparr[i+1].size()-1][1]-ytuparr[i][ytuparr[i].size()-1][1])*dx;
                break;
            }
            else
            {
                for (size_t j = 0; j < ytuparr[0].size(); j++)
                {
                    ytup_intrp.push_back({ytuparr[i][j][0]+(ytuparr[i+1][j][0]-ytuparr[i][j][0])*dx,ytuparr[i][j][1]+(ytuparr[i+1][j][1]-ytuparr[i][j][1])*dx});
                }
            }
            break;
        }
    }
    return linterp(y0, ytup_intrp, xfind);
}

double fzero(double (*f)(double,double,double), double x0, double x1, double ww, double qq, double tol)
{
    double xm=0.0;
    if ((*f)(x0,ww,qq)*(*f)(x1,ww,qq)<0)
    {
        while (abs(x1-x0)>tol)
        {
            xm = (x0+x1)/2.; 
            if ((*f)(xm,ww,qq)==0.0) { return xm; }
            else if ((*f)(xm,ww,qq)*(*f)(x0,ww,qq)<0)
            {
                x1 = xm;
            } else
            {
                x0 = xm;
            }
        }    
    }
    return xm;
}

double sspa_elf(double w, double q)
{
    return (w-(q*q)/2.)/w*linterp(w-(q*q)/2.,ene_elf[0]);
}

double spa_dispers(double w0, double w, double q)
{
    return w-sqrt(w0*w0+(1./3.)*pow((q*pow(0.75*PI,1./3.)*pow(w0,2./3.)),2)+pow(q,4)/4.);
}
double spa_elf(double w, double q)
{
    double omz = fzero(&spa_dispers,0.0,50.,w,q);
    return linterp(omz,ene_elf[0])/(1.+(PI*q*q)/(6.*pow(0.75*PI,1./3.)*pow(omz,2./3.)));
}

double fpa_elf(double w, double q)
{
    double integ = 0;
    return integ;
}

double dft_elf(double w, double q)
{
    return linterp2d(q,w,q_arr,ene_elf,false,false);
}

double elfq(double q, double om, int dq)
{
    if (qdep == 2)
    {
        if (dq==1) { return sspa_elf(om,q)/q; }
        else if (dq==2) { return sspa_elf(om,q)/(q*q); }
        else { return sspa_elf(om,q); }
    } else if (qdep == 1)
    {
        if (dq==1) { return spa_elf(om,q)/q; }
        else if (dq==2) { return spa_elf(om,q)/(q*q); }
        else { return spa_elf(om,q); }
    } else if (qdep == 0)
    {
        if (dq==1) { return dft_elf(om,q)/q; }
        else if (dq==2) { return dft_elf(om,q)/(q*q); }
        else { return dft_elf(om,q); }
    } else if (qdep == 3)
    {
        if (dq==1) { return fpa_elf(om,q)/q; }
        else if (dq==2) { return fpa_elf(om,q)/(q*q); }
        else { return fpa_elf(om,q); }
    }
    else {
        cout << "Wrong option in qdep" << endl;
        exit(1);
        return 0.0;
    }
}

double qIntFun(double om, double ie, int dummy)
{
    double tp = ie;
    double qm, qp, c;
    if (ins)
    {
        c = pow((1.+(tp - eg)/(137.*137.)),2)/(1.+(tp - eg)/(2.*137.*137.))*1/(PI*(tp - eg));

        qm = sqrt((tp - eg)*(2.+(tp - eg)/(cc*cc)))-sqrt((tp-eg-om)*(2.+(tp-eg-om)/(cc*cc)));
        qp = sqrt((tp - eg)*(2.+(tp - eg)/(cc*cc)))+sqrt((tp-eg-om)*(2.+(tp-eg-om)/(cc*cc)));
    }
    else
    {
        c = pow((1.+(tp)/(137.*137.)),2)/(1.+tp/(2.*137.*137.))*1/(PI*tp);
    
        qm = sqrt(tp*(2.+tp/(cc*cc)))-sqrt((tp-om)*(2.+(tp-om)/(cc*cc)));
        qp = sqrt(tp*(2.+tp/(cc*cc)))+sqrt((tp-om)*(2.+(tp-om)/(cc*cc)));
    }

    vector<array<double,2> > qint = int_inelastic_ene(&elfq,qm,qp,qintgrid,om);
    return qint[0][1]*c;
}

double jdos(double e, double de, double r)
{
    return (de*sqrt(e*(de+e))+2.*e*sqrt(e*(de+e))-pow(de,2.)*asinh(sqrt(e/de)))/(de*sqrt(ef*(de+ef))+2.*ef*sqrt(ef*(de+ef))-pow(de,2.)*asinh(sqrt(ef/de)))-r;
}

void printVector(vector<double> &arr)
{
    for (size_t i = 0; i < arr.size(); i++)
    {
        cout << setprecision(17) << arr[i] << endl;
    }
}
void printVector(vector<array<double,2> > &arr)
{
    for (size_t i = 0; i < arr.size(); i++)
    {
        cout << arr[i][0] << " " << arr[i][1] << endl;
    }
}
void printVector(vector<vector<array<double,2> > > &arr)
{
    for (size_t i = 0; i < arr.size(); i++)
    {
        for (size_t j = 0; j < arr[0].size(); j++)
        {
            cout << setprecision(17) << arr[i][j][0] << " " << arr[i][j][1] << endl;
        }
        cout << "----------------------------------" << endl;
    }
}

void marker(string num)
{ cout << "HERE_" << num << endl; }
void print(string text)
{ cout << text << endl; }
void print(int num)
{ cout << num << endl; }
void print(double num)
{ cout << num << endl; }
void print(string text, int num)
{ cout << text <<  num << endl; }
void print(string text, double num)
{ cout << text <<  num << endl; }
void printProgress()
{
    print("# PROGRESS:");
    print("# 0%       25%       50%       75%      100%");
    print("# |_________|_________|_________|_________|");
    cout << "#  ";
    cout.flush();
}
int printStars(int progr, int is, int size)
{
    if (progr < is*40/size)
    {
        cout << "*";
        cout.flush();
        progr = is*40/size;
    }
    return progr;
}
double rand01()
{
    uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(mt);
}
string checkName(string name)
{
    int n = 1;
    string new_name = name;
    ifstream f;
    f.open(new_name);
    while (f)
    {
        f.close();
        new_name = name+"_"+to_string(n);
        f.open(new_name);
        n++;
    }
    f.close();
    return new_name;
}
