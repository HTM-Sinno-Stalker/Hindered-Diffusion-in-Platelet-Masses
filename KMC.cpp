#include <iostream>
#include <cmath>
#include <fstream>
#include <unordered_map>
#include <set>
#include <vector>
#include <random>
#include <string>
#include <chrono>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std::chrono;
using namespace std;

typedef double T;

// ***********************************************************************************
// ************************** SECTION 1: INITIALIZATION ******************************
// ***********************************************************************************

// 1.1: INITIALIZE STRUCTURES *****************************************************

T dblmax = std::numeric_limits<T>::max(); // The maximum possible number in c++. Here, this value is functionally infinity

struct particle // A structure that defines information about a particle
{
    int type;      // For if you have many kinds of particles in circulation
    int center[3]; // Position of the particle in x, y, z coords
    T dp;          // Diameter of particle

    T motion[4]; //  Motion rates in +/- x, y, z directions
    T tau[4];    // Time in which each motion above would occue

    particle() // constructor, is called to initialize the structure 
    {
        type = 0;
        center[0] = center[1] = center[2] = 0;
        dp = 0;

        for (int i = 0; i < 4; ++i)
        {
            motion[i] = 0.;
        }

        for (int i = 0; i < 4; ++i)
        {
            tau[i] = dblmax;
        }
    }
};

struct event // A structure which is used to store the events in the sorted event list
{
    int idx; // Location (see unordered_map)
    int type; // motion with some direction, or a reaction 
    T tau; // time for event to occur

    event(int idx_, int type_, T tau_) : idx(idx_), type(type_), tau(tau_)
    {
    }
};

struct orderEvent // Orders events according to increasing tau values (next reaction method)
{
    bool operator()(const event& e1, const event& e2) const
    {
        return e1.tau < e2.tau;
    }
};

// 1.2: INITIALIZE NUMBER VALUES *****************************************************

int SiZe = 0;             // The number of simulations you'd like to run

vector <int> IDrun;       // The ID number on the file output of the simulation
int ident = 0;

vector <int> NumberParticles;        // Number of particles in the system 
int nparticles = 0;

vector <int> FINITEVOL; // Whether or not the particles have a finite volume (1) or not (0)
int FiniteVol = 0;
vector<vector<int>> neigh; // The neighbor list, for finite particles

vector <T> cellThreshold; // For CellPose runs. The % threshold beyond which a node is binarized as a cellspace (i.e. impassable)
T CellThreshold = 0;

vector<T> diameters; // Particle diameters in Base-Units
T dia = 0;

vector <int> SEED;      // The random # seeds for the simulations
int seed = 0;

vector <string> FileOpen;  // For CellPose runs. List of domain files to open in turn for each simulation in a run.
T ImSize = 0;             // The size of the CellPose image domain. Standard is 12.5 Base-Units (for a 1250 x 1250 pixel image)
T subIS = 0;                // the size of subdomain. 
vector<int> bumpx;          // for subdomains. what the starting value is along x.
vector<int> bumpy;          // for subdomains. what the starting value is along y.
int startx;
int starty;

vector<vector<T>> ClotDomain;

T totalTimeNRM = 0.;    // The current simulation time of a run
T simulationTime = 0;   // The desired final simulation time of a run
T dtOutput = 0;          // The simulation time elapsed between recordings of particle positions
vector<T> TIME;         // Vector of time stamps corrsponding to each particle record
vector<T> ExeTime;      // The real time it takes each simulation to complete
int ITER = 0;           // counter of KMC iterations for file-writing purposes
int neighlimit = 0;

vector<particle> particleList; // Instantaneous list of particles
vector<vector<particle>> particleListOverall; // List of particles throughout time

T diffusion[4]; // Diffusive rate of motion in +/- x, y, z directions

std::mt19937 engine(0); // Random number generator
std::uniform_real_distribution<T> distribution(0., 1.);  // set random number range to between 0, 1 

std::unordered_map<long long int, int> occupation; // Keeps track of particle occupation as function of lattice spacing. Returns the index of the particle in particleList.

std::set<event, orderEvent> eventList; // Running list of ordered events for KMC

// 1.3: INITIALIZE FUNCTIONS *****************************************************

T randn(); // Generates random number between 0 and 1
bool Wall(int cx, int cy, int cz); // This function tells you if you have crossed a "wall" (boundary) in the domain
void CheckWithinDomain(int current); // Defines rates such that a move where a particle leaves the domain is not allowed

bool Overlap(int current, int cx, int cy, int cz); // For finite particles; tells you if a particle's position overlaps another
void NeighList(); // Constructs a list of each particle's neighbors in order to make the overlap check more efficient


void iniKMCLattice();// Used to initalize the LKMC lattice with given initial particle distribution
void setRates();// This function sets the rate of each event for each particle (for initialization or finite particles) 
void setRates(int current, int type);// This function sets the rate of each event for a given particle that moves
void runKMC();// This function runs one iteration of the LKMC algorithm
void writeReadMe();
void writeOutputFiles();// Write out the velocity field and particle positions for visualization


// ***********************************************************************************
// ************************** SECTION 2: MAIN FUNCTION ******************************
// ***********************************************************************************

int main(int argc, char* argv[])
{

    // 2.1: SET RUN ARRAYS *****************************************************
    SiZe = 4;

    IDrun = { 0 };

    ImSize = 1250; // number of pixels in image (10mmmm x 10 mmmm), or (10000nm x 10000 nm)
    subIS = 250; // subdomain image size
    bumpx = { 0,125,250,375,500,625,750,875,1000};
    bumpy = { 0,125,250,375,500,625,750,875,1000};

    FINITEVOL = { 0, 1 }; // 0=tracer 1=finite 

    NumberParticles = { 1000 };

    cellThreshold = {20,25}; // Between 0 and 255. 

    SEED = { 163748,283745,827465, 283745, 827465, 827465 }; // Six digits please

    diameters = { 8 };  // nm.

    dia = 1;// ImSize / 10000 * 8;// diameters[mm];
    CellThreshold = 12;// cellThreshold[mm];
    
    //string flowsname = "291_2x_flows.txt";

    FileOpen = { "flowfield_295.txt"};

    simulationTime = 0.003;

    dtOutput = 0.000005;

    ClotDomain.clear();
    ClotDomain.resize(ImSize, vector<T>(ImSize));
    ifstream ClotFile{ FileOpen[0]};
    T temp = 0.0;
    for (int i = 0; i < ImSize; ++i)
    {
        for (int j = 0; j < ImSize; ++j)
        {
            ClotFile >> temp;
            if (temp > CellThreshold)
            {
                temp = 1; // wall space
            }
            else
            {
                temp = 0; // pore space
            }
            ClotDomain[j][ImSize-i-1] = temp;
        }
    }

    // set diffusivity in nodes^2/sec
    for (int j = 0; j < 4; ++j)
    {
        //diffusion[j] = 10000. * (ImSize / 10000.) * (ImSize / 10000.);
        diffusion[j] = 1e8 *(1250 / 10000.)* (1250 / 10000.);
    }

    // 2.2: SIMULATING EACH RUN *****************************************************
    for (int mm = 3; mm < SiZe; ++mm)
    {
        printf("starting up...\n");
        auto start = high_resolution_clock::now(); // start clock of simulation runtime

        // pull simulation values from arrays
        ident = IDrun[0]+mm; //IDrun[mm];
        FiniteVol = 0;// FINITEVOL[mm];
        nparticles = 100;// NumberParticles[mm];
        seed = 564389;//SEED[mm];
        startx = bumpx[mm]; 
        starty = bumpy[mm]; 
        engine.seed(seed);

        // zero out arrays
        totalTimeNRM = 0.;
        particleList.clear();
        particleListOverall.clear();
        TIME.clear();

        neighlimit = 100 ^ 2;// (2 * diffusivity[0] / (h * h) * dtOutput)* (2 * diffusivity[0] / (h * h) * dtOutput);

        iniKMCLattice();// Initialize the lattice with particles
        printf("initial configuration set\n");
        setRates();// Create a rate database of all possible events and set rates and putative times for each possible event
        printf("rates set\n");
        ITER = 1;
        particleListOverall.push_back(particleList);
        TIME.push_back(0);
        writeReadMe();
        writeOutputFiles();

        // Run KMC
        while (totalTimeNRM <= simulationTime)
        {
            runKMC();

            if (totalTimeNRM - TIME.back() >= dtOutput) // if the time since last recording particles = dtOutput
            {
                particleListOverall.push_back(particleList); // record particle positions
                TIME.push_back(totalTimeNRM); // record time
                //system("cls");
                if (ITER%100 ==0)
                {
                    writeOutputFiles();
                    printf(" %f percent complete. \n", (totalTimeNRM + simulationTime * mm) / (simulationTime * SiZe) * 100.0);
                }
                
               /* if (mm > 0)
                {
                    for (int g = 0; g < mm; ++g)
                    {
                        printf(" Test %i results available. Execution time was %f mmms. \n", IDrun[0]+g, ExeTime[g]);
                    }
                }*/
                ++ITER;
                if (FiniteVol == 1)
                {
                    NeighList();
                }
            }
        }
        particleListOverall.push_back(particleList);
        writeOutputFiles();

        auto stop = high_resolution_clock::now(); // stop clock of simulation runtime
        auto duration = duration_cast<minutes>(stop - start); // calculate duration of simulation 
        ExeTime.push_back(T(duration.count())); // Record simulation execution time for this run
          system("cls");
          printf(" %f percent complete\n", 100.0 * (mm + 1) / SiZe);
          for (int g = 0; g < mm + 1; ++g)
          {
              printf(" Test %i results available. Execution time was %f mmms. \n", IDrun[0]+g, ExeTime[g]);
          }

    }

    //std::ofstream Timefile("ExeTime.txt"); // Write simulation execution time file 
    //for (int j = 0; j < SiZe; ++j)
    //{
    //    Timefile << std::to_string(IDrun[j]) << " took " << ExeTime[j] << " minutes to run. " << "\n";
    //}
    //Timefile.close();
    return 0;
}

// ******************************************************************************************
// ************************** SECTION 3:FUNCTION DEFININITIONS ******************************
// ******************************************************************************************

T randn()
{
    return distribution(engine);
}

bool Wall(int cx, int cy, int cz) // Checks if move trips accross wall (true=yes)
{

    // how many domain spaces (DS) has the partcile traversed in each direction?
    int DSx = floor(T(cx) / (subIS));
    int DSy = floor(T(cy) / (subIS));
    // what's the reflection factor (RF) of each dimension? 0 for nonreflected, 1 for reflected
    int RFx = DSx % 2;
    int RFy = DSy % 2;
    // now let's find the particle position in the reference cell. 
    int cxref = cx - DSx * subIS;
    int cyref = cy - DSy * subIS;
    if (RFx > 0)
    {
        cxref = subIS - cxref - 1;
    }
    if (RFy > 0)
    {
        cyref = subIS - cyref - 1;
    }
    if (ClotDomain[int(cxref + startx)][int(cyref + starty)]> 0)
    {
        return 1;
    }
    return 0;
}

void NeighList() // Construct Neighbor List
{
    for (int j = 0; j < particleList.size(); ++j)
    {
        for (int i = 0; i < (particleList.size()); ++i)
        {
            neigh[j][i] = 0;
            neigh[i][j] = 0;

            int rx = (particleList[i].center[0] - particleList[j].center[0]);
            int ry = (particleList[i].center[1] - particleList[j].center[1]);
            int rz = (particleList[i].center[2] - particleList[j].center[2]);
            int rr = rx * rx + ry * ry + rz * rz; // distance between two particles 
            if (rr <= neighlimit)
            {
                neigh[j][i] = 1; // YES a neighbor
                neigh[i][j] = 1;
            }
        }
    }
}

bool Overlap(int current, int cx, int cy, int cz) // checks if particle overlaps any other particles
{
    if (FiniteVol == 1)
    {
        for (int ijk = 0; ijk < particleList.size(); ++ijk)
        {
            if (neigh[current][ijk] == 1 && ijk != current)
            {
                particle ptemp = particleList[ijk];
                T dist = sqrt((cx - ptemp.center[0]) * (cx - ptemp.center[0]) + (cy - ptemp.center[1]) * (cy - ptemp.center[1]) + (cz - ptemp.center[2]) * (cz - ptemp.center[2]));
                if (dist < dia)
                {
                    return true;
                }
                else
                {
                    continue;
                }
            }
            else
            {
                continue;
            }

        }
        return false;
    }
    else
    {
        return false;
    }
}

void CheckWithinDomain(int current) // This function ensures that all the rates are defined such that a move where a particle leaves the domain is not allowed
{
    particle pt = particleList[current];

    int i, j, k;
    i = pt.center[0];
    j = pt.center[1];
    k = pt.center[2];

    if (Overlap(current, i + 1, j, k) || Wall(i + 1, j, k))
    {
        auto iT = eventList.find(event(current, 0, pt.tau[0]));
        if (iT != eventList.end())
        {
            eventList.erase(iT);
        }
        pt.motion[0] = 0.;
        pt.tau[0] = dblmax;
    }
    if (Overlap(current, i - 1, j, k) || Wall(i - 1, j, k))
    {
        auto iT = eventList.find(event(current, 1, pt.tau[1]));
        if (iT != eventList.end())
        {
            eventList.erase(iT);
        }
        pt.motion[1] = 0.;
        pt.tau[1] = dblmax;
    }
    if (Overlap(current, i, j + 1, k) || Wall(i, j + 1, k))
    {
        auto iT = eventList.find(event(current, 2, pt.tau[2]));
        if (iT != eventList.end())
        {
            eventList.erase(iT);
        }
        pt.motion[2] = 0.;
        pt.tau[2] = dblmax;
    }
    if (Overlap(current, i, j - 1, k) || Wall(i, j - 1, k))
    {
        auto iT = eventList.find(event(current, 3, pt.tau[3]));
        if (iT != eventList.end())
        {
            eventList.erase(iT);
        }
        pt.motion[3] = 0.;
        pt.tau[3] = dblmax;
    }
   /* if (Overlap(current, i, j, k + 1) || Wall(i, j, k + 1))
    {
        auto iT = eventList.find(event(current, 4, pt.tau[4]));
        if (iT != eventList.end())
        {
            eventList.erase(iT);
        }
        pt.motion[4] = 0.;
        pt.tau[4] = dblmax;
    }
    if (Overlap(current, i, j, k - 1) || Wall(i, j, k - 1))
    {
        auto iT = eventList.find(event(current, 5, pt.tau[5]));
        if (iT != eventList.end())
        {
            eventList.erase(iT);
        }
        pt.motion[5] = 0.;
        pt.tau[5] = dblmax;
    }*/

    particleList[current] = pt;
}

void setRates(int current, int type) // This function sets the rate of each event for a given particle that moves
{
    particle pt = particleList[current];

    for (int j = 0; j < 4; ++j)
    {
        pt.motion[j] = diffusion[j];
        auto iT = eventList.find(event(current, j, pt.tau[j]));
        if (iT != eventList.end())
        {
            eventList.erase(iT);
        }
        pt.tau[j] = dblmax;
        T r = randn();
        T dtau = totalTimeNRM - log(r) / pt.motion[j];
        eventList.insert(event(current, j, dtau));
        pt.tau[j] = dtau;
    }
    particleList[current] = pt;
    CheckWithinDomain(current);
    pt = particleList[current];
}

void setRates() // This function sets the rate of each event for each particle
{
    int current;

    particle pt;

    for (current = 0; current < particleList.size(); ++current)
    {

        pt = particleList[current];

        for (int j = 0; j < 4; ++j)
        {
            pt.motion[j] = diffusion[j];

            T dtau = totalTimeNRM - log(randn()) / pt.motion[j];

            auto iT = eventList.find(event(current, j, pt.tau[j]));
            if (iT != eventList.end())
            {
                eventList.erase(iT);
            }
            pt.tau[j] = dtau;
            eventList.insert(event(current, j, dtau));
        }

        particleList[current] = pt;
        CheckWithinDomain(current);
        pt = particleList[current];
    }
}

void iniKMCLattice() // Used to initalize the LKMC lattice with given initial particle distribution
{
    neigh.resize(nparticles, vector<int>(nparticles));
    int quarter = floor(nparticles / 4);
    int xbump = 0;
    int ybump = 0;
    int i = 0;
    while (i < nparticles)
    {
       particle pt;
       pt.center[0] = int(randn() * subIS);
       pt.center[1] = int(randn() * subIS);
        pt.center[2] = 1;
       /* if (i < quarter * 4)
        {
            ybump = floor((i + 1) / (quarter * 2));
            xbump = floor((i + 1) / quarter) - 2 * floor((i + 1) / (quarter * 2));
            pt.center[0] = int(randn() * ImSize / 2) + ImSize * xbump / 2;
            pt.center[1] = int(randn() * ImSize / 2) + ImSize * ybump / 2;
            pt.center[2] = 1;
        }
        else
        {
            pt.center[0] = int(randn() * ImSize);
            pt.center[1] = int(randn() * ImSize);
            pt.center[2] = 1;

        }*/
        if (FiniteVol == 1)
        {
            for (int j = 0; j < (particleList.size()); ++j)
            {
                neigh[i][j] = 0;
                neigh[j][i] = 0;

                int rx = (pt.center[0] - particleList[j].center[0]);
                int ry = (pt.center[1] - particleList[j].center[1]);
                int rz = (pt.center[2] - particleList[j].center[2]);
                int rr = rx * rx + ry * ry + rz * rz; // distance between two particles 
                if (rr <= neighlimit)
                {
                    neigh[j][i] = 1; // YES a neighbor
                    neigh[i][j] = 1;
                }
            }
            if (!Wall(pt.center[0], pt.center[1], pt.center[2]) || !Overlap(i, pt.center[0], pt.center[1], pt.center[2]))
            {
                pt.dp = dia;
                particleList.push_back(pt);
                ++i;
            }
        }
        else
        {
            if (!Wall(pt.center[0], pt.center[1], pt.center[2]))
            {
                pt.dp = dia;
                particleList.push_back(pt);
                ++i;
            }
        }
    }

}

void runKMC() // one iteration of KMC
{
    int current;
    int j;
    T oldTime = totalTimeNRM;
    auto executedEvent = eventList.begin(); // choose fastest event to execute
    current = executedEvent->idx;
    j = executedEvent->type;
    totalTimeNRM = executedEvent->tau;
    eventList.erase(executedEvent);
    particle pt = particleList[current];
    pt.center[j / 2] = pt.center[j / 2] + (int)pow(-1, j); //move particle by one node in the direction specified by j
    if (FiniteVol == 1)
    {
        particleList[current] = pt; // update particle list
        setRates(current, j);
        int part = current;
        for (int current = 0; current < particleList.size(); ++current)
        {

            if (neigh[part][current] == 1 && part != current)
            {
                setRates(current, j);
            }
            else
            {
                continue;
            }
        }
    }
    else
    {
        particleList[current] = pt; // update particle list
        setRates(current, j); // update new rates for system
    }
}

void writeReadMe() // Write out particle positions for visualization / analysis
{
    // first write the ReadMe
    std::ofstream fileReadMe("ReadMe_C" + std::to_string(ident) + ".txt");

    // write ReadMe file
    fileReadMe << "ID:" << ident << "\n";
    fileReadMe << "diffusion runs\n";
    fileReadMe << "number particles:" << nparticles << "\n";
    fileReadMe << "pixels are (8nm)^3/pixel.\n";
    fileReadMe << "particle diameter:" << dia << " nodes \n";
    fileReadMe << "particle diffusivity:" << diffusion[0] << " nodes^2/s \n";
    fileReadMe << "RNG seed:" << seed << " s^-1 \n";
    fileReadMe.close();
}

void writeOutputFiles() // Write out particle positions for visualization / analysis
{
    std::ofstream file("ForOvito_C" + std::to_string(ident) + ".txt");
    std::ofstream ffile("Pos_C" + std::to_string(ident) + ".csv");
    std::ofstream fffile("Time_C" + std::to_string(ident) + ".csv");
    for (int j = 0; j < ITER; ++j)
    {
        file << nparticles << " \n";
        file << j << "th iteration, each pixel =" << 10000 / ImSize << " squared.\n";
        for (int i = 0; i < particleListOverall[j].size(); ++i)
        {
            particle pt = particleListOverall[j][i];
            file << "H  ";
            file << pt.center[0] << "   ";
            file << pt.center[1] << "   ";
            file << pt.center[2] << endl;
            ffile << pt.center[0] << ",";
            ffile << pt.center[1] << ",";
            ffile << pt.center[2] << endl;
        }
        fffile << TIME[j] << ",";
    }
    file.close();
    ffile.close();
    fffile.close();
}