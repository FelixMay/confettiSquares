//---------------------------------------------------------------------------

#ifndef ForestH
#define ForestH

#include "Tree.h"
#include "Square.h"
#include "randomc.h"
//#include "randoma\randoma.h"
#include "stocc.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <sstream>
#include <cmath>

//---------------------------------------------------------------------------
class CPara
{
public:
	int    Tmax;    //number of time steps (one time step = 5 years)
	double pImmi;     //probability of immigration of a novel species during an recruitment event
	double pMortRec;  //mortality rate for neutral model
	double pLDD;      //probability of global dispersal in the plot
	double meanDisp;     //mean dispersal distance in meters
	double sdDisp;    //standard deviation if dispersal distance
	double alpha;     //parameter of dispersal kernel (Clark et al. 1999 Eq. 5b
   double muDisp;    // parameter of log-normal dispersal kernel
   double sigmaDisp; // parameter of log-normal dispersal kernel

   int Scenario;   //1 ... neutral mode
                   //2 ... mortality and recruitment rates from data
                   //3 ... mortality and recruitment rates + correlation from data
                   //4 ... mortality and recruitment rates + local density dependence from data
	CPara()
	{
	   Tmax = 10;
		pImmi = 4.0;
		pMortRec = 0.14;
		pLDD  = 0.1;

		meanDisp = 30;
		sdDisp = 30;

		alpha = 2.0*meanDisp/sqrt(Pi);
		sigmaDisp = sqrt(log(1.0 + (sdDisp*sdDisp)/(meanDisp*meanDisp)));
		muDisp = log(meanDisp) - 0.5 * sigmaDisp*sigmaDisp;

		Scenario = 1;
   }

	CPara(int    tmax,
         double pimmi,
         double pmortrec,
         double pldd,
         double mdisp,
         double sddisp,
         int scena
        )
   {
      Tmax = tmax;
      pImmi = pimmi;
      pMortRec = pmortrec;
      pLDD = pldd;
      meanDisp = mdisp;
      sdDisp = sddisp;

      alpha = 2.0*meanDisp/sqrt(Pi);
      sigmaDisp = sqrt(log(1.0 + (sdDisp*sdDisp)/(meanDisp*meanDisp)));
	   muDisp = log(meanDisp) - 0.5 * sigmaDisp*sigmaDisp;

      Scenario = scena;
   }

	~CPara(){};

	void GetDispPars(){
	   alpha = 2.0*meanDisp/sqrt(Pi);

      sigmaDisp = sqrt(log(1.0 + (sdDisp*sdDisp)/(meanDisp*meanDisp)));
		muDisp = log(meanDisp) - 0.5 * sigmaDisp*sigmaDisp;
	}
};

//---------------------------------------------------------------------------
class CForest
{
public:
	CPara* Pars;
	bool spatial_output;
	bool time_output;

	double Xmax;
	double Ymax;

	int64_t MaxTreeID;
	int nTreesPlot;

	int SpecMax;
	std::map<int,int> SpecAbund;  // map with first --> key, second --> abund

	//Raster
	int SquareSize;  //side length of one square

	int XSquares;  //number of squares in x-direction
	int YSquares;  //number of squares in y-direction

	int nSquares;

	std::vector<CSquare*> SquareVec;  //vector of pointers
	CSquare ***SquareGrid;       //matrix of pointers

	//random number generators
	CRandomMersenne* RandGen1;
	StochasticLib1* RandGen2;

	//index variables
	int isim;
	int irep;

	//Input Data
	std::vector<int> nRecruitObs;
	std::vector<int> nDeathObs;

	std::vector<double> pRecruitObs; //recruitment rates
	std::vector<double> pDeathObs; //mortality rates
	//std::vector<int> nTreesObs;

	int DensBin;
	int DensCl_Min;                 //smallest and largest density class
	int DensCl_Max;

	std::map<int, std::vector<int> > nRecruitObsDens; //numbers of recruits in density classes
	std::map<int, std::vector<int> > nDeathObsDens; //numbers of mortality events in density classes
	std::map<int, std::vector<double> > pDeathObsDens; //mortality rates in density classes
	std::map<int, std::vector<double> > pRecruitObsDens; //recruitment rates in density classes

	//Files
	std::fstream AbundFile;
	std::fstream DivFile;
	std::fstream SAD_File;
	std::fstream SquareFile;

	static const int MaxSAD = 17;
	int SAD[MaxSAD];    //Species abundance distributions as octave curve 2^0 - 2^16

	//Private functions
	inline double Distance(double x1, double y1, double x2, double y2);
	inline void BoundIntGrid(int& xx, int& yy, int Xmax, int Ymax);

	//inline void BoundGrid(int& xx, int& yy, double& xb, double& yb);
	inline void PeriodBound(double& xx, double& yy);
	void GetNewXY(double &x1, double &y1);
	int GetSpecID_mother(double x0, double y0, bool nearest = false);

	CForest(int seed, int square_size, bool steps_out);
	~CForest();

	void FileOpen(std::string label);
	void Init();

	void ClearForest();
	void WriteTrees(int step,int sim, int rep, std::string label);
	void WriteSquares(int step,int sim, int rep, std::string label);
	void WriteOutput(int step, int sim, int rep);
	double GetShannon();

	int GetSAD();

	int GetRandSquare(); //sample squares with probability of their relative abundance

	void Recruit_Mort_OneStep();
	void UpdateSquares();

	//Neutral dynamics
	void mort_neutral_random_tree(int nmort);

	void recruit_neutral_random_mother_fw(int nrecruit);
	void recruit_neutral_bw(int nrecruit, bool dens_dep);

	//sampling from data
	void recruit_data_bw(int nrecruit, bool replace_sq, bool dens_dep);
	void mort_data(int nmort, bool replace_sq, bool dens_dep);
	void mort_data_cluster(int nmort, bool replace_sq, bool dens_dep);
	void mort_data_cluster2(int nmort, bool replace_sq, bool dens_dep);

   void Mort_Dens_Data3 (int nmort); //density dependence from data
   void Mort_Dens_Data4a (int nmort);
   void Mort_Dens_Data4b (int nmort);

   void Rec_Dens_Data3 (int nmort);
   void Rec_Dens_Data4 (int nmort);

	void OneRun(int isim, int irep, std::string label);

	std::string IntToString(int i)
	{
	  std::ostringstream os;
	  os<<i;
	  return os.str();
	}
};
//---------------------------------------------------------------------------
#endif
