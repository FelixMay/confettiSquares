//---------------------------------------------------------------------------

#include "Forest.h"
#include <time.h>
#include <string>
using namespace std;
//---------------------------------------------------------------------------

int StringToInt(std::string S)
{
	int result;
	std::istringstream is(S);
	is>>result;
	return result;
}

int main(int argc, char* argv[])
{
   int SquareSize = 20;
   int NRep = 30;  //number replicates
   bool StepsOut = true;

	string SimFileName = "Para";
	string FileLabel = "1";

	if (argc == 2){
		SimFileName = argv[1];
	}

	if (argc == 3){
		SimFileName = argv[1];
		FileLabel = argv[2];
	}

	if (argc == 4){
		SimFileName = argv[1];
		FileLabel = argv[2];
		NRep = StringToInt(argv[3]);
	}

	if (argc == 5){
		SimFileName = argv[1];
		FileLabel = argv[2];
		NRep = StringToInt(argv[3]);
		SquareSize = StringToInt(argv[4]);
	}

	SimFileName = "InOut/"+SimFileName + FileLabel +".txt";

	ifstream InFile;
	InFile.open(SimFileName.c_str());

	cout<<"Input-File:\t"<<SimFileName<<endl;
	cout<<"Sim-Label:\t"<<FileLabel<<endl;
	cout<<"Size of squares:\t"<<SquareSize<<endl;
	time_t start, end;

	start = time(0);

	string line1;
	getline(InFile,line1);

	CPara* pPara = new CPara();

	int seed = (int) start;
	//int seed = 99;

	CForest* pForest = new CForest(seed, SquareSize, StepsOut);
	pForest->FileOpen(FileLabel);

	int isim;

	if (InFile.good()) {

		//read first parameter set
		InFile>>isim;
		InFile>>pPara->Tmax;
		InFile>>pPara->pImmi;
		InFile>>pPara->pMortRec;
		InFile>>pPara->pLDD;
		InFile>>pPara->meanDisp;
		InFile>>pPara->sdDisp;
		InFile>>pPara->Scenario;
		pPara->GetDispPars();

		while (InFile.good()) {

			cout<<"Sim "<<isim<<endl;
			pForest->Pars = pPara;

			//run simulations
			for (int irep=1; irep <= NRep; ++irep) {

				cout<<"  Rep "<<irep<<endl;
				pForest->OneRun(isim, irep, FileLabel);
				pForest->ClearForest();
			}  // end irep

			//try to read new parameter set
			InFile>>isim;
         InFile>>pPara->Tmax;
         InFile>>pPara->pImmi;
         InFile>>pPara->pMortRec;
         InFile>>pPara->pLDD;
         InFile>>pPara->meanDisp;
         InFile>>pPara->sdDisp;
         InFile>>pPara->Scenario;
         pPara->GetDispPars();
      }
	}
	else cout<<"Error SimFile"<<endl;

	delete pForest;
	delete pPara;

	end = time(0);

	cout<<"\nRuntime: "<<end - start<<" seconds"<<endl;

	cin.ignore();

	return 0;
}
//---------------------------------------------------------------------------
