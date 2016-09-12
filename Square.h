//---------------------------------------------------------------------------

#ifndef SquareH
#define SquareH

#include <list>
#include "Tree.h"

//---------------------------------------------------------------------------
const double Pi = 3.14159265358979323846;

class CTree;

class CSquare
{
public:
	static int sqSize;  //side length of the squares

	//index coordinates
	int iX;
	int iY;

	//coordinates in meters
	double X;
	double Y;

	int nTrees;

	int nTrees_t0;
	int nSpec_t0;
	int nRecruits_t1;
	int nDeath_t1;

   double pDeath; //specific mortality rate for each square
	double pRec; //specific mortality rate for each square

	std::vector<CTree*> TreeList; //list of trees in the square
	std::vector<CTree*> RecruitList;
	std::map<int,int> SpecMap;

	CSquare();
	CSquare(int ix, int iy);
	~CSquare();

	void InitSquare(int ix, int iy);
};

typedef std::vector<CSquare*>::iterator SquareIterV;

//---------------------------------------------------------------------------
#endif
