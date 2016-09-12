#include "Square.h"

// ---------------------------------------------------------------------------
int CSquare::sqSize = 20.0;

// ---------------------------------------------------------------------------
CSquare::CSquare()
{
	X = 0;
	Y = 0;

	nTrees_t0 = 0;
	nSpec_t0 = 0;
   nRecruits_t1 = 0;
   nDeath_t1 = 0;
};

CSquare::CSquare(int ix, int iy)
{
	iX = ix;
   iY = iy;

	X = (double) iX*sqSize;
	Y = (double) iY*sqSize;

	nTrees_t0 = 0;
   nRecruits_t1 = 0;
   nDeath_t1 = 0;
};

// ---------------------------------------------------------------------------
CSquare::~CSquare()
{
   TreeList.clear();
   RecruitList.clear();
	SpecMap.clear();
};

// ---------------------------------------------------------------------------
void CSquare::InitSquare(int ix, int iy)
{
	iX = ix;
   iY = iy;

	X = (double) iX*sqSize;
	Y = (double) iY*sqSize;

   nTrees_t0 = 0;
   nSpec_t0 = 0;
   nRecruits_t1 = 0;
   nDeath_t1 = 0;

   TreeList.clear();
   RecruitList.clear();
   SpecMap.clear();
};
