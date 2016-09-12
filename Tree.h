//---------------------------------------------------------------------------

#ifndef TreeH
#define TreeH

#include <math.h>
#include <vector>
#include <map>

//---------------------------------------------------------------------------
class CTree
{
public:
	long int TreeID;

	double X;
	double Y;

	unsigned int SpecID;
	//bool dead;

	CTree(long int id, double x, double y, int spec)
	{
		TreeID = id;
		X = x;
		Y = y;
		SpecID = spec;
		//dead = false;
	};

	~CTree(){};
};

typedef std::vector<CTree*>::iterator TreeIterV;
//typedef std::list<CTree*>::iterator TreeIterL;

//---------------------------------------------------------------------------
#endif
