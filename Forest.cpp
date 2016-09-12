// ---------------------------------------------------------------------------
#include "Forest.h"
#include <time.h>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iomanip>
#include <cstdlib>
using namespace std;

// ---------------------------------------------------------------------------
CForest::CForest(int seed, int square_size, bool steps_out) {

   time_output = steps_out;

   Xmax = 1000.0;
   Ymax = 500.0;

   // Grid
   SquareSize = square_size;
   CSquare::sqSize = square_size;
   //SquareSize = CSquare::sqSize; // square cell size

   XSquares = (int) Xmax / SquareSize;
   YSquares = (int) Ymax / SquareSize;

   nSquares = XSquares*YSquares;

   CSquare* pSquare;

   SquareGrid = new (CSquare**[XSquares]);
   for (int iX = 0; iX < XSquares; iX++){
      SquareGrid[iX] = new (CSquare*[YSquares]);
      for (int iY = 0; iY < YSquares; iY++){
         pSquare = new CSquare(iX,iY);
         SquareVec.push_back(pSquare);
         SquareGrid[iX][iY] = pSquare;
      }
   }

   // Random number generators
   //int seed = (int) time(0);            // random seed
   // RandGen1 = new CRandomMersenneA(seed);

   srand (seed);
   RandGen1 = new CRandomMersenne(seed);
   RandGen2 = new StochasticLib1(seed);

   // Input-Output
   //isim = 1;
   //irep = 1;

   ifstream InputFile;

   //read input data
   if (SquareSize == 10.0){
      InputFile.open("InOut/Data/SquaresInput10.txt");
      DensBin = 10;
   }
   if (SquareSize == 20.0){
      InputFile.open("InOut/Data/SquaresInput20.txt");
      DensBin = 40;
   }
   if (SquareSize == 50.0){
      InputFile.open("InOut/Data/SquaresInput50.txt");
      DensBin = 100;
   }
   if (SquareSize == 100.0){
      InputFile.open("InOut/Data/SquaresInput100.txt");
      DensBin = 200;
   }

   //string line, sq_label, year;
   //int x,y, ntrees, nspec, nrec, ndead, dens_class;

   string line, sq_label;
   int ntrees, nrec, ndead, dens_class;

   DensCl_Min = 100;
   DensCl_Max = 0;

   if (InputFile.good()){

     getline(InputFile,line);

      do {
         //InputFile>>sq_label;
         //InputFile>>x;
         //InputFile>>y;
         InputFile>>ntrees;
         //InputFile>>nspec;
         InputFile>>ndead;
         InputFile>>nrec;
         //InputFile>>year;

         if (!InputFile.eof()){
            //nTreesObs.push_back(ntrees);
            nRecruitObs.push_back(nrec);
            nDeathObs.push_back(ndead);

            dens_class = (int) floor(ntrees/DensBin);
            if (dens_class < DensCl_Min) DensCl_Min = dens_class;
            if (dens_class > DensCl_Max) DensCl_Max = dens_class;

            nRecruitObsDens[dens_class].push_back(nrec);
            nDeathObsDens[dens_class].push_back(ndead);

            double prec = static_cast<double>(nrec)/ntrees;
            double pmort = static_cast<double>(ndead)/ntrees;

            pRecruitObs.push_back(prec);
            pDeathObs.push_back(pmort);

            pRecruitObsDens[dens_class].push_back(prec);
            pDeathObsDens[dens_class].push_back(pmort);
         }
      } while(!InputFile.eof());

      InputFile.close();
   }
   else cout<<"Error InFile!"<<endl;
}

// ---------------------------------------------------------------------------
CForest::~CForest() {

  SpecAbund.clear();

   for (SquareIterV isquare = SquareVec.begin();
                    isquare != SquareVec.end(); isquare++){
      delete (*isquare);
   }

   for (int ix = 0; ix < XSquares; ix++)
      delete[] SquareGrid[ix];
   delete[] SquareGrid;

   delete RandGen1;
   delete RandGen2;

   //nTreesObs.clear();
   nRecruitObs.clear();
   nDeathObs.clear();
   SquareFile.close();
}

// ---------------------------------------------------------------------------
void CForest::FileOpen(string label) {
   string FileName;
   string FileNameEnd = label + ".csv";

   FileName = "InOut/Diversity" + FileNameEnd;
   DivFile.open(FileName.c_str(), ios::in); {
      if (DivFile.good()) {
         DivFile.close();
         DivFile.open(FileName.c_str(), ios::out | ios::app);
      }
      else {
         DivFile.clear();
         DivFile.open(FileName.c_str(), ios::out);
         DivFile << "SimNr; RepNr; Step; NTrees; NSpec; Shannon" << endl;
      }
   }

   /*
   FileName = "InOut\\Abund" + FileNameEnd;
   AbundFile.open(FileName.c_str(), ios::in); {
      if (AbundFile.good()) {
         AbundFile.close();
         AbundFile.open(FileName.c_str(), ios::out | ios::app);
      }
      else {
         AbundFile.clear();
         AbundFile.open(FileName.c_str(), ios::out);
      }
   }
   */

   /*
   FileName = "InOut\\SAD" + FileNameEnd;
   SAD_File.open(FileName.c_str(), ios::in); {
      if (SAD_File.good()) {
         SAD_File.close();
         SAD_File.open(FileName.c_str(), ios::out | ios::app);
      }
      else {
         SAD_File.clear();
         SAD_File.open(FileName.c_str(), ios::out);

         SAD_File << "A1; A2_3; A4_7; A8_15; A16_31; A32_63; A64_127; A128_255; "
                  << "A256_511; A512_1023; A1024_2047; A2048-4095; A4096-8191; "
                  << "A8192_16383; A16384_32767; A32768_65535; A65536-Inf"
                  << endl;
         // SAD_File<<"A1; A2; A3_4; A5_8; A9_16; A17_32; A33_64; A65_128;"
         // <<"A129_256; A257_512; A513_1024; A1025_2048; A2049-Inf;"<<endl;
      }
   }
   */
}

// ---------------------------------------------------------------------------
void CForest::ClearForest() {

   for (SquareIterV isquare = SquareVec.begin();
                    isquare != SquareVec.end(); ++isquare){
      for (TreeIterV itree = (*isquare)->TreeList.begin();
                     itree != (*isquare)->TreeList.end(); ++itree)
         delete (*itree);

      (*isquare)->TreeList.clear();
      (*isquare)->RecruitList.clear();
   }

   SpecAbund.clear();
}

// ---------------------------------------------------------------------------
inline double CForest::Distance(double x1, double y1, double x2, double y2) {
   return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

// ---------------------------------------------------------------------------
inline void CForest::BoundIntGrid(int& xx, int& yy, int Xmax, int Ymax) {
   xx = xx % Xmax;
   if (xx < 0)
      xx = Xmax + xx;

   yy = yy % Ymax;
   if (yy < 0)
      yy = Ymax + yy;
}

// ---------------------------------------------------------------------------
inline void CForest::PeriodBound(double& xx, double& yy) {
   xx = Xmax * (xx / Xmax - floor(xx / Xmax));
   yy = Ymax * (yy / Ymax - floor(yy / Ymax));
}

// ---------------------------------------------------------------------------
void CForest::Init() {

   CTree *pTree;

   int iX, iY;

   // Init Raster Squares
   for (iX = 0; iX < XSquares; iX++)
      for (iY = 0; iY < YSquares; iY++)
         SquareGrid[iX][iY]->InitSquare(iX, iY);

   // Init Trees
   double x, y;
   unsigned int SpecID;

   SpecMax = 0;
   MaxTreeID = 0;

   //NTrees = Pars->NTrees;

   // random init
   /*
   for (int i = 0; i < NTrees; i++) {
      x = RandGen1->Random() * Xmax;
      y = RandGen1->Random() * Ymax;
      SpecID = GetRandSpec();

      pTree1 = new CTree(x, y, SpecID, Pars->r_max);
      TreeList.push_back(pTree1);

      iX1 = (int) floor(x / CellSize);
      iY1 = (int) floor(y / CellSize);
      Grid[iX1][iY1].TreeList.push_back(pTree1);

      ++SpecAbund[pTree1->SpecID];
   }
   */

   // initialization from input-file

   ifstream BCI_File("InOut/Data/bci1995_1cm_20140212.txt");

   string line, tag, SpecStr;
   double dbh;

   nTreesPlot = 0;

   if (BCI_File.good()){

     getline(BCI_File,line);

      do {
         BCI_File>>tag;
         BCI_File>>x;
         BCI_File>>y;
         BCI_File>>dbh;
         BCI_File>>SpecStr;
         BCI_File>>SpecID;

         if (!BCI_File.eof()){
            if ((x<Xmax) && (y<Ymax)){

               MaxTreeID++;
               pTree = new CTree(MaxTreeID,x,y,SpecID);
               ++nTreesPlot;
               ++SpecAbund[pTree->SpecID];

               iX = (int) floor(pTree->X/SquareSize);
               iY = (int) floor(pTree->Y/SquareSize);

               SquareGrid[iX][iY]->TreeList.push_back(pTree);
               ++SquareGrid[iX][iY]->SpecMap[SpecID];
               //++SquareGrid[iX][iY]->nTrees_t0;
            }

            if (SpecID>SpecMax) SpecMax = SpecID;
         }
      } while(!BCI_File.eof());

      BCI_File.close();
   }
}

// ---------------------------------------------------------------------------
void CForest::GetNewXY(double &x1, double &y1) {

   double x0 = x1;
   double y0 = y1;

   double r, r_dist, r_angle;

   if (Pars->sdDisp<0.0001){
      r = RandGen1->Random();

      //r_dist = sqrt(-Pars->alpha*Pars->alpha*(log(1.0-r)));     //Gaussian dispersal kernel following Clark et al. 1999
                                                                //this kernel produces the correct mean, but has a unimodal shape?
      r_dist = sqrt(Pi/2.0) * std::abs(RandGen2->Normal(0.0,Pars->meanDisp)); //Gaussian kernel with the same mean,
                                                                         //but different normalization constant
                                                                         //maximum density at distance zero
      //r_dist = r*2.0*Pars->mDisp;               //uniform dispersal
   }
   else r_dist = exp(RandGen2->Normal(Pars->muDisp, Pars->sigmaDisp));  //log-normal dispersal kernel

   r_angle = RandGen1->Random() * 2.0 * Pi;

   x1 = x0 + cos(r_angle) * r_dist;
   y1 = y0 + sin(r_angle) * r_dist;

   PeriodBound(x1, y1);
}

// ---------------------------------------------------------------------------
int CForest::GetSpecID_mother(double x0, double y0, bool nearest) {

   int ntrees{0}, iX_mother, iY_mother, SpecID_mother{0};
   double r, r_dist, r_angle, x1,y1;

   int ntrial{0};
   int max_trial = 20;

   do {

      if (Pars->sdDisp<0.0001){
         r_dist = sqrt(Pi/2.0) * std::abs(RandGen2->Normal(0.0,Pars->meanDisp));
      }
      else r_dist = exp(RandGen2->Normal(Pars->muDisp, Pars->sigmaDisp));  //log-normal dispersal kernel

      r_angle = RandGen1->Random() * 2.0 * Pi;

      x1 = x0 + cos(r_angle) * r_dist;
      y1 = y0 + sin(r_angle) * r_dist;

      PeriodBound(x1, y1);

      iX_mother = floor(x1/SquareSize);
      iY_mother = floor(y1/SquareSize);

      //random tree in cell
      ntrees = SquareGrid[iX_mother][iY_mother]->TreeList.size();
      ++ntrial;
   } while ((ntrees == 0) && (ntrial < max_trial));

   if (ntrial < max_trial){

      if (nearest){//search nearest tree

         double dmin = 2.0*SquareSize;
         double d1;

         CTree *itree_next = SquareGrid[iX_mother][iY_mother]->TreeList[0];

//         for (TreeIterV itree = SquareGrid[iX_mother][iY_mother]->TreeList.begin();
//                        itree != SquareGrid[iX_mother][iY_mother]->TreeList.end();itree++)
         for (const auto itree: SquareGrid[iX_mother][iY_mother]->TreeList)
         {
            d1 = Distance(x1,y1,itree->X,itree->Y);
            if (d1 < dmin){
               dmin = d1;
               itree_next = itree;
            }
         }
         SpecID_mother = itree_next->SpecID;
      }
      else {
         //random tree from cell
         int irand = RandGen1->IRandom(0,ntrees-1);
         SpecID_mother = SquareGrid[iX_mother][iY_mother]->TreeList[irand]->SpecID;
      }
   }

   else { //long distance dispersal

      CSquare *pSquare_mother = nullptr;

      do {
         pSquare_mother = SquareVec[RandGen1->IRandom(0,nSquares-1)];
         ntrees = pSquare_mother->TreeList.size();
      } while (ntrees == 0);

      //random tree in Square
      SpecID_mother = pSquare_mother->TreeList[RandGen1->IRandom(0,ntrees-1)]->SpecID;
   }

   return(SpecID_mother);
}

//---------------------------------------------------------------------------
void CForest::mort_neutral_random_tree(int nmort)
{
   int countMort{0};
   int isquare{0};
   int itree{0};

   CSquare *pSquare = nullptr;
   CTree *pTree = nullptr;

   do {

      //choose random square with probabilities of relative abundances
      isquare = GetRandSquare();
      pSquare = SquareVec[isquare];

      //random tree in square
      itree = RandGen1->IRandom(0,pSquare->TreeList.size()-1);
      pTree = pSquare->TreeList[itree];

      //update square
      --pSquare->SpecMap[pTree->SpecID];
      if (pSquare->SpecMap[pTree->SpecID] == 0)
         pSquare->SpecMap.erase(pTree->SpecID);
      ++pSquare->nDeath_t1;
      --pSquare->nTrees;

      //update forest
      --SpecAbund[pTree->SpecID];
      if (SpecAbund[pTree->SpecID] == 0)
         SpecAbund.erase(pTree->SpecID);

      --nTreesPlot;
      pSquare->TreeList.erase(pSquare->TreeList.begin()+itree);
      delete pTree;

      ++countMort;
      if (countMort == nmort) break; //finishes FOR loop

   } while (countMort < nmort);
}

//---------------------------------------------------------------------------
void CForest::recruit_neutral_random_mother_fw(int nrecruit)
{
   int countRecPlot{0};
   int isquare{0};
   int itree{0};
   int SpecID_mother{0};

   double x_new{0.0}, y_new{0.0};
   double rand1;

   CSquare *pSquare = nullptr;
   CTree *pTree1 = nullptr;
   CTree *pTree2 = nullptr;

   int iX{0};
   int iY{0};

   do {
      //choose random square with probabilities of relative abundances
      isquare = GetRandSquare();
      pSquare = SquareVec[isquare];

      //random tree in square
      itree = RandGen1->IRandom(0,pSquare->TreeList.size()-1);
      pTree1 = pSquare->TreeList[itree];

		rand1 = RandGen1->Random();

		//immigration of novel species
		if (rand1 < Pars->pImmi){
			SpecMax++;
			SpecID_mother = SpecMax;

			//random position in plot
			x_new = RandGen1->Random()*Xmax;
			y_new = RandGen1->Random()*Ymax;
		}

      else {
      //long-distance dispersal in the plot
         if ((rand1 >= Pars->pImmi) & (rand1 < (Pars->pImmi+Pars->pLDD))){

            //random position in plot
            x_new = RandGen1->Random()*Xmax;
            y_new = RandGen1->Random()*Ymax;

            //random tree in Square
            SpecID_mother = pTree1->SpecID;
         }

         else {
            //local dispersal
            x_new = pTree1->X;
            y_new = pTree1->Y;

            GetNewXY(x_new,y_new);

            SpecID_mother = pTree1->SpecID;
         }
      } // if not immigration

      MaxTreeID++;
      pTree2 = new CTree(MaxTreeID,x_new,y_new,SpecID_mother);

      iX = floor(pTree2->X/SquareSize);
      iY = floor(pTree2->Y/SquareSize);

      pSquare = SquareGrid[iX][iY];
      pSquare->RecruitList.push_back(pTree2);

      ++pSquare->nRecruits_t1;
      ++countRecPlot;

   } while (countRecPlot < nrecruit);
}

//---------------------------------------------------------------------------
void CForest::recruit_neutral_bw(int nrecruit, bool dens_dep)
{
   CTree* pTree = nullptr;
   CSquare* pSquare = nullptr;
   CSquare* pSquare_mother = nullptr;

   //int iX, iY
   int SpecID_mother, irand;
   double x_new, y_new, rand1;

   int countRecPlot{0};
   int isquare{0};
   int sq_ntrees{0};

   do {
      //choose random square with or without density dependence
      if (dens_dep == true) isquare = GetRandSquare();
      else                  isquare = RandGen1->IRandom(0,nSquares-1);
      pSquare = SquareVec[isquare];

      //random position in cell
      x_new = (pSquare->iX + RandGen1->Random() )*SquareSize;
      y_new = (pSquare->iY + RandGen1->Random() )*SquareSize;

      //find mother tree
      rand1 = RandGen1->Random();

      //immigration of novel species
      if (rand1 < Pars->pImmi){
         SpecMax++;
         SpecID_mother = SpecMax;
      }

      else {
      //long-distance dispersal in the plot
         if ((rand1 >= Pars->pImmi) & (rand1 < (Pars->pImmi+Pars->pLDD))){
            //random square with trees
            do {
               irand = RandGen1->IRandom(0,nSquares-1);
               pSquare_mother = SquareVec[irand];
               sq_ntrees = pSquare_mother->TreeList.size();
            } while (sq_ntrees == 0);

            //random tree in Square
            SpecID_mother = pSquare_mother->TreeList[RandGen1->IRandom(0,sq_ntrees-1)]->SpecID;
         }

         else {
            //local dispersal from one of the eight neighbors
            SpecID_mother = GetSpecID_mother(x_new,y_new);
         }
      } // if not immigration

      ++MaxTreeID;
      pTree = new CTree(MaxTreeID,x_new,y_new,SpecID_mother);
      pSquare->RecruitList.push_back(pTree);
      ++pSquare->nRecruits_t1;

      ++countRecPlot;

   } while (countRecPlot < nrecruit);
}

//---------------------------------------------------------------------------
void CForest::recruit_data_bw(int nrecruit, bool replace_sq, bool dens_dep)
{
   CTree* pTree = nullptr;
   CSquare* pSquare = nullptr;
   CSquare* pSquare_mother = nullptr;

   //int iX, iY
   int SpecID_mother, irand;
   double x_new, y_new, rand1;

   int countRecPlot{0}, countRecSq{0};
   int isquare{0};
   int sq_ntrees{0};
   int nRecSq{0.0};

   int dens_cl{0};
   int nsq_dens{0};

   int nData = nRecruitObs.size();

   if (!replace_sq) random_shuffle(SquareVec.begin(),SquareVec.end());  //random permutation of square order

   do {
      //choose random square without density dependence
      if (replace_sq) isquare = RandGen1->IRandom(0,nSquares-1);
      pSquare = SquareVec[isquare];

      if (dens_dep){
         dens_cl = floor(pSquare->TreeList.size()/DensBin);

         if (dens_cl < DensCl_Min) dens_cl = DensCl_Min;
         if (dens_cl > DensCl_Max) dens_cl = DensCl_Max;

         nsq_dens = nRecruitObsDens[dens_cl].size();
         nRecSq = nRecruitObsDens[dens_cl][RandGen1->IRandom(0,nsq_dens-1)];
      }
      else {
         nRecSq = nRecruitObs[RandGen1->IRandom(0,nData-1)];
      }

      countRecSq = 0;

      do {
         //random position in cell
         x_new = (pSquare->iX + RandGen1->Random() )*SquareSize;
         y_new = (pSquare->iY + RandGen1->Random() )*SquareSize;

         //find mother tree
         rand1 = RandGen1->Random();

         //immigration of novel species
         if (rand1 < Pars->pImmi){
            SpecMax++;
            SpecID_mother = SpecMax;
         }

         else {
         //long-distance dispersal in the plot
            if ((rand1 >= Pars->pImmi) & (rand1 < (Pars->pImmi+Pars->pLDD))){
               //random square with trees
               do {
                  irand = RandGen1->IRandom(0,nSquares-1);
                  pSquare_mother = SquareVec[irand];
                  sq_ntrees = pSquare_mother->TreeList.size();
               } while (sq_ntrees == 0);

               //random tree in Square
               SpecID_mother = pSquare_mother->TreeList[RandGen1->IRandom(0,sq_ntrees-1)]->SpecID;
            }

            else {
               //local dispersal from one of the eight neighbors
               SpecID_mother = GetSpecID_mother(x_new,y_new);
            }
         } // if not immigration

         ++MaxTreeID;
         pTree = new CTree(MaxTreeID,x_new,y_new,SpecID_mother);
         pSquare->RecruitList.push_back(pTree);
         ++pSquare->nRecruits_t1;

         ++countRecSq;
         ++countRecPlot;

      } while ((countRecSq < nRecSq) && (countRecPlot < nrecruit));

      if (!replace_sq){
         ++isquare;
         if (isquare == nSquares){
            random_shuffle(SquareVec.begin(),SquareVec.end());
            isquare = 0;
         }
      }

   } while (countRecPlot < nrecruit);
}

//---------------------------------------------------------------------------
void CForest::mort_data(int nmort, bool replace_sq, bool dens_dep)
{
   CSquare *pSquare = nullptr;
   CTree *pTree = nullptr;

   int countMortPlot{0};
   int isquare{0};
   int sq_ntrees{0};
   int itree{0};
   int dens_cl{0};
   int nsq_dens{0};

   double pMortSq{0};

   int nData = pDeathObs.size();

    // run over all squares and assign potential mortality
   for (auto &square: SquareVec){

      if (dens_dep){
         dens_cl = (int) floor(square->TreeList.size()/DensBin);

         if (dens_cl < DensCl_Min) dens_cl = DensCl_Min;
         if (dens_cl > DensCl_Max) dens_cl = DensCl_Max;

         nsq_dens = pDeathObsDens[dens_cl].size();
         pMortSq = pDeathObsDens[dens_cl][RandGen1->IRandom(0,nsq_dens-1)];
      }
      else {
         pMortSq = pDeathObs[RandGen1->IRandom(0,nSquares-1)];
      }

      square->pDeath = pMortSq;
   }

   if (!replace_sq) random_shuffle(SquareVec.begin(),SquareVec.end());  //random permutation of square order

    // kill trees in random order
   do {

      //if (replace_sq) isquare = RandGen1->IRandom(0,nSquares-1);
      isquare = GetRandSquare();
      pSquare = SquareVec[isquare];

      sq_ntrees = pSquare->TreeList.size();
      if (sq_ntrees > 0){
         itree = RandGen1->IRandom(0,sq_ntrees-1); //random tree from square
         pTree = pSquare->TreeList[itree];

         //check if tree dies
         if (RandGen1->Random() < pSquare->pDeath){

            //update local abundances
            pSquare->SpecMap[pTree->SpecID]--;
            if (pSquare->SpecMap[pTree->SpecID] == 0)
               pSquare->SpecMap.erase(pTree->SpecID);

            ++pSquare->nDeath_t1;
            --pSquare->nTrees;

            //update forest abundances
            SpecAbund[pTree->SpecID]--;
            if (SpecAbund[pTree->SpecID] == 0)
               SpecAbund.erase(pTree->SpecID);

            --nTreesPlot;
            ++countMortPlot;

            //remove tree
            pSquare->TreeList.erase(pSquare->TreeList.begin()+itree);
            delete pTree;
         }
      } // if sq_nTrees >0

      if (!replace_sq){
         ++isquare;
         if (isquare == nSquares){
            random_shuffle(SquareVec.begin(),SquareVec.end());
            isquare = 0;
         }
      }

   } while (countMortPlot < nmort);
}

//---------------------------------------------------------------------------
void CForest::mort_data_cluster(int nmort, bool replace_sq, bool dens_dep)
{
   CSquare *pSquare = nullptr;
   CTree *pTree = nullptr;

   int countMortPlot{0};
   int isquare{0};
   int sq_ntrees{0};
   int itree{0};
   int dens_cl{0};
   int nsq_dens{0};

   double pMortSq{0.0};

   int nData = nDeathObs.size();

    // run over all squares and assign potential mortality
   for (auto square: SquareVec){

      if (dens_dep){
         dens_cl = (int) floor(square->TreeList.size()/DensBin);

         if (dens_cl < DensCl_Min) dens_cl = DensCl_Min;
         if (dens_cl > DensCl_Max) dens_cl = DensCl_Max;

         nsq_dens = pDeathObsDens[dens_cl].size();
         pMortSq = pDeathObsDens[dens_cl][RandGen1->IRandom(0,nsq_dens-1)];

      }
      else {
         pMortSq = pDeathObs[RandGen1->IRandom(0,nSquares-1)];
      }

      square->pDeath = pMortSq;
   }

   if (!replace_sq) random_shuffle(SquareVec.begin(),SquareVec.end());  //random permutation of square order

    // kill trees in random order
   do {

      if (replace_sq) isquare = RandGen1->IRandom(0,nSquares-1);
      pSquare = SquareVec[isquare];

      sq_ntrees = pSquare->TreeList.size();
      if (sq_ntrees > 0){

         //check all trees from square
         for (TreeIterV itree = pSquare->TreeList.begin();
                        itree != pSquare->TreeList.end();){

            pTree = (*itree);

            //check if tree dies
            if (RandGen1->Random() < pSquare->pDeath){
               //update local abundances
               pSquare->SpecMap[pTree->SpecID]--;
               if (pSquare->SpecMap[pTree->SpecID] == 0)
                  pSquare->SpecMap.erase(pTree->SpecID);

               pSquare->nDeath_t1++;

               //update forest abundances
               SpecAbund[pTree->SpecID]--;
               if (SpecAbund[pTree->SpecID] == 0)
                  SpecAbund.erase(pTree->SpecID);

               --nTreesPlot;
               ++countMortPlot;

               //remove tree
               itree = pSquare->TreeList.erase(itree);
               delete pTree;
               if (countMortPlot == nmort)
                  break;
            } // if death
            else {
               ++itree;
            }
         } //for all trees in
      } // if sq_nTrees >0

      if (!replace_sq){
         ++isquare;
         if (isquare == nSquares){
            random_shuffle(SquareVec.begin(),SquareVec.end());
            isquare = 0;
         }
      }

   } while (countMortPlot < nmort);
}

//---------------------------------------------------------------------------
void CForest::mort_data_cluster2(int nmort, bool replace_sq, bool dens_dep)
{
   CSquare *pSquare = nullptr;
   CTree *pTree = nullptr;

   int countMortPlot{0};
   int countMortSq{0};
   int isquare{0};
   int sq_ntrees{0};
   int itree{0};
   int dens_cl{0};
   int nsq_dens{0};

   int nMortSq{0};

   if (!replace_sq) random_shuffle(SquareVec.begin(), SquareVec.end());  //random permutation of square order

    // kill trees in random order
   do {

      isquare = GetRandSquare();
      //if (replace_sq) isquare = RandGen1->IRandom(0,nSquares-1);
      pSquare = SquareVec[isquare];

      if (dens_dep){
         dens_cl = (int) floor(pSquare->TreeList.size()/DensBin);

         if (dens_cl < DensCl_Min) dens_cl = DensCl_Min;
         if (dens_cl > DensCl_Max) dens_cl = DensCl_Max;

         nsq_dens = nDeathObsDens[dens_cl].size();
         nMortSq = nDeathObsDens[dens_cl][RandGen1->IRandom(0,nsq_dens-1)];

      }
      else {
         nMortSq = nDeathObs[RandGen1->IRandom(0,nSquares-1)];
      }

      countMortSq = 0;

      while (((countMortSq < nMortSq) && (pSquare->TreeList.size()>0)) && (countMortPlot < nmort)){

         //random tree from square
         sq_ntrees = pSquare->TreeList.size();
         itree = RandGen1->IRandom(0,sq_ntrees-1);
         pTree = pSquare->TreeList[itree];

         //update local abundances
         --pSquare->SpecMap[pTree->SpecID];
         if (pSquare->SpecMap[pTree->SpecID] == 0)
            pSquare->SpecMap.erase(pTree->SpecID);

         ++pSquare->nDeath_t1;

         //update forest abundances
         -- SpecAbund[pTree->SpecID];
         if (SpecAbund[pTree->SpecID] == 0)
            SpecAbund.erase(pTree->SpecID);

         --nTreesPlot;
         ++countMortPlot;
         ++countMortSq;

         //remove tree
         pSquare->TreeList.erase(pSquare->TreeList.begin() + itree);
         delete pTree;

      } // while ((countMortSq < nMortSq) && (countMortPlot < nmort)){

      if (!replace_sq){
         ++isquare;
         if (isquare == nSquares){
            random_shuffle(SquareVec.begin(),SquareVec.end());
            isquare = 0;
         }
      }

   } while (countMortPlot < nmort);
}


//---------------------------------------------------------------------------
void CForest::Rec_Dens_Data3 (int nrecruit)
{
   CTree* pTree;
   CSquare* pSquare;
   CSquare* pSquare_mother;

   //int iX, iY;
   int SpecID_mother, irand;
   double x_new, y_new, rand1;

   int countRecPlot = 0;
   int countRecSq = 0;
   int isquare = 0;
   int nRecSq, sq_ntrees;

   int dens_cl, nsq_dens;

   random_shuffle(SquareVec.begin(),SquareVec.end());  //random permutation of square order

   do {

      pSquare = SquareVec[isquare];

      //dens_cl = (int) floor(pSquare->nTrees_t0/DensBin);
      dens_cl = (int) floor(pSquare->TreeList.size()/DensBin);

      if (dens_cl < DensCl_Min) dens_cl = DensCl_Min;
      if (dens_cl > DensCl_Max) dens_cl = DensCl_Max;

      nsq_dens = nRecruitObsDens[dens_cl].size();
      nRecSq = nRecruitObsDens[dens_cl][RandGen1->IRandom(0,nsq_dens-1)];

      countRecSq = 0;

      do {

         //random position in cell
         x_new = (pSquare->iX + RandGen1->Random() )*SquareSize;
         y_new = (pSquare->iY + RandGen1->Random() )*SquareSize;

         //find mother tree
         rand1 = RandGen1->Random();

         //immigration of novel species
         if (rand1 < Pars->pImmi){
            SpecMax++;
            SpecID_mother = SpecMax;
         }

         else {
         //long-distance dispersal in the plot
            if ((rand1 >= Pars->pImmi) & (rand1 < (Pars->pImmi+Pars->pLDD))){
               //random square with trees
               do {
                  irand = RandGen1->IRandom(0,nSquares-1);
                  pSquare_mother = SquareVec[irand];
                  sq_ntrees = pSquare_mother->TreeList.size();
               } while (sq_ntrees == 0);

               //random tree in Square
               SpecID_mother = pSquare_mother->TreeList[RandGen1->IRandom(0,sq_ntrees-1)]->SpecID;
            }

            else {
               //local dispersal from one of the eight neighbors
               SpecID_mother = GetSpecID_mother(x_new,y_new);
            }
         } // if not immigration

         MaxTreeID++;
         pTree = new CTree(MaxTreeID,x_new,y_new,SpecID_mother);
         pSquare->RecruitList.push_back(pTree);
         pSquare->nRecruits_t1++;

         countRecSq++;
         countRecPlot++;

      } while ((countRecSq < nRecSq) && (countRecPlot < nrecruit));

      isquare++;
      if (isquare == nSquares){
         random_shuffle(SquareVec.begin(),SquareVec.end());
         isquare = 0;
      }

   } while (countRecPlot < nrecruit);
}


//---------------------------------------------------------------------------
void CForest::Mort_Dens_Data3 (int nmort)
{
   CSquare *pSquare;
   CTree *pTree;

   int countMortPlot = 0;
   int countMortSq = 0;
   int nMortSq;
   int isquare = 0;

   int itree, sq_ntrees;
   int dens_cl, nsq_dens;

   random_shuffle(SquareVec.begin(),SquareVec.end());  //random permutation of square order

   do {

      pSquare = SquareVec[isquare];

      //dens_cl = (int) floor(pSquare->nTrees_t0/DensBin);
      dens_cl = (int) floor(pSquare->TreeList.size()/DensBin);

      if (dens_cl < DensCl_Min) dens_cl = DensCl_Min;
      if (dens_cl > DensCl_Max) dens_cl = DensCl_Max;

      nsq_dens = nDeathObsDens[dens_cl].size();
      nMortSq = nDeathObsDens[dens_cl][RandGen1->IRandom(0,nsq_dens-1)];
      sq_ntrees = pSquare->TreeList.size();

      if (nMortSq<pSquare->TreeList.size()){

      countMortSq = 0;
      //sq_ntrees = SquareVec[isquare]->TreeList.size();

      while (((countMortSq < nMortSq) && (pSquare->TreeList.size()>0)) && (countMortPlot < nmort)){

         //random tree from square
         sq_ntrees = pSquare->TreeList.size();
         itree = RandGen1->IRandom(0,sq_ntrees-1);
         pTree = pSquare->TreeList[itree];

         //update local abundances
         pSquare->SpecMap[pTree->SpecID]--;
         if (pSquare->SpecMap[pTree->SpecID] == 0)
            pSquare->SpecMap.erase(pTree->SpecID);

         pSquare->nDeath_t1++;

         //update forest abundances
         SpecAbund[pTree->SpecID]--;
         if (SpecAbund[pTree->SpecID] == 0)
            SpecAbund.erase(pTree->SpecID);

         nTreesPlot--;
         countMortSq++;
         countMortPlot++;

         //remove tree
         pSquare->TreeList.erase(pSquare->TreeList.begin()+itree);
         delete pTree;

         //kill and remove tree
//         if (pTree->dead == false){
//            pTree->dead = true;
//
//            pSquare->nTrees_t0--;
//            pSquare->nDeath_t1++;
//
//            countMortSq++;
//            countMortPlot++;
//         } // if tree not dead

      } // while ((countMortSq < nMortSq) && (countMortPlot < nmort)){

      } // if nMortSq < nTrees in square

      isquare++;
      if (isquare == nSquares){
        random_shuffle(SquareVec.begin(),SquareVec.end());
        isquare = 0;
      }

   } while (countMortPlot < nmort);
}

//---------------------------------------------------------------------------
void CForest::Rec_Dens_Data4 (int nrecruit)
{
   CTree* pTree;
   CSquare* pSquare;
   CSquare* pSquare_mother;

   //int iX, iY;
   int SpecID_mother, irand;
   double x_new, y_new, rand1;

   int countRecPlot = 0;
   //int countRecSq = 0;
   int isquare = 0;
   //int nRecSq;
   int sq_ntrees;

   int dens_cl, nsq_dens;

   double pRecSq;

   //random_shuffle(SquareVec.begin(),SquareVec.end());  //random permutation of square order

   // run over all squares and assign recruitment rate
   for (SquareIterV isq = SquareVec.begin(); isq != SquareVec.end();isq++){
      pSquare = *isq;

      dens_cl = (int) floor(pSquare->TreeList.size()/DensBin);

      if (dens_cl < DensCl_Min) dens_cl = DensCl_Min;
      if (dens_cl > DensCl_Max) dens_cl = DensCl_Max;

      nsq_dens = pRecruitObsDens[dens_cl].size();
      pRecSq = pRecruitObsDens[dens_cl][RandGen1->IRandom(0,nsq_dens-1)];

      pSquare->pRec = pRecSq;
   }

   //generate recruits in random order
   while (countRecPlot < nrecruit){

      //randomly sample from quadrats with replacement
      isquare = RandGen1->IRandom(0,nSquares-1);
      pSquare = SquareVec[isquare];

      if (RandGen1->Random() < pSquare->pRec){

         //random position in cell
         x_new = (pSquare->iX + RandGen1->Random() )*SquareSize;
         y_new = (pSquare->iY + RandGen1->Random() )*SquareSize;

         //find mother tree
         rand1 = RandGen1->Random();

         //immigration of novel species
         if (rand1 < Pars->pImmi){
            SpecMax++;
            SpecID_mother = SpecMax;
         }

         else {
         //long-distance dispersal in the plot
            if ((rand1 >= Pars->pImmi) & (rand1 < (Pars->pImmi+Pars->pLDD))){
               //random square with trees
               do {
                  irand = RandGen1->IRandom(0,nSquares-1);
                  pSquare_mother = SquareVec[irand];
                  sq_ntrees = pSquare_mother->TreeList.size();
               } while (sq_ntrees == 0);

               //random tree in Square
               SpecID_mother = pSquare_mother->TreeList[RandGen1->IRandom(0,sq_ntrees-1)]->SpecID;
            }

            else {
               //local dispersal from one of the eight neighbors
               SpecID_mother = GetSpecID_mother(x_new,y_new);
            }
         } // if not immigration

         MaxTreeID++;
         pTree = new CTree(MaxTreeID,x_new,y_new,SpecID_mother);
         pSquare->RecruitList.push_back(pTree);
         pSquare->nRecruits_t1++;

         countRecPlot++;
      } //if (RandGen1->Random() < pSquare->pRec)

//      isquare++;
//      if (isquare == nSquares){
//         random_shuffle(SquareVec.begin(),SquareVec.end());
//         isquare = 0;
//      }
   } //while (countRecPlot < nrecruit)
}

//---------------------------------------------------------------------------
void CForest::Mort_Dens_Data4a (int nmort)
{
   CSquare *pSquare;
   CTree *pTree;

   int countMortPlot = 0;
   //int nMortSq;
   int isquare = 0;

   int itree, sq_ntrees;
   int dens_cl, nsq_dens;

   double pMortSq;

   //int sum_mort_pot = 0;

   // run over all squares and assign potential mortality
   for (SquareIterV isq = SquareVec.begin(); isq != SquareVec.end();isq++){
      pSquare = *isq;

      dens_cl = (int) floor(pSquare->TreeList.size()/DensBin);

      if (dens_cl < DensCl_Min) dens_cl = DensCl_Min;
      if (dens_cl > DensCl_Max) dens_cl = DensCl_Max;

      nsq_dens = pDeathObsDens[dens_cl].size();
      //nMortSq = nDeathObsDens[dens_cl][RandGen1->IRandom(0,nsq_dens-1)];
      pMortSq = pDeathObsDens[dens_cl][RandGen1->IRandom(0,nsq_dens-1)];

      //pSquare->nDeath_pot = nMortSq;
      //sum_mort_pot += nMortSq;

      //pSquare->pDeath = (double) nMortSq / pSquare->TreeList.size();
      pSquare->pDeath = pMortSq;
   }

//   int nMortPlot;
//   if (sum_mort_pot<nmort)
//      nMortPlot = sum_mort_pot;
//   else
//      nMortPlot = nmort;

   random_shuffle(SquareVec.begin(),SquareVec.end());  //random permutation of square order

   // kill trees in random order
   //while (countMortPlot < nmort){

   do {
       //randomly sample from quadrats with replacement
      //isquare = RandGen1->IRandom(0,nSquares-1);
      pSquare = SquareVec[isquare];

      //random tree from square
      sq_ntrees = pSquare->TreeList.size();
      if (sq_ntrees > 0){
         itree = RandGen1->IRandom(0,sq_ntrees-1);
         pTree = pSquare->TreeList[itree];

         //check if tree dies
         if (RandGen1->Random() < pSquare->pDeath){
            //update local abundances
            pSquare->SpecMap[pTree->SpecID]--;
            if (pSquare->SpecMap[pTree->SpecID] == 0)
               pSquare->SpecMap.erase(pTree->SpecID);

            pSquare->nDeath_t1++;

            //update forest abundances
            SpecAbund[pTree->SpecID]--;
            if (SpecAbund[pTree->SpecID] == 0)
               SpecAbund.erase(pTree->SpecID);

            nTreesPlot--;
            countMortPlot++;

            //remove tree
            pSquare->TreeList.erase(pSquare->TreeList.begin()+itree);
            delete pTree;
         }
      } // if sq_nTrees >0

      ++isquare;
      if (isquare == nSquares){
        //random_shuffle(SquareVec.begin(),SquareVec.end());
        isquare = 0;
      }

   //} //while (countMortPlot < nMortPlot){

   } while (countMortPlot < nmort);

}

//---------------------------------------------------------------------------
void CForest::Mort_Dens_Data4b (int nmort)
{
   CSquare *pSquare;
   CTree *pTree;

   int countMortPlot = 0;

   int isquare = 0;

   int itree, sq_ntrees;
   int dens_cl, nsq_dens;

   double pMortSq;


   // run over all squares and assign potential mortality
   for (SquareIterV isq = SquareVec.begin(); isq != SquareVec.end();isq++){
      pSquare = *isq;

      dens_cl = (int) floor(pSquare->TreeList.size()/DensBin);

      if (dens_cl < DensCl_Min) dens_cl = DensCl_Min;
      if (dens_cl > DensCl_Max) dens_cl = DensCl_Max;

      nsq_dens = pDeathObsDens[dens_cl].size();
      pMortSq = pDeathObsDens[dens_cl][RandGen1->IRandom(0,nsq_dens-1)];

      pSquare->pDeath = pMortSq;
   }

   random_shuffle(SquareVec.begin(),SquareVec.end());  //random permutation of square order
   //while (countMortPlot < nmort){
   do {

       //randomly sample from quadrats with replacement
      //isquare = RandGen1->IRandom(0,nSquares-1);
      pSquare = SquareVec[isquare];

      sq_ntrees = pSquare->TreeList.size();
      if (sq_ntrees > 0){

         //check all trees from square
         for (TreeIterV itree = pSquare->TreeList.begin();
                        itree != pSquare->TreeList.end();){

            pTree = (*itree);

            //check if tree dies
            if (RandGen1->Random() < pSquare->pDeath){
               //update local abundances
               pSquare->SpecMap[pTree->SpecID]--;
               if (pSquare->SpecMap[pTree->SpecID] == 0)
                  pSquare->SpecMap.erase(pTree->SpecID);

               pSquare->nDeath_t1++;

               //update forest abundances
               SpecAbund[pTree->SpecID]--;
               if (SpecAbund[pTree->SpecID] == 0)
                  SpecAbund.erase(pTree->SpecID);

               nTreesPlot--;
               countMortPlot++;

               //remove tree
               itree = pSquare->TreeList.erase(itree);
               delete pTree;
            } // if death
            else {
               ++itree;
            }
         } //for all trees in quadrat
      } // if sq_nTrees >0

      isquare++;
      if (isquare == nSquares){
        //random_shuffle(SquareVec.begin(),SquareVec.end());
        isquare = 0;
      }

   } while (countMortPlot < nmort);

   //} //while (countMortPlot < nMortPlot){
}

//---------------------------------------------------------------------------
void CForest::Recruit_Mort_OneStep()
{
   CSquare *pSquare;

   //Set nRecruits and nDeath to zero
   //Set number of trees to length of local tree list
   for (SquareIterV isq = SquareVec.begin(); isq != SquareVec.end();isq++){
      pSquare = *isq;

      pSquare->nTrees_t0 = pSquare->TreeList.size();
      pSquare->nTrees = pSquare->nTrees_t0;

      pSquare->nSpec_t0 = pSquare->SpecMap.size();

      pSquare->nRecruits_t1 = 0;
      pSquare->nDeath_t1 = 0;
   }

   double meanMortRec = nTreesPlot*Pars->pMortRec;    //mean number of mortality-recruitment events in 5 years
   int nMortRec = RandGen2->Poisson(meanMortRec); //number of mortality-recruitment for current time step

   switch (Pars->Scenario){
      case 1:
         recruit_neutral_random_mother_fw(nMortRec);
         mort_neutral_random_tree(nMortRec);
         break;
      case 2:
         recruit_data_bw(nMortRec,false,false);
         mort_neutral_random_tree(nMortRec);
         break;
      case 3:
         recruit_data_bw(nMortRec,false,true);
         mort_neutral_random_tree(nMortRec);
         break;
      case 4:
         recruit_neutral_random_mother_fw(nMortRec);
         mort_data_cluster2(nMortRec,false,false);
         break;
      case 5:
         recruit_data_bw(nMortRec,false,false);
         mort_data_cluster2(nMortRec,false,false);
         break;
      case 6:
         recruit_data_bw(nMortRec,false,true);
         mort_data_cluster2(nMortRec,false,false);
         break;
      case 7:
         recruit_neutral_random_mother_fw(nMortRec);
         mort_data(nMortRec,true,false);
         break;
      case 8:
         recruit_data_bw(nMortRec,false,false);
         mort_data(nMortRec,true,false);
         break;
      case 9:
         recruit_data_bw(nMortRec,false,true);
         mort_data(nMortRec,true,false);
         break;
      default:
         recruit_neutral_random_mother_fw(nMortRec);
         mort_neutral_random_tree(nMortRec);
   }
}

//---------------------------------------------------------------------------
void CForest::UpdateSquares(){

   CSquare *pSquare = nullptr;
   CTree *pTree = nullptr;

   //Change recruits to trees and delete dead trees
   for (SquareIterV isq = SquareVec.begin(); isq != SquareVec.end(); ++isq){
      pSquare = *isq;

      //update local species abundances with recruits
      for (TreeIterV itree =  pSquare->RecruitList.begin();
                     itree !=  pSquare->RecruitList.end(); itree++){
         ++pSquare->SpecMap[(*itree)->SpecID];
         ++SpecAbund[(*itree)->SpecID];
      }

      //change recruits to trees
      pSquare->TreeList.insert(pSquare->TreeList.end(),
                               pSquare->RecruitList.begin(),
                               pSquare->RecruitList.end());
      nTreesPlot += pSquare->RecruitList.size();
      pSquare->RecruitList.clear();
   } // for squares
}

//---------------------------------------------------------------------------
//neutral model
/*
void CForest::Recruit_Mort_M0()
{
   //TreeIterV itree;
   CTree* pTree1;
   CTree* pTree2;
   CSquare* pSquare;

   int iX, iY, SpecID_mother;
   double x_new, y_new, x_mother, y_mother, rand1;

   //Set nRecruits and nDeath to zero
   for (SquareIterV isq = SquareVec.begin(); isq != SquareVec.end();isq++){
      pSquare = *isq;
      pSquare->nRecruits = 0;
      pSquare->nDeath = 0;
   }

   double meanMortRec = NTrees*Pars->pMortRec;    //mean number of mortality-recruitment events in 5 years
   int nMortRec = RandGen2->Poisson(meanMortRec); //number of mortality-recruitment for current time step

   //simulate mortality & recruitment events until nMortRec is reached

   //1. Recruitment
   int countRec = 0;
   int isquare = 0;

   random_shuffle(SquareVec.begin(),SquareVec.end());  //random permutation of square order

   do {

      pSquare = SquareVec[isquare];

      //draw number of recruits from

      for (TreeIterV itree = pSquare->TreeList.begin();
                     itree != pSquare->TreeList.end();itree++){
         if (RandGen1->Random() < Pars->pMortRec){
            pTree1 = (*itree);

            //generate recruit
            rand1 = RandGen1->Random();

            //(a) immigration of novel species
            if (rand1 < Pars->pImmi){
               SpecMax++;
               SpecID_mother = SpecMax;

               //random position
               x_new = RandGen1->Random()*Xmax;
               y_new = RandGen1->Random()*Ymax;

               iX = (int) floor(pTree1->X/SquareSize);
               iY = (int) floor(pTree1->Y/SquareSize);

            }
            else {
               SpecID_mother = pTree1->SpecID;

               //(b) long-distance dispersal in the plot
               if ((rand1 >= Pars->pImmi) & (rand1 < (Pars->pImmi+Pars->pLDD))){
                  //random position
                  x_new = RandGen1->Random()*Xmax;
                  y_new = RandGen1->Random()*Ymax;

                  iX = (int) floor(pTree1->X/SquareSize);
                  iY = (int) floor(pTree1->Y/SquareSize);
               } // if global dispersal
               else {

                  // (c) local dispersal
                  x_mother = pTree1->X;
                  y_mother = pTree1->Y;
                  x_new = x_mother;
                  y_new = y_mother;

                  //generate new position for recruit
                  GetNewXY(x_new,y_new);

                  iX = (int) floor(pTree1->X/SquareSize);
                  iY = (int) floor(pTree1->Y/SquareSize);
               } //local dispersal
            } //else

            pTree2 = new CTree(x_new,y_new,SpecID_mother);
            SquareGrid[iX][iY]->RecruitList.push_back(pTree2);

            countRec++;
            if (countRec == nMortRec) break; //finishes FOR loop

         } //if recruitment event
      } // for trees in square

      isquare++;
      if (isquare == nSquares) isquare = 0;
   } while (countRec < nMortRec);


   //2. Mortality
   random_shuffle(SquareVec.begin(),SquareVec.end());  //random permutation of square order

   int countMort = 0;
   isquare = 0;

   do {

      pSquare = SquareVec[isquare];

      for (TreeIterV itree = pSquare->TreeList.begin();
                     itree != pSquare->TreeList.end();){
         if (RandGen1->Random() < Pars->pMortRec){
            pTree1 = (*itree);

            //update square
            pSquare->SpecMap[pTree1->SpecID]--;
            if (pSquare->SpecMap[pTree1->SpecID] == 0)
               pSquare->SpecMap.erase(pTree1->SpecID);

            pSquare->TreeList.erase(itree);
            pSquare->nDeath++;

            //update forest
            SpecAbund[pTree1->SpecID]--;
            if (SpecAbund[pTree1->SpecID] == 0)
               SpecAbund.erase(pTree1->SpecID);
            NTrees--;
            delete pTree1;

            countMort++;
            if (countMort == nMortRec) break; //finishes FOR loop
         } //mortality event
         else { itree++;}
      }  // for trees in square
      isquare++;
      if (isquare == nSquares) isquare = 0;
   } while (countMort < nMortRec);

   //Change recruits to trees
   for (SquareIterV isq = SquareVec.begin(); isq != SquareVec.end();isq++){
      pSquare = *isq;

      for (TreeIterV itree =  pSquare->RecruitList.begin();
                     itree !=  pSquare->RecruitList.end(); itree++){
         pSquare->SpecMap[(*itree)->SpecID]++;
         SpecAbund[(*itree)->SpecID]++;
      }

      pSquare->nRecruits =  pSquare->RecruitList.size();
      NTrees = NTrees + pSquare->nRecruits;

      pSquare->TreeList.insert(pSquare->TreeList.end(),
                               pSquare->RecruitList.begin(),
                               pSquare->RecruitList.end());
      pSquare->RecruitList.clear();
   }
};
*/

//---------------------------------------------------------------------------
//Model with uncorrelated recruitment and death events from data
/*
void CForest::Recruit_Mort_M1a()
{
   //TreeIterV itree;
   CTree* pTree;
   CSquare* pSquare;

   int iX, iY, SpecID_mother, irand;
   double x_new, y_new, rand1;

   //Set nRecruits and nDeath to zero
   for (SquareIterV isq = SquareVec.begin(); isq != SquareVec.end();isq++){
      pSquare = *isq;
      pSquare->nRecruits = 0;
      pSquare->nDeath = 0;
   }

   double meanMortRec = NTrees*Pars->pMortRec;    //mean number of mortality-recruitment events in 5 years
   int nMortRec = RandGen2->Poisson(meanMortRec); //number of mortality-recruitment for current time step

   //simulate mortality & recruitment events until nMortRec is reached

   //1. Recruitment
   int countRecPlot = 0;
   int countRecSq = 0;
   int isquare = 0;
   int nRecSq, sq_ntrees;

   random_shuffle(SquareVec.begin(),SquareVec.end());  //random permutation of square order

   do {

      pSquare = SquareVec[isquare];

      nRecSq = nRecruitObs[RandGen1->IRandom(0,nSquares-1)];
      countRecSq = 0;

      do {

         //random position in cell
         x_new = (pSquare->iX + RandGen1->Random() )*SquareSize;
         y_new = (pSquare->iY + RandGen1->Random() )*SquareSize;

         //find mother tree
         rand1 = RandGen1->Random();

         //immigration of novel species
         if (rand1 < Pars->pImmi){
            SpecMax++;
            SpecID_mother = SpecMax;
         }

         else {
         //long-distance dispersal in the plot
            if ((rand1 >= Pars->pImmi) & (rand1 < (Pars->pImmi+Pars->pLDD))){
               //random square with trees
               do {
                  irand = RandGen1->IRandom(0,nSquares-1);
                  pSquare = SquareVec[irand];
                  sq_ntrees = pSquare->TreeList.size();
               } while (sq_ntrees == 0);

               //random tree in Square
               SpecID_mother = pSquare->TreeList[RandGen1->IRandom(0,sq_ntrees-1)]->SpecID;
            }

            else {
               //local dispersal from one of the eight neighbors
               SpecID_mother = GetSpecID_mother(x_new,y_new);
            }
         } // if not immigration

         pTree = new CTree(x_new,y_new,SpecID_mother);
         pSquare->RecruitList.push_back(pTree);
         pSquare->nRecruits++;
         countRecSq++;
         countRecPlot++;

      } while ((countRecSq < nRecSq) && (countRecPlot < nMortRec));

      isquare++;
      if (isquare == nSquares) isquare = 0;

   } while (countRecPlot < nMortRec);

   //2. Mortality
   random_shuffle(SquareVec.begin(),SquareVec.end());  //random permutation if square order

   int countMortPlot = 0;
   int countMortSq = 0;
   int nMortSq;
   isquare = 0;

   int itree;

   do {

      pSquare = SquareVec[isquare];

      nMortSq = nDeathObs[RandGen1->IRandom(0,nSquares-1)];
      countMortSq = 0;

      while ((countMortSq < nMortSq) && (countMortPlot < nMortRec) && (pSquare->TreeList.size()>0)){
         //random tree from square
         sq_ntrees = pSquare->TreeList.size();
         itree = RandGen1->IRandom(0,sq_ntrees-1);
         pTree = pSquare->TreeList[itree];

         //update square
         pSquare->SpecMap[pTree->SpecID]--;
         if (pSquare->SpecMap[pTree->SpecID] == 0)
            pSquare->SpecMap.erase(pTree->SpecID);

         pSquare->TreeList.erase(pSquare->TreeList.begin()+itree);
         pSquare->nDeath++;

         //update forest
         SpecAbund[pTree->SpecID]--;
         if (SpecAbund[pTree->SpecID] == 0)
            SpecAbund.erase(pTree->SpecID);
         NTrees--;

         delete pTree;

         countMortSq++;
         countMortPlot++;
      }

      isquare++;
      if (isquare == nSquares) isquare = 0;

   } while (countMortPlot < nMortRec);


   //3. Change recruits to trees
   for (SquareIterV isq = SquareVec.begin(); isq != SquareVec.end();isq++){
      pSquare = *isq;

      for (TreeIterV itree =  pSquare->RecruitList.begin();
                     itree !=  pSquare->RecruitList.end(); itree++){
         pSquare->SpecMap[(*itree)->SpecID]++;
         SpecAbund[(*itree)->SpecID]++;
      }

      pSquare->nRecruits =  pSquare->RecruitList.size();
      NTrees = NTrees + pSquare->nRecruits;

      pSquare->TreeList.insert(pSquare->TreeList.end(),
                               pSquare->RecruitList.begin(),
                               pSquare->RecruitList.end());
      pSquare->RecruitList.clear();

   } //for squares
};
*/

// ---------------------------------------------------------------------------
void CForest::OneRun(int isim, int irep, string label)
{
   Init();

//   if (time_output==true){
//      WriteOutput(0,isim,irep);
//      WriteSquares(0,isim,irep);
//   }

   for (int tstep = 1; tstep <= (Pars->Tmax + 5); tstep++){

      Recruit_Mort_OneStep();
      UpdateSquares();
      //cout << "     Step " << tstep << endl;

      if (time_output==true){
         //if (tstep == 1) WriteSquares(tstep,isim,irep);
         //if ((tstep >= 0) && (tstep % 100 == 0))
         if (tstep >= Pars->Tmax){
            cout << "     Step " << tstep << endl;
            WriteOutput(tstep, isim, irep);
            WriteSquares(tstep, isim, irep, label);
            //WriteTrees(tstep, isim, irep, label);
         }
      }
   }

   if (!time_output){
      WriteSquares(Pars->Tmax + 5,isim,irep,label);
      WriteOutput(Pars->Tmax + 5,isim,irep);
   }
}

// ---------------------------------------------------------------------------
void CForest::WriteOutput(int step, int sim = 1, int rep = 1) {

   // Community data  ----------------------------------------------
   double Shannon = GetShannon();
   int NSpec = GetSAD();

   DivFile << sim  << "; "
           << rep  << "; "
           << step << "; "
           << nTreesPlot<< "; "
           << NSpec << "; "
           << Shannon
           << endl;

   /*
   for (int i = 0; i < (MaxSAD-1); ++i)
      SAD_File << SAD[i] << "; ";
   SAD_File << SAD[MaxSAD-1] << endl;
   */


   // Species data ------------------------------------------------------------
   /*
   map<int, int>::iterator spec_it;

   for (spec_it = SpecAbund.begin(); spec_it!=SpecAbund.end(); ++spec_it)
      AbundFile << spec_it->second <<"; ";
   AbundFile<<endl;
   */
}

// ---------------------------------------------------------------------------
void CForest::WriteTrees(int step, int sim=1, int rep=1, string label="1") {

   //string FileName = "InOut\\TreesOut" + IntToString(tstep) + ".csv";

   string FileName = "InOut/TreesOut"    + label + "_"
                                          + IntToString(sim) + "_"
                                          + IntToString(rep) + "_"
                                          + IntToString(step) + ".csv";

   ofstream OutFile1(FileName.c_str());
   OutFile1 << "Nr;"
            << "X;"
            << "Y;"
            << "SpecID"
            << endl;

   CTree* pTree;

   TreeIterV itree;

   int i=0;
   for (int iX = 0; iX < XSquares; iX++){
      for (int iY = 0; iY < YSquares; iY++){

         for (TreeIterV itree = SquareGrid[iX][iY]->TreeList.begin();
                        itree != SquareGrid[iX][iY]->TreeList.end(); itree++){

            pTree = (*itree);

            OutFile1 << pTree->TreeID << ";"
                     << pTree->X << ";"
                     << pTree->Y << ";"
                     << pTree->SpecID
                     << endl;
            i++;

         } // for trees in square


      } // iY
   } //iX

   OutFile1.close();
}

// ---------------------------------------------------------------------------
void CForest::WriteSquares(int step, int sim=1, int rep=1, string label="1") {

   string FileName = "InOut/SquaresOut"  + label + "_"
                                          + IntToString(sim) + "_"
                                          + IntToString(rep) + "_"
                                          + IntToString(step) + ".csv";

   ofstream OutFile1(FileName.c_str());
   OutFile1 //<< "SquareID;"
            //<< "X;"
            //<< "Y;"
            << "nTrees;"
            << "nSpecies;"
            << "nDeath;"
            << "nRecruits"
            //<< "pDeath;"
            //<< "pRec"
            << endl;

   //int iSq = 0;
   for (int iX = 0; iX < XSquares; iX++){
      for (int iY = 0; iY < YSquares; iY++){
         OutFile1 //<< iSq << ";"
                  //<< SquareGrid[iX][iY]->X << ";"
                  //<< SquareGrid[iX][iY]->Y << ";"
                  << SquareGrid[iX][iY]->nTrees_t0 << ";"
                  << SquareGrid[iX][iY]->nSpec_t0 << ";"
                  << SquareGrid[iX][iY]->nDeath_t1 << ";"
                  << SquareGrid[iX][iY]->nRecruits_t1
                  //<< SquareGrid[iX][iY]->pDeath << ";"
                  //<< SquareGrid[iX][iY]->pRec
                  << endl;
         //iSq++;
      }
   }
   OutFile1.close();
}

// ---------------------------------------------------------------------------
double CForest::GetShannon() {

   map<int,int>::iterator spec_it;

   double shannon = 0, relabund;
   double sum = 0.0;

   int ntrees = 0;
   for (spec_it = SpecAbund.begin(); spec_it != SpecAbund.end(); ++spec_it) {
      ntrees = ntrees + spec_it->second;
   }

   for (spec_it = SpecAbund.begin(); spec_it != SpecAbund.end(); ++spec_it) {
      if (spec_it->second > 0) {
         relabund = (double) spec_it->second / ntrees;
         shannon += log(relabund) * relabund;
         sum += relabund;
      }
   }

   return(-shannon);
}

// ---------------------------------------------------------------------------
int CForest::GetSAD() {

   map<int,int>::iterator spec_it;

   double Abund;
   int log2Abund;
   int nspec = 0;


   for (int i = 0; i < MaxSAD; ++i)
      SAD[i] = 0;

   for (spec_it = SpecAbund.begin(); spec_it != SpecAbund.end(); ++spec_it) {
      Abund = spec_it->second;
      if (Abund > 0) {
         ++nspec;
         log2Abund = (int) floor(log(Abund) / log(2.0));
         // log2Abund = (int) ceil(log(Abund)/log(2.0));
         if (log2Abund >= MaxSAD)
            ++SAD[MaxSAD - 1];
         else
            ++SAD[log2Abund];
      }
   }

   return(nspec);
}

// ---------------------------------------------------------------------------
int CForest::GetRandSquare() {

	int isq = 0;
	double cum_prob = static_cast<double>(SquareVec[0]->nTrees)/nTreesPlot;

	double r1 = RandGen1->Random();
	bool choose = false;

	while (choose == false) {
		if (r1 <= cum_prob)
			choose = true;
		else {
			++isq;
			cum_prob += static_cast<double>(SquareVec[isq]->nTrees)/nTreesPlot;
      }
	}

	return(isq);
}
