#ifndef __MUSETELETRACKER__
#define __MUSETELETRACKER__

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TH1I.h"
#include "TH2I.h"
#include <iostream>

#include "lumigemtree.h"
#include "teletracktree.h"
#include "museadctree.h"
#include "musetdc1190tree.h"
#include<muserawtree.h>


class MUSEteleTracker:public Plugin
{
 private:
  LumiGEM    *clusters;
  TeleTracks *teletracks;


// testbeamanalysistree *tdc1190;


 public:
  MUSEteleTracker(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~MUSEteleTracker();
  // add funtions with return value Long_t here:
  
  Long_t defineHistograms();
  Long_t startup();
  Long_t process();
  Long_t finalize();

  Long_t process_ethan();
  Long_t finalize_ethan();
 // Methods for efficiency (in addition to normal functions above!)
  Long_t startup_efficiency();
  Long_t histos_efficiency(); //make histograms
  //  Long_t get_clusters(); // get cluster info
  Long_t process_efficiency(); // filter and fill histograms
  Long_t finalize_efficiency(); //get the average effifiency


  //Some APVs did not work on DS GEM at DEcember, 2017 data. This flag turn DS GEM OFF by taking in account to define the tracks. 
   bool DSOFF;

  virtual Long_t cmdline(char * cmd);

  ClassDef(MUSEteleTracker,1);

  //Definitions Ishara
  //int *totaltracks,*falsetracksUS,*falsetracks4th,*falsetracksMS;
  int event=0;
  int totaltracks=0;
  int falsetracksUS=0;
  int falsetracks4th=0;
  int falsetracksMS=0;

  int truehitsUS=0;
  int truehits4th=0;
  int truehitsMS=0;

  double threshold_residue=50.0;

  int hits0=0;
  int hits1=0;
  int hits2=0;

  int gem2realhits[3]={0,0,0};

  int truetracksUS=0;
  int truetracks4th=0;
  int truetracksMS=0;

  float pre_eff_US=0;
  float pre_eff_4th=0;
  float pre_eff_MS=0;

  int Eff_Den=0;
  int US_count=0;
  int Fourth_count=0;
  int MS_count=0;

  double pre_overall_eff=0;

  int inefficiency_gem2=0;
  

  TH1D *multiplicityUS,*multiplicity4TH,*multiplicityMS,*multiplicityDS,*projectedX;
  TH2D *mult2dUS_MS,*mult2dUS_DS,*mult2dMS_DS;
  TH1D *trackmultiplicityUS,*trackmultiplicityMS,*trackmultiplicityDS;


  TH2D *effgemUS, *effgem4TH, *effgemMS;

  TH1D *xresidual;
  TH1D *yresidual;
  //  TH1D *USX_eff,*MSX_eff,*DSX_eff;
  //TH1D *USY_eff,*MSY_eff,*DSY_eff;

  TH1D *tracksprojectedgemUSX_onecut,*tracksprojectedgemUSX_any;
  TH1D *tracksprojectedgemMSX_onecut,*tracksprojectedgemMSX_any;
  TH1D *tracksprojectedgemDSX_onecut,*tracksprojectedgemDSX_any;

  TH1D *tracksprojectedgemUSY_onecut,*tracksprojectedgemUSY_any;
  TH1D *tracksprojectedgemMSY_onecut,*tracksprojectedgemMSY_any;
  TH1D *tracksprojectedgemDSY_onecut,*tracksprojectedgemDSY_any;

  
  
  
  TH1D * totresx[4]; 
  TH1D * totresy[4];
  TH1D * threehitresx[4];
  TH1D * threehitresy[4];
  TH1D * fourhitresx[4];
  TH1D * fourhitresy[4];

  TH1D *totChi2x, *totChi2y, *threehitChi2x, *threehitChi2y, *fourhitChi2x, *fourhitChi2y;
  
  TH1D * firsthitresx[4];
  TH1D * firsthitresy[4];
  TH1D *firsthitchi2x, *firsthitchi2y;

  TH1D * firsthitresxcut[4];
  TH1D * firsthitresycut[4];

  TH1D * threechi2cutxyx[4];
  TH1D * threechi2cutxyy[4];
  TH1D * fourchi2cutxyx[4];
  TH1D * fourchi2cutxyy[4];
  
  TH1D * totxdist[4];
  TH1D * totydist[4];
  TH1D * multixdist[3][4];
  TH1D * multiydist[3][4];
  
  TH2D * tot2dmap[4];
  TH2D * multi2dmap[3][4];

  TH1D *totmx, *totmy, *totZX, *totZY;
  TH1D * multimx[4];
  TH1D * multimy[4];
  TH1D * multiZX[4];
  TH1D * multiZY[4];

  TH1D * largestremovedx[4][4];
  TH1D * largestremovedy[4][4];
  TH1D * largestremovedmx[4];
  TH1D * largestremovedmy[4];
  TH1D * largestremovedZX[4];
  TH1D * largestremovedZY[4];
  
  TH1D * resremovedx[4];
  TH1D * resremovedy[4];
  TH1D * resweightedx[4];
  TH1D * resweightedy[4];
  TH2D * h2resmapx[4];
  TH2D * h2resmapy[4];

  TH2D * totbest[4];
  TH2D * multibest[3][4];

  TH2D *h2resmapsx,*h2resmapsy,*singleGEMclust;


};


void getLeastSquaresLine(std::vector<double> x, std::vector<double> y, std::vector<double> z, double &slopezy, double &slopezx,
 double &bzy, double &bzx);

#endif
