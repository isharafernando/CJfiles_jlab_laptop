#ifndef __MUSETELETRACKER__
#define __MUSETELETRACKER__

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TH1I.h"
#include "TH2I.h"
#include <iostream>

#include "GEMhittree.h"
#include "teletracktree.h"
#include "Scinttree.h"
#include "Scinthittree.h"
#include<rawtree.h>
#include "Geometry.h"


class MUSEteleTracker:public Plugin
{
 private:
  LumiGEM    *clusters;
  TeleTracks *teletracks;
  ScintHits  *Hits;


// testbeamanalysistree *tdc1190;


 public:
  MUSEteleTracker(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~MUSEteleTracker();
  geometry::geomHandler handle;
  // add funtions with return value Long_t here:
  
  Long_t defineHistograms();
  Long_t startup();
  Long_t process();
  Long_t finalize();

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



  TH1D *multiplicityUS,*multiplicity4TH,*multiplicityMS,*multiplicityDS;
  TH2D *mult2dUS_MS,*mult2dUS_DS,*mult2dMS_DS;
  TH1D *trackmultiplicityUS,*trackmultiplicityMS,*trackmultiplicityDS;

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

  
  
  TH2D * GEMXvsBHX[4];
  TH2D * GEMYvsBHY[4];
  TH2D * GEMXvsBHY[4];
  TH2D * GEMYvsBHX[4];
  TH2D * GEMXvsTotBHbar[4];
  TH2D * GEMYvsTotBHbar[4];

  TH1D * GEMXvsBHbar[16][4];
  TH1D * GEMYvsBHbar[16][4];

  TH2D * projBHXvsBHX[4];
  TH2D * projBHYvsBHY[4];

  TH1D * xres4mmBHcut[4][2][4];
  TH1D * yres4mmBHcut[4][2][4];
  TH1D * xres8mmBHcut[4][2][4];
  TH1D * yres8mmBHcut[4][2][4];

  
  TH1D * totresx[4]; 
  TH1D * totresy[4];
  TH1D * threehitresx[4];
  TH1D * threehitresy[4];
  TH1D * fourhitresx[4];
  TH1D * fourhitresy[4];

  TH1D * totBHresx[4][2];
  TH1D * totBHresy[4][2];
  TH1D * totBHplaneresx[4][2][3];
  TH1D * totBHplaneresy[4][2][3];
  TH1D * totBHparticleresx[4][2][3][3];
  TH1D * totBHparticleresy[4][2][3][3];

  TH1D * fourBHresx[4][2];
  TH1D * fourBHresy[4][2];
  TH1D * fourBHplaneresx[4][2][3];
  TH1D * fourBHplaneresy[4][2][3];
  TH1D * fourBHparticleresx[4][2][3][3];
  TH1D * fourBHparticleresy[4][2][3][3];

  TH1D * threeBHresx[4][2];
  TH1D * threeBHresy[4][2];
  TH1D * threeBHplaneresx[4][2][3];
  TH1D * threeBHplaneresy[4][2][3];
  TH1D * threeBHparticleresx[4][2][3][3];
  TH1D * threeBHparticleresy[4][2][3][3];



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


/***** The following definitions are for GEM-BH calculations (added by Ishara: July 2020) ***************/

  Long_t startup_GEMBH_efficiency();
  Long_t histos_GEMBH_efficiency(); //make histograms
  Long_t process_GEMBH_efficiency(); // filter and fill histograms
  Long_t finalize_GEMBH_efficiency(); //get the average effifiency


  double threshold_residue=10.0;

  float pre_eff_US=0;
  float pre_eff_4th=0;
  float pre_eff_MS=0;


};


void getLeastSquaresLine(std::vector<double> x, std::vector<double> y, std::vector<double> z, double &slopezy, double &slopezx,
 double &bzy, double &bzx);

#endif
