#include "TH1.h"
#include "TGraph.h"
#include "TAttLine.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TGraphEditor.h"
#include "TNamed.h"
#include "TObject.h"
#include "TStyle.h"
#include "TAttFill.h"
#include "TAttMarker.h"
#include <iostream>
#include <fstream>
#include "stdio.h"
#include "stdlib.h"
#include <cmath>
#include "TPaveLabel.h"
#include "TApplication.h"

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

using namespace std;

void getavgeff(int run1, int mode=1) //mode: 1=DAQevent, 2=MIDAS standalone, 3=MIDAS)
{
   
    if (mode == 1) {
        exit;
    }
    else if (mode == 2) {
        exit;
    }
    else if (mode == 3) {
        TFile* file1 = TFile::Open(Form("./trackedfiles/GEM_tracked0%d.root",run1));
        TTree *MUSEteleTracks = (TTree*)file1->Get("MUSEteleTracks");
       // MUSEteleTracks->AddFriend("lumigemcooked",Form("./cookedfiles/cookedGEMrun0%d.root",run1));
    }
    else {
        printf("Mode not given");
        exit;
    }
   
    /////JUst any cluster hit on GEMs ///////////////////////////////
    

    TH1 *tracksUSall;
    gDirectory->GetObject("Efficiency/US/trackmultiplicityUS",tracksUSall);

    TH1 *tracksUSproj;
    gDirectory->GetObject("MUSEteleTracker/Number of possible clusters US",tracksUSproj);
  
    Double_t USX_n = tracksUSproj->GetEntries();
    cout<< "GEM US projected tracks: " << USX_n <<endl;
    Double_t USX_d = tracksUSall->GetEntries();
    cout<<"GEM US all tracks: " << USX_d << endl;
    Double_t USX_eff=USX_n/USX_d;
   
    Double_t USerr=1/USX_d*(sqrt(USX_n*(1-USX_eff)));
    cout<<" GEM US Efficiency    : "<<USX_eff<<" +/- "<<USerr<<endl;
   
   
}


/*
 Useful functions from MUSEteleTracker
 
 
 
*/


