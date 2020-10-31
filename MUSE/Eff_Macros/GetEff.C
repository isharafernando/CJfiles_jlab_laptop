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

void GetEff(int run_num)
{
   
    TFile *file1 = new TFile(Form("/home/ishara/MUSE/trackedfiles/GEM_BH_tracked0%d.root",run_num));
    
    TCanvas *C1 = new TCanvas("C1",Form("Effciency Plot"), 2000,2000);

    C1->cd(1);
    TH2F *USEfficiency = new TH2F("USEfficiency","Efficiency on US",50,-50.0,50.0,50,-50.0,50.0);
   
    TH2 *USnumerator; gDirectory->GetObject("Efficiency/US/ActualProjectedGoodClusters_GEM_US",USnumerator);
    TH2 *USdenominator; gDirectory->GetObject("Efficiency/US/tracksprojectedgem_active_US",USdenominator);
    
    USEfficiency->Divide(USnumerator,USdenominator,1,1,"");
    USEfficiency->GetXaxis()->SetTitle("X Position (mm)");
    USEfficiency->GetYaxis()->SetTitle("Y Position (mm)");
    USEfficiency->Draw("colz");
    
    Double_t USX_n = USnumerator->GetEntries();
    Double_t USX_d = USdenominator->GetEntries();
    Double_t USX_eff=USX_n/USX_d;
    
    Double_t USerr=1/USX_d*(sqrt(USX_n*(1-USX_eff)));
    cout << " GEM US Efficiency: " << USX_eff << " +/-" << USerr << endl;

    C1->SaveAs("EfficiencyUS.pdf");
return;
   
}


/*
 Useful functions from MUSEteleTracker
 
 
 
*/


