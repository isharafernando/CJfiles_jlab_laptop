
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TChain.h>

using namespace std;

void gem_efficiency()
{
    gSystem->Load("/usr/local/lib/libmusetree.so");
    TFile *file1 = new TFile(Form("/home/ishara/MUSE/trackedfiles/GEM_BH_tracked0%d.root",6290));
    
    //TTree *test1 = (TTree*)file1->Get("MUSEteleTracks");
    TH1F*h1 = (TH1F*)file1->Get("Efficiency/US/multiplicityUS");
    TH1F*h2 = (TH1F*)file1->Get("Efficiency/4TH/multiplicity4TH");
    TH1F*h3 = (TH1F*)file1->Get("Efficiency/MS/multiplicityMS");
    TH1F*h4 = (TH1F*)file1->Get("Efficiency/US/max_charge_clusters_USGEM");
    TH1F*h5 = (TH1F*)file1->Get("Efficiency/4TH/max_charge_clusters_4THGEM");
    TH1F*h6 = (TH1F*)file1->Get("Efficiency/MS/max_charge_clusters_MSGEM");
    
    
    TCanvas *C1 = new TCanvas("C1",Form("Effciency Plots_minADC200"), 3000,1000);
    C1->Divide(3,1);
    
    C1->cd(1);
    TH2F *USany = new TH2F("USany","Efficiencies on US (cluster>=1, cuts)",50,-50.0,50.0,50,-50.0,50.0);
    TH2F *Ftany = new TH2F("4THany","Efficiencies on 4TH (cluster>=1, cuts)",50,-50.0,50.0,50,-50.0,50.0);
    TH2F *MSany = new TH2F("MSany","Efficiencies on MS (cluster>=1, cuts)",50,-50.0,50.0,50,-50.0,50.0);
    
    TH2 *tracksprojectedgemUS_oneormore_cut; gDirectory->GetObject("Efficiency/US/tracksprojectedgemUS_oneormore_cut",tracksprojectedgemUS_oneormore_cut);
    TH2 *tracksprojectedgemUS_any; gDirectory->GetObject("Efficiency/US/tracksprojectedgemUS_any",tracksprojectedgemUS_any);
    
    USany->Divide(tracksprojectedgemUS_oneormore_cut,tracksprojectedgemUS_any,1,1,"");
    //USany->SetMaximum(1);
    USany->GetXaxis()->SetTitle("Hori. Position (mm)");
    USany->GetYaxis()->SetTitle("Vert. Position (mm)");
    USany->Draw("colz");
    
    Double_t USX_n = tracksprojectedgemUS_oneormore_cut->GetEntries();
    Double_t USX_d = tracksprojectedgemUS_any->GetEntries();
    Double_t USX_eff=USX_n/USX_d;
    
    Double_t USerr=1/USX_d*(sqrt(USX_n*(1-USX_eff)));
    cout << " GEM US Efficiency: " << USX_eff << " +/-" << USerr << endl;
    
    C1->cd(2);
    TH2 *tracksprojectedgem4TH_oneormore_cut; gDirectory->GetObject("Efficiency/4TH/tracksprojectedgem4TH_oneormore_cut",tracksprojectedgem4TH_oneormore_cut);
    TH2 *tracksprojectedgem4TH_any; gDirectory->GetObject("Efficiency/4TH/tracksprojectedgem4TH_any",tracksprojectedgem4TH_any);
    
    Ftany->Divide(tracksprojectedgem4TH_oneormore_cut,tracksprojectedgem4TH_any,1,1,"");
    //MSany->SetMaximum(1);
    Ftany->GetXaxis()->SetTitle("Hori. Position (mm)");
    Ftany->GetYaxis()->SetTitle("Vert. Position (mm)");
    Ftany->Draw("colz");
    
    Double_t FtX_n = tracksprojectedgem4TH_oneormore_cut->GetEntries();
    Double_t FtX_d = tracksprojectedgem4TH_any->GetEntries();
    Double_t FtX_eff=FtX_n/FtX_d;
    
    Double_t Fterr=1/FtX_d*(sqrt(FtX_n*(1-FtX_eff)));
    cout << " GEM 4TH Efficiency: "<< FtX_eff << " +/-" << Fterr << endl;
    
    C1->cd(3);
    TH2 *tracksprojectedgemMS_oneormore_cut; gDirectory->GetObject("Efficiency/MS/tracksprojectedgemMS_oneormore_cut",tracksprojectedgemMS_oneormore_cut);
    TH2 *tracksprojectedgemMS_any; gDirectory->GetObject("Efficiency/MS/tracksprojectedgemMS_any",tracksprojectedgemMS_any);
    
    MSany->Divide(tracksprojectedgemMS_oneormore_cut,tracksprojectedgemMS_any,1,1,"");
    //MSany->SetMaximum(1);
    MSany->GetXaxis()->SetTitle("Hori. Position (mm)");
    MSany->GetYaxis()->SetTitle("Vert. Position (mm)");
    MSany->Draw("colz");
    
    Double_t MSX_n = tracksprojectedgemMS_oneormore_cut->GetEntries();
    Double_t MSX_d = tracksprojectedgemMS_any->GetEntries();
    Double_t MSX_eff=MSX_n/MSX_d;
    
    Double_t MSerr=1/MSX_d*(sqrt(MSX_n*(1-MSX_eff)));
    cout << " GEM MS Efficiency: " << MSX_eff << " +/- " << MSerr << endl;
    
    C1->SaveAs("Efficiency plots_minPeakAdc200_6290.pdf");
    
    /*TCanvas *canv = new TCanvas("c1",Form("Multiplicity Plots"), 3000,1000);
     //gStyle->SetOptStat(kFALSE);
     canv->Divide(3,1);
     canv->cd(1);
     h1->Draw("hist");
     //TH1F *h1 = new TH1F("h1","h1",100,-200,200);
     //test1->Draw("teletracks.tracks.x0 >> h1");
     canv->cd(2);
     h2->Draw("hist");
     canv->cd(3);
     h3->Draw("hist");
     canv->SaveAs("Multiplicity plots.pdf");
     TCanvas *canv1 = new TCanvas("c1",Form("Max_Charge Cluster Distribution"), 3000,1000);
     canv1->Divide(3,1);
     canv1->cd(1);
     h4->Draw("colz");
     canv1->cd(2);
     h5->Draw("colz");
     canv1->cd(3);
     h6->Draw("colz");
     canv1->SaveAs("Max Charge Cluster Distribusion.pdf");*/
    
    
    
    return;
    
}
