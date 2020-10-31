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

void getavgeff_Anusha(int run1, int mode=1) //mode: 1=DAQevent, 2=MIDAS standalone, 3=MIDAS)
{

    //gStyle->SetOptStat("0111");


    // int run1 = 2038;
    // int run2 = 1755;
    // int run3 = 1756;
    //int run4 = 1757;
    //int run5 = 1758;
    //int run6 = 1759;
    //int run7 = 1760;


    if (mode == 1) {
        TFile* file1 = TFile::Open(Form("/home/anusha/work/muse/muse/Analysis/PSI2017/data/ASCII/cookedfiles/GEMtracks0%d.root",run1));
        TTree *MUSEteleTracks = (TTree*)file1->Get("MUSEteleTracks");
        MUSEteleTracks->AddFriend("lumigemcooked",Form("/home/anusha/work/muse/muse/Analysis/PSI2017/data/ASCII/cookedfiles/GEMclusters0%d.root",run1));
    }
    else if (mode == 2) {
        TFile* file1 = TFile::Open(Form("/home/anusha/work/muse/muse/Analysis/PSI2017/data/midas_standalone/cookedfiles/GEMtracks00%d.root",run1));
        TTree *MUSEteleTracks = (TTree*)file1->Get("MUSEteleTracks");
        MUSEteleTracks->AddFriend("lumigemcooked",Form("/home/anusha/work/muse/muse/Analysis/PSI2017/data/midas_standalone/cookedfiles/GEMclusters00%d.root",run1));
    }
    else if (mode == 3) {
        TFile* file1 = TFile::Open(Form("./trackedfiles/GEM_tracked0%d.root",run1));
        TTree *MUSEteleTracks = (TTree*)file1->Get("MUSEteleTracks");
        MUSEteleTracks->AddFriend("lumigemcooked",Form("/cookedfiles/cookedGEMrun0%d.root",run1));
    }
    else {
        printf("Mode not given");
        exit;
    }

    /////JUst any cluster hit on GEMs //////////////////////////////
  
    TCanvas *c1 = new TCanvas("c1","All Cluster Positions ",1800,600);
   /*
    c1->Divide(3,2);
    c1->cd(1);
    TH2F *USall = new TH2F("USall","All Cluster Positions on US",50,-50.0,50.0,50,-50.0,50.0);
    MUSEteleTracks->Draw("(xl*0.4-50):(yl*0.4-50) >> USall","GEMid==3","colz");
    USall->GetXaxis()->SetTitle("Hori. Position (mm)");
    USall->GetYaxis()->SetTitle("Vert. Position (mm)");

    c1->cd(2);
    TH2F *MSall = new TH2F("MSall","All Cluster Positions on MS",50,-50.0,50.0,50,-50.0,50.0);
    MUSEteleTracks->Draw("xl*0.4-50:yl*0.4-50 >> MSall","GEMid==4","colz");
    MSall->GetXaxis()->SetTitle("Hori. Position (mm)");
    MSall->GetYaxis()->SetTitle("Vert. Position (mm)");

    c1->cd(3);
    TH2F *DSall = new TH2F("DSall","All Cluster Positions on DS",50,-50.0,50.0,50,-50.0,50.0);
    MUSEteleTracks->Draw("xl*0.4-50:yl*0.4-50 >> DSall","GEMid==5","colz");
    DSall->GetXaxis()->SetTitle("Hori. Position (mm)");
    DSall->GetYaxis()->SetTitle("Vert. Position (mm)");
    */

    ///////////////// GEM efficiencies (The one cluster charge cluster is selected on the first two GEMs and at least one track on the third GEM within the visinity of cut)//////////////////

    //  TCanvas *c5 = new TCanvas("c5","Efficiencies on the GEMs -  projected (in vicinity cut,cluster =1/any)",1200,400);
    //c5->Divide(3,1);
    //c5->cd(1);
    c1->cd(4);
    TH2F *USany = new TH2F("USany","Efficiencies on US (cluster>=1, cuts)",50,-50.0,50.0,50,-50.0,50.0);
    TH2F *MSany = new TH2F("MSany","Efficiencies on MS (cluster>=1, cuts)",50,-50.0,50.0,50,-50.0,50.0);
    TH2F *DSany = new TH2F("DSany","Efficiencies on DS (cluster>=1, cuts)",50,-50.0,50.0,50,-50.0,50.0);

    TH2 *tracksprojectedgemUS_oneormore_cut; gDirectory->GetObject("TeleTracks/US/tracksprojectedgemUS_oneormore_cut",tracksprojectedgemUS_oneormore_cut);
    TH2 *tracksprojectedgemUS_any; gDirectory->GetObject("TeleTracks/US/tracksprojectedgemUS_any",tracksprojectedgemUS_any);

    USany->Divide(tracksprojectedgemUS_oneormore_cut,tracksprojectedgemUS_any,1,1,"");
    //USany->SetMaximum(1);
    USany->GetXaxis()->SetTitle("Hori. Position (mm)");
    USany->GetYaxis()->SetTitle("Vert. Position (mm)");
    USany->Draw("colz");

    Double_t USX_n = tracksprojectedgemUS_oneormore_cut->GetEntries();
    Double_t USX_d = tracksprojectedgemUS_any->GetEntries();
    Double_t USX_eff=USX_n/USX_d;

    Double_t USerr=1/USX_d*(sqrt(USX_n*(1-USX_eff)));
    cout<<" GEM US Efficiency    : "<<USX_eff<<" +/- "<<USerr<<endl;

    // TLatex *us= new TLatex();
    // us->DrawLatex(10,-40,"#color[1]{Eff= 95.16 }""%""");

    c1->cd(5);
    TH2 *tracksprojectedgemMS_oneormore_cut; gDirectory->GetObject("TeleTracks/MS/tracksprojectedgemMS_oneormore_cut",tracksprojectedgemMS_oneormore_cut);
    TH2 *tracksprojectedgemMS_any; gDirectory->GetObject("TeleTracks/MS/tracksprojectedgemMS_any",tracksprojectedgemMS_any);

    MSany->Divide(tracksprojectedgemMS_oneormore_cut,tracksprojectedgemMS_any,1,1,"");
    //MSany->SetMaximum(1);
    MSany->GetXaxis()->SetTitle("Hori. Position (mm)");
    MSany->GetYaxis()->SetTitle("Vert. Position (mm)");
    MSany->Draw("colz");

    Double_t MSX_n = tracksprojectedgemMS_oneormore_cut->GetEntries();
    Double_t MSX_d = tracksprojectedgemMS_any->GetEntries();
    Double_t MSX_eff=MSX_n/MSX_d;

    Double_t MSerr=1/MSX_d*(sqrt(MSX_n*(1-MSX_eff)));
    cout<<" GEM MS Efficiency    : "<<MSX_eff<<" +/- "<<MSerr<<endl;


    //TLatex *ms= new TLatex();
    // ms->DrawLatex(10,-40,"#color[1]{Eff= 95.03 }""%""");

    c1->cd(6);

    TH2 *tracksprojectedgemDS_oneormore_cut; gDirectory->GetObject("TeleTracks/DS/tracksprojectedgemDS_oneormore_cut",tracksprojectedgemDS_oneormore_cut);
    TH2 *tracksprojectedgemDS_any; gDirectory->GetObject("TeleTracks/DS/tracksprojectedgemDS_any",tracksprojectedgemDS_any);

    DSany->Divide(tracksprojectedgemDS_oneormore_cut,tracksprojectedgemDS_any,1,1,"");
    //DSany->SetMaximum(1);
    DSany->GetXaxis()->SetTitle("Hori. Position (mm)");
    DSany->GetYaxis()->SetTitle("Vert. Position (mm)");
    DSany->Draw("colz");

    Double_t DSX_n = tracksprojectedgemDS_oneormore_cut->GetEntries();
    Double_t DSX_d = tracksprojectedgemDS_any->GetEntries();
    Double_t DSX_eff=DSX_n/DSX_d;

    Double_t DSerr=1/DSX_d*(sqrt(DSX_n*(1-DSX_eff)));
    cout<<" GEM DS Efficiency    : "<<DSX_eff<<" +/- "<<DSerr<<endl;

}
