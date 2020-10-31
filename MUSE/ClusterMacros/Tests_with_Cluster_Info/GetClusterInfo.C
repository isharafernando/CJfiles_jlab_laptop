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

void GetClusterInfo(int run_num)
{
   
    // This code extracts the cluster information

    TFile *file1 = new TFile(Form("/home/ishara/MUSE/cookedfiles/cookedGEMrun0%d.root",run_num));
    
    TCanvas *C1 = new TCanvas("C1",Form("Hitmap US"), 3000,3000);
    C1->cd(1);
    TH1 *hitmap_gem_US;
    gDirectory->GetObject("US/hitmap_gem_0",hitmap_gem_US);
    hitmap_gem_US->GetEntries();
    hitmap_gem_US->GetXaxis()->SetTitle("Horizontal Channels");
    hitmap_gem_US->GetYaxis()->SetTitle("Vertical Channels");
    hitmap_gem_US->Draw("colz");
    C1->SaveAs("Hitmap_US.pdf");


   TCanvas *C2 = new TCanvas("C2",Form("Clustermap US"), 3000,3000);
    C2->cd(2);
    TH1 *clustermap_gem_US;
    gDirectory->GetObject("US/clustermap_gem_0",clustermap_gem_US);
    clustermap_gem_US->GetEntries();
    clustermap_gem_US->GetXaxis()->SetTitle("Horizontal Axis");
    clustermap_gem_US->GetYaxis()->SetTitle("Vertical Axis");
    clustermap_gem_US->Draw("colz");
    C2->SaveAs("Clustermap_US.pdf");	


   TCanvas *C3 = new TCanvas("C3",Form("Cluster multiplicity US"), 3000,2000);
    C3->cd(3);
    TH1 *Cluster_Mult_US;
    gDirectory->GetObject("US/cluster_multiplicity_gem_0",Cluster_Mult_US);
    Cluster_Mult_US->GetEntries();
    Cluster_Mult_US->GetXaxis()->SetTitle("Counts");
    Cluster_Mult_US->GetYaxis()->SetTitle("Multiplicity");
    Cluster_Mult_US->Draw("colz");
    C3->SaveAs("Cluster_Mult_US.pdf");	


    TCanvas *C4 = new TCanvas("C4",Form("Cluster Charge Distribution US"), 3000,2000);
    C4->cd(4);
    TH1 *Cluster_Charge_US;
    gDirectory->GetObject("US/cluster_charge_gem_0",Cluster_Charge_US);
    Cluster_Charge_US->GetEntries();
    Cluster_Charge_US->GetXaxis()->SetTitle("ADC");
    Cluster_Charge_US->GetYaxis()->SetTitle("Counts");
    Cluster_Charge_US->Draw("colz");
    Cluster_Charge_US->Fit("gaus","R","",15500,19000);
    gStyle->SetOptFit(0111);
    C4->SaveAs("2DCluster_Charge_US.pdf");	


    TCanvas *C5 = new TCanvas("C5",Form("Cluster Weight Distribution US"), 3000,2000);
    C5->cd(5);
    TH1 *Cluster_Weight_US;
    gDirectory->GetObject("US/cluster_weight_gem_0",Cluster_Weight_US);
    Cluster_Weight_US->GetEntries();
    Cluster_Weight_US->GetXaxis()->SetTitle("ADC");
    Cluster_Weight_US->GetYaxis()->SetTitle("Counts");
    Cluster_Weight_US->Draw("colz");
    //Cluster_Weight_US->Fit("gaus","R","",15500,19000);
    //gStyle->SetOptFit(0111);
    C5->SaveAs("2DCluster_Weight_US.pdf");	


   


return;
   
}


/*
 Useful functions from MUSEteleTracker
 
 
 
*/


