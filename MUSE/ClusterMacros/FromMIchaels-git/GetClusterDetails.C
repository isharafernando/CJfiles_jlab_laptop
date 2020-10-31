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

void GetClusterDetails(int run_num)
{
   
    TFile *file1 = new TFile(Form("/home/ishara/MUSE/cookedfiles/cookedGEMrun0%d.root",run_num));
    
    TCanvas *C1 = new TCanvas("C1",Form("ADC values of x"), 3000,2000);
    C1->cd(1);
    TH1 *ADCSpec_gem_x;
    gDirectory->GetObject("US/adcSpec_gem_0_x",ADCSpec_gem_x);
    ADCSpec_gem_x->GetEntries();
    ADCSpec_gem_x->GetXaxis()->SetTitle("APV strip number");
    ADCSpec_gem_x->GetYaxis()->SetTitle("ADC value");
    ADCSpec_gem_x->Draw("colz");
    C1->SaveAs("Michael_v_ClusterInfo_US_x.pdf");


   TCanvas *C2 = new TCanvas("C2",Form("ADC values of y"), 3000,2000);
    C2->cd(2);
    TH1 *ADCSpec_gem_y;
    gDirectory->GetObject("US/adcSpec_gem_0_y",ADCSpec_gem_y);
    ADCSpec_gem_y->GetEntries();
    ADCSpec_gem_y->GetXaxis()->SetTitle("APV strip number");
    ADCSpec_gem_y->GetYaxis()->SetTitle("ADC value");
    ADCSpec_gem_y->Draw("colz");
    C2->SaveAs("Michael_v_ClusterInfo_US_y.pdf");	


   TCanvas *C3 = new TCanvas("C3",Form("ADC values of x after cmode subtraction"), 3000,2000);
    C3->cd(3);
    TH1 *ADCSpec_gem_x_cmode;
    gDirectory->GetObject("US/adcSpec_gem_0_x_cmode",ADCSpec_gem_x_cmode);
    ADCSpec_gem_x_cmode->GetEntries();
    ADCSpec_gem_x_cmode->GetXaxis()->SetTitle("APV strip number");
    ADCSpec_gem_x_cmode->GetYaxis()->SetTitle("ADC value");
    ADCSpec_gem_x_cmode->Draw("colz");
    C3->SaveAs("Michael_v_ClusterInfo_US_x_cmode.pdf");	


   TCanvas *C4 = new TCanvas("C4",Form("ADC values of y after cmode subtraction"), 3000,2000);
    C4->cd(4);
    TH1 *ADCSpec_gem_y_cmode;
    gDirectory->GetObject("US/adcSpec_gem_0_y_cmode",ADCSpec_gem_y_cmode);
    ADCSpec_gem_y_cmode->GetEntries();
    ADCSpec_gem_y_cmode->GetXaxis()->SetTitle("APV strip number");
    ADCSpec_gem_y_cmode->GetYaxis()->SetTitle("ADC value");
    ADCSpec_gem_y_cmode->Draw("colz");
    C4->SaveAs("Michael_v_ClusterInfo_US_y_cmode.pdf");	


   TCanvas *C5 = new TCanvas("C5",Form("ADC values of x after pedestal subtraction"), 3000,2000);
    C5->cd(5);
    TH1 *ADCSpec_gem_x_pedestal;
    gDirectory->GetObject("US/adcSpec_gem_0_x_ped",ADCSpec_gem_x_pedestal);
    ADCSpec_gem_x_pedestal->GetEntries();
    ADCSpec_gem_x_pedestal->GetXaxis()->SetTitle("APV strip number");
    ADCSpec_gem_x_pedestal->GetYaxis()->SetTitle("ADC value");
    ADCSpec_gem_x_pedestal->Draw("colz");
    C5->SaveAs("Michael_v_ClusterInfo_US_x_pedestal.pdf");	


   TCanvas *C6 = new TCanvas("C6",Form("ADC values of y after pedestal subtraction"), 3000,2000);
    C6->cd(6);
    TH1 *ADCSpec_gem_y_pedestal;
    gDirectory->GetObject("US/adcSpec_gem_0_y_ped",ADCSpec_gem_y_pedestal);
    ADCSpec_gem_y_pedestal->GetEntries();
    ADCSpec_gem_y_pedestal->GetXaxis()->SetTitle("APV strip number");
    ADCSpec_gem_y_pedestal->GetYaxis()->SetTitle("ADC value");
    ADCSpec_gem_y_pedestal->Draw("colz");
    C6->SaveAs("Michael_v_ClusterInfo_US_y_pedestal.pdf");		


    TCanvas *C7 = new TCanvas("C7",Form("ADC values of x after ped. & cmode"), 3000,2000);
    C7->cd(7);
    TH1 *ADCSpec_gem_x_ped_cmode;
    gDirectory->GetObject("US/adcSpec_gem_0_x_ped_cmode",ADCSpec_gem_x_ped_cmode);
    ADCSpec_gem_x_ped_cmode->GetEntries();
    ADCSpec_gem_x_ped_cmode->GetXaxis()->SetTitle("APV strip number");
    ADCSpec_gem_x_ped_cmode->GetYaxis()->SetTitle("ADC value");
    ADCSpec_gem_x_ped_cmode->Draw("colz");
    C7->SaveAs("Michael_v_US_x_ped_cmode.pdf");



    TCanvas *C8 = new TCanvas("C8",Form("ADC values of y after ped. & cmode"), 3000,2000);
    C8->cd(8);
    TH1 *ADCSpec_gem_y_ped_cmode;
    gDirectory->GetObject("US/adcSpec_gem_0_y_ped_cmode",ADCSpec_gem_y_ped_cmode);
    ADCSpec_gem_y_ped_cmode->GetEntries();
    ADCSpec_gem_y_ped_cmode->GetXaxis()->SetTitle("APV strip number");
    ADCSpec_gem_y_ped_cmode->GetYaxis()->SetTitle("ADC value");
    ADCSpec_gem_y_ped_cmode->Draw("colz");
    C8->SaveAs("Michael_v_US_y_ped_cmode.pdf");


    TCanvas *C9 = new TCanvas("C9",Form("Hit Map"), 3000,2000);
    C9->cd(9);
    TH1 *HITMAP;
    gDirectory->GetObject("US/hitmap_gem_0",HITMAP);
    HITMAP->GetEntries();
    HITMAP->GetXaxis()->SetTitle("Horizontal");
    HITMAP->GetYaxis()->SetTitle("Vertical");
    HITMAP->Draw("colz");
    C9->SaveAs("Michael_v_US_hitmap.pdf");


    TCanvas *C10 = new TCanvas("C10",Form("ADC wieghted Ped Cmode X"), 3000,2000);
    C10->cd(10);
    TH1 *ADCwieghtedPedCmodeX;
    gDirectory->GetObject("US/adcweightSpec_gem_0_x_ped_cmode",ADCwieghtedPedCmodeX);
    ADCwieghtedPedCmodeX->GetEntries();
    ADCwieghtedPedCmodeX->GetXaxis()->SetTitle("APV strip number");
    ADCwieghtedPedCmodeX->GetYaxis()->SetTitle("ADC weight");
    ADCwieghtedPedCmodeX->Draw("colz");
    C10->SaveAs("ADCwieghtedPedCmodeX_US.pdf");


    TCanvas *C11 = new TCanvas("C11",Form("ADC wieghted Ped Cmode Y"), 3000,2000);
    C11->cd(11);
    TH1 *ADCwieghtedPedCmodeY;
    gDirectory->GetObject("US/adcweightSpec_gem_0_y_ped_cmode",ADCwieghtedPedCmodeY);
    ADCwieghtedPedCmodeY->GetEntries();
    ADCwieghtedPedCmodeY->GetXaxis()->SetTitle("APV strip number");
    ADCwieghtedPedCmodeY->GetYaxis()->SetTitle("ADC weight");
    ADCwieghtedPedCmodeY->Draw("colz");
    C11->SaveAs("ADCwieghtedPedCmodeY_US.pdf");

return;
   
}


/*
 Useful functions from MUSEteleTracker
 
 
 
*/


