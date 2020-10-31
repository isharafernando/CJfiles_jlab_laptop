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
    C1->SaveAs("ClusterInfo_US_x.pdf");


   TCanvas *C2 = new TCanvas("C2",Form("ADC values of y"), 3000,2000);
    C2->cd(2);
    TH1 *ADCSpec_gem_y;
    gDirectory->GetObject("US/adcSpec_gem_0_y",ADCSpec_gem_y);
    ADCSpec_gem_y->GetEntries();
    ADCSpec_gem_y->GetXaxis()->SetTitle("APV strip number");
    ADCSpec_gem_y->GetYaxis()->SetTitle("ADC value");
    ADCSpec_gem_y->Draw("colz");
    C2->SaveAs("ClusterInfo_US_y.pdf");	


   TCanvas *C3 = new TCanvas("C3",Form("ADC values of x after cmode subtraction"), 3000,2000);
    C3->cd(3);
    TH1 *ADCSpec_gem_x_cmode;
    gDirectory->GetObject("US/adcSpec_gem_0_x_cmode",ADCSpec_gem_x_cmode);
    ADCSpec_gem_x_cmode->GetEntries();
    ADCSpec_gem_x_cmode->GetXaxis()->SetTitle("APV strip number");
    ADCSpec_gem_x_cmode->GetYaxis()->SetTitle("ADC value");
    ADCSpec_gem_x_cmode->Draw("colz");
    C3->SaveAs("ClusterInfo_US_x_cmode.pdf");	


   TCanvas *C4 = new TCanvas("C4",Form("ADC values of y after cmode subtraction"), 3000,2000);
    C4->cd(4);
    TH1 *ADCSpec_gem_y_cmode;
    gDirectory->GetObject("US/adcSpec_gem_0_y_cmode",ADCSpec_gem_y_cmode);
    ADCSpec_gem_y_cmode->GetEntries();
    ADCSpec_gem_y_cmode->GetXaxis()->SetTitle("APV strip number");
    ADCSpec_gem_y_cmode->GetYaxis()->SetTitle("ADC value");
    ADCSpec_gem_y_cmode->Draw("colz");
    C4->SaveAs("ClusterInfo_US_y_cmode.pdf");	


   TCanvas *C5 = new TCanvas("C5",Form("US APV 0 Spectrum"), 3000,2000);
    C5->cd(5);
    TH1 *APV0;
    gDirectory->GetObject("US/adcSpec_gem_0_apv_0",APV0);
    APV0->GetEntries();
    APV0->GetXaxis()->SetTitle("APV strip number");
    APV0->GetYaxis()->SetTitle("ADC value");
    APV0->Draw("colz");
    C5->SaveAs("APV0_US_x.pdf");	


   TCanvas *C6 = new TCanvas("C6",Form("US APV 1 Spectrum"), 3000,2000);
    C6->cd(6);
    TH1 *APV1;
    gDirectory->GetObject("US/adcSpec_gem_0_apv_1",APV1);
    APV1->GetEntries();
    APV1->GetXaxis()->SetTitle("APV strip number");
    APV1->GetYaxis()->SetTitle("ADC value");
    APV1->Draw("colz");
    C6->SaveAs("APV1_US_x.pdf");	

   TCanvas *C7 = new TCanvas("C7",Form("US APV 2 Spectrum"), 3000,2000);
    C7->cd(7);
    TH1 *APV2;
    gDirectory->GetObject("US/adcSpec_gem_0_apv_2",APV2);
    APV2->GetEntries();
    APV2->GetXaxis()->SetTitle("APV strip number");
    APV2->GetYaxis()->SetTitle("ADC value");
    APV2->Draw("colz");
    C7->SaveAs("APV2_US_x.pdf");


   TCanvas *C8 = new TCanvas("C8",Form("US APV 3 Spectrum"), 3000,2000);
    C8->cd(8);
    TH1 *APV3;
    gDirectory->GetObject("US/adcSpec_gem_0_apv_3",APV3);
    APV3->GetEntries();
    APV3->GetXaxis()->SetTitle("APV strip number");
    APV3->GetYaxis()->SetTitle("ADC value");
    APV3->Draw("colz");
    C8->SaveAs("APV3_US_x.pdf");	


   TCanvas *C9 = new TCanvas("C9",Form("ADC values of x after pedestal subtraction"), 3000,2000);
    C9->cd(9);
    TH1 *ADCSpec_gem_x_pedestal;
    gDirectory->GetObject("US/adcSpec_gem_0_x_pedestal",ADCSpec_gem_x_pedestal);
    ADCSpec_gem_x_pedestal->GetEntries();
    ADCSpec_gem_x_pedestal->GetXaxis()->SetTitle("APV strip number");
    ADCSpec_gem_x_pedestal->GetYaxis()->SetTitle("ADC value");
    ADCSpec_gem_x_pedestal->Draw("colz");
    C9->SaveAs("ClusterInfo_US_x_pedestal.pdf");	


   TCanvas *C10 = new TCanvas("C10",Form("ADC values of y after pedestal subtraction"), 3000,2000);
    C10->cd(10);
    TH1 *ADCSpec_gem_y_pedestal;
    gDirectory->GetObject("US/adcSpec_gem_0_y_pedestal",ADCSpec_gem_y_pedestal);
    ADCSpec_gem_y_pedestal->GetEntries();
    ADCSpec_gem_y_pedestal->GetXaxis()->SetTitle("APV strip number");
    ADCSpec_gem_y_pedestal->GetYaxis()->SetTitle("ADC value");
    ADCSpec_gem_y_pedestal->Draw("colz");
    C10->SaveAs("ClusterInfo_US_y_pedestal.pdf");		


    TCanvas *C11 = new TCanvas("C11",Form("ADC values of x after ped. & cmode"), 3000,2000);
    C11->cd(11);
    TH1 *ADCSpec_gem_x_ped_cmode;
    gDirectory->GetObject("US/adcSpec_gem_0_x_ped_cmode",ADCSpec_gem_x_ped_cmode);
    ADCSpec_gem_x_ped_cmode->GetEntries();
    ADCSpec_gem_x_ped_cmode->GetXaxis()->SetTitle("APV strip number");
    ADCSpec_gem_x_ped_cmode->GetYaxis()->SetTitle("ADC value");
    ADCSpec_gem_x_ped_cmode->Draw("colz");
    C11->SaveAs("US_x_ped_cmode.pdf");



    TCanvas *C12 = new TCanvas("C12",Form("ADC values of y after ped. & cmode"), 3000,2000);
    C12->cd(12);
    TH1 *ADCSpec_gem_y_ped_cmode;
    gDirectory->GetObject("US/adcSpec_gem_0_y_ped_cmode",ADCSpec_gem_y_ped_cmode);
    ADCSpec_gem_y_ped_cmode->GetEntries();
    ADCSpec_gem_y_ped_cmode->GetXaxis()->SetTitle("APV strip number");
    ADCSpec_gem_y_ped_cmode->GetYaxis()->SetTitle("ADC value");
    ADCSpec_gem_y_ped_cmode->Draw("colz");
    C12->SaveAs("US_y_ped_cmode.pdf");

return;
   
}


/*
 Useful functions from MUSEteleTracker
 
 
 
*/


