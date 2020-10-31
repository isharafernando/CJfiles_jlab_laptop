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

void GetClusterInfo()
{
   
    // This code extracts the cluster information

    TFile *file1 = new TFile(Form("/home/ishara/MUSE/cookedfiles/cookedGEMrun0%d_median.root",6290));
    
    TCanvas *C1 = new TCanvas("C1",Form("ADC values of x"), 3000,2000);
    C1->cd(1);
    TH1 *ADCSpec_gem_x;
    gDirectory->GetObject("US/adcSpec_gem_0_x",ADCSpec_gem_x);
    ADCSpec_gem_x->GetEntries();
    ADCSpec_gem_x->GetXaxis()->SetTitle("APV strip number");
    ADCSpec_gem_x->GetYaxis()->SetTitle("ADC value");
    ADCSpec_gem_x->Draw("colz");
    C1->SaveAs("Median_Method_ClusterInfo_US_x.pdf");


   TCanvas *C2 = new TCanvas("C2",Form("ADC values of y"), 3000,2000);
    C2->cd(2);
    TH1 *ADCSpec_gem_y;
    gDirectory->GetObject("US/adcSpec_gem_0_y",ADCSpec_gem_y);
    ADCSpec_gem_y->GetEntries();
    ADCSpec_gem_y->GetXaxis()->SetTitle("APV strip number");
    ADCSpec_gem_y->GetYaxis()->SetTitle("ADC value");
    ADCSpec_gem_y->Draw("colz");
    C2->SaveAs("Median_Method_ClusterInfo_US_y.pdf");	


   TCanvas *C3 = new TCanvas("C3",Form("ADC values of x after cmode subtraction"), 3000,2000);
    C3->cd(3);
    TH1 *ADCSpec_gem_x_cmode;
    gDirectory->GetObject("US/adcSpec_gem_0_x_cmode",ADCSpec_gem_x_cmode);
    ADCSpec_gem_x_cmode->GetEntries();
    ADCSpec_gem_x_cmode->GetXaxis()->SetTitle("APV strip number");
    ADCSpec_gem_x_cmode->GetYaxis()->SetTitle("ADC value");
    ADCSpec_gem_x_cmode->Draw("colz");
    C3->SaveAs("Median_Method_ClusterInfo_US_x_cmode.pdf");	


   


return;
   
}


/*
 Useful functions from MUSEteleTracker
 
 
 
*/


