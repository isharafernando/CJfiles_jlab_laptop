#include <MUSEteleTracker.h>
#include<iostream>
#include<cmath>
#include "lumigemtree.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1.h"
#include "TH2.h"


Long_t MUSEteleTracker::startup_efficiency()
{
    
    
    return Plugin::ok;
    
    
    ///////////////////////////////////////////
}
Long_t MUSEteleTracker::histos_efficiency()
{
    effgemUS=dH2("effgemUS","Efficiency on GEM US",50,-50.0,50.0,50,-50.0,50.0);
    effgem4TH=dH2("effgem4TH","Efficiency on GEM MS",50,-50.0,50.0,50,-50.0,50.0);
    effgemMS=dH2("effgemMS","Efficiency on GEM DS",50,-50.0,50.0,50,-50.0,50.0);

    //projectedX=dH1("Efficiency/test/tracksprojectedgemX_any", "Tracks Projected on GEM X-Any cluster",50,-50.0,50.0,50);
     //projectedX=dH1("Efficiency/test/tracksprojectedgemX_any","Multiplicity US GEM",10,-50.0,50.0);
    
    
    multiplicityUS=dH1("Efficiency/US/multiplicityUS","Multiplicity US GEM",11,-0.5,10.5);
    multiplicity4TH=dH1("Efficiency/4TH/multiplicity4TH","Multiplicity 4TH GEM",11,-0.5,10.5);
    multiplicityMS=dH1("Efficiency/MS/multiplicityMS","Multiplicity MS GEM",11,-0.5,10.5);
    //multiplicityDS=dH1("Efficiency/DS/multiplicityDS","Multiplicity DS GEM",11,-0.5,10.5);
    
    ///2d multiplicity
    //mult2dUS_MS=dH2("Efficiency/mult2dUS_MS","Multiplicity US Vs MS",9,-0.5,8.5,9,-0.5,8.5);
    //mult2dMS_DS=dH2("Efficiency/mult2dMS_DS","Multiplicity MS Vs DS",9,-0.5,8.5,9,-0.5,8.5);
    //mult2dUS_DS=dH2("Efficiency/mult2dUS_DS","Multiplicity US Vs DS",9,-0.5,8.5,9,-0.5,8.5);
    
    //trackmultiplicityUS=dH1("Efficiency/US/trackmultiplicityUS","Track Multiplicity US GEM",11,-0.5,10.5);
    //trackmultiplicityMS=dH1("Efficiency/MS/trackmultiplicityMS","Track Multiplicity MS GEM",11,-0.5,10.5);
    //trackmultiplicityDS=dH1("Efficiency/DS/trackmultiplicityDS","Track Multiplicity DS GEM",11,-0.5,10.5);
    
    //singleGEMclust=dH2("Efficiency/US/singleGEMclust","Single GEM clusters",50,-50,50,50,-50,50);
    
    
    //USX_eff=dH1("USX_eff","O GEM X Efficiency",50,-50.0,50.0);
    //MSX_eff=dH1("MSX_eff","1 GEM X Efficiency",50,-50.0,50.0);
    //DSX_eff=dH1("DSX_eff","2 GEM X Efficiency",50,-50.0,50.0);
    
    //USY_eff=dH1("USY_eff","O GEM Y Efficiency",50,-50.0,50.0);
    //MSY_eff=dH1("MSY_eff","1 GEM Y Efficiency",50,-50.0,50.0);
    //DSY_eff=dH1("DSY_eff","2 GEM Y Efficiency",50,-50.0,50.0);
    
    return 0;
    
}

Long_t MUSEteleTracker::process_efficiency()
{
    // SHOULD USE THESE DISTANCES FOR THE DATA AFTER December, 2016
    const double zgem[4] =
    { 1840.7, 1902.7, 1964.7, 2026.7 };
    
    
    //const double zgempos[8]= { 1846.2, 1946.2, 2046.2, 2146.2, 1840.7, 1902.7, 1964.7, 2026.7};
    //const double zgempos[4]= { -536.0, -474.0 , -412.0, -350.0};
    const double zgempos[3]= {-350.0,-412.0,-474.0};
    int zgemd=0;
    
    const char* gemdet;
    
    // loop over both telescopes:
    // loop over all possible combinations of clusters:
    double xe1[GEM_NUM]={-10000}, ye1[GEM_NUM]={-10000}, dx1={-10000}, dy1={-10000};
    double xe[GEM_NUM]={-10000}, ye[GEM_NUM]={-10000}, xe_max[GEM_NUM]={-10000}, ye_max[GEM_NUM]={-10000};
    double dx_max[GEM_NUM]={-10000},dy_max[GEM_NUM]={-10000}, dxe={-10000}, dye={-10000};
    double x[GEM_NUM]={-10000.0},y[GEM_NUM]={-10000.0},dx={-10000.0},dy={-10000.0};
    double xmaxchg[GEM_NUM]={-10000};
    double ymaxchg[GEM_NUM]={-10000};
    bool anyhit=true;
    bool nohit=true;
    int gemd;
    
    double cluster_gem_charge[20];
    int gemanyhit[3] = { 0, 0, 0};
    int tracks[3] = { 0, 0, 0 };
    
    //This For loop collect the number of clusters on each GEM and x,y coordinates of max_charge cluster
    for (int gems=0; gems<GEM_NUM; gems++)
    {
        double maxcharge=0;
        bool test=false;
        for (unsigned int g1=0; g1<clusters->hits.size(); g1++)
        {
            if (clusters->hits[g1].GEMid!=gems) continue;
            gemanyhit[gems]++;
            cluster_gem_charge[gemanyhit[gems]]=clusters->hits[g1].charge;
            if(cluster_gem_charge[gemanyhit[gems]]>maxcharge) {
                maxcharge =cluster_gem_charge[gemanyhit[gems]];
                test=true;
                xmaxchg[1] = clusters->hits[g1].xl*0.4-50.;
                ymaxchg[1] = clusters->hits[g1].yl*0.4-50.;
            }
            else test=false;
        };
        if (gems==0)  multiplicityUS->Fill(gemanyhit[gems]);
        if (gems==1)  multiplicity4TH->Fill(gemanyhit[gems]);
        if (gems==2)  multiplicityMS->Fill(gemanyhit[gems]);

        // mult2dUS_MS->Fill(gemanyhit[3],gemanyhit[4]);//test
        // mult2dUS_DS->Fill(gemanyhit[3],gemanyhit[5]);
        // mult2dMS_DS->Fill(gemanyhit[4],gemanyhit[5]);
        
        if ((gems==0)&&test) H2(xmaxchg[1],ymaxchg[1], Form("Efficiency/US/max_charge_clusters_USGEM"), Form("Maximum Charge Cluster positions on US GEM"),50,-50.0,50.0,50,-50.0,50.0);
        if ((gems==1)&&test) H2(xmaxchg[1],ymaxchg[1], Form("Efficiency/4TH/max_charge_clusters_4THGEM"), Form("Maximum Charge Cluster positions on 4TH GEM"),50,-50.0,50.0,50,-50.0,50.0);
        if ((gems==2)&&test) H2(xmaxchg[1],ymaxchg[1], Form("Efficiency/MS/max_charge_clusters_MSGEM"), Form("Maximum Charge Cluster positions on MS GEM"),50,-50.0,50.0,50,-50.0,50.0);
        //if ((gems==9)&&test) H2(xmaxchg[1],ymaxchg[1], Form("Efficiency/DS/max_charge_clusters_DSGEM"), Form("Maximum Charge Cluster positions on DS GEM"),50,-50.0,50.0,50,-50.0,50.0);
    };
    
    // loop over all the GEMs and get the efficiency related histograms/////////////////////////////
    // Assume GEM 1,2,3,4 and each time it goto one GEM, say GEM 1, then loop over all the other GEMs, GEM, 2,3,4. If Gem 2 does not have any clusters found, the loop continues to search clusters on the GEM 3. If the GEM 3 has a cluster, then the loop continue to search any clusters on the rest GEM 4.
    
    int NG=3;
    int gem0hits[3]={0,0,0};
    int gem1hits[3]={0,0,0};
    int gem2hits[3]={0,0,0};
    
    for (int gem0=0; gem0<NG; gem0++)  // loop 1 over the first GEM
    {
        
        //Here we take the number of clusters on GEM0
        for (unsigned int g00=0; g00<clusters->hits.size(); g00++)
        {
            if ((clusters->hits[g00].GEMid!=gem0))  continue;
            gem0hits[gem0]++;  // total # of hits on the first gem
        };
        
        //Here we switch to the next GEM if we don't get clusters
        if (gem0hits[gem0]==0) continue;
        
        //Here is the case if we have clusters >= 1
        double cluster_charge_gem0[gem0hits[gem0]];
        
        if (gem0hits[gem0]>=1)
        {
            bool good_cluster_gem0=false;
            double max_charge_gem0=0;
            int gem0_cluster_number=-1;
            for (unsigned int g0=0; g0<clusters->hits.size(); g0++)
            {
               if ((clusters->hits[g0].GEMid!=gem0)) continue;
               cluster_charge_gem0[g0]=clusters->hits[g0].charge;
                if(cluster_charge_gem0[g0]>max_charge_gem0)
                {
                    good_cluster_gem0=true;
                    max_charge_gem0=cluster_charge_gem0[g0];
                    gem0_cluster_number=g0;
                };
            };
            
            if(good_cluster_gem0)
            {
              xe[1] = clusters->hits[gem0_cluster_number].xl*0.4-50.;
              ye[1] = clusters->hits[gem0_cluster_number].yl*0.4-50.;
                
              
              // Then start the next loop over next 2 GEMs only if there is a good cluster on first GEM
              
              for (int gem1=1; gem1<NG; gem1++)
              {
                  for (unsigned int g10=0; g10<clusters->hits.size(); g10++)
                  {
                      if ((clusters->hits[g10].GEMid!=gem1))  continue;
                      gem1hits[gem1]++;  // total # of hits on the second GEM
                  };
                  
                  if (gem1hits[gem1]==0) continue;
                  
                  double cluster_charge_gem1[gem1hits[gem1]];
                  if (gem1hits[gem1]>=1)
                  {
                      bool good_cluster_gem1=false;
                      double max_charge_gem1=0;
                      int gem1_cluster_number=-1;
                      for (unsigned int g1=0; g1<clusters->hits.size(); g1++)
                      {
                          if ((clusters->hits[g1].GEMid!=gem1)) continue;
                          cluster_charge_gem1[g1]=clusters->hits[g1].charge;
                          if(cluster_charge_gem1[g1]>max_charge_gem1)
                          {
                              good_cluster_gem1=true;
                              max_charge_gem1=cluster_charge_gem1[g1];
                              gem1_cluster_number=g1;
                          };
                      };
                      
                      if(good_cluster_gem0)
                      {
                          xe[2] = clusters->hits[gem1_cluster_number].xl*0.4-50.;
                          ye[2] = clusters->hits[gem1_cluster_number].yl*0.4-50.;
                          
                          double slopexe=(xe[2]-xe[1])/(zgempos[gem1]-zgempos[gem0]);
                          double slopeye=(ye[2]-ye[1])/(zgempos[gem1]-zgempos[gem0]);
                          
                          if ((gem0==0)&&(gem1==1))
                          {
                              //std::cout << "MS" << std::endl;
                              zgemd=2;
                              gemdet="MS";
                          };
                          
                          if ((gem0==0)&&(gem1==2))
                          {
                              //std::cout << "4TH" << std::endl;
                              zgemd=1;
                              gemdet="4TH";
                          };
                          
                          if ((gem0==1)&&(gem1==2))
                          {
                              //std::cout << "US" << std::endl;
                              zgemd=2;
                              gemdet="US";
                          };
                          
                          //Project the tracks on GEM2 (thrid gem)
                          double X3 = xe[1] + slopexe*(zgempos[zgemd]-zgempos[gem0]);
                          double Y3 = ye[1] + slopeye*(zgempos[zgemd]-zgempos[gem0]);
                          
			  //printf("Projected coordinates are = (%d,%d) \n",X3,Y3);
                          
                          H2(X3,Y3, Form("Efficiency/%s/tracksprojectedgem%s_any", gemdet,gemdet), Form("Tracks Projected on %s GEM-Any cluster", gemdet),50,-50.0,50.0,50,-50.0,50.0);
                          
                          //projectedX=dH1("Efficiency/test/tracksprojectedgemX_any","Multiplicity US GEM");
                          //projectedX->Fill(X3);
                          
                          H1(Y3, Form("Efficiency/%s/tracksprojectedgem%sY_any", gemdet, gemdet), Form("Tracks Projected on %s GEM Y-Any cluster", gemdet),50,-50.0,50.0);
                        
                          
                          

                      } // If GEM1 has a good cluster with max charge
                  }; //If GEM1 has clusters >=
                  
              }; // Loop for GEM1
                
            }; // If GEM0 has a good cluster with max charge
        }; // If GEM0 has clusters >=1
            
            
    };//Loop for GEM0
        
    
    
    
    
    return 0;
    //  return Plugin::ok;
};

Long_t MUSEteleTracker::finalize_efficiency()
{
    
    //Double USX_n = tracksprojectedgemUSX_onecut->GetEntries();
    //Double USX_d = tracksprojectedgemUSX_any->GetEntries();
    //Double USX_eff=USX_n/USX_d;
    
    //printf("US X efficiency %d %3.2lf\n",USX_n,USX_eff);
    
    //USX_eff->Divide(tracksprojectedgemUSX_onecut,tracksprojectedgemUSX_any,1,1,"");
    //MSX_eff->Divide(tracksprojectedgemMSX_onecut,tracksprojectedgemMSX_any,1,1,"");
    //DSX_eff->Divide(tracksprojectedgemDSX_onecut,tracksprojectedgemDSX_any,1,1,"");
    
    //USY_eff->Divide(tracksprojectedgemUSY_onecut,tracksprojectedgemUSY_any,1,1,"");
    //MSY_eff->Divide(tracksprojectedgemMSY_onecut,tracksprojectedgemMSY_any,1,1,"");
    //DSY_eff->Divide(tracksprojectedgemDSY_onecut,tracksprojectedgemDSY_any,1,1,"");
    
    return Plugin::ok;
}

