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
    
    //multiplicityUS=dH1("Efficiency/US/multiplicityUS","Multiplicity US GEM",11,-0.5,10.5);
    //multiplicity4TH=dH1("Efficiency/4TH/multiplicity4TH","Multiplicity 4TH GEM",11,-0.5,10.5);
    //multiplicityMS=dH1("Efficiency/MS/multiplicityMS","Multiplicity MS GEM",11,-0.5,10.5);
    //multiplicityDS=dH1("Efficiency/DS/multiplicityDS","Multiplicity DS GEM",11,-0.5,10.5);
    
    
    return 0;
    
}

Long_t MUSEteleTracker::process_efficiency()
{
    // Distances from PSI
    //const double zgempos[3]= {-350.0,-412.0,-474.0};
    const double zgempos[3]= {-536.0,-474.0,-412.0};
    int zgemd=0;
    
    const char* gemdet;
    
    // loop over all possible combinations of clusters:
    double xe1[GEM_NUM]={-10000}, ye1[GEM_NUM]={-10000}, dx1={-10000}, dy1={-10000};
    double xe[GEM_NUM]={-10000}, ye[GEM_NUM]={-10000}, xe_max[GEM_NUM]={-10000}, ye_max[GEM_NUM]={-10000};
    double dx_max[GEM_NUM]={-10000},dy_max[GEM_NUM]={-10000}, dxe={-10000}, dye={-10000};
    double x[GEM_NUM]={-10000.0},y[GEM_NUM]={-10000.0},dx={-10000.0},dy={-10000.0};
    double xmaxchg[GEM_NUM]={-10000};
    double ymaxchg[GEM_NUM]={-10000};
    bool anyhit=false;
    bool nohit=false;
    int gemd;
    

    
    double cluster_gem_charge[20];
    int gemanyhit[3] = { 0, 0, 0};
    int tracks[3] = { 0, 0, 0 };
    
    //This For loop collect the number of clusters on each GEM and x,y coordinates of max_charge cluster
    for (int gems=0; gems<GEM_NUM; gems++)
    {
        double maxcharge=0;
        bool test=false;
        //cluster_gem_charge[gemanyhit[gems]]=0;
        for (unsigned int g1=0; g1<clusters->hits.size(); g1++)
        {
            if (clusters->hits[g1].GEMid!=gems) continue;
            gemanyhit[gems]++;
            cluster_gem_charge[gemanyhit[gems]]=clusters->hits[g1].quality;
	    //printf("cluster charge is %d\n",cluster_gem_charge[gemanyhit[gems]]);
            if(cluster_gem_charge[gemanyhit[gems]]>maxcharge) {
	        //printf("found max cluster charge\n");
                maxcharge =cluster_gem_charge[gemanyhit[gems]];
                test=true;
                xmaxchg[1] = clusters->hits[g1].xl*0.4-50.;
                ymaxchg[1] = clusters->hits[g1].yl*0.4-50.;
            }
            else test=false;
        };
        //if (gems==0)  multiplicityUS->Fill(gemanyhit[gems]);
        //if (gems==1)  multiplicity4TH->Fill(gemanyhit[gems]);
        //if (gems==2)  multiplicityMS->Fill(gemanyhit[gems]);

        
        // These plots show a cluster distribution among the GEMs for max charge (quaity)
        if ((gems==0)&&test) H2(xmaxchg[1],ymaxchg[1], Form("Efficiency/US/max_charge_clusters_USGEM"), Form("Maximum Charge Cluster positions on US GEM"),50,-50.0,50.0,50,-50.0,50.0);
        if ((gems==1)&&test) H2(xmaxchg[1],ymaxchg[1], Form("Efficiency/4TH/max_charge_clusters_4THGEM"), Form("Maximum Charge Cluster positions on 4TH GEM"),50,-50.0,50.0,50,-50.0,50.0);
        if ((gems==2)&&test) H2(xmaxchg[1],ymaxchg[1], Form("Efficiency/MS/max_charge_clusters_MSGEM"), Form("Maximum Charge Cluster positions on MS GEM"),50,-50.0,50.0,50,-50.0,50.0);


        
    };
    
    // This print statement will let you to see how many clusters recorded in each GEM
    //printf("clusters US=%d ; 4th=%d ; MS=%d; DS=%d ; \n",gemanyhit[0], gemanyhit[1], gemanyhit[2], gemanyhit[3]);
    
    // loop over all the GEMs and get the efficiency related histograms/////////////////////////////
    
    int NG=3;  // GEM_NUM or Number of GEMs
    // The following arrays are for each search loop (for-loop) over the GEMs
    int gem0hits[3]={0,0,0};
    int gem1hits[3]={0,0,0};
    int gem2hits[3]={0,0,0};



int gemmulti[3] = { 0, 0, 0 };
    for (unsigned int g1=0; g1<clusters->hits.size(); g1++)
    {
        // If the GEMid is different than the cluster's gem ID, then g1 gets incremented by 1 unit
        if (clusters->hits[g1].GEMid!=0) continue;
        gemmulti[0]++;
        for (unsigned int g2=0; g2<clusters->hits.size(); g2++)
        {
            if (clusters->hits[g2].GEMid!=1) continue;
            gemmulti[1]++;
            for (unsigned int g3=0; g3<clusters->hits.size(); g3++)
            {
                if (clusters->hits[g3].GEMid!=2) continue;
                gemmulti[2]++;
                
            };
        };
    };
    

 H1(gemmulti[0], Form("Efficiency/US/Number of possible clusters US"), Form("Number of possible clusters - US GEM"),
       11,-0.5,10.5);
    H1(gemmulti[1], Form("Efficiency/4TH/Number of possible clusters 4TH"), Form("Number of possible clusters - 4TH GEM"),
       11,-0.5,10.5);
    H1(gemmulti[2], Form("Efficiency/MS/Number of possible clusters MI"), Form("Number of possible clusters - MI GEM"),
       11,-0.5,10.5);
    
    
    
    for (int gem0=0; gem0<NG; gem0++)  // loop 1 over the first GEM
    {
        
        //Here we take the number of clusters on GEM0
        for (unsigned int g00=0; g00<clusters->hits.size(); g00++)
        {
            if ((clusters->hits[g00].GEMid!=gem0))  continue;
            gem0hits[gem0]++;  // total # of hits on the first gem
            hits0=hits0+1;
        };
        
        //Here we switch to the next GEM if we don't get clusters
        if (gem0hits[gem0]==0) continue;
        
        //Here is the case if we have clusters >= 1
        int cluster_charge_gem0[gem0hits[gem0]];
	int test_charge;
        //printf("checking for good cluster on GEM0\n");
        if (gem0hits[gem0]>=1)
        {
	    //printf("found one or more clusters on GEM0\n");
            bool good_cluster_gem0=false;
            double max_charge_gem0=0;
            int gem0_cluster_number=-1;
            for (unsigned int g0=0; g0<clusters->hits.size(); g0++)
            {
               if ((clusters->hits[g0].GEMid!=gem0)) continue;
               cluster_charge_gem0[g0]=0;
               cluster_charge_gem0[g0]=clusters->hits[g0].quality;
        
                if(cluster_charge_gem0[g0]>max_charge_gem0)
                {
                    good_cluster_gem0=true;
                    max_charge_gem0=cluster_charge_gem0[g0];
                    gem0_cluster_number=g0;
		         //printf("found the max charge cluseter on GEM0 \n");
                };
            };
            
            if(good_cluster_gem0)
            {
              //printf("good cluster found on GEM0\n");
              xe[1] = clusters->hits[gem0_cluster_number].xl*0.4-50.;
              ye[1] = clusters->hits[gem0_cluster_number].yl*0.4-50.;
              //printf("GEM0 cluster position is  (%5.2f,%5.2f)\n",xe[1],ye[1]);
                
              
              // Then start the next loop over next 2 GEMs only if there is a good cluster on first GEM
              
              for (int gem1=gem0+1; gem1<NG; gem1++)
              {
                  for (unsigned int g10=0; g10<clusters->hits.size(); g10++)
                  {
                      if ((clusters->hits[g10].GEMid!=gem1))  continue;
                      gem1hits[gem1]++;  // total # of hits on the second GEM
 		      
                  };
                  
                  if (gem1hits[gem1]==0) continue;
                  
                  int cluster_charge_gem1[gem1hits[gem1]];
                  if (gem1hits[gem1]>=1)
                  {
                      bool good_cluster_gem1=false;
                      double max_charge_gem1=0;
                      int gem1_cluster_number=-1;
                      for (unsigned int g1=0; g1<clusters->hits.size(); g1++)
                      {
                          if ((clusters->hits[g1].GEMid!=gem1)) continue;
                          cluster_charge_gem1[g1]=clusters->hits[g1].quality;
                          if(cluster_charge_gem1[g1]>max_charge_gem1)
                          {
                              good_cluster_gem1=true;
                              max_charge_gem1=cluster_charge_gem1[g1];
                              gem1_cluster_number=g1;
                          };
                      };
                      
                      if(good_cluster_gem1)
                      {
                          //printf("good cluster found on GEM1\n");
                          
                          xe[2] = clusters->hits[gem1_cluster_number].xl*0.4-50.;
                          ye[2] = clusters->hits[gem1_cluster_number].yl*0.4-50.;
                          //printf("GEM1 cluster position is  (%5.2f,%5.2f)\n",xe[2],ye[2]);
                          
                          double slopexe=(xe[2]-xe[1])/(zgempos[gem1]-zgempos[gem0]);
                          double slopeye=(ye[2]-ye[1])/(zgempos[gem1]-zgempos[gem0]);
                          
                          // Here we collect the number of events which can create tracks
                          // Also we can identify the 3rd GEM


                          if ((gem0==0)&&(gem1==1))
                          {
                              //std::cout << "MS" << std::endl;
                              //printf("MS GEM found\n");
                              zgemd=2;
                              gemdet="MS";
                              MS_count=MS_count+1;
                          };
                          
                          if ((gem0==0)&&(gem1==2))
                          {
                              //printf("4TH GEM found\n");
                              //std::cout << "4TH" << std::endl;
                              zgemd=1;
                              gemdet="4TH";
                              Fourth_count=Fourth_count+1;
                          };
                          
                          if ((gem0==1)&&(gem1==2))
                          {
                              //printf("US GEM found\n");
                              //std::cout << "US" << std::endl;
                              zgemd=0;
                              gemdet="US";
                              US_count=US_count+1;
                          };
                          
                          //Project the tracks on GEM2 (third gem)
                          double X3 = xe[2] + slopexe*(zgempos[zgemd]-zgempos[gem0]);
                          double Y3 = ye[2] + slopeye*(zgempos[zgemd]-zgempos[gem0]);
                          
                          //printf("x,y slopes are = (%5.2f,%5.2f) \n",slopexe,slopeye);
                          
                          H2(X3,Y3, Form("Efficiency/%s/tracksprojectedgem_%s", gemdet,gemdet), Form("Distribution of Tracks Projected on %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);
                          
                           H1(X3, Form("Efficiency/%s/tracksprojectedgem_%s_X", gemdet, gemdet), Form("X-Distribution of Tracks Projected on %s GEM ", gemdet),50,-50.0,50.0);
                          
                          H1(Y3, Form("Efficiency/%s/tracksprojectedgem_%s_Y", gemdet, gemdet), Form("Y-distribution of Tracks Projected on %s GEM ", gemdet),50,-50.0,50.0);
                        
                          
                          //Checking for hits on 3rd GEM
                          for (unsigned int g20=0; g20<clusters->hits.size(); g20++)
                          {
                              if ((clusters->hits[g20].GEMid!=zgemd))  continue;
                              gem2hits[zgemd]++;  // total # of hits on the third GEM
			      
                              
                              xe1[3] = clusters->hits[g20].xl*0.4-50.;
                              ye1[3] = clusters->hits[g20].yl*0.4-50.;
                              
                              
                              dx1=xe1[3]-X3; //vertical residue of any cluseter on GEM2
                              dy1=ye1[3]-Y3; //hosrizontal residue of any cluseter on GEM2
                              
                              //H1(dx1, Form("Efficiency/%s/xresiduagem%s", gemdet, gemdet), Form("Vert. Residua on %s GEM (clusters >=1)", gemdet),100, -50.0, 50.0);
                              //H1(dy1, Form("Efficiency/%s/yresiduagem%s", gemdet, gemdet), Form("Hori. Residua on %s GEM (clusters >=1)", gemdet),100, -50.0, 50.0);
                              
                          };
                          
                          if(gem2hits[zgemd]==0)
                          {
                              
                              //H2(X3,Y3, Form("Efficiency/%s/tracksprojectedgem%s_nohit", gemdet, gemdet), Form("No Clusters On %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);
                          };
                          
                          int cluster_charge_gem2[gem2hits[zgemd]];
                          if(gem2hits[zgemd]>=1)
                          {
                              bool good_cluster_gem2=false;
                              double max_charge_gem2=0;
                              int gem2_cluster_number=-1;
                              for (unsigned int g2=0; g2<clusters->hits.size(); g2++)
                              {
                                  if ((clusters->hits[g2].GEMid!=zgemd)) continue;
                                  cluster_charge_gem2[g2]=clusters->hits[g2].quality;
                                  if(cluster_charge_gem2[g2]>max_charge_gem2)
                                  {
                                      good_cluster_gem2=true;
                                      max_charge_gem2=cluster_charge_gem2[g2];
                                      gem2_cluster_number=g2;
                                  };
                              };
                              
                              if(good_cluster_gem2)
                              { 
				  //H2(X3,Y3, Form("Efficiency/%s/tracksprojectedgem_%s", gemdet,gemdet), Form("Distribution of Tracks Projected on %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);
                                  //printf("good cluster found on %d GEM2\n",zgemd);
                                  xe[3] = clusters->hits[gem2_cluster_number].xl*0.4-50.;
                                  ye[3] = clusters->hits[gem2_cluster_number].yl*0.4-50.;
                                  //printf("GEM2 cluster position is  (%5.2f,%5.2f)\n",xe[3],ye[3]);
                                 
                                  dx=xe[3]-X3; //vertical residue of max.quality cluster
                                  dy=ye[3]-Y3; //hosrizontal residue of max.quality cluster

				  //H2(xe[3],ye[3], Form("Efficiency/%s/ActualGoodClusters_GEM_%s", gemdet,gemdet), Form("Actual clusters on %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);
                          
                           H1(xe[3], Form("Efficiency/%s/ActualGoodClusters_GEM_%s_X", gemdet, gemdet), Form("X-Distribution of Actual clusters on %s GEM ", gemdet),50,-50.0,50.0);
                          
                          H1(ye[3], Form("Efficiency/%s/ActualGoodClusters_GEM_%s_Y", gemdet, gemdet), Form("Y-distribution of Actual clsuters on %s GEM ", gemdet),50,-50.0,50.0);

                                  // Here we apply threshold condition for the residue
                                  if((fabs(dx)<threshold_residue) && (fabs(dy)<threshold_residue))
                                  {
                                      if(zgemd==0)
                                      { truehitsUS=truehitsUS+1;};
                                      if(zgemd==1)
                                      { truehits4th=truehits4th+1;};
                                      if(zgemd==2)
                                      { truehitsMS=truehitsMS+1;};

				     //H2(xe[3],ye[3], Form("Efficiency/%s/ActualGoodClusters_GEM_%s", gemdet,gemdet), Form("Actual clusters on %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);
			 	    // Here one has to fill the original projected location if there is a hit close to that projected location in order to just obtain the efficiency of the interested GEM element
				   H2(X3,Y3, Form("Efficiency/%s/ActualGoodClusters_GEM_%s", gemdet,gemdet), Form("Actual clusters on %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);
                                  };
                                  
                                  //printf("GEM2 residue (x,y) = (%5.2f,%5.2f)\n",dx,dy);
                         	                          
                                  H1(dx, Form("Efficiency/%s/x_residue%s", gemdet, gemdet), Form("Verticle residue on %s GEM", gemdet),50,-50.0,50.0);
                          
                                  H1(dy, Form("Efficiency/%s/y_residue%s", gemdet, gemdet), Form("Horizontal residue on %s GEM", gemdet),50,-50.0,50.0);
                                  
                                  
                                  
                                  
                              }; //If GEM2 has a good cluster with max charge
                          }; //If GEM2 has clusters >
                          
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
    
    auto effUS=dH2(TString::Format("Efficiency/%s/Efficiency_Plot_GEM_%s","US","US"),"Efficency Plot US GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto didUS=dH2(TString::Format("Efficiency/%s/ActualGoodClusters_GEM_%s","US","US"),"Actual clusters on US GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto shouldUS=dH2(TString::Format("Efficiency/%s/tracksprojectedgem_%s","US","US"),"Distribution of Tracks Projected on US GEM",50,-50.0,50.0,50,-50.0,50.0);
    effUS->Divide(didUS,shouldUS);

    auto eff4TH=dH2(TString::Format("Efficiency/%s/Efficiency_Plot_GEM_%s","4TH","4TH"),"Efficency Plot 4TH GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto did4TH=dH2(TString::Format("Efficiency/%s/ActualGoodClusters_GEM_%s","4TH","4TH"),"Actual clusters on 4TH GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto should4TH=dH2(TString::Format("Efficiency/%s/tracksprojectedgem_%s","4TH","4TH"),"Distribution of Tracks Projected on 4TH GEM",50,-50.0,50.0,50,-50.0,50.0);
    eff4TH->Divide(did4TH,should4TH);

    auto effMS=dH2(TString::Format("Efficiency/%s/Efficiency_Plot_GEM_%s","MS","MS"),"Efficency Plot MS GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto didMS=dH2(TString::Format("Efficiency/%s/ActualGoodClusters_GEM_%s","MS","MS"),"Actual clusters on MS GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto shouldMS=dH2(TString::Format("Efficiency/%s/tracksprojectedgem_%s","MS","MS"),"Distribution of Tracks Projected on MS GEM",50,-50.0,50.0,50,-50.0,50.0);
    effMS->Divide(didMS,shouldMS);

    double testUSnumer=didUS->GetEntries();
    double testUSdenom=shouldUS->GetEntries();
    pre_eff_US=testUSnumer/testUSdenom*100;
    //pre_eff_US=truehitsUS/(double)US_count*100;
    
    double test4THnumer=did4TH->GetEntries();
    double test4THdenom=should4TH->GetEntries();
    pre_eff_4th=test4THnumer/test4THdenom*100;
    //pre_eff_4th=truehits4th/(double)Fourth_count*100;
    
    double testMSnumer=didMS->GetEntries();
    double testMSdenom=shouldMS->GetEntries();
    pre_eff_MS=testMSnumer/testMSdenom*100;
    //pre_eff_MS=truehitsMS/(double)MS_count*100;

    printf("US GEM efficiency is %5.2f%\n",pre_eff_US);
    printf("4TH GEM efficiency is %5.2f%\n",pre_eff_4th);
    printf("MS GEM efficiency is %5.2f%\n",pre_eff_MS);

    printf("GEM_NUM = %d\n",GEM_NUM);

    
    return Plugin::ok;
}

