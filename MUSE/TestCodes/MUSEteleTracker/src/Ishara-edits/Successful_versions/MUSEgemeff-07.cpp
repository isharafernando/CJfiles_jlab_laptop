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

    
    multiplicityUS=dH1("Efficiency/US/multiplicityUS","Multiplicity US GEM",11,-0.5,10.5);
    multiplicity4TH=dH1("Efficiency/4TH/multiplicity4TH","Multiplicity 4TH GEM",11,-0.5,10.5);
    multiplicityMS=dH1("Efficiency/MS/multiplicityMS","Multiplicity MS GEM",11,-0.5,10.5);
    //multiplicityDS=dH1("Efficiency/DS/multiplicityDS","Multiplicity DS GEM",11,-0.5,10.5);
    
    
    return 0;
    
}

Long_t MUSEteleTracker::process_efficiency()
{
    // Distances from PSI
    const double zgempos[3]= {-350.0,-412.0,-474.0};
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
    
    // define this in .h file > inefficiency_gem2 = 0
    
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
	        printf("found max cluster charge\n");
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

        
        
        if ((gems==0)&&test) H2(xmaxchg[1],ymaxchg[1], Form("Efficiency/US/max_charge_clusters_USGEM"), Form("Maximum Charge Cluster positions on US GEM"),50,-50.0,50.0,50,-50.0,50.0);
        if ((gems==1)&&test) H2(xmaxchg[1],ymaxchg[1], Form("Efficiency/4TH/max_charge_clusters_4THGEM"), Form("Maximum Charge Cluster positions on 4TH GEM"),50,-50.0,50.0,50,-50.0,50.0);
        if ((gems==2)&&test) H2(xmaxchg[1],ymaxchg[1], Form("Efficiency/MS/max_charge_clusters_MSGEM"), Form("Maximum Charge Cluster positions on MS GEM"),50,-50.0,50.0,50,-50.0,50.0);
        
    };

    printf("clusters US=%d ; 4th=%d ; MS=%d; DS=%d ; \n",gemanyhit[0], gemanyhit[1], gemanyhit[2], gemanyhit[3]);
    
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
	       //test_charge=0;
               cluster_charge_gem0[g0]=clusters->hits[g0].quality;
	       //test_charge=clusters->hits[g0].quality;
               //printf("checking max charge on GEM0 \n");
		//printf("cluster charge is %d\n",cluster_charge_gem0[g0]);
		//printf("test charge is %d\n",test_charge);
                if(cluster_charge_gem0[g0]>max_charge_gem0)
                {
                    good_cluster_gem0=true;
                    max_charge_gem0=cluster_charge_gem0[g0];
                    gem0_cluster_number=g0;
		    printf("found the max charge cluseter on GEM0 \n");
                };
            };
            
            if(good_cluster_gem0)
            {
              printf("good cluster found on GEM0\n");
              xe[1] = clusters->hits[gem0_cluster_number].xl*0.4-50.;
              ye[1] = clusters->hits[gem0_cluster_number].yl*0.4-50.;
              printf("GEM0 cluster position is  (%5.2f,%5.2f)\n",xe[1],ye[1]);
                
              
              // Then start the next loop over next 2 GEMs only if there is a good cluster on first GEM
              
              for (int gem1=gem0+1; gem1<NG; gem1++)
              {
                  for (unsigned int g10=0; g10<clusters->hits.size(); g10++)
                  {
                      if ((clusters->hits[g10].GEMid!=gem1))  continue;
                      gem1hits[gem1]++;  // total # of hits on the second GEM
 		      //hits1=hits1+1;
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
			  printf("good cluster found on GEM1\n");
                          
                          xe[2] = clusters->hits[gem1_cluster_number].xl*0.4-50.;
                          ye[2] = clusters->hits[gem1_cluster_number].yl*0.4-50.;
                          printf("GEM1 cluster position is  (%5.2f,%5.2f)\n",xe[2],ye[2]);
                          
                          double slopexe=(xe[2]-xe[1])/(zgempos[gem1]-zgempos[gem0]);
                          double slopeye=(ye[2]-ye[1])/(zgempos[gem1]-zgempos[gem0]);
                          
                          if ((gem0==0)&&(gem1==1))
                          {
                              //std::cout << "MS" << std::endl;
			      printf("MS GEM found\n");
                              zgemd=2;
                              gemdet="MS";
 			      MS_count=MS_count+1;
                          };
                          
                          if ((gem0==0)&&(gem1==2))
                          {
			      printf("4TH GEM found\n");
                              //std::cout << "4TH" << std::endl;
                              zgemd=1;
                              gemdet="4TH";
			      Fourth_count=Fourth_count+1;
                          };
                          
                          if ((gem0==1)&&(gem1==2))
                          {
			      printf("US GEM found\n");
                              //std::cout << "US" << std::endl;
                              zgemd=0;
                              gemdet="US";
			      US_count=US_count+1;
                          };
                          
                          //Project the tracks on GEM2 (third gem)
                          double X3 = xe[2] + slopexe*(zgempos[zgemd]-zgempos[gem0]);
                          double Y3 = ye[2] + slopeye*(zgempos[zgemd]-zgempos[gem0]);
                          
			  printf("x,y slopes are = (%5.2f,%5.2f) \n",slopexe,slopeye);
                          
                          H2(X3,Y3, Form("Efficiency/%s/tracksprojectedgem%s_any", gemdet,gemdet), Form("Tracks Projected on %s GEM-Any cluster", gemdet),50,-50.0,50.0,50,-50.0,50.0);
                          
                           H1(X3, Form("Efficiency/%s/tracksprojectedgem%sX_any", gemdet, gemdet), Form("Tracks Projected on %s GEM X-Any cluster", gemdet),50,-50.0,50.0);
                          
                          H1(Y3, Form("Efficiency/%s/tracksprojectedgem%sY_any", gemdet, gemdet), Form("Tracks Projected on %s GEM Y-Any cluster", gemdet),50,-50.0,50.0);
                        
                          
                          bool gem2hit=true;
                          //Checking for hits on 3rd GEM
                          for (unsigned int g20=0; g20<clusters->hits.size(); g20++)
                          {
                              if ((clusters->hits[g20].GEMid!=zgemd))  continue;
                              gem2hits[zgemd]++;  // total # of hits on the third GEM
			      //hits2=hits2+1;
                              
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
				 printf("good cluster found on %d GEM2\n",zgemd);
                                 Eff_Den=Eff_Den+1;
				 gem2realhits[zgemd]=gem2realhits[zgemd]+1;
				// if(zgemd==0)
				//{ truehitsUS=truehitsUS+1;};
    				// if(zgemd==1)
				//{ truehits4th=truehits4th+1;};
				// if(zgemd==2)
				//{ truehitsMS=truehitsMS+1;};
                                  xe[3] = clusters->hits[gem2_cluster_number].xl*0.4-50.;
                                  ye[3] = clusters->hits[gem2_cluster_number].yl*0.4-50.;
                                  printf("GEM2 cluster position is  (%5.2f,%5.2f)\n",xe[3],ye[3]);
                                 
                                  dx=xe[3]-X3; //vertical residue
                                  dy=ye[3]-Y3; //hosrizontal residue

				  if((fabs(dx)<threshold_residue) && (fabs(dy)<threshold_residue))
				  {
				   if(zgemd==0)
				   { truehitsUS=truehitsUS+1;};
    				   if(zgemd==1)
				   { truehits4th=truehits4th+1;};
				   if(zgemd==2)
				   { truehitsMS=truehitsMS+1;};
				  };
                                  
				  printf("GEM2 residue (x,y) = (%5.2f,%5.2f)\n",dx,dy);
                         	                          
                           	  H1(dx, Form("Efficiency/%s/x_residue%s", gemdet, gemdet), Form("Verticle residue on %s GEM", gemdet),50,-50.0,50.0);
                          
                                  H1(dy, Form("Efficiency/%s/y_residue%s", gemdet, gemdet), Form("Horizontal residue on %s GEM", gemdet),50,-50.0,50.0);
                                  //H2(X3,Y3, Form("Efficiency/%s/tracksprojectedgem%s_oneormore_cut", gemdet, gemdet), Form("Tracks Projected on %s GEM--cluster >=1, Inside Vicinity", gemdet),50,-50.0,50.0,50,-50.0,50.0);
                                  //H2(xe[3],ye[3], Form("Efficiency/%s/tracksdetectedgem%s_oneormore_cut", gemdet, gemdet), Form("Cluster positions on %s GEM--cluster >=1, Inside Vicinity", gemdet),50,-50.0,50.0,50,-50.0,50.0);
                                  
                                  
                                  
                              }; //If GEM2 has a good cluster with max charge
                          }; //If GEM2 has clusters >
                          
                      } // If GEM1 has a good cluster with max charge
                  }; //If GEM1 has clusters >=
                  
              }; // Loop for GEM1
                
            }; // If GEM0 has a good cluster with max charge
        }; // If GEM0 has clusters >=1
            
            
    };//Loop for GEM0
        
    //printf("gem0hits is %d, %d, %d \n",gem0hits[0],gem0hits[1],gem0hits[2]);
    //printf("gem1hits is %d, %d, %d \n",gem1hits[0],gem1hits[1],gem1hits[2]);
    //printf("gem2hits is %d, %d, %d \n",gem2hits[0],gem2hits[1],gem2hits[2]);
    
    
    
    return 0;
    //  return Plugin::ok;
};

Long_t MUSEteleTracker::finalize_efficiency()
{
    
    //Double USX_n = tracksprojectedgemUSX_onecut->GetEntries();
    //Double USX_d = tracksprojectedgemUSX_any->GetEntries();
    //Double USX_eff=USX_n/USX_d;
    //pre_eff_US=US_count/Eff_Den;
    //pre_eff_4th=Fourth_count/Eff_Den;
    //pre_eff_MS=MS_count/Eff_Den;
    
    //printf("Number of good tracks %d\n",Eff_Den);

    //printf("US GEM projected hits as 3rd is %d\n",US_count);
    //printf("4TH GEM projected hits as 3rd is %d\n",Fourth_count);
    //printf("MS GEM projected hits as 3rd is %d\n",MS_count);
     
    //truehitsUS=gem2realhits[0];
    //truehits4th=gem2realhits[1];
    //truehitsMS=gem2realhits[2];

    //printf("US GEM real hits as 3rd is %d\n",truehitsUS);
    //printf("4TH GEM real hits as 3rd is %d\n",truehits4th);
    //printf("MS GEM real hits as 3rd is %d\n",truehitsMS);


    pre_eff_US=truehitsUS/(double)US_count*100;
    pre_eff_4th=truehits4th/(double)Fourth_count*100;
    pre_eff_MS=truehitsMS/(double)MS_count*100;

    printf("US GEM efficiency is %5.2f%\n",pre_eff_US);
    printf("4TH GEM efficiency is %5.2f%\n",pre_eff_4th);
    printf("MS GEM efficiency is %5.2f%\n",pre_eff_MS);
    //printf("US X efficiency %d %3.2lf\n",USX_n,USX_eff);
    //printf("Inefficiency of GEM2 is = %d \n",inefficiency_gem2);
    //printf("GEM0 hits = %d,  GEM1 hits = %d,  GEM2 hits = %d \n",hits0,hits1,hits2);
    
    //USX_eff->Divide(tracksprojectedgemUSX_onecut,tracksprojectedgemUSX_any,1,1,"");
    //MSX_eff->Divide(tracksprojectedgemMSX_onecut,tracksprojectedgemMSX_any,1,1,"");
    //DSX_eff->Divide(tracksprojectedgemDSX_onecut,tracksprojectedgemDSX_any,1,1,"");
    
    //USY_eff->Divide(tracksprojectedgemUSY_onecut,tracksprojectedgemUSY_any,1,1,"");
    //MSY_eff->Divide(tracksprojectedgemMSY_onecut,tracksprojectedgemMSY_any,1,1,"");
    //DSY_eff->Divide(tracksprojectedgemDSY_onecut,tracksprojectedgemDSY_any,1,1,"");
    
    return Plugin::ok;
}

