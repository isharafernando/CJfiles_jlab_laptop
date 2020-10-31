/* GEM efficiency (integrated with BH) Version by Ishara Fernando */
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
    
    handle.geom();
    // get input branch with GEM clusters:
    clusters=NULL;
    getBranchObject("LumiGEMhits", (TObject**)&clusters);
    if (clusters==NULL)
    {
        printf(" Cannot find branch >LumiGEMhits< in input ROOT file - trying output branch\n");
        getOutBranchObject("LumiGEMhits",(TObject**)&clusters);
        if(clusters==NULL)
        {
            printf("Couldn't find any clusters in any ROOT file :(\n");
            return -1;
        }
    };
    printf(" Scinthits (Hits) @%p\n", Hits);
    
    
    // get input branch for BH
    Hits=NULL;
    getBranchObject("ScintHits", (TObject**)&Hits);
    if (Hits==NULL)
    {
        printf(" Cannot find branch >Scinthits< in input ROOT file - trying output branch\n");
        getOutBranchObject("ScintHits",(TObject**)&Hits);
        if(Hits==NULL)
        {
            printf("Couldn't find any Hits in any ROOT file :(\n");
            return -1;
        }
    };
    printf(" Scinthits (Hits) @%p\n", Hits);

    
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
    //double zgempos[3]={0,0,0};

    //for(unsigned int k=0; k<3; k++){
//		  double loc[3]={0,0,0};
//		  double mas[3]={0,0,0};
//		  
//		  handle.GEM[k]->LocalToMaster(loc,mas);
//
//		  zgempos[k] = mas[2]*10;
//		  
//	  }
    
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
    
    
    /* BH section start */
    
    double alignment[16] = {3.0,1.25,3.8,0.75,3.5,1.75,1.6,2.45,1.2,.25,.1,0,4.5,3.0,2.5,2.5};
    int particletype = 0;
    int plane2hits = 0;
    int plane3hits = 0;
    std::string realornoise[] = {"Real Hits", "Noise Hits"};
    std::string particlenames[] = {"Pion", "Electron", "Muon"};
    
    std::vector< std::vector<int> > BHhits(Hits->hits.size(), std::vector<int>(2));
    std::vector< std::vector<double> > BHhitxyz(Hits->hits.size(), std::vector<double>(3));

    
    int BHPlaneCount[4]={0,0,0,0};
    int Plane2=0;
    int Plane2bar=0;
    double Plane2crds[3]={0,0,0};
    int Plane3=0;
    int Plane3bar=0;
    double Plane3crds[3]={0,0,0};
    for(size_t i = 0; i < Hits->hits.size(); i++)
    {
        double RF = Hits->hits[i].rf;
        int plane, side, bar;
        
        BH_internal_to_logic(Hits->hits[i].id,&plane,&bar,&side);
        
        double loc[3] = {0,0,0};
        double mas[3] = {0,0,0};
        handle.BH[plane][bar]->LocalToMaster(loc,mas);
        
        BHhits[i][0] = plane;
        BHhits[i][1] = bar;
        
        BHhitxyz[i][0] = mas[0]*10;
        BHhitxyz[i][1] = mas[1]*10;
        BHhitxyz[i][2] = mas[2]*10;

        if(plane==2)
        { 
 	Plane2=plane;
        Plane2bar = bar;
        Plane2crds[0] = mas[0]*10;
        Plane2crds[1] = mas[1]*10;
        Plane2crds[2] = mas[2]*10;
        BHPlaneCount[plane]++;
        }; 

       if(plane==3)
        { 
	Plane3=plane;
        Plane3bar = bar;
        Plane3crds[0] = mas[0]*10;
        Plane3crds[1] = mas[1]*10;
        Plane3crds[2] = mas[2]*10;
        BHPlaneCount[plane]++;
        };      
       
    };
    
    //printf("BH plane 2 hits = %d ; and plane 3 hits = %d\n",BHPlaneCount[2],BHPlaneCount[3]);

    if(BHPlaneCount[2]==1 and BHPlaneCount[3]==1)
    {
     //printf("BH plane 2 hits = %d ; and plane 3 hits = %d\n",BHPlaneCount[2],BHPlaneCount[3]);
     //printf("BH plane %d bar %d : hit coordinates are  (%7.4f,%7.4f,%7.4f)\n",Plane2,Plane2bar,Plane2crds[0],Plane2crds[1],Plane2crds[2]);
     //printf("BH plane %d bar %d : hit coordinates are  (%7.4f,%7.4f,%7.4f)\n",Plane3,Plane3bar,Plane3crds[0],Plane3crds[1],Plane3crds[2]);
  
    /* BH section finished */

    
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
                
              if((good_cluster_gem0)&&(fabs(xe[1])<50)&&(fabs(ye[1])<50))
                {
              
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
                          
                          if((good_cluster_gem1)&&(fabs(xe[2])<50)&&(fabs(ye[2])<50))
                          {
                          
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

                            
                           if((fabs(X3)<50)&&(fabs(Y3)<50))
                              {

                          
			  //Project the tracks on BH planes
			  double BHplane2X = xe[1] + slopexe*(Plane2crds[2]-zgempos[gem0]);
                          double BHplane2Y = ye[1] + slopeye*(Plane2crds[2]-zgempos[gem0]);

                          double BHplane3X = xe[1] + slopexe*(Plane3crds[2]-zgempos[gem0]);
                          double BHplane3Y = ye[1] + slopeye*(Plane3crds[2]-zgempos[gem0]); 
                          
			  //printf("Projected BH plane %d : hit coordinates are  (%7.4f,%7.4f,%7.4f)\n",Plane2,BHplane2X,BHplane2Y,Plane2crds[2]);
     			  //printf("Projected BH plane %d : hit coordinates are  (%7.4f,%7.4f,%7.4f)\n",Plane3,BHplane3X,BHplane3Y,Plane3crds[2]);

                           double BHdx=Plane3crds[0]-BHplane3X;
			   double BHdy=Plane2crds[1]-BHplane2Y;

                           H1(BHdx, Form("Efficiency/%s/BH_residuals_of_gem_%s_X", gemdet, gemdet), Form("Residual X-Distribution of Tracks Projected on BH to evaluate %s GEM ", gemdet),50,-50.0,50.0);
                          
                           H1(BHdy, Form("Efficiency/%s/BH_residuals_of_gem_%s_Y", gemdet, gemdet), Form("Residual Y-distribution of Tracks Projected on BH to evaluate %s GEM ", gemdet),50,-50.0,50.0);


			   if((fabs(BHdx)<4.0000) && (fabs(BHdy)<4.0000))
			   {
                            
                           //printf("Filteration for BH\n");
                           //printf("Plane 2 x coordinates of the bar = %7.4f; and projected %7.4f\n",Plane2crds[0],BHplane2X);
			   //printf("Plane 3 x coordinates of the bar = %7.4f, and projected %7.4f\n",Plane3crds[0],BHplane3X);
                           //printf("Plane 2 y coordinates of the bar = %7.4f; and projected %7.4f\n",Plane2crds[1],BHplane2Y);
			   //printf("Plane 3 y coordinates of the bar = %7.4f, and projected %7.4f\n",Plane3crds[1],BHplane3Y);
                           //printf("BH dx = %7.4f\n",fabs(BHdx));
			   //printf("BH dy = %7.4f\n",fabs(BHdy));

                          
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
                          
                          H1(xe[3], Form("Efficiency/%s/ActualClusters_GEM_%s_X", gemdet, gemdet), Form("X-Distribution of Actual clusters on %s GEM ", gemdet),50,-50.0,50.0);
                          
                          H1(ye[3], Form("Efficiency/%s/ActualClusters_GEM_%s_Y", gemdet, gemdet), Form("Y-distribution of Actual clsuters on %s GEM ", gemdet),50,-50.0,50.0);
                          H2(xe[3],ye[3], Form("Efficiency/%s/ActualClusters_GEM_%s", gemdet,gemdet), Form("Actual clusters on %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);

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
				   H2(xe[3],ye[3], Form("Efficiency/%s/ActualGoodClusters_GEM_%s", gemdet,gemdet), Form("Actual Good clusters on %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);
				   H2(X3,Y3, Form("Efficiency/%s/ActualProjectedGoodClusters_GEM_%s", gemdet,gemdet), Form("Actual Projected Good clusters on %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);
                                  };
                                  
                                  //printf("GEM2 residue (x,y) = (%5.2f,%5.2f)\n",dx,dy);
                         	                          
                                  H1(dx, Form("Efficiency/%s/x_residue%s", gemdet, gemdet), Form("Verticle residue on %s GEM", gemdet),50,-50.0,50.0);
                          
                                  H1(dy, Form("Efficiency/%s/y_residue%s", gemdet, gemdet), Form("Horizontal residue on %s GEM", gemdet),50,-50.0,50.0);
                                  
                                  
                                  
                                  
                              }; //If GEM2 has a good cluster with max charge
                          }; //If GEM2 has clusters >
                          
			}; //BH filtering
                          };// Projected cluster on GEM2 in the active area < 50mm
                     };// GEM1 good cluster in the active area < 50mm

                      }; // If GEM1 has a good cluster with max charge
                  }; //If GEM1 has clusters >=
                  
              }; // Loop for GEM1
                
                
            };// GEM0 good cluster in the active area < 50mm
            
            }; // If GEM0 has a good cluster with max charge
        }; // If GEM0 has clusters >=1
            
            
    };//Loop for GEM0
    
    }; //BH hits 1 on each plane condition closing
    
    return 0;
    //  return Plugin::ok;
};

Long_t MUSEteleTracker::finalize_efficiency()
{
    
    auto effUS=dH2(TString::Format("Efficiency/%s/Efficiency_Plot_GEM_%s","US","US"),"Efficency Plot US GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto didUS=dH2(TString::Format("Efficiency/%s/ActualProjectedGoodClusters_GEM_%s","US","US"),"Actual Projected Good clusters on US GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto shouldUS=dH2(TString::Format("Efficiency/%s/tracksprojectedgem_%s","US","US"),"Distribution of Tracks Projected on US GEM",50,-50.0,50.0,50,-50.0,50.0);
    effUS->Divide(didUS,shouldUS);

    auto eff4TH=dH2(TString::Format("Efficiency/%s/Efficiency_Plot_GEM_%s","4TH","4TH"),"Efficency Plot 4TH GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto did4TH=dH2(TString::Format("Efficiency/%s/ActualProjectedGoodClusters_GEM_%s","4TH","4TH"),"Actual Projected Good clusters on 4TH GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto should4TH=dH2(TString::Format("Efficiency/%s/tracksprojectedgem_%s","4TH","4TH"),"Distribution of Tracks Projected on 4TH GEM",50,-50.0,50.0,50,-50.0,50.0);
    eff4TH->Divide(did4TH,should4TH);

    auto effMS=dH2(TString::Format("Efficiency/%s/Efficiency_Plot_GEM_%s","MS","MS"),"Efficency Plot MS GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto didMS=dH2(TString::Format("Efficiency/%s/ActualProjectedGoodClusters_GEM_%s","MS","MS"),"Actual Projected Good clusters on MS GEM",50,-50.0,50.0,50,-50.0,50.0);
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

