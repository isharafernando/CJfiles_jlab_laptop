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
      
    return 0;
    
}

Long_t MUSEteleTracker::process_efficiency()
{
    int NG=3;  // GEM_NUM or Number of GEMs

    const double zgempos[3]= {-536.0,-474.0,-412.0};
    const char* gemdet;
    
    // loop over all possible combinations of clusters:
    double xe[GEM_NUM]={-10000}, ye[GEM_NUM]={-10000}, xe_max[GEM_NUM]={-10000}, ye_max[GEM_NUM]={-10000};
    double dx_max[GEM_NUM]={-10000},dy_max[GEM_NUM]={-10000}, dxe={-10000}, dye={-10000};
    double x[GEM_NUM]={-10000.0},y[GEM_NUM]={-10000.0},dx={-10000.0},dy={-10000.0};
    double xmaxchg[GEM_NUM]={-10000};
    double ymaxchg[GEM_NUM]={-10000};
      
    
    /* BH filtering 1 definitions start */
    
    double alignment[16] = {3.0,1.25,3.8,0.75,3.5,1.75,1.6,2.45,1.2,.25,.1,0,4.5,3.0,2.5,2.5};
    int plane2hits = 0;
    int plane3hits = 0;
    
    
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
      
    /* BH filtering 1 definitions finished */

    

    
        
      
    /* Cluster Multiplicity before BH filtering 1 */
    

    int gemmultibefore[3] = { 0, 0, 0 };
    for (unsigned int g1=0; g1<clusters->hits.size(); g1++)
    {
        // If the GEMid is different than the cluster's gem ID, then g1 gets incremented by 1 unit
        if (clusters->hits[g1].GEMid!=0) continue;
        gemmultibefore[0]++;
        for (unsigned int g2=0; g2<clusters->hits.size(); g2++)
        {
            if (clusters->hits[g2].GEMid!=1) continue;
            gemmultibefore[1]++;
            for (unsigned int g3=0; g3<clusters->hits.size(); g3++)
            {
                if (clusters->hits[g3].GEMid!=2) continue;
                gemmultibefore[2]++;
                
            };
        };
    };
    

    H1(gemmultibefore[0], Form("Efficiency/US/Cluster_Multiplicity_US"), Form("Cluster Multiplicity of US GEM before BH cluster filtering"),
       11,-0.5,10.5);
    H1(gemmultibefore[1], Form("Efficiency/4TH/Cluster_Multiplicity_4TH"), Form("Cluster Multiplicity of 4TH GEM before BH cluster filtering"),
       11,-0.5,10.5);
    H1(gemmultibefore[2], Form("Efficiency/MS/Cluster_Multiplicity_MS"), Form("Cluster Multiplicity of MS GEM before BH cluster filtering"),
       11,-0.5,10.5);
    
    
    // BH Cluster filtering start

    if(BHPlaneCount[2]==1 and BHPlaneCount[3]==1)
    {

       /* Cluster Multiplicity after BH filtering 1 */

     int gemmultiafter[3] = { 0, 0, 0 };
    for (unsigned int g1=0; g1<clusters->hits.size(); g1++)
    {
        if (clusters->hits[g1].GEMid!=0) continue;
        gemmultiafter[0]++;
        for (unsigned int g2=0; g2<clusters->hits.size(); g2++)
        {
            if (clusters->hits[g2].GEMid!=1) continue;
            gemmultiafter[1]++;
            for (unsigned int g3=0; g3<clusters->hits.size(); g3++)
            {
                if (clusters->hits[g3].GEMid!=2) continue;
                gemmultiafter[2]++;
                
            };
        };
    };
    

    H1(gemmultiafter[0], Form("Efficiency/US/BH_filtered_Cluster_Multiplicity_US"), Form("Cluster Multiplicity of US GEM after BH cluster filtering"),
       11,-0.5,10.5);
    H1(gemmultiafter[1], Form("Efficiency/4TH//BH_filteredCluster_Multiplicity_4TH"), Form("Cluster Multiplicity of 4TH GEM after BH cluster filtering"),
       11,-0.5,10.5);
    H1(gemmultiafter[2], Form("Efficiency/MS//BH_filteredCluster_Multiplicity_MS"), Form("Cluster Multiplicity of MS GEM after BH cluster filtering"),
       11,-0.5,10.5);


    
    
    // Loop over the GEM in interest
    int gem0hits[3]={0,0,0};

      for (int gem0=0; gem0<NG; gem0++)  
    {
        
        int gem1hits[3]={0,0,0};
        double xg[2]={0.0,0.0};
        double yg[2]={0.0,0.0};
        int gemID[2]={0,0};
        for (int gem1=0; gem1<NG; gem1++)
              {

                if (gem1!=gem0)
                {
                  for (unsigned int g10=0; g10<clusters->hits.size(); g10++)
                  {
                      if ((clusters->hits[g10].GEMid!=gem1))  continue;
                      gem1hits[gem1]++;  // total # of hits on the 2nd and 3rd GEM
          
                  };
                  
                  if (gem1hits[gem1]==0) continue;

                  int cluster_charge_gem1[gem1hits[gem1]];

                  // Checking one or more clusters found

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
                      }; // Looking for a good cluster (loop over clusters)

                      if(good_cluster_gem1)
                      {
                        xg[gem1] = clusters->hits[gem1_cluster_number].xl*0.4-50.;
                        yg[gem1] = clusters->hits[gem1_cluster_number].yl*0.4-50.;

                        gemID[gem1]=gem1;

                      };

                  }; // If one or more clusters found in jth GEM element

                }; // Closing the if condition for not to loop over again the GEM in interest


              }; // Loop over other two GEMs

              //Calculating the slope only if the clusters are in the active area
              if ((fabs(xg[gemID[1]])<50)&&(fabs(yg[gemID[1]])<50)&&(fabs(xg[gemID[2]])<50)&&(fabs(yg[gemID[2]])<50))
              {
                double slopexg=(xg[2]-xg[1])/(zgempos[gemID[2]]-zgempos[gemID[1]]);
                double slopeyg=(yg[2]-yg[1])/(zgempos[gemID[2]]-zgempos[gemID[1]]);

                //printf("Slope x = %5.2f \n",slopexg);
                //printf("Slope y = %5.2f \n",slopeyg);
                printf("zgemposes gem ID %d = %5.2f\n",gemID[1],zgempos[gemID[1]]);
                printf("zgemposes gem ID %d = %5.2f\n",gemID[2],zgempos[gemID[2]]);


                // Looking for a gool cluster on GEM0
                int cluster_charge_gem0[gem0hits[gem0]];
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
                      printf("found the max charge cluseter on GEM0 \n");
                    };
                };


                }; // If condition for checking the GEM1 clusters are in the active area

    }; // Loop over the GEM in interest


    }; //BH cluster filtering finish

    return 0;
    //  return Plugin::ok;
};

Long_t MUSEteleTracker::finalize_efficiency()
{
    
    /*
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
    */

    
    return Plugin::ok;
}


