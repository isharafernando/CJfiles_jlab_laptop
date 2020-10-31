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
    int NG=GEM_NUM;  // GEM_NUM or Number of GEMs

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
          //Plane2crds[0] = mas[0]*10;
          //Plane2crds[1] = mas[1]*10;
          //Plane2crds[2] = mas[2]*10;
          BHPlaneCount[plane]++;
        }; 

       if(plane==3)
        { 
	        Plane3=plane;
          Plane3bar = bar;
          //Plane3crds[0] = mas[0]*10;
          //Plane3crds[1] = mas[1]*10;
          //Plane3crds[2] = mas[2]*10;
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

  //bool need_BH_filtering=true;

  if(BHPlaneCount[2]==1 and BHPlaneCount[3]==1)
  {



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
          Plane2crds[0] = BHhitxyz[i][0];
          Plane2crds[1] = BHhitxyz[i][1];
          Plane2crds[2] = BHhitxyz[i][2];
          
        }; 

       if(plane==3)
        { 
          Plane3=plane;
          Plane3bar = bar;
          Plane3crds[0] = BHhitxyz[i][0];
          Plane3crds[1] = BHhitxyz[i][1];
          Plane3crds[2] = BHhitxyz[i][2];
          
        };      
       
    };


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
    //int gem0hits[3]={0,0,0};

    for (int gem0=0; gem0<NG; gem0++)  
    {   
      if (gem0==0)
      {gemdet="US";};
      if (gem0==1)
      {gemdet="4TH";};
      if (gem0==2)
      {gemdet="MS";};
      
        
      int gem1hits[3]={0,0,0};
      double xg[2]={0.0000,0.0000};
      double yg[2]={0.0000,0.0000};
      int gemID[2]={0,0};
      //printf("GEM0 ID = %d\n",gem0);
      int gem1count=0; //count over GEM1
      double slopexg=0.0;
      double slopeyg=0.0;


      bool GoToNextGEM0=false;
      if (GoToNextGEM0==true) continue;
      

      for (int gem1=0; ((gem1<NG) and (GoToNextGEM0==false)); gem1++)
      {


        if (gem1==gem0) continue; // Not considering the GEM0
        
         //printf("GEM1 ID = %d\n",gem1);
        for (unsigned int g10=0; g10<clusters->hits.size(); g10++)
        {
         if ((clusters->hits[g10].GEMid!=gem1))  continue;
          gem1hits[gem1]++;  // total # of hits on the 2nd and 3rd GEM
          
        };
                  
        
        if (gem1hits[gem1]==0) continue; 
        // This condition will eleminate the case where there are no good clusters found for any GEM1


        // Checking one or more clusters found

        
          //printf("GEM1 ID %d has hits = %d\n",gem1,gem1hits[gem1]);
                    
        double cluster_charge_gem1=0.0;
        bool good_cluster_gem1=false;
        double max_charge_gem1=0.0;
        int gem1_cluster_number=-1;
        // Looking for a good cluster (loop over clusters)          
        for (unsigned int g1=0; g1<clusters->hits.size(); g1++)
        {
          if ((clusters->hits[g1].GEMid!=gem1)) continue;
          cluster_charge_gem1=clusters->hits[g1].quality;
          //printf("cluster %d charge is = %5.2f \n",cluster_charge_gem1);
                          
          if(cluster_charge_gem1>max_charge_gem1)
          {
            good_cluster_gem1=true;
            max_charge_gem1=cluster_charge_gem1;
            gem1_cluster_number=g1;
          };
                          
        }; // Looking for a good cluster (loop over clusters)
                      
        if(good_cluster_gem1==false) 
        {
          continue;
          GoToNextGEM0=true;
        }; 
          
        //printf("GEM1 ID %d 's good cluster is cluster # %d \n",gem1,gem1_cluster_number); 
        //xg[gem1count] = clusters->hits[gem1_cluster_number].xl*0.4-50.;
        //yg[gem1count] = clusters->hits[gem1_cluster_number].yl*0.4-50.;

        double tempx = clusters->hits[gem1_cluster_number].xl*0.4-50.;
        double tempy = clusters->hits[gem1_cluster_number].yl*0.4-50.;

        if ((fabs(tempx)>50.0)or(fabs(tempy)>50.0))
        {
          continue;
          GoToNextGEM0=true;
        };

        xg[gem1count]=tempx;
        yg[gem1count]=tempy;


        gemID[gem1count]=gem1;

        //printf("GEM ID count  %d \n",gem1count);
        //printf("GEMID %d has a good cluster at (%5.2f,%5.2f)\n",gemID[gem1count],xg[gem1],yg[gem1]);
        //printf("GEMID %d has a good cluster at (%5.2f,%5.2f)\n",gem1count,xg[gem1count],yg[gem1count]);

        gem1count=gem1count+1;

                        
        //printf("GEM ID count  %d \n",gem1count);
        //printf("GEMID %d has a good cluster \n",gemID[gem1count]);      
      
      }; // Loop over other two GEMs

      // The following condition will skip the GEM0 to next one if we don't get two good clsuters
      if (gemID[0]==gemID[1]) continue;

      //printf("GEM in interest spot 1 %s\n",gemdet);
      // Here I add an aditional cut to filter out the cases where the cluster are found outside the active area

      //if ((fabs(xg[gemID[0]])>50.0)or(fabs(yg[gemID[0]])>50.0)or(fabs(xg[gemID[1]])>50.0)or(fabs(yg[gemID[1]])>50.0)) continue;
      //printf("zgemposes gem ID %d = %5.2f\n",gemID[0],zgempos[gemID[0]]);
      //printf("zgemposes gem ID %d = %5.2f\n",gemID[1],zgempos[gemID[1]]);
      //printf("GEM in interest spot 1 %s\n",gemdet);

      slopexg=(xg[1]-xg[0])/(zgempos[gemID[1]]-zgempos[gemID[0]]);
      slopeyg=(yg[1]-yg[0])/(zgempos[gemID[1]]-zgempos[gemID[0]]);

      //printf("GEM IDs %d, %d \n",gemID[0],gemID[1]);

      //printf("Slope x = %5.2f \n",slopexg);
      //printf("Slope y = %5.2f \n",slopeyg);

      //printf("GEM in interest %s\n",gemdet);

      //Project the tracks on BH planes
                    
      double BHplane2X = xg[1] + slopexg*(Plane2crds[2]-zgempos[gemID[1]]);
      double BHplane2Y = yg[1] + slopeyg*(Plane2crds[2]-zgempos[gemID[1]]);

      double BHplane3X = xg[1] + slopexg*(Plane3crds[2]-zgempos[gemID[1]]);
      double BHplane3Y = yg[1] + slopeyg*(Plane3crds[2]-zgempos[gemID[1]]);  

      //printf("Projected BH plane %d : hit coordinates are  (%7.4f,%7.4f),(%7.4f,%7.4f)\n",Plane2,BHplane2X,BHplane2Y,Plane2crds[0],Plane2crds[1]);
      //printf("Projected BH plane %d : hit coordinates are  (%7.4f,%7.4f),(%7.4f,%7.4f)\n",Plane3,BHplane3X,BHplane3Y,Plane3crds[0],Plane3crds[1]);

      //printf("Plane 2 z-coords = %5.2f\n", Plane2crds[2]);
      //printf("Plane 3 z-coords = %5.2f\n", Plane3crds[2]);

      //printf("Plane 2 (x,y) are (%7.4f,%7.4f)\n", Plane2crds[0],Plane2crds[1]);
      //printf("Plane 2 projected (x,y) are (%7.4f,%7.4f)\n", BHplane2X,BHplane2Y);
      //printf("Plane 3 (x,y) are (%7.4f,%7.4f)\n", Plane3crds[0],Plane3crds[1]);
      //printf("Plane 3 projected (x,y) are (%7.4f,%7.4f)\n", BHplane3X,BHplane3Y);



      double BHdx=Plane3crds[0]-BHplane3X-0.000;
      double BHdy=Plane2crds[1]-BHplane2Y;

      //printf("BH x residue %5.4f \n",BHdx);
      //printf("BH y residue %5.4f \n",BHdy );


      H1(BHdx, Form("Efficiency/%s/BH_residuals_of_gem_%s_X", gemdet, gemdet), Form("Residual X-Distribution of Tracks Projected on BH to evaluate %s GEM ", gemdet),50,-50.0,50.0);
                          
      H1(BHdy, Form("Efficiency/%s/BH_residuals_of_gem_%s_Y", gemdet, gemdet), Form("Residual Y-distribution of Tracks Projected on BH to evaluate %s GEM ", gemdet),50,-50.0,50.0);



      if((fabs(BHdx)>4.0) or (fabs(BHdy)>4.0)) continue;
      
      H1(BHdx, Form("Efficiency/%s/BH_residuals_of_gem_%s_X_after_cut", gemdet, gemdet), Form("Residual X-Distribution of Tracks Projected on BH after cut to evaluate %s GEM ", gemdet),50,-50.0,50.0);
                          
      H1(BHdy, Form("Efficiency/%s/BH_residuals_of_gem_%s_Y_after_cut", gemdet, gemdet), Form("Residual Y-distribution of Tracks Projected on BH after cut to evaluate %s GEM ", gemdet),50,-50.0,50.0);

      // Checking whether the projected tracks to BH are valid
      // If the projected tracks to BH are not valid, then go to the nexg GEM0
  
      //Project the tracks on GEM0
      double X3 = xg[0] + slopexg*(zgempos[gem0]-zgempos[gemID[0]]);
      double Y3 = yg[0] + slopeyg*(zgempos[gem0]-zgempos[gemID[0]]);


      H2(X3,Y3, Form("Efficiency/%s/tracksprojectedgem_%s", gemdet,gemdet), Form("Distribution of Tracks Projected on %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);

      if((fabs(X3)>50.0) or (fabs(Y3)>50.0)) continue;

      H2(X3,Y3, Form("Efficiency/%s/tracksprojectedgem_active_%s", gemdet,gemdet), Form("Distribution of Tracks Projected inside active area on %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);
                          
      //H1(X3, Form("Efficiency/%s/tracksprojectedgem_%s_X", gemdet, gemdet), Form("X-Distribution of Tracks Projected on %s GEM ", gemdet),50,-50.0,50.0);
                          
      //H1(Y3, Form("Efficiency/%s/tracksprojectedgem_%s_Y", gemdet, gemdet), Form("Y-distribution of Tracks Projected on %s GEM ", gemdet),50,-50.0,50.0);
                  
      // Looking for a gool cluster on GEM0
      
      int gem0hits[3]={0,0,0};
      for (unsigned int g0=0; g0<clusters->hits.size(); g0++)
      {
        if ((clusters->hits[g0].GEMid!=gem0)) continue;
        gem0hits[gem0]++;  // total number of hits on the GEM0

      };      

      if(gem0hits[gem0]==0)
      {
        H2(X3,Y3, Form("Efficiency/%s/No_hits_gem_%s", gemdet,gemdet), Form("No clusters on %s GEM after projecting track", gemdet),50,-50.0,50.0,50,-50.0,50.0);
        continue;
      };
        
      //If we have one or more clusters
      //double cluster_charge_gem0;
      bool good_cluster_gem0=false;
      double max_charge_gem0=0;
      int gem0_cluster_number=-1;

      for (unsigned int g0=0; g0<clusters->hits.size(); g0++)
      {
        if ((clusters->hits[g0].GEMid!=gem0)) continue;
        double cluster_charge_gem0=0.0;
        cluster_charge_gem0=clusters->hits[g0].quality;
        if(cluster_charge_gem0>max_charge_gem0)
        {
          good_cluster_gem0=true;
          max_charge_gem0=cluster_charge_gem0;
          gem0_cluster_number=g0;
          //printf("found the max charge cluseter on GEM0 \n");
        };
      };
      

      if(!good_cluster_gem0) continue;

      xe[gem0] = clusters->hits[gem0_cluster_number].xl*0.4-50.;
      ye[gem0] = clusters->hits[gem0_cluster_number].yl*0.4-50.;

      dx=xe[gem0]-X3; //vertical residue 
      dy=ye[gem0]-Y3; //hosrizontal residue

      //printf("residues are (%5.2f,%5.2f)\n",dx,dy);

      // && (fabs(xe[gem0])>50.0) && (fabs(ye[gem0])>50.0)

      if((fabs(dx)>5.0) and (fabs(dy)>5.0) and (fabs(xe[gem0])>50.0) and (fabs(ye[gem0])>50.0))
      {
        H2(X3,Y3, Form("Efficiency/%s/No_good_hits_gem_%s", gemdet,gemdet), Form("No good clusters on %s GEM after projecting track", gemdet),50,-50.0,50.0,50,-50.0,50.0);
          continue;
      };
        

       H2(xe[gem0],ye[gem0], Form("Efficiency/%s/ActualGoodClusters_GEM_%s", gemdet,gemdet), Form("Actual Good clusters on %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);
       

       H2(xe[gem0],Plane2crds[0], Form("Efficiency/%s/BHPlane2_corr_X_GEM_%s", gemdet,gemdet), Form("Correlation of BH Plane 2 X  Vs %s GEM X", gemdet),50,-50.0,50.0,50,-50.0,50.0);
       
       H2(xe[gem0],Plane2crds[1], Form("Efficiency/%s/BHPlane2_corr_XY_GEM_%s", gemdet,gemdet), Form("Correlation of BH Plane 2 Y and %s GEM X", gemdet),50,-50.0,50.0,50,-50.0,50.0);
       
       H2(xe[gem0],Plane3crds[0], Form("Efficiency/%s/BHPlane3_corr_X_GEM_%s", gemdet,gemdet), Form("Correlation of BH Plane 3 X and %s GEM X", gemdet),50,-50.0,50.0,50,-50.0,50.0);
       
       H2(xe[gem0],Plane3crds[1], Form("Efficiency/%s/BHPlane3_corr_XY_GEM_%s", gemdet,gemdet), Form("Correlation of BH Plane 3 Y and %s GEM X", gemdet),50,-50.0,50.0,50,-50.0,50.0);
       

       H2(ye[gem0],Plane2crds[0], Form("Efficiency/%s/BHPlane2_corr_YX_GEM_%s", gemdet,gemdet), Form("Correlation of BH Plane 2 X and %s GEM Y", gemdet),50,-50.0,50.0,50,-50.0,50.0);
       
       H2(ye[gem0],Plane2crds[1], Form("Efficiency/%s/BHPlane2_corr_Y_GEM_%s", gemdet,gemdet), Form("Correlation of BH Plane 2 Y  and %s GEM Y", gemdet),50,-50.0,50.0,50,-50.0,50.0);
       
       H2(ye[gem0],Plane3crds[0], Form("Efficiency/%s/BHPlane3_corr_YX_GEM_%s", gemdet,gemdet), Form("Correlation of BH Plane 3 X and %s GEM Y", gemdet),50,-50.0,50.0,50,-50.0,50.0);
       
       H2(ye[gem0],Plane3crds[1], Form("Efficiency/%s/BHPlane3_corr_Y_GEM_%s", gemdet,gemdet), Form("Correlation of BH Plane 3 Y and %s GEM Y", gemdet),50,-50.0,50.0,50,-50.0,50.0);
       


       H2(X3,Y3, Form("Efficiency/%s/ActualProjectedGoodClusters_GEM_%s", gemdet,gemdet), Form("Actual Projected Good clusters on %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);
      
       H1(dx, Form("Efficiency/%s/x_residue%s", gemdet, gemdet), Form("Verticle residue on %s GEM", gemdet),200,-50.0,50.0);
                          
       H1(dy, Form("Efficiency/%s/y_residue%s", gemdet, gemdet), Form("Horizontal residue on %s GEM", gemdet),200,-50.0,50.0);

      
      



     }; // Loop over the GEM in interest
  }; //BH cluster filtering finish

  return 0;
  //  return Plugin::ok;
};

Long_t MUSEteleTracker::finalize_efficiency()
{
    
    
    auto effUS=dH2(TString::Format("Efficiency/%s/Efficiency_Plot_GEM_%s","US","US"),"Efficency Plot US GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto didUS=dH2(TString::Format("Efficiency/%s/ActualProjectedGoodClusters_GEM_%s","US","US"),"Actual Projected Good clusters on US GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto shouldUS=dH2(TString::Format("Efficiency/%s/tracksprojectedgem_active_%s","US","US"),"Distribution of Tracks Projected inside active area on US GEM",50,-50.0,50.0,50,-50.0,50.0);
    effUS->Divide(didUS,shouldUS);

    auto eff4TH=dH2(TString::Format("Efficiency/%s/Efficiency_Plot_GEM_%s","4TH","4TH"),"Efficency Plot 4TH GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto did4TH=dH2(TString::Format("Efficiency/%s/ActualProjectedGoodClusters_GEM_%s","4TH","4TH"),"Actual Projected Good clusters on 4TH GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto should4TH=dH2(TString::Format("Efficiency/%s/tracksprojectedgem_active_%s","4TH","4TH"),"Distribution of Tracks Projected inside active area on 4TH GEM",50,-50.0,50.0,50,-50.0,50.0);
    eff4TH->Divide(did4TH,should4TH);

    auto effMS=dH2(TString::Format("Efficiency/%s/Efficiency_Plot_GEM_%s","MS","MS"),"Efficency Plot MS GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto didMS=dH2(TString::Format("Efficiency/%s/ActualProjectedGoodClusters_GEM_%s","MS","MS"),"Actual Projected Good clusters on MS GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto shouldMS=dH2(TString::Format("Efficiency/%s/tracksprojectedgem_active_%s","MS","MS"),"Distribution of Tracks Projected inside active area on MS GEM",50,-50.0,50.0,50,-50.0,50.0);
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

