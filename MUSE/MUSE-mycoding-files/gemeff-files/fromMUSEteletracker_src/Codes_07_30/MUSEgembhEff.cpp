/* GEM efficiency (integrated with BH) Version by Ishara Fernando */
#include <MUSEteleTracker.h>
#include<iostream>
#include<cmath>
#include "GEMhittree.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1.h"
#include "TH2.h"


#include<iostream>
#include<cmath>
#include<vector>
#include<numeric>
#include<algorithm>
#include<array>

#include <fstream>
#include <sstream>
#include <string>
#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"






Long_t MUSEteleTracker::startup_GEMBH_efficiency()
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
    printf(" BH_Hits (Hits) @%p\n", Hits);
    
    
    // get input branch for BH
    Hits=NULL;
    getBranchObject("BH_Hits", (TObject**)&Hits);
    if (Hits==NULL)
    {
        printf(" Cannot find branch >Scinthits< in input ROOT file - trying output branch\n");
        getOutBranchObject("BH_Hits",(TObject**)&Hits);
        if(Hits==NULL)
        {
            printf("Couldn't find any Hits in any ROOT file :(\n");
            return -1;
        }
    };
    printf(" BH_Hits (Hits) @%p\n", Hits);

    
    return Plugin::ok;
    
    
    ///////////////////////////////////////////
}
Long_t MUSEteleTracker::histos_GEMBH_efficiency()
{
      
    return 0;
    
}

Long_t MUSEteleTracker::process_GEMBH_efficiency()
{
    int NG=GEM_NUM;  // GEM_NUM or Number of GEMs

    //const double zgempos[3]= {-536.0,-474.0,-412.0};
    std::vector<double> zgempos = {0, 0, 0};

	  for(unsigned int k=0; k<3; k++){
		  double loc[3]={0,0,0};
		  double mas[3]={0,0,0};
		  
		  handle.GEM[k]->LocalToMaster(loc,mas);

		  zgempos[k] = mas[2]*10;
		  
	  }
    const char* gemdet;
    const char* gemdets;

    
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
        //double RF = Hits->hits[i].rf;
        //int plane, side, bar;
        //BH_internal_to_logic(Hits->hits[i].id,&plane,&bar,&side);
        int plane = Hits->hits[i].wall_id;
        int bar = Hits->hits[i].bar_id;
        
        double loc[3] = {0,0,0};
        double mas[3] = {0,0,0};
        handle.BH[plane][bar]->LocalToMaster(loc,mas);
        
      
        if(plane==2)
        { 
 	        Plane2=plane;
          Plane2bar = bar;
          Plane2crds[0] = mas[0]*10;
          Plane2crds[1] = mas[1]*10;
          Plane2crds[2] = mas[2]*10;
          H1(Plane2crds[0], Form("GEMBHEfficiency/BH/Plane2/BH_Plane2_x_positions"), Form("Cluster x positions of BH Plane 2"),50,-50.0,50.0);
          H1(Plane2crds[1], Form("GEMBHEfficiency/BH/Plane2/BH_Plane2_y_positions"), Form("Cluster y positions of BH Plane 2"),50,-50.0,50.0);
          H1(Plane2crds[2], Form("GEMBHEfficiency/BH/Plane2/BH_Plane2_z_positions"), Form("Cluster z positions of BH Plane 2"),50,-650.0,-600.0);
          H1(bar, Form("GEMBHEfficiency/BH/Plane2/BH_Plane2_bar_number"), Form("Bar number of BH Plane 2"),16,-0.5,16.5);
          BHPlaneCount[plane]++;
        }; 

       if(plane==3)
        { 
	        Plane3=plane;
          Plane3bar = bar;
          Plane3crds[0] = mas[0]*10;
          Plane3crds[1] = mas[1]*10;
          Plane3crds[2] = mas[2]*10;
           H1(Plane3crds[0], Form("GEMBHEfficiency/BH/Plane3/BH_Plane3_x_positions"), Form("Cluster x positions of BH Plane 3"),50,-50.0,50.0);
           H1(Plane3crds[1], Form("GEMBHEfficiency/BH/Plane3/BH_Plane3_y_positions"), Form("Cluster y positions of BH Plane 3"),50,-50.0,50.0);
           H1(Plane3crds[2], Form("GEMBHEfficiency/BH/Plane3/BH_Plane3_z_positions"), Form("Cluster z positions of BH Plane 3"),50,-650.0,600.0);
           H1(bar, Form("GEMBHEfficiency/BH/Plane3/BH_Plane3_bar_number"), Form("Bar number of BH Plane 3"),15,-0.5,14.5);
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



    




    /*
    //This For loop collect the number of clusters on each GEM and x,y coordinates of max_charge cluster
    int gemmulti3gems[3] = { 0, 0, 0 }; 
    double cluster_gem_charge[20];
    //double xmaxchg[GEM_NUM]={-10000};
    //double ymaxchg[GEM_NUM]={-10000};
    for (int gems=0; gems<GEM_NUM; gems++)
    {
        double maxcharge=0;
        bool test=false;
        //cluster_gem_charge[gemanyhit[gems]]=0;
        for (unsigned int g1=0; g1<clusters->hits.size(); g1++)
        {
            if (clusters->hits[g1].GEMid!=gems) continue;
            gemmulti3gems[gems]++;
            cluster_gem_charge[gemmulti3gems[gems]]=clusters->hits[g1].quality;
      //printf("cluster charge is %d\n",cluster_gem_charge[gemanyhit[gems]]);
            if(cluster_gem_charge[gemmulti3gems[gems]]>maxcharge) 
            {
                //printf("found max cluster charge\n");
                maxcharge =cluster_gem_charge[gemmulti3gems[gems]];
                test=true;
                
                xmaxchg[gems] = clusters->hits[g1].xl*0.4-50.;
                ymaxchg[gems] = clusters->hits[g1].yl*0.4-50.;


            }
            else test=false;
        };
         
        
    };
    */
    



    

   // Multiplicities in BH planes
  H1(BHPlaneCount[2], Form("GEMBHEfficiency/BH/Plane2/BH_Plane2_Cluster_Multiplicity"), Form("Cluster Multiplicity of BH Plane 2"),11,-0.5,10.5);
  H1(BHPlaneCount[3], Form("GEMBHEfficiency/BH/Plane3/BH_Plane3_Cluster_Multiplicity"), Form("Cluster Multiplicity of BH Plane 3"),11,-0.5,10.5);

  //Multiplicities in GEM planes
  H1(gemmultibefore[0], Form("GEMBHEfficiency/US/Cluster_Multiplicity_US"), Form("Cluster Multiplicity of US GEM before BH cluster filtering"),11,-0.5,10.5);
  H1(gemmultibefore[1], Form("GEMBHEfficiency/4TH/Cluster_Multiplicity_4TH"), Form("Cluster Multiplicity of 4TH GEM before BH cluster filtering"),11,-0.5,10.5);
  H1(gemmultibefore[2], Form("GEMBHEfficiency/MS/Cluster_Multiplicity_MS"), Form("Cluster Multiplicity of MS GEM before BH cluster filtering"),11,-0.5,10.5);

  
    
  // BH 1 hit filetring filtering start

  //bool need_BH_filtering=true;

  if(BHPlaneCount[2]>0 and BHPlaneCount[3]>0)
  {


    int BHPlaneCountAC[4]={0,0,0,0};
    int Plane2AC=0;
    int Plane2barAC=0;
    double Plane2crdsAC[3]={0.0,0.0,0.0};
    int Plane3AC=0;
    int Plane3barAC=0;
    double Plane3crdsAC[3]={0.0,0.0,0.0};
    for(size_t i = 0; i < Hits->hits.size(); i++)
    {
        //int plane, side, bar;       
        //BH_internal_to_logic(Hits->hits[i].id,&plane,&bar,&side);
        int plane = Hits->hits[i].wall_id;
        int bar = Hits->hits[i].bar_id;
        
        double loc[3] = {0,0,0};
        double mas[3] = {0,0,0};
        handle.BH[plane][bar]->LocalToMaster(loc,mas);
        

        if(plane==2)
        { 
          Plane2AC=plane;
          Plane2barAC = bar;
          Plane2crdsAC[0] = mas[0]*10;
          Plane2crdsAC[1] = mas[1]*10;
          Plane2crdsAC[2] = mas[2]*10;
          H1(Plane2crdsAC[0], Form("GEMBHEfficiency/BH/Plane2/BH_Plane2_1hit_x_positions"), Form("Cluster x positions of BH Plane 2 1hit"),50,-50.0,50.0);
          H1(Plane2crdsAC[1], Form("GEMBHEfficiency/BH/Plane2/BH_Plane2_1hit_y_positions"), Form("Cluster y positions of BH Plane 2 1hit"),50,-50.0,50.0);
          H1(Plane2crdsAC[2], Form("GEMBHEfficiency/BH/Plane2/BH_Plane2_1hit_z_positions"), Form("Cluster z positions of BH Plane 2 1hit"),50,-650.0,-600.0);
          H1(Plane2barAC, Form("GEMBHEfficiency/BH/Plane2/BH_Plane2_1hit_bar_number"), Form("Bar number of BH Plane 2 1hit"),16,-0.5,16.5);
          BHPlaneCountAC[plane]++;
          
        }; 

       if(plane==3)
        { 
          Plane3AC=plane;
          Plane3barAC = bar;
          Plane3crdsAC[0] = mas[0]*10;
          Plane3crdsAC[1] = mas[1]*10;
          Plane3crdsAC[2] = mas[2]*10;
          H1(Plane3crdsAC[0], Form("GEMBHEfficiency/BH/Plane3/BH_Plane3_1hit_x_positions"), Form("Cluster x positions of BH Plane 3 1hit"),50,-50.0,50.0);
          H1(Plane3crdsAC[1], Form("GEMBHEfficiency/BH/Plane3/BH_Plane3_1hit_y_positions"), Form("Cluster y positions of BH Plane 3 1hit"),50,-50.0,50.0);
          H1(Plane3crdsAC[2], Form("GEMBHEfficiency/BH/Plane3/BH_Plane3_1hit_z_positions"), Form("Cluster z positions of BH Plane 3 1hit"),50,-650.0,-600.0);
          H1(Plane3barAC, Form("GEMBHEfficiency/BH/Plane3/BH_Plane3_1hit_bar_number"), Form("Bar number of BH Plane 3 1hit"),15,-0.5,14.5);
          BHPlaneCountAC[plane]++;
          
        };      
       
    };

    //printf("(%d,%d)\n",BHPlaneCountAC[2],BHPlaneCountAC[3]);

    // BH multiplicities to cross check wehter there are only 1 hit per plane
    H1(BHPlaneCountAC[2], Form("GEMBHEfficiency/BH/Plane2/BH_Plane2_1Hitluster_Multiplicity"), Form("Cluster Multiplicity of BH Plane 2 with 1 hit"),11,-0.5,10.5);
    H1(BHPlaneCountAC[3], Form("GEMBHEfficiency/BH/Plane3/BH_Plane3_1Hit_Cluster_Multiplicity"), Form("Cluster Multiplicity of BH Plane 3 with 1 hit"),11,-0.5,10.5);



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
    
     // GEM multiplicities after BH 1-hit cut
    H1(gemmultiafter[0], Form("GEMBHEfficiency/US/BH_filtered_Cluster_Multiplicity_US"), Form("Cluster Multiplicity of US GEM after BH cluster filtering"),11,-0.5,10.5);
    H1(gemmultiafter[1], Form("GEMBHEfficiency/4TH//BH_filteredCluster_Multiplicity_4TH"), Form("Cluster Multiplicity of 4TH GEM after BH cluster filtering"),11,-0.5,10.5);
    H1(gemmultiafter[2], Form("GEMBHEfficiency/MS//BH_filteredCluster_Multiplicity_MS"), Form("Cluster Multiplicity of MS GEM after BH cluster filtering"),11,-0.5,10.5);


    // Loop over the GEM in interest
    int gemmultiaftertrack[3] = { 0, 0, 0 };
    for (int gem0=0; gem0<NG; gem0++)  
    {

      if (gem0==0){gemdet="US";};
      if (gem0==1){gemdet="4TH";};
      if (gem0==2){gemdet="MS";};

      // Loop over the other two GEMs
      bool GoToNextGEM0=false;
      int gem1hits[3]={0,0,0};
      double xg[2]={0.0000,0.0000};
      double yg[2]={0.0000,0.0000};
      double zg[2]={0.0000,0.0000};
      int gemID[2]={0,0};
      int gem1count=0; //count over GEM1
      double slopexg=0.0;
      double slopeyg=0.0;
      for (int gem1=0; ((gem1<NG) and (GoToNextGEM0==false) and (gem1count<2)); gem1++)
      {
        if (gem1==gem0) continue; // Not considering the GEM
        //printf("GEM1 IDs = (%d,%d)\n",gem0,gem1);
        for (unsigned int g10=0; g10<clusters->hits.size(); g10++)
        {
         if ((clusters->hits[g10].GEMid!=gem1))  continue;
          gem1hits[gem1]++;  // total # of hits on the 2nd and 3rd GEM
          
        }; // Loop to get number of clusters in each GEM among 2 GEMs
          //printf("Clusters on GEM %d is = %d\n",gem1,gem1hits[gem1]);
        if (gem1hits[gem1]==0) 
        {
          GoToNextGEM0=true;
          continue; 
        };
        // This condition will eleminate the cases where there are no good clusters found for any GEM1
        //printf("Clusters on GEM %d is = %d\n",gem1,gem1hits[gem1]);

        double cluster_charge_gem1=0.0;
        bool good_cluster_gem1=false;
        double max_charge_gem1=0.0;
        int gem1_cluster_number=-1;
        // Looking for a good cluster (loop over clusters)          
        for (unsigned int g1=0; g1<clusters->hits.size(); g1++)
        {
          if ((clusters->hits[g1].GEMid!=gem1)) continue;
          cluster_charge_gem1=clusters->hits[g1].charge;
          //cluster_charge_gem1=clusters->hits[g1].quality;
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
          GoToNextGEM0=true;
          continue;
        };

        double tempx = clusters->hits[gem1_cluster_number].xl;
        double tempy = clusters->hits[gem1_cluster_number].yl;
        double tempz = zgempos[gem1];

        //printf("cluster position is = (%5.2f,%5.2f,%5.2f) \n",tempx,tempy,tempz);

        //Checking whether the cluster is in the active area if the GEM plane
        if ((fabs(tempx)>50.0) and (fabs(tempy)>50.0))
        {
          GoToNextGEM0=true;
          continue;
        };


        xg[gem1count]=tempx;
        yg[gem1count]=tempy;
        zg[gem1count]=tempz;
        gemID[gem1count]=gem1;

        gem1count=gem1count+1;

      }; // Loop closes for the otsher two GEMs

      // The following condition will skip the GEM0 to next one if we don't get two good clsuters
      //if (gemID[0]==gemID[1]) continue;
      if (gemID[0]!=gemID[1])
      {

      //slopexg=(xg[1]-xg[0])/(zgempos[gemID[1]]-zgempos[gemID[0]]);
      //slopeyg=(yg[1]-yg[0])/(zgempos[gemID[1]]-zgempos[gemID[0]]);

      slopexg=(xg[1]-xg[0])/(zg[1]-zg[0]);
      slopeyg=(yg[1]-yg[0])/(zg[1]-zg[0]);

     
      for (unsigned int g1=0; g1<clusters->hits.size(); g1++)
      {
        // If the GEMid is different than the cluster's gem ID, then g1 gets incremented by 1 unit
        if (clusters->hits[g1].GEMid!=gem0) continue;
        gemmultiaftertrack[gem0]++;
      };

      // GEM multiplicities after generating tracks (before projecting it to anywhere)
      H1(gemmultiaftertrack[gem0], Form("GEMBHEfficiency/%s/Cluster_Multiplicity_%s_After_tracks_by_other2",gemdet,gemdet), Form("Cluster Multiplicity of %s GEM After tracks by other2",gemdet),11,-0.5,10.5);


      //Project the tracks on BH planes
      // Notet that plane 2 has only y coordinates varying (x is fixed to 5mm): Horizontal ,where as plane 3 is verticle
                    
      //double BHplane2X = xg[1] + slopexg*(Plane2crds[2]-zgempos[gemID[1]]);
      //double BHplane2Y = yg[1] + slopeyg*(Plane2crds[2]-zgempos[gemID[1]]);

      //double BHplane3X = xg[1] + slopexg*(Plane3crds[2]-zgempos[gemID[1]]);
      //double BHplane3Y = yg[1] + slopeyg*(Plane3crds[2]-zgempos[gemID[1]]); 

      double BHplane2X = xg[1] + slopexg*(Plane2crds[2]-zg[1]);
      double BHplane2Y = yg[1] + slopeyg*(Plane2crds[2]-zg[1]);

      double BHplane3X = xg[1] + slopexg*(Plane3crds[2]-zg[1]);
      double BHplane3Y = yg[1] + slopeyg*(Plane3crds[2]-zg[1]); 



      H1(BHplane2Y, Form("GEMBHEfficiency/%s/BH_plane2Y_projected_%s", gemdet, gemdet), Form("BH_plane2Y_projected %s GEM", gemdet),14,-50.0,50.0);
      H1(Plane2crdsAC[1], Form("GEMBHEfficiency/%s/BH_plane_2Y_actual_%s", gemdet, gemdet), Form("BH_plane_2Y_actual %s GEM", gemdet),14,-50.0,50.0); 


      for(size_t i = 0; i < Hits->hits.size(); i++)
      {

        //int plane, side, bar; 
        //BH_internal_to_logic(Hits->hits[i].id,&plane,&bar,&side);
        int plane = Hits->hits[i].wall_id;
        int bar = Hits->hits[i].bar_id;
        
        double loc[3] = {0,0,0};
        double mas[3] = {0,0,0};
        handle.BH[plane][bar]->LocalToMaster(loc,mas);
        

        if(plane==2)
        { 
          Plane2AC=plane;
          Plane2barAC = bar;
          H1(Plane2barAC, Form("GEMBHEfficiency/BH/Plane2/BH_Plane2_1hit_bar_after_projecting_track"), Form("Bar number of BH Plane 2 1hit after projecting track"),16,-0.5,16.5);
          
          
        }; 

       if(plane==3)
        { 
          Plane3AC=plane;
          Plane3barAC = bar;
          H1(Plane3barAC, Form("GEMBHEfficiency/BH/Plane3/BH_Plane3_1hit_bar_after_projecting_track"), Form("Bar number of BH Plane 3 1hit after projecting track"),15,-0.5,14.5);

        };
          
      
      };
       

       //printf("BH plane 2 projected %5.4f \n",BHplane2Y);
       //printf("BH plane 2 actual %5.4f \n",Plane2crds[1]);

      //Note that there should be a negative sign for x y for GEM coordinates since there is a sign flip with BH (MUSE) coordinate frame
      //double BHdx=Plane3crds[0]-(-BHplane3X)-11.000;
      // negative sign is removed because the coorrdinate frames are corrected
      double BHdx=Plane3crds[0]-(BHplane3X)-0.000;
      double BHdy=Plane2crds[1]-(BHplane2Y)-0.000;

      // Coordinate frames between BH and GEMs are corrected by Michael in the version migrated to github.com from git@mit.com

      //BH correlations
      //Her we use a negative sign because x and y should be flipped to get the MUSE coordinate frame
       H2(BHplane2X,Plane2crds[0], Form("GEMBHEfficiency/%s/BHPlane2_corr_X%s", gemdet,gemdet), Form("Correlation of BH Plane 2: ActualX vs ProjectedX : %s GEM", gemdet),50,-50.0,50.0,13,-50.0,50.0);
       
       H2(BHplane2X,Plane2crds[1], Form("GEMBHEfficiency/%s/BHPlane2_corr_XY%s", gemdet,gemdet), Form("Correlation of BH Plane 2 : ActualY vs ProjectedX: %s GEM", gemdet),50,-50.0,50.0,13,-50.0,50.0);

       H2(BHplane2Y,Plane2crds[0], Form("GEMBHEfficiency/%s/BHPlane2_corr_YX%s", gemdet,gemdet), Form("Correlation of BH Plane 2 X : ActualX vs ProjectedY: %s GEM", gemdet),50,-50.0,50.0,13,-50.0,50.0);
       
       H2(BHplane2Y,Plane2crds[1], Form("GEMBHEfficiency/%s/BHPlane2_corr_Y%s", gemdet,gemdet), Form("Correlation of BH Plane 2 Y  : ActualY vs ProjectedY: %s GEM", gemdet),50,-50.0,50.0,13,-50.0,50.0);
  

       
       H2(BHplane3X,Plane3crds[0], Form("GEMBHEfficiency/%s/BHPlane3_corr_X%s", gemdet,gemdet), Form("Correlation of BH Plane 3 X : ActualX vs ProjectedX: %s GEM", gemdet),50,-50.0,50.0,16,-50.0,50.0);
       
       H2(BHplane3X,Plane3crds[1], Form("GEMBHEfficiency/%s/BHPlane3_corr_XY%s", gemdet,gemdet), Form("Correlation of BH Plane 3 Y : ActualY vs ProjectedX: %s GEM", gemdet),50,-50.0,50.0,16,-50.0,50.0);
            
       H2(BHplane3Y,Plane3crds[0], Form("GEMBHEfficiency/%s/BHPlane3_corr_YX%s", gemdet,gemdet), Form("Correlation of BH Plane 3 X : ActualX vs ProjectedY: %s GEM", gemdet),50,-50.0,50.0,16,-50.0,50.0);
       
       H2(BHplane3Y,Plane3crds[1], Form("GEMBHEfficiency/%s/BHPlane3_corr_Y%s", gemdet,gemdet), Form("Correlation of BH Plane 3 Y : ActualY vs ProjectedY vs: %s GEM", gemdet),50,-50.0,50.0,16,-50.0,50.0);
       

      H1(BHdx, Form("GEMBHEfficiency/%s/BH_residuals_of_gem_%s_X", gemdet, gemdet), Form("Residual X-Distribution of Tracks Projected on BH to evaluate %s GEM ", gemdet),14,-50.0,50.0);
                          
      H1(BHdy, Form("GEMBHEfficiency/%s/BH_residuals_of_gem_%s_Y", gemdet, gemdet), Form("Residual Y-distribution of Tracks Projected on BH to evaluate %s GEM ", gemdet),14,-50.0,50.0);


          H1(Plane2crds[1], Form("GEMBHEfficiency/BH/Plane2/BH_Plane2_y_positions_After_tracks"), Form("Cluster y positions of BH Plane 2 after projecting tracks"),50,-50.0,50.0);
          H1(Plane3crds[0], Form("GEMBHEfficiency/BH/Plane3/BH_Plane2_x_positions_After_tracks"), Form("Cluster x positions of BH Plane 3 after projecting tracks"),50,-50.0,50.0);



        if((fabs(BHdx)<10.0) and (fabs(BHdy)<10.0)) 
        {
          H1(BHdx, Form("GEMBHEfficiency/%s/BH_residuals_of_gem_%s_X_after_cut", gemdet, gemdet), Form("Residual X-Distribution of Tracks Projected on BH after cut to evaluate %s GEM ", gemdet),50,-50.0,50.0);
                          
          H1(BHdy, Form("GEMBHEfficiency/%s/BH_residuals_of_gem_%s_Y_after_cut", gemdet, gemdet), Form("Residual Y-distribution of Tracks Projected on BH after cut to evaluate %s GEM ", gemdet),50,-50.0,50.0);


          

          //Project the tracks on GEM0
          double X3 = xg[1] + slopexg*(zgempos[gem0]-zgempos[gemID[1]]);
          double Y3 = yg[1] + slopeyg*(zgempos[gem0]-zgempos[gemID[1]]);


          H2(X3,Y3, Form("GEMBHEfficiency/%s/tracksprojectedgem_%s", gemdet,gemdet), Form("Distribution of Tracks Projected on %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);

          if((fabs(X3)>50.0) and (fabs(Y3)>50.0)) continue;

          H2(X3,Y3, Form("GEMBHEfficiency/%s/tracksprojectedgem_active_%s", gemdet,gemdet), Form("Distribution of Tracks Projected inside active area on %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);

          int gem0hits[3]={0,0,0};
          for (unsigned int g0=0; g0<clusters->hits.size(); g0++)
          {
            if ((clusters->hits[g0].GEMid!=gem0)) continue;
            gem0hits[gem0]++;  // total number of hits on the GEM0

          };      

         if(gem0hits[gem0]==0)
          {
            H2(X3,Y3, Form("GEMBHEfficiency/%s/No_hits_gem_%s", gemdet,gemdet), Form("No clusters on %s GEM after projecting track", gemdet),50,-50.0,50.0,50,-50.0,50.0);
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
            cluster_charge_gem0=clusters->hits[g0].charge;
            //cluster_charge_gem0=clusters->hits[g0].quality;
            if(cluster_charge_gem0>max_charge_gem0)
            {
              good_cluster_gem0=true;
              max_charge_gem0=cluster_charge_gem0;
              gem0_cluster_number=g0;
              //printf("found the max charge cluseter on GEM0 \n");
            };
          };
      

          if(!good_cluster_gem0) continue;

           xe[gem0] = clusters->hits[gem0_cluster_number].xl;
           ye[gem0] = clusters->hits[gem0_cluster_number].yl;

           dx=xe[gem0]-X3; //vertical residue 
           dy=ye[gem0]-Y3; //hosrizontal residue


            //if((fabs(dx)>5.0) or (fabs(dy)>5.0) or (fabs(xe[gem0])>50.0) or (fabs(ye[gem0])>50.0))
            //{
            //  H2(X3,Y3, Form("GEMBHEfficiency/%s/No_good_hits_gem_%s", gemdet,gemdet), Form("No good clusters on %s GEM after projecting track", gemdet),50,-50.0,50.0,50,-50.0,50.0);
            //  continue;
            //};
            //(fabs(dx)<10.0) and (fabs(dy)<10.0) and 
            if((fabs(dx)<10.0) and (fabs(dy)<10.0) and(fabs(xe[gem0])<50.0) and (fabs(ye[gem0])<50.0))
          {
            H2(xe[gem0],ye[gem0], Form("GEMBHEfficiency/%s/ActualGoodClusters_GEM_%s", gemdet,gemdet), Form("Actual Good clusters on %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);

            H2(X3,Y3, Form("GEMBHEfficiency/%s/ActualProjectedGoodClusters_GEM_%s", gemdet,gemdet), Form("Actual Projected Good clusters on %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);
      
            H1(dx, Form("GEMBHEfficiency/%s/x_residue%s", gemdet, gemdet), Form("X residue on %s GEM", gemdet),200,-50.0,50.0);
                          
            H1(dy, Form("GEMBHEfficiency/%s/y_residue%s", gemdet, gemdet), Form("Y residue on %s GEM", gemdet),200,-50.0,50.0);
           };

           H2(X3,Y3, Form("GEMBHEfficiency/%s/No_good_hits_gem_%s", gemdet,gemdet), Form("No good clusters on %s GEM after projecting track", gemdet),50,-50.0,50.0,50,-50.0,50.0);
            
         }; // Bh threshold cut  if condition

      }; // If condition for valid track from other two GEMs

     }; // Loop closes for the GEM in interest

     // GEM multiplicities after generating tracks (before projecting it to anywhere)
      //H1(gemmultiaftertrack[0], Form("GEMBHEfficiency/US/Cluster_Multiplicity_US_After_tracks_by_other2"), Form("Cluster Multiplicity of US GEM After tracks by other2"),11,-0.5,10.5);
      //H1(gemmultiaftertrack[1], Form("GEMBHEfficiency/4TH//Cluster_Multiplicity_4TH_After_tracks_by_other2"), Form("Cluster Multiplicity of 4TH GEM After tracks by other2"),11,-0.5,10.5);
      //H1(gemmultiaftertrack[2], Form("GEMBHEfficiency/MS//Cluster_Multiplicity_MS_After_tracks_by_other2"), Form("Cluster Multiplicity of MS GEM After tracks by other2"),11,-0.5,10.5);   
    
   }; //BH 1-hit filtering finish

  return 0;
  //  return Plugin::ok;
};

Long_t MUSEteleTracker::finalize_GEMBH_efficiency()
{
    
    
    auto effUS=dH2(TString::Format("GEMBHEfficiency/%s/Efficiency_Plot_GEM_%s","US","US"),"Efficency Plot US GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto didUS=dH2(TString::Format("GEMBHEfficiency/%s/ActualProjectedGoodClusters_GEM_%s","US","US"),"Actual Projected Good clusters on US GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto shouldUS=dH2(TString::Format("GEMBHEfficiency/%s/tracksprojectedgem_active_%s","US","US"),"Distribution of Tracks Projected inside active area on US GEM",50,-50.0,50.0,50,-50.0,50.0);
    effUS->Divide(didUS,shouldUS);

    auto eff4TH=dH2(TString::Format("GEMBHEfficiency/%s/Efficiency_Plot_GEM_%s","4TH","4TH"),"Efficency Plot 4TH GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto did4TH=dH2(TString::Format("GEMBHEfficiency/%s/ActualProjectedGoodClusters_GEM_%s","4TH","4TH"),"Actual Projected Good clusters on 4TH GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto should4TH=dH2(TString::Format("GEMBHEfficiency/%s/tracksprojectedgem_active_%s","4TH","4TH"),"Distribution of Tracks Projected inside active area on 4TH GEM",50,-50.0,50.0,50,-50.0,50.0);
    eff4TH->Divide(did4TH,should4TH);

    auto effMS=dH2(TString::Format("GEMBHEfficiency/%s/Efficiency_Plot_GEM_%s","MS","MS"),"Efficency Plot MS GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto didMS=dH2(TString::Format("GEMBHEfficiency/%s/ActualProjectedGoodClusters_GEM_%s","MS","MS"),"Actual Projected Good clusters on MS GEM",50,-50.0,50.0,50,-50.0,50.0);
    auto shouldMS=dH2(TString::Format("GEMBHEfficiency/%s/tracksprojectedgem_active_%s","MS","MS"),"Distribution of Tracks Projected inside active area on MS GEM",50,-50.0,50.0,50,-50.0,50.0);
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

