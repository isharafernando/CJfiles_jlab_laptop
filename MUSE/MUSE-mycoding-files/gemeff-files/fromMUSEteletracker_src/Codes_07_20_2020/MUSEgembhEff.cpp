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
        printf(" Cannot find branch >BH_Hits< in input ROOT file - trying output branch\n");
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
Long_t MUSEteleTracker::histos_efficiency()
{
      
    return 0;
    
}

Long_t MUSEteleTracker::process_efficiency()
{
    int NG=GEM_NUM;  // GEM_NUM or Number of GEMs

    const double zgempos[3]= {-536.0,-474.0,-412.0};
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
        int plane, side, bar;
        
        BH_internal_to_logic(Hits->hits[i].id,&plane,&bar,&side);
        
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
          H1(Plane2crds[0], Form("Efficiency/BH/Plane2/BH_Plane2_x_positions"), Form("Cluster x positions of BH Plane 2"),50,-50.0,50.0);
          H1(Plane2crds[1], Form("Efficiency/BH/Plane2/BH_Plane2_y_positions"), Form("Cluster y positions of BH Plane 2"),50,-50.0,50.0);
          H1(Plane2crds[2], Form("Efficiency/BH/Plane2/BH_Plane2_z_positions"), Form("Cluster z positions of BH Plane 2"),50,-650.0,-600.0);
          H1(bar, Form("Efficiency/BH/Plane2/BH_Plane2_bar_number"), Form("Bar number of BH Plane 2"),16,-0.5,16.5);
          BHPlaneCount[plane]++;
        }; 

       if(plane==3)
        { 
	        Plane3=plane;
          Plane3bar = bar;
          Plane3crds[0] = mas[0]*10;
          Plane3crds[1] = mas[1]*10;
          Plane3crds[2] = mas[2]*10;
           H1(Plane3crds[0], Form("Efficiency/BH/Plane3/BH_Plane3_x_positions"), Form("Cluster x positions of BH Plane 3"),50,-50.0,50.0);
           H1(Plane3crds[1], Form("Efficiency/BH/Plane3/BH_Plane3_y_positions"), Form("Cluster y positions of BH Plane 3"),50,-50.0,50.0);
           H1(Plane3crds[2], Form("Efficiency/BH/Plane3/BH_Plane3_z_positions"), Form("Cluster z positions of BH Plane 3"),50,-650.0,600.0);
           H1(bar, Form("Efficiency/BH/Plane3/BH_Plane3_bar_number"), Form("Bar number of BH Plane 3"),15,-0.5,14.5);
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
  H1(BHPlaneCount[2], Form("Efficiency/BH/Plane2/BH_Plane2_Cluster_Multiplicity"), Form("Cluster Multiplicity of BH Plane 2"),11,-0.5,10.5);
  H1(BHPlaneCount[3], Form("Efficiency/BH/Plane3/BH_Plane3_Cluster_Multiplicity"), Form("Cluster Multiplicity of BH Plane 3"),11,-0.5,10.5);

  //Multiplicities in GEM planes
  H1(gemmultibefore[0], Form("Efficiency/US/Cluster_Multiplicity_US"), Form("Cluster Multiplicity of US GEM before BH cluster filtering"),11,-0.5,10.5);
  H1(gemmultibefore[1], Form("Efficiency/4TH/Cluster_Multiplicity_4TH"), Form("Cluster Multiplicity of 4TH GEM before BH cluster filtering"),11,-0.5,10.5);
  H1(gemmultibefore[2], Form("Efficiency/MS/Cluster_Multiplicity_MS"), Form("Cluster Multiplicity of MS GEM before BH cluster filtering"),11,-0.5,10.5);

  
    
  // BH 1 hit filetring filtering start

  //bool need_BH_filtering=true;

  if(BHPlaneCount[2]==1 and BHPlaneCount[3]==1)
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
        int plane, side, bar;
        
        BH_internal_to_logic(Hits->hits[i].id,&plane,&bar,&side);
        
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
          H1(Plane2crdsAC[0], Form("Efficiency/BH/Plane2/BH_Plane2_1hit_x_positions"), Form("Cluster x positions of BH Plane 2 1hit"),50,-50.0,50.0);
          H1(Plane2crdsAC[1], Form("Efficiency/BH/Plane2/BH_Plane2_1hit_y_positions"), Form("Cluster y positions of BH Plane 2 1hit"),50,-50.0,50.0);
          H1(Plane2crdsAC[2], Form("Efficiency/BH/Plane2/BH_Plane2_1hit_z_positions"), Form("Cluster z positions of BH Plane 2 1hit"),50,-650.0,-600.0);
          H1(Plane2barAC, Form("Efficiency/BH/Plane2/BH_Plane2_1hit_bar_number"), Form("Bar number of BH Plane 2 1hit"),16,-0.5,16.5);
          BHPlaneCountAC[plane]++;
          
        }; 

       if(plane==3)
        { 
          Plane3AC=plane;
          Plane3barAC = bar;
          Plane3crdsAC[0] = mas[0]*10;
          Plane3crdsAC[1] = mas[1]*10;
          Plane3crdsAC[2] = mas[2]*10;
          H1(Plane3crdsAC[0], Form("Efficiency/BH/Plane3/BH_Plane3_1hit_x_positions"), Form("Cluster x positions of BH Plane 3 1hit"),50,-50.0,50.0);
          H1(Plane3crdsAC[1], Form("Efficiency/BH/Plane3/BH_Plane3_1hit_y_positions"), Form("Cluster y positions of BH Plane 3 1hit"),50,-50.0,50.0);
          H1(Plane3crdsAC[2], Form("Efficiency/BH/Plane3/BH_Plane3_1hit_z_positions"), Form("Cluster z positions of BH Plane 3 1hit"),50,-650.0,-600.0);
          H1(Plane3barAC, Form("Efficiency/BH/Plane3/BH_Plane3_1hit_bar_number"), Form("Bar number of BH Plane 3 1hit"),15,-0.5,14.5);
          BHPlaneCountAC[plane]++;
          
        };      
       
    };

    //printf("(%d,%d)\n",BHPlaneCountAC[2],BHPlaneCountAC[3]);

    // BH multiplicities to cross check wehter there are only 1 hit per plane
    H1(BHPlaneCountAC[2], Form("Efficiency/BH/Plane2/BH_Plane2_1Hitluster_Multiplicity"), Form("Cluster Multiplicity of BH Plane 2 with 1 hit"),11,-0.5,10.5);
    H1(BHPlaneCountAC[3], Form("Efficiency/BH/Plane3/BH_Plane3_1Hit_Cluster_Multiplicity"), Form("Cluster Multiplicity of BH Plane 3 with 1 hit"),11,-0.5,10.5);



    /* Cluster Multiplicity after BH filtering 1 */

    int gemmultiafter[3] = { 0, 0, 0 };
    for (unsigned int g1=0; g1<clusters->hits.size(); g1++)
    {
	    // bool gem1hit=true;
      	int gemhittot[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	  	for (unsigned int g1=0; g1<clusters->hits.size(); g1++) //get the # of clusters on first GEM in each event
	  	{
	      //		if ((clusters->hits[g1].GEMid!=gem) && ((clusters->hits[g1].ampl<400)))  continue;
	  		if ((clusters->hits[g1].GEMid!=gem))  continue;
	      	gemhittot[gem]++;  // total # of hits on the first gem
	  	};
	  //  if (gemhittot[gem]>1) gem1hit=false;
	  	if (gemhittot[gem]==0) continue; // If nothing on the first gem, then loop back continues to search gem hits on the second gem
	   //  if the first GEM has some hits, then gemhittot[gem]>0 and continue the following////////////////////
	  	double cluster_gem1_charge[gemhittot[gem]];
	  	int gemhittot1omore[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	  
	  //xxxxxxxx if ((gemhittot[gem]==1)&&(gem1hit)) { //If the first GEM has only one cluster
	    if ((gemhittot[gem]==1)) 
	    { //If the first GEM has only one cluster
	       //if (gemhittot[gem]>=1) { // If the first GEM has any number of clusters, clusters >=1 

	     	bool found_gem1_gdcluster=false;
	     	double maxcharge_gem1=0;
	     	int clustnmgem=-1;
			//printf("pass  gem1 %d %d\n",gem,gemhittot[gem] );	 
	     	for (unsigned int g1=0; g1<clusters->hits.size(); g1++)
	     	{
		    //  if ((clusters->hits[g1].GEMid!=gem)&& ((clusters->hits[g1].ampl<400))) continue;
	     		if ((clusters->hits[g1].GEMid!=gem)) continue;
	     		gemhittot1omore[gem]++;
	     		cluster_gem1_charge[gemhittot1omore[gem]]=clusters->hits[g1].charge;

		 		if (cluster_gem1_charge[gemhittot1omore[gem]]>maxcharge_gem1)///select the cluster with the maximum charge on the first GEM
		   		{
				  	found_gem1_gdcluster=true;
				  	maxcharge_gem1=cluster_gem1_charge[gemhittot1omore[gem]];
				  	clustnmgem=g1;
				}
			};

	      	if (found_gem1_gdcluster) xe[1] = clusters->hits[clustnmgem].xl*0.4-50.; ye[1] = clusters->hits[clustnmgem].yl*0.4-50.; // get x/y coordinates for the max. charge cluster on the first GEM at loop 1
	      	//loop 1 absolutely does not end here
	      	//Ethan
	      //////////////************loop 1 ends and loop 2 starts ***********************//////////////
	      	if ((found_gem1_gdcluster)&&((fabs(xe[1])<40) && (fabs(ye[1])<40)))
	      	{
		  		for (int gemt=gem+1; gemt<10; gemt++) // loop 2 start looping on gem 2 if a "good" cluster found on gem 1
		  		{ 
		  			int gem1hittot[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, };
		      		for (unsigned int g2=0; g2<clusters->hits.size(); g2++) //get the # of clusters on the second GEM if the first GEM have selected the maximum cluster
		      		{ 
			  //	if ((clusters->hits[g2].GEMid!=gemt)&& ((clusters->hits[g2].ampl<400))) continue;	
			  			if ((clusters->hits[g2].GEMid!=gemt)) continue; //Get the total # of clusters on GEM 2
			  			gem1hittot[gemt]++;  // total number of clusters on GEM 2
					};
		      //std::cout << "gem1 hit tot " << gem1hittot[gemt] << std::endl;
		      		if (gem1hittot[gemt]==0) continue; // If nothing on the second gem, then the loop 2 continues to search gem clusters on the third GEM. Still we are inside loop 1 with "good" GEM clusters on the GEM 1
		    ///////////////////////////////
		            double cluster_gem2_charge[gem1hittot[gemt]];
			        int gem1hittot1omore[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
				    bool found_gem2_gdcluster=false;
				      //if (gem1hittot[gemt]>=1) { // clusters >=1 on GEM 2 //If the second GEM has any number of clusters, cluster >=1
				    if (gem1hittot[gemt]==1) 
				    {// clusters ==1 on GEM 2 //If the second GEM has only one cluster
					//	{
					   	double maxcharge_gemt=0;
				      	int clustnmgemt=-1;
					//std::cout << "pass gem 2 " << gemt << " " << gem1hittot[gemt] << std::endl; 
				      	for (unsigned int g2=0; g2<clusters->hits.size(); g2++)
				      	{
				      		if ((clusters->hits[g2].GEMid!=gemt)) continue;
				      		gem1hittot1omore[gemt]++;
				      		cluster_gem2_charge[gem1hittot1omore[gemt]]=clusters->hits[g2].charge;
						    if (cluster_gem2_charge[gem1hittot1omore[gemt]]>maxcharge_gemt)///select the cluster with the maximum charge
						    {
						    	found_gem2_gdcluster=true;
						    	maxcharge_gemt=cluster_gem2_charge[gem1hittot1omore[gemt]];
						    	clustnmgemt=g2;
						    }
						};
						if (found_gem2_gdcluster)
						{ 
							xe[2] = clusters->hits[clustnmgemt].xl*0.4-50.; 
							ye[2] = clusters->hits[clustnmgemt].yl*0.4-50.; // get x/y coordinates for the max. charge cluster on the second GEM at loop 2
						}
			////////////////////////////////************loop 2 ends and loop 3 starts ***********************///////////////////
			
						if ((found_gem2_gdcluster)&&((fabs(xe[2])<40) && (fabs(ye[2])<40)))
						{
			    			for (int gem3=gemt+1; gem3<10; gem3++) // loop 3 start looping on gem 3 if a "good" cluster found on gem 2
			    			{ 
			    				int gem3hittot[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
								for (unsigned int g3=0; g3<clusters->hits.size(); g3++) //get the # of clusters on the  GEM 3 if the first 2 GEMs  have selected the maximum cluster
								{ 
								    //	if ((clusters->hits[g2].GEMid!=gemt)&& ((clusters->hits[g2].ampl<400))) continue;	
								    if ((clusters->hits[g3].GEMid!=gem3)) continue; //Get the total # of clusters on GEM 2
								    gem3hittot[gem3]++;  // total number of clusters on GEM 3
								};
								if (gem3hittot[gem3]==0) continue; // If nothing on the third gem, then the loop 3 continues to search gem clusters on the other GEMs on the same loop (GEM 3 and GEM 4). Still we are inside loop 1 and 2 with "good" GEM clusters on the GEM 1 and GEM 2
								///////////////////////////////
								double cluster_gem3_charge[gem3hittot[gem3]];
								int gem3hittot1omore[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
								bool found_gem3_gdcluster=false;
								//if (gem3hittot[gem3]>=1) { // clusters >=1 on GEM 3 //If the third GEM has any number of clusters, cluster >=1
								if (gem3hittot[gem3]==1) 
								{// clusters ==1 on GEM 3 //If the third GEM has only one cluster
					  //	{
									double maxcharge_gem3=0;
									int clustnmgem3=-1;
								  //std::cout << "pass gem3 " << gem3 << " " << gem3hittot[gem3] << std::endl; 
									for (unsigned int g3=0; g3<clusters->hits.size(); g3++)
									{
										if ((clusters->hits[g3].GEMid!=gem3)) continue;
										gem3hittot1omore[gem3]++;
										cluster_gem3_charge[gem3hittot1omore[gem3]]=clusters->hits[g3].charge;
								      if (cluster_gem3_charge[gem3hittot1omore[gem3]]>maxcharge_gem3)///select the cluster with the maximum charge
								      {
								      	found_gem3_gdcluster=true;
								      	maxcharge_gem3=cluster_gem3_charge[gem3hittot1omore[gem3]];
								      	clustnmgem3=g3;
								      }
								    };
								  	if (found_gem3_gdcluster)
								  	{ 
								  		xe[3] = clusters->hits[clustnmgem3].xl*0.4-50.; 
								  		ye[3] = clusters->hits[clustnmgem3].yl*0.4-50.; // get x/y coordinates for the max. charge cluster on the third GEM at loop 3
				  					}
									//  if ((found_gem3_gdcluster)&&((fabs(xe[3])<40) && (fabs(ye[3])<40))) {
									////////////////////////////////************loop 3 ends and loop 4 starts ***********************///////////////////
									// if ((gem1hittot[gemt]>0))

									if ((found_gem3_gdcluster)&&((fabs(xe[3])<40) && (fabs(ye[3])<40)))  // If the maximum charge cluster found on the third GEM at loop 3
									{
										double slopexe = (xe[3]-xe[1])/(zgempos[gem3]-zgempos[gem]);
										double slopeye = (ye[3]-ye[1])/(zgempos[gem3]-zgempos[gem]);

										if ((gem==6)&&(gemt==7)&&(gem3==8))
										{
									  //std::cout << "DS" << std::endl;
											zgemd=7;
											gemdet="DS";
										};
										if ((gem==6)&&(gemt==8)&&(gem3==9))
										{
										//std::cout << "4th" << std::endl;
											zgemd=5;
											gemdet="4TH";
										};
										if ((gem==6)&&(gemt==7)&&(gem3==9))
										{
										//std::cout << "MS" << std::endl;
											zgemd=6;
											gemdet="MS";
										};
										if ((gem==7)&&(gemt==8)&&(gem3==9))
										{
									  //std::cout << "US" << std::endl;
											zgemd=4;
											gemdet="US";
										};
									      //Project tracks on the fourth GEM (if all  both first, second and third GEMs have only maximum charge cluster selected)
									  	double X3 = xe[1] + slopexe*(zgempos[zgemd]-zgempos[gem]);
									  	double Y3 = ye[1] + slopeye*(zgempos[zgemd]-zgempos[gem]); 
									    if ((fabs(X3)<40) && (fabs(Y3)<40))  // this includes the clusters >=1 as well as "no hits at all" on the third GEM
									    {
									    	H2(X3,Y3, Form("Efficiency/%s/tracksprojectedgem%s_any", gemdet,gemdet), Form("Tracks Projected on %s GEM-Any cluster", gemdet),50,-50.0,50.0,50,-50.0,50.0);

									      	H1(X3, Form("Efficiency/%s/tracksprojectedgem%sX_any", gemdet, gemdet), Form("Tracks Projected on %s GEM X-Any cluster", gemdet),50,-50.0,50.0);
									      	H1(Y3, Form("Efficiency/%s/tracksprojectedgem%sY_any", gemdet, gemdet), Form("Tracks Projected on %s GEM Y-Any cluster", gemdet),50,-50.0,50.0);
									    };
				      
				      
										int gem4hittot[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
										bool gem4hit=true;
										if (anyhit) // look for any number of hits on the fourth GEM
										{
											for (unsigned int g4=0; g4<clusters->hits.size(); g4++)
											{
										  //	if ((clusters->hits[g3].GEMid!=zgemd)&& ((clusters->hits[g3].ampl<400))) continue;
												if (clusters->hits[g4].GEMid!=zgemd) continue;
										  		gem4hittot[zgemd]++;// search clusters on fourth gem (for g4)
										  
												// printf("how many hits %d %d %d %d %d %d\n",gem,gemhittot[gem],gemt,gem1hittot[gemt],zgemd,gem2hittot[zgemd]);
												xe1[4] = clusters->hits[g4].xl*0.4-50.; ye1[4] = clusters->hits[g4].yl*0.4-50.;
												dx1=xe1[4]-X3; // vertical residual
												dy1=ye1[4]-Y3; // horizontal residual					   
												if ((fabs(xe1[4])<40) && (fabs(ye1[4])<40))
												{
													H1(dx1, Form("Efficiency/%s/xresiduagem%s", gemdet, gemdet), Form("Vert. Residua on %s GEM (clusters >=1)", gemdet),100, -50.0, 50.0);
													H1(dy1, Form("Efficiency/%s/yresiduagem%s", gemdet, gemdet), Form("Hori. Residua on %s GEM (clusters >=1)", gemdet),100, -50.0, 50.0);
													H2(xe1[4],ye1[4], Form("Efficiency/%s/tracksdetectedgem%s_any", gemdet, gemdet), Form("Cluster positions on %s GEM-Any cluster (clusters >=1)", gemdet),50,-50.0,50.0,50,-50.0,50.0);
													/////////use for mapping check ///////////
													H2(xe1[4],dx1, Form("Efficiency/%s/vert_res_vs_vert_detec_clust%s_any", gemdet, gemdet), Form("Vert. residual Vs Vert. cluster positions on %s GEM-Any cluster (clusters >=1)", gemdet),50,-50.0,50.0,50,-25.0,25.0);
													H2(xe1[4],dy1, Form("Efficiency/%s/hor_res_vs_vert_detec_clust%s_any", gemdet, gemdet), Form("Hor. residual Vs vert. cluster positions on %s GEM-Any cluster (clusters >=1)", gemdet),50,-50.0,50.0,50,-25.0,25.0);

													H2(ye1[4],dy1, Form("Efficiency/%s/hor_res_vs_hor_detec_clust%s_any", gemdet, gemdet), Form("Hor. residual Vs hor. cluster positions on %s GEM-Any cluster (clusters >=1)", gemdet),50,-50.0,50.0,50,-25.0,25.0);

													H2(ye1[4],dx1, Form("Efficiency/%s/vert_res_vs_hor_detec_clust%s_any", gemdet, gemdet), Form("Vert. residual Vs hor. cluster positions on %s GEM-Any cluster (clusters >=1)", gemdet),50,-50.0,50.0,50,-25.0,25.0);

													H2(X3,dx1, Form("Efficiency/%s/vert_res_vs_vert_proj_clust%s_any", gemdet, gemdet), Form("Vert. residual Vs vert. proj. cluster positions on %s GEM-Any cluster (clusters >=1)", gemdet),50,-50.0,50.0,50,-25.0,25.0);
													H2(X3,dy1, Form("Efficiency/%s/hor_res_vs_vert_proj_clust%s_any", gemdet, gemdet), Form("Hor. residual Vs vert. proj. cluster positions on %s GEM-Any cluster (clusters >=1)", gemdet),50,-50.0,50.0,50,-25.0,25.0);

													H2(Y3,dy1, Form("Efficiency/%s/hor_res_vs_hor_proj_clust%s_any", gemdet, gemdet), Form("Hor. residual Vs hor. proj. cluster positions on %s GEM-Any cluster (clusters >=1)", gemdet),50,-50.0,50.0,50,-25.0,25.0);

													H2(Y3,dx1, Form("Efficiency/%s/vert_res_vs_hor_proj_clust%s_any", gemdet, gemdet), Form("Vert. residual Vs hor. proj. cluster positions on %s GEM-Any cluster (clusters >=1)", gemdet),50,-50.0,50.0,50,-25.0,25.0);
												/////////////////////////////

					      						};
					    					};	// g4 loop for 	anyhit		
					  					};//anyhit

				    ////////////////////

					  					H1(gem4hittot[zgemd], Form("Efficiency/%s/multiplicity%s_any", gemdet, gemdet), Form("Multiplicity - Any cluster at %s GEM", gemdet),11,-0.5,10.5);

			  ////////////////
									   	if (gem4hittot[zgemd]==0) // no clusters at all on the 4th GEM
									   	{				
									   		if ((fabs(X3)<40) && (fabs(Y3)<40)) H2(X3,Y3, Form("Efficiency/%s/tracksprojectedgem%s_nohit", gemdet, gemdet), Form("No Clusters On %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);
									   	};
				   ////////////////


									   double dxe1omore[gem4hittot[zgemd]],dye1omore[gem4hittot[zgemd]];
									   int gem4hittot1omore[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

									   double cluster_gem4_charge[gem4hittot[zgemd]];
									   bool found_gem4_gdcluster=false;
								       /////////////////////
								       /*
								       for (unsigned int g4=0; g4<clusters->hits.size(); g4++)
									 {
										      //	if ((clusters->hits[g3].GEMid!=zgemd)&& ((clusters->hits[g3].ampl<400))) continue;
									   if (clusters->hits[g4].GEMid!=zgemd) continue;
									   gem4hittot[zgemd]++;// search clusters on fourth gem (for g4)
					*/
									   /////////////

									    if (gem4hittot[zgemd]>=1)// clusters >=1
									    {
											double maxcharge_zgemd=0;
											int clustnmzgemd=-1;
											for (unsigned int g4=0; g4<clusters->hits.size(); g4++)
											{
												if ((clusters->hits[g4].GEMid!=zgemd)) continue;
												gem4hittot1omore[zgemd]++;
												xe[4] = clusters->hits[g4].xl*0.4-50.; ye[4] = clusters->hits[g4].yl*0.4-50.;
												dxe1omore[gem4hittot1omore[zgemd]]=xe[4]-X3;
												dye1omore[gem4hittot1omore[zgemd]]=ye[4]-Y3;

											////////////////////////////////////////////////////////////////////
												cluster_gem4_charge[gem4hittot1omore[zgemd]]=clusters->hits[g4].charge;
												if (cluster_gem4_charge[gem4hittot1omore[zgemd]]>maxcharge_zgemd)///select the cluster with the maximum charge
												{
													found_gem4_gdcluster=true;
													maxcharge_zgemd=cluster_gem4_charge[gem4hittot1omore[zgemd]];
													clustnmzgemd=g4;
												};
											////////////////////////////////////////////////////////////////////////				    			 
											};
				  
											if (found_gem4_gdcluster)
											{
												xe_max[4] = clusters->hits[clustnmzgemd].xl*0.4-50.; 
												ye_max[4] = clusters->hits[clustnmzgemd].yl*0.4-50.;
												dx_max[4] =xe_max[4]-X3;
												dy_max[4] =ye_max[4]-Y3;

											  //  if ((fabs(dx_max[3])<=10) && (fabs(dy_max[3])<=10))
											  //	{
											  //printf("pass gem3 %d %d %d %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf\n",gem, gemt,zgemd,ye[1],xe[1],ye[2],xe[2],xe_max[3],ye_max[3] );
											  // H2(ye[1],xe[1], Form("cluster_hit_map_USGEM"), Form("cluster hit map US GEM"),50,-50.0,50.0,50,-50.0,50.0);
											  //H2(ye[2],xe[2], Form("cluster_hit_map_MSGEM"), Form("cluster hit map MS GEM"),50,-50.0,50.0,50,-50.0,50.0);
											  // H2(ye_max[3],xe_max[3], Form("cluster_hit_map_DSGEM"), Form("cluster hit map DS GEM"),50,-50.0,50.0,50,-50.0,50.0);	    

											  //	}
											};
				  							bool good=false;
				  							for (int clust=0; clust<gem4hittot[zgemd]; clust++)
											{
												if (((fabs(dxe1omore[clust])<=10) && (fabs(dye1omore[clust])<=10))) 
												{ 
													good=true;
												//	tracks[zgemd]++;
												} 
												else good=false;      
											};

	  										//  printf("pass  gem1 %d %d \n",zgemd,tracks[zgemd]);

											if ((good) && ((fabs(X3)<40) && (fabs(Y3)<40)))
											{
											// tracks[zgemd]++;
												H2(X3,Y3, Form("Efficiency/%s/tracksprojectedgem%s_oneormore_cut", gemdet, gemdet), Form("Tracks Projected on %s GEM--cluster >=1, Inside Vicinity", gemdet),50,-50.0,50.0,50,-50.0,50.0);
											//   H1(X3, Form("tracksprojectedgem%sX_oneormore_cut", gemdet), Form("Tracks Projected on %s GEM X-Inside Vicinity", gemdet),50,-50.0,50.0);
											//   H1(Y3, Form("tracksprojectedgem%sY_oneormore_cut", gemdet), Form("Tracks Projected on %s GEM X-Inside Vicinity", gemdet),50,-50.0,50.0);	
											};

											//	if ((fabs(X3)<50) && (fabs(Y3)<50))
											//	  {
											//	    H2(Y3,X3, Form("tracksprojectedgem%s_onecut", gemdet), Form("Tracks Projected on %s GEM-Inside Vicinity", gemdet),50,-50.0,50.0,50,-50.0,50.0);
											//  };

											if ((good) && ((fabs(xe[4])<40) && (fabs(ye[4])<40)))
											{
												H2(xe[4],ye[4], Form("Efficiency/%s/tracksdetectedgem%s_oneormore_cut", gemdet, gemdet), Form("Cluster positions on %s GEM--cluster >=1, Inside Vicinity", gemdet),50,-50.0,50.0,50,-50.0,50.0);
											};

											if ((!good) && ((fabs(xe[4])<40) && (fabs(ye[4])<40)))
											{

												H2(xe[4],ye[4], Form("Efficiency/%s/tracksdetectedgem%s_oneormore_cutout", gemdet, gemdet), Form("Cluster positions on %s GEM--cluster >=1, Outside Vicinity", gemdet),50,-50.0,50.0,50,-50.0,50.0);
											//   printf("pass  gem1 %d %5.2lf %5.2lf \n",zgemd,ye[3],xe[3] );
											};

											if ((fabs(X3)<40) && (fabs(Y3)<40)) H2(X3,Y3, Form("Efficiency/%s/tracksprojectedgem%s_oneormore", gemdet, gemdet), Form("Tracks Projected on %s GEM-cluster >=1", gemdet),50,-50.0,50.0,50,-50.0,50.0);
											if ((!good) && ((fabs(X3)<40) && (fabs(Y3)<40))) H2(X3,Y3, Form("Efficiency/%s/tracksprojectedgem%s_oneormore_cutout", gemdet, gemdet), Form("Tracks Projected on %s GEM-cluster >=1, Outside Vicinity", gemdet),50,-50.0,50.0,50,-50.0,50.0);

										};//cluster >=1

			      ////////// Check if gem 4 has found only 1 cluster  /////

				     					if (gem4hittot[zgemd]>1) gem4hit=false;
										if ((gem4hittot[zgemd]!=0)&&(gem4hit)) //clusters =1
										{
											for (unsigned int g4=0; g4<clusters->hits.size(); g4++)
											{
										  // if ((clusters->hits[g3].GEMid!=zgemd)&& ((clusters->hits[g3].ampl<400))) continue;
												if ((clusters->hits[g4].GEMid!=zgemd)) continue;
										  //	    printf("get x/y  gem3 %d %d\n",zgemd,gem2hittot[zgemd]);
												xe[4] = clusters->hits[g4].xl*0.4-50.; ye[4] = clusters->hits[g4].yl*0.4-50.;
												dxe=xe[4]-X3;
												dye=ye[4]-Y3;
												if ((fabs(X3)<40) && (fabs(Y3)<40)) H2(X3,Y3, Form("Efficiency/%s/tracksprojectedgem%s_one", gemdet, gemdet), Form("Tracks Projected on %s GEM-One cluster", gemdet),50,-50.0,50.0,50,-50.0,50.0);
												if ((fabs(xe[4])<40) && (fabs(ye[4])<40))
												{
													H1(dxe, Form("Efficiency/%s/xresiduagem%s_one", gemdet, gemdet), Form("Vert. Residua on %s GEM", gemdet),100, -50.0, 50.0);
													H1(dye, Form("Efficiency/%s/yresiduagem%s_one", gemdet, gemdet), Form("Hori. Residua on %s GEM", gemdet),100, -50.0, 50.0);
													H2(xe[4],ye[4], Form("Efficiency/%s/tracksdetectedgem%s_one", gemdet, gemdet), Form("Cluster positions on %s GEM-One cluster", gemdet),50,-50.0,50.0,50,-50.0,50.0);

												}; // inside the gem fabs(x and y)<40

										      //	    if (((fabs(dxe)<=10) && (fabs(dye)<=10))&& ((clusters->hits[g3].ampl>=400) || (clusters->hits[g3].ampl<=1500)))
												if (((fabs(dxe)<=10) && (fabs(dye)<=10)))
												{
													if ((fabs(xe[4])<40) && (fabs(ye[4])<40)) H2(xe[3],ye[3], Form("Efficiency/%s/tracksdetectedgem%s_onecut", gemdet, gemdet), Form("Cluster positions on %s GEM-Inside Vicinity", gemdet),50,-50.0,50.0,50,-50.0,50.0);
													if ((fabs(X3)<40) && (fabs(Y3)<40))
													{
														H2(X3,Y3, Form("Efficiency/%s/tracksprojectedgem%s_onecut", gemdet, gemdet), Form("Tracks Projected on %s GEM-Inside Vicinity", gemdet),50,-50.0,50.0,50,-50.0,50.0);
														H1(X3, Form("Efficiency/%s/tracksprojectedgem%sX_onecut", gemdet,gemdet), Form("Tracks Projected on %s GEM X-Inside Vicinity", gemdet),50,-50.0,50.0);
														H1(Y3, Form("Efficiency/%s/tracksprojectedgem%sY_onecut", gemdet,gemdet), Form("Tracks Projected on %s GEM X-Inside Vicinity", gemdet),50,-50.0,50.0);

													};
												};
												if ((fabs(dxe)>10) || (fabs(dye)>10))
												{
													if ((fabs(xe[4])<40) && (fabs(ye[4])<40)) H2(xe[4],ye[4], Form("Efficiency/%s/tracksdetectedgem%s_onecutout", gemdet, gemdet), Form("Cluster positions on %s GEM-Outside Vicinity", gemdet),50,-50.0,50.0,50,-50.0,50.0);
													if ((fabs(Y3)<50) && (fabs(X3)<50)) H2(X3,Y3, Form("Efficiency/%s/tracksprojectedgem%s_onecutout", gemdet, gemdet), Form("Tracks Projected on %s GEM-Outside Vicinity", gemdet),50,-50.0,50.0,50,-50.0,50.0);
												};
				    						}; // cluster =1, g4 loop
				  						};// if gem3hit , clusters =1
			      // }; // g4 loop
			    					};//if (found_gem3_gdcluster=true) a maximum charge cluster found in second GEM. THese are inside loop 4


			  // if((gem==3)&&(gemt==4))mult2dUS_MS->Fill(gemhittot[gem],gem1hittot[gemt]);

			  //  mult2dUS_MS->Fill(gemhittot[0],gem1hittot[1]);//test
			  //  mult2dUS_DS->Fill(gemhittot[0],gem1hittot[2]);
			  // mult2dMS_DS->Fill(gemhittot[1],gem1hittot[2]);

								};// // clusters >=1 on the third  GEM , loop 3  
			      			}; //start looping gem 3, loop 3
		  //printf ("'''''test %d %d %d %d \n",gem,gemhittot[gem],gemt, gem1hittot[gemt]); 
		  // if ((gem==3)&&(gemt==4)) mult2dUS_MS->Fill(gemhittot[gem],gem1hittot[gemt]);  
		  //  mult2dUS_MS->Fill(gemhittot[3],gem1hittot[4]);
		  //mult2dUS_DS->Fill(gemhittot[3],gem1hittot[5]);
		  //mult2dMS_DS->Fill(gemhittot[4],gem1hittot[5]);
			  			};//if (found_gem2_gdcluster=true) a maximum charge cluster found in the second GEM
		      		};/// clusters >=1 on the second  GEM   
		    	}; //start looping gem 2, loop 2
			};//	if (found_gem1_gdcluster=true) a maximum charge cluster found in the first GEM
	    };// // clusters >=1 on the first GEM   
	}; //start looping gem 1, loop 1


      //I just did
      ////////////////////////Track MUltiplicity ////////////
      //   int gemanyhit[6] = { 0, 0, 0, 0, 0, 0 };
      //   int tracks[6] = { 0, 0, 0, 0, 0, 0 };
      // double minchi2fit=1000;
     
      for (unsigned int g1=0; g1<clusters->hits.size(); g1++)
      {
        // If the GEMid is different than the cluster's gem ID, then g1 gets incremented by 1 unit
        if (clusters->hits[g1].GEMid!=gem0) continue;
        gemmultiaftertrack[gem0]++;
      };

      // GEM multiplicities after generating tracks (before projecting it to anywhere)
      H1(gemmultiaftertrack[gem0], Form("Efficiency/%s/Cluster_Multiplicity_%s_After_tracks_by_other2",gemdet,gemdet), Form("Cluster Multiplicity of %s GEM After tracks by other2",gemdet),11,-0.5,10.5);


      //Project the tracks on BH planes
                    
      double BHplane2X = xg[1] + slopexg*(Plane2crds[2]-zgempos[gemID[1]]);
      double BHplane2Y = yg[1] + slopeyg*(Plane2crds[2]-zgempos[gemID[1]]);

      double BHplane3X = xg[1] + slopexg*(Plane3crds[2]-zgempos[gemID[1]]);
      double BHplane3Y = yg[1] + slopeyg*(Plane3crds[2]-zgempos[gemID[1]]); 


      H1(-BHplane2Y, Form("Efficiency/%s/BH_plane2Y_projected_%s", gemdet, gemdet), Form("BH_plane2Y_projected %s GEM", gemdet),14,-50.0,50.0);
      H1(Plane2crds[1], Form("Efficiency/%s/BH_plane_2Y_actual_%s", gemdet, gemdet), Form("BH_plane_2Y_actual %s GEM", gemdet),14,-50.0,50.0); 


      for(size_t i = 0; i < Hits->hits.size(); i++)
      {

        int plane, side, bar;
        
        BH_internal_to_logic(Hits->hits[i].id,&plane,&bar,&side);
        
        double loc[3] = {0,0,0};
        double mas[3] = {0,0,0};
        handle.BH[plane][bar]->LocalToMaster(loc,mas);
        

        if(plane==2)
        { 
          Plane2AC=plane;
          Plane2barAC = bar;
          H1(Plane2barAC, Form("Efficiency/BH/Plane2/BH_Plane2_1hit_bar_after_projecting_track"), Form("Bar number of BH Plane 2 1hit after projecting track"),16,-0.5,16.5);
          
          
        }; 

       if(plane==3)
        { 
          Plane3AC=plane;
          Plane3barAC = bar;
          H1(Plane3barAC, Form("Efficiency/BH/Plane3/BH_Plane3_1hit_bar_after_projecting_track"), Form("Bar number of BH Plane 3 1hit after projecting track"),15,-0.5,14.5);

        };
          
      
      };
       

       //printf("BH plane 2 projected %5.4f \n",BHplane2Y);
       //printf("BH plane 2 actual %5.4f \n",Plane2crds[1]);

      //Note that there should be a negative sign for x y for GEM coordinates since there is a sign flip with BH (MUSE) coordinate frame
      double BHdx=Plane3crds[0]-(-BHplane3X)-11.000;
      double BHdy=Plane2crds[1]-(-BHplane2Y)-0.000;

      //BH correlations
      //Her we use a negative sign because x and y should be flipped to get the MUSE coordinate frame
       H2(-BHplane2X,Plane2crds[0], Form("Efficiency/%s/BHPlane2_corr_X%s", gemdet,gemdet), Form("Correlation of BH Plane 2: ActualX vs ProjectedX : %s GEM", gemdet),50,-50.0,50.0,13,-50.0,50.0);
       
       H2(-BHplane2X,Plane2crds[1], Form("Efficiency/%s/BHPlane2_corr_XY%s", gemdet,gemdet), Form("Correlation of BH Plane 2 : ActualY vs ProjectedX: %s GEM", gemdet),50,-50.0,50.0,13,-50.0,50.0);

       H2(-BHplane2Y,Plane2crds[0], Form("Efficiency/%s/BHPlane2_corr_YX%s", gemdet,gemdet), Form("Correlation of BH Plane 2 X : ActualX vs ProjectedY: %s GEM", gemdet),50,-50.0,50.0,13,-50.0,50.0);
       
       H2(-BHplane2Y,Plane2crds[1], Form("Efficiency/%s/BHPlane2_corr_Y%s", gemdet,gemdet), Form("Correlation of BH Plane 2 Y  : ActualY vs ProjectedY: %s GEM", gemdet),50,-50.0,50.0,13,-50.0,50.0);
  

       
       H2(-BHplane3X,Plane3crds[0], Form("Efficiency/%s/BHPlane3_corr_X%s", gemdet,gemdet), Form("Correlation of BH Plane 3 X : ActualX vs ProjectedX: %s GEM", gemdet),50,-50.0,50.0,16,-50.0,50.0);
       
       H2(-BHplane3X,Plane3crds[1], Form("Efficiency/%s/BHPlane3_corr_XY%s", gemdet,gemdet), Form("Correlation of BH Plane 3 Y : ActualY vs ProjectedX: %s GEM", gemdet),50,-50.0,50.0,16,-50.0,50.0);
            
       H2(-BHplane3Y,Plane3crds[0], Form("Efficiency/%s/BHPlane3_corr_YX%s", gemdet,gemdet), Form("Correlation of BH Plane 3 X : ActualX vs ProjectedY: %s GEM", gemdet),50,-50.0,50.0,16,-50.0,50.0);
       
       H2(-BHplane3Y,Plane3crds[1], Form("Efficiency/%s/BHPlane3_corr_Y%s", gemdet,gemdet), Form("Correlation of BH Plane 3 Y : ActualY vs ProjectedY vs: %s GEM", gemdet),50,-50.0,50.0,16,-50.0,50.0);
       

      H1(BHdx, Form("Efficiency/%s/BH_residuals_of_gem_%s_X", gemdet, gemdet), Form("Residual X-Distribution of Tracks Projected on BH to evaluate %s GEM ", gemdet),14,-50.0,50.0);
                          
      H1(BHdy, Form("Efficiency/%s/BH_residuals_of_gem_%s_Y", gemdet, gemdet), Form("Residual Y-distribution of Tracks Projected on BH to evaluate %s GEM ", gemdet),14,-50.0,50.0);


          H1(Plane2crds[1], Form("Efficiency/BH/Plane2/BH_Plane2_y_positions_After_tracks"), Form("Cluster y positions of BH Plane 2 after projecting tracks"),50,-50.0,50.0);
          H1(Plane3crds[0], Form("Efficiency/BH/Plane3/BH_Plane2_x_positions_After_tracks"), Form("Cluster x positions of BH Plane 3 after projecting tracks"),50,-50.0,50.0);



        if((fabs(BHdx)<8.0) and (fabs(BHdy)<8.0)) 
        {
          H1(BHdx, Form("Efficiency/%s/BH_residuals_of_gem_%s_X_after_cut", gemdet, gemdet), Form("Residual X-Distribution of Tracks Projected on BH after cut to evaluate %s GEM ", gemdet),50,-50.0,50.0);
                          
          H1(BHdy, Form("Efficiency/%s/BH_residuals_of_gem_%s_Y_after_cut", gemdet, gemdet), Form("Residual Y-distribution of Tracks Projected on BH after cut to evaluate %s GEM ", gemdet),50,-50.0,50.0);


          

          //Project the tracks on GEM0
          double X3 = xg[1] + slopexg*(zgempos[gem0]-zgempos[gemID[1]]);
          double Y3 = yg[1] + slopeyg*(zgempos[gem0]-zgempos[gemID[1]]);


          H2(X3,Y3, Form("Efficiency/%s/tracksprojectedgem_%s", gemdet,gemdet), Form("Distribution of Tracks Projected on %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);

          if((fabs(X3)>50.0) or (fabs(Y3)>50.0)) continue;

          H2(X3,Y3, Form("Efficiency/%s/tracksprojectedgem_active_%s", gemdet,gemdet), Form("Distribution of Tracks Projected inside active area on %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);

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


            //if((fabs(dx)>5.0) or (fabs(dy)>5.0) or (fabs(xe[gem0])>50.0) or (fabs(ye[gem0])>50.0))
            //{
            //  H2(X3,Y3, Form("Efficiency/%s/No_good_hits_gem_%s", gemdet,gemdet), Form("No good clusters on %s GEM after projecting track", gemdet),50,-50.0,50.0,50,-50.0,50.0);
            //  continue;
            //};
            //(fabs(dx)<10.0) and (fabs(dy)<10.0) and 
            if((fabs(xe[gem0])<50.0) and (fabs(ye[gem0])<50.0))
          {
            H2(xe[gem0],ye[gem0], Form("Efficiency/%s/ActualGoodClusters_GEM_%s", gemdet,gemdet), Form("Actual Good clusters on %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);

            H2(X3,Y3, Form("Efficiency/%s/ActualProjectedGoodClusters_GEM_%s", gemdet,gemdet), Form("Actual Projected Good clusters on %s GEM", gemdet),50,-50.0,50.0,50,-50.0,50.0);
      
            H1(dx, Form("Efficiency/%s/x_residue%s", gemdet, gemdet), Form("X residue on %s GEM", gemdet),200,-50.0,50.0);
                          
            H1(dy, Form("Efficiency/%s/y_residue%s", gemdet, gemdet), Form("Y residue on %s GEM", gemdet),200,-50.0,50.0);
           };

           H2(X3,Y3, Form("Efficiency/%s/No_good_hits_gem_%s", gemdet,gemdet), Form("No good clusters on %s GEM after projecting track", gemdet),50,-50.0,50.0,50,-50.0,50.0);
            
         }; // Bh threshold cut  if condition

      }; // If condition for valid track from other two GEMs

     }; // Loop closes for the GEM in interest

     // GEM multiplicities after generating tracks (before projecting it to anywhere)
      //H1(gemmultiaftertrack[0], Form("Efficiency/US/Cluster_Multiplicity_US_After_tracks_by_other2"), Form("Cluster Multiplicity of US GEM After tracks by other2"),11,-0.5,10.5);
      //H1(gemmultiaftertrack[1], Form("Efficiency/4TH//Cluster_Multiplicity_4TH_After_tracks_by_other2"), Form("Cluster Multiplicity of 4TH GEM After tracks by other2"),11,-0.5,10.5);
      //H1(gemmultiaftertrack[2], Form("Efficiency/MS//Cluster_Multiplicity_MS_After_tracks_by_other2"), Form("Cluster Multiplicity of MS GEM After tracks by other2"),11,-0.5,10.5);   
    
   }; //BH 1-hit filtering finish

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

