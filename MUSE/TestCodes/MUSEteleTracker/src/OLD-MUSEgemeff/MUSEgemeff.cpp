#include <MUSEteleTracker.h>
#include<iostream>
#include<cmath>
#include "lumigemtree.h"
#include "TGraph.h"
#include "TF1.h"


Long_t MUSEteleTracker::startup_efficiency()
{ 
// try to get the SiPM TDC output tree:


//     tdc1190=NULL;
//      getBranchObject("v1190",(TObject **) &tdc1190);
//      if (tdc1190==NULL)
//        {
//          debug(0,"Could not find v1190 branch\n");
//        }
//      else printf("*** tdc1190 is %p\n", tdc1190);

  // For now, let's give a warning and leave early
  //debug(0,"MUSEteleTracker plugin is unavailable as it made use of testbeamanalysistree (no longer available)\n");
  //return -1;

	return Plugin::ok;


 ///////////////////////////////////////////
}
Long_t MUSEteleTracker::histos_efficiency()
{
  //effgemUS=dH2("effgemUS","Efficiency on GEM US",50,-50.0,50.0,50,-50.0,50.0);
  //effgemMS=dH2("effgemMS","Efficiency on GEM MS",50,-50.0,50.0,50,-50.0,50.0);
  //effgemDS=dH2("effgemDS","Efficiency on GEM DS",50,-50.0,50.0,50,-50.0,50.0); 


	multiplicityUS=dH1("Efficiency/US/multiplicityUS","Multiplicity US GEM",11,-0.5,10.5);
	multiplicity4TH=dH1("Efficiency/4TH/multiplicity4TH","Multiplicity 4TH GEM",11,-0.5,10.5);
	multiplicityMS=dH1("Efficiency/MS/multiplicityMS","Multiplicity MS GEM",11,-0.5,10.5);
	multiplicityDS=dH1("Efficiency/DS/multiplicityDS","Multiplicity DS GEM",11,-0.5,10.5);

  ///2d multiplicity
	mult2dUS_MS=dH2("Efficiency/mult2dUS_MS","Multiplicity US Vs MS",9,-0.5,8.5,9,-0.5,8.5);
	mult2dMS_DS=dH2("Efficiency/mult2dMS_DS","Multiplicity MS Vs DS",9,-0.5,8.5,9,-0.5,8.5);
	mult2dUS_DS=dH2("Efficiency/mult2dUS_DS","Multiplicity US Vs DS",9,-0.5,8.5,9,-0.5,8.5);

	trackmultiplicityUS=dH1("Efficiency/US/trackmultiplicityUS","Track Multiplicity US GEM",11,-0.5,10.5);
	trackmultiplicityMS=dH1("Efficiency/MS/trackmultiplicityMS","Track Multiplicity MS GEM",11,-0.5,10.5);
	trackmultiplicityDS=dH1("Efficiency/DS/trackmultiplicityDS","Track Multiplicity DS GEM",11,-0.5,10.5);

	singleGEMclust=dH2("Efficiency/US/singleGEMclust","Single GEM clusters",50,-50,50,50,-50,50);


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
    { 1840.7, 1902.7, 1964.7, 2026.7 }; // OLYMPUS, right sector GEMs // use at PSI // changed the GEM detector spaces at PSI, alomost close to the MUSE final setup. distance are in "mm" // SHOULD USE THESE DISTANCES FOR THE DATA AFTER AUGUST, 2016
  // 1840.7, 1940.7, 2040.7  }; // OLYMPUS right sector GEMs, use at PSI, distance are in "mm" SHOULD USE THE FOLLOWING DISTANCES FOR THE DATA AFTER AUGUST, 2015 TO 2016
  // { 1846.2, 2180.3, 2600.7,    // OLYMPUS, left sector GEMs // USE at t=0 in = left  for MUSE
    //  1840.7, 1920.7, 2000.7  }; // OLYMPUS right sector GEMs // use at PSI // distance are in "mm" 

      //1840.7, 2176.5, 2597.2  }; // OLYMPUS, right sector GEMs (US, MS, DS) // distance are in "mm" // SHOULD USE THESE DISTANCES FOR THE DATA BEFORE AUGUST, 2015

    const double zgempos[8]= { 1846.2, 1946.2, 2046.2, 2146.2, 1840.7, 1902.7, 1964.7, 2026.7};
    int zgemd=0;

    const char* gemdet;

  //std::string gemdet[] = {"US","DS","MI","4t"};
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
    int gemanyhit[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    int tracks[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    for (int gems=GEM_NUM; gems<GEM_NUM+4; gems++)
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
    			xmaxchg[1] = clusters->hits[g1].xl*0.4-50.; ymaxchg[1] = clusters->hits[g1].yl*0.4-50.;
    		}
    		else test=false;
    	};
    	if (gems==6)  multiplicityUS->Fill(gemanyhit[gems]);
    	if (gems==7)  multiplicity4TH->Fill(gemanyhit[gems]);
    	if (gems==8)  multiplicityMS->Fill(gemanyhit[gems]);
    	if (gems==9)  multiplicityDS->Fill(gemanyhit[gems]);

	  // mult2dUS_MS->Fill(gemanyhit[3],gemanyhit[4]);//test
	  // mult2dUS_DS->Fill(gemanyhit[3],gemanyhit[5]);
	  // mult2dMS_DS->Fill(gemanyhit[4],gemanyhit[5]);

    	if ((gems==6)&&test) H2(xmaxchg[1],ymaxchg[1], Form("Efficiency/US/max_charge_clusters_USGEM"), Form("Maximum Charge Cluster positions on US GEM"),50,-50.0,50.0,50,-50.0,50.0);
    	if ((gems==7)&&test) H2(xmaxchg[1],ymaxchg[1], Form("Efficiency/4TH/max_charge_clusters_4THGEM"), Form("Maximum Charge Cluster positions on 4TH GEM"),50,-50.0,50.0,50,-50.0,50.0);
    	if ((gems==8)&&test) H2(xmaxchg[1],ymaxchg[1], Form("Efficiency/MS/max_charge_clusters_MSGEM"), Form("Maximum Charge Cluster positions on MS GEM"),50,-50.0,50.0,50,-50.0,50.0);
    	if ((gems==9)&&test) H2(xmaxchg[1],ymaxchg[1], Form("Efficiency/DS/max_charge_clusters_DSGEM"), Form("Maximum Charge Cluster positions on DS GEM"),50,-50.0,50.0,50,-50.0,50.0);	   
    };

	// loop over all the GEMs and get the efficiency related histograms/////////////////////////////
	// Assume GEM 1,2,3,4 and each time it goto one GEM, say GEM 1, then loop over all the other GEMs, GEM, 2,3,4. If Gem 2 does not have any clusters found, the loop continues to search clusters on the GEM 3. If the GEM 3 has a cluster, then the loop continue to search any clusters on the rest GEM 4.

    for (int gem=GEM_NUM; gem<GEM_NUM+4; gem++)  // loop 1 over the first GEM
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
				  							for (unsigned int clust=0; clust<gem4hittot[zgemd]; clust++)
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


      /*I just did
      ////////////////////////Track MUltiplicity ////////////
      //   int gemanyhit[6] = { 0, 0, 0, 0, 0, 0 };
      //   int tracks[6] = { 0, 0, 0, 0, 0, 0 };
      // double minchi2fit=1000;
     
      for (int gem=t*3+0; gem<t*3+3; gem++)  // loop over the first GEM
	{
	    // bool gem1hit=true;
	  int gemhittot[6] = { 0, 0, 0, 0, 0, 0 };
	  for (unsigned int g1=0; g1<clusters->hits.size(); g1++) //get the # of clusters on first GEM in each event
	    {
	      //		if ((clusters->hits[g1].GEMid!=gem) && ((clusters->hits[g1].ampl<400)))  continue;
	      if ((clusters->hits[g1].GEMid!=gem))  continue;
	      gemhittot[gem]++;
	    };
	 
	  //  if (gemhittot[gem]>1) gem1hit=false;
	  if (gemhittot[gem]==0) continue;  
	  if (gemhittot[gem]>=1)// clusters >=1 on the first GEM
	    {	
	      for (unsigned int g1=0; g1<clusters->hits.size(); g1++)
		{
		  if (clusters->hits[g1].GEMid!=gem) continue;		
		  x[0] = clusters->hits[g1].xl*0.4-50.; y[0] = clusters->hits[g1].yl*0.4-50.;

		  if (fabs(x[0]<40)&& fabs(y[0]<40))
		    {
		      //  printf("gem 1 after %d %d %5.2lf %5.2lf \n",g1,gem,x[0],y[0]);
		      for (int gemt=gem+1; gemt<6; gemt++) // start looping on gem 2
			{ 
			  int gem1hittot[6] = { 0, 0, 0, 0, 0, 0 };
			  for (unsigned int g2=0; g2<clusters->hits.size(); g2++)//get the # of clusters on the second GEM if the first GEM have selected the a cluster
			    { 
			//	if ((clusters->hits[g2].GEMid!=gemt)&& ((clusters->hits[g2].ampl<400))) continue;	
			      if ((clusters->hits[g2].GEMid!=gemt)) continue; //Get the total # of clusters on GEM 2
			      gem1hittot[gemt]++;
			    };
			  if (gem1hittot[gemt]==0) continue;
			  if (gem1hittot[gemt]>=1)// clusters >=1 on GEM 2
			    {
			      for (unsigned int g2=0; g2<clusters->hits.size(); g2++)
				{
				  if ((clusters->hits[g2].GEMid!=gemt)) continue;
				  x[1] = clusters->hits[g2].xl*0.4-50.; y[1] = clusters->hits[g2].yl*0.4-50.;
				  if (fabs(x[1]<40)&& fabs(y[1]<40))
				    {
				      // printf("gem 2 after %d %d %d %5.2lf %5.2lf \n",g2,gem, gemt,x[1],y[1]);

				      double slopex = (x[1]-x[0])/(zgempos[gemt]-zgempos[gem]);
				      double slopey = (y[1]-y[0])/(zgempos[gemt]-zgempos[gem]);			
				      if ((gem==3)&&(gemt==4))
					{
					  zgemd=5;
					  gemdet="DS";
					};
				      if ((gem==3)&&(gemt==5))
					{
					  zgemd=4;
					  gemdet="MS";
					};
				      if ((gem==4)&&(gemt==5))
					{
					  zgemd=3;
					  gemdet="US";
					};
				      double X3 = x[0] + slopex*(zgempos[zgemd]-zgempos[gem]);
				      double Y3 = y[0] + slopey*(zgempos[zgemd]-zgempos[gem]); 
				      int gem2hittot[6] = { 0, 0, 0, 0, 0, 0 };
				      for (unsigned int g3=0; g3<clusters->hits.size(); g3++)
					{
					  if (clusters->hits[g3].GEMid!=zgemd) continue;
					  gem2hittot[zgemd]++;// search clusters on third gem (for g3)
					  x[2] = clusters->hits[g3].xl*0.4-50.; y[2] = clusters->hits[g3].yl*0.4-50.;
					  dx=x[2]-X3;
					  dy=y[2]-Y3;
					  //	  printf("gem 3 after %d %d %5.2lf %5.2lf \n",g3,clusters->hits[g3].GEMid,x[2],y[2]);
					  //printf("tracks %d %d %d %d %d %d %f %f %f %f %f %f \n",gem,g1,gemt,g2,zgemd,g3,x[0],y[0],x[1],y[1],x[2],y[2]);
					}

				      if (fabs(x[2]<40)&& fabs(y[2]<40))
					{
					  tracks[zgemd]++;
					  ///////////////////// do straight line fit through all 3 GEMs
					  /* 
					  TGraph *gtmp  = new TGraph(3, x, y);
					  TF1    *fpol1 = new TF1("fpol1", "pol1", 0, 2); // do polynomial fit among 3 points
					  gtmp->Fit(fpol1, "QN", "");
					  double b = fpol1->GetParameter(0);// get parameters of fpol1(x) = m*x + b
					  double m = fpol1->GetParameter(1);
					  Double_t chi2fit = fpol1->GetChisquare();
					  delete   fpol1;  // clean up
					  delete   gtmp;
					  if (chi2fit<minchi2fit) minchi2fit=chi2fit;
					  //printf("fit param:  %d, %d, %d, %d, %f, %f, %f\n",gem,gemt,zgemd,tracks[zgemd],m, b,chi2fit );
					  // printf("tracks %d %d %d %d %f %f %f %f %f %f %f\n",gem,gemt,zgemd,tracks[zgemd],x[0],y[0],x[1],y[1],x[2],y[2],chi2fit);
					  ////////////////////
					 / */
					  // if (zgemd==3)  trackmultiplicityUS->Fill(tracks[zgemd]);
					  //if (zgemd==4)  trackmultiplicityMS->Fill(tracks[zgemd]);
					  //if (zgemd==5)  trackmultiplicityDS->Fill(tracks[zgemd]);
					  // printf("tracks %d %d %d %f %f %f %f %f %f %d \n",gem,gemt,zgemd,x[0],y[0],x[1],y[1],x[2],y[2],tracks[zgemd]);
      /* I just did
					}
				   
				    }
				}
			    }//clusters >=1 on GEM 2 
			} //start looping on gem 2
		    }
		}
	    }//clusters >=1 on the first GEM
	}// loop over the first GEM

      */ ///  I just did
      // printf("min chi 2: %f\n",minchi2fit);
      ////////////////////////////////
      //////////////////////////////////////////////////////////////


    //////////////////////////////////////////get rid of cosmic clusters
  /*
      int k=0;
      double x[3]={-10000.0},y[3]={-10000.0};
      for (int gem=k*3+0; gem<k*3+2; gem++)  // loop over the first GEM
	{
	    // bool gem1hit=true;
	  int gemhittot[6] = { 0, 0, 0, 0, 0, 0 };
	  for (unsigned int g1=0; g1<clusters->hits.size(); g1++) //get the # of clusters on first GEM in each event
	    {
	      //		if ((clusters->hits[g1].GEMid!=gem) && ((clusters->hits[g1].ampl<400)))  continue;
	      if ((clusters->hits[g1].GEMid!=gem))  continue;
	      gemhittot[gem]++;
	    };
 printf("how many hits %d %d \n",gem,gemhittot[gem]);
 if (gem==0) multiplicityUS->Fill(gemhittot[gem]);
 
	  //  if (gemhittot[gem]>1) gem1hit=false;
	  // printf("how many %d %d\n",gem,gemhittot[gem]);
	  if (gemhittot[gem]==0) {  // No hits on the first GEM
	    // printf("No hits on the first GEM %d \n",gem);
	    for (int gemt=gem+1; gemt<6; gemt++) // start looping on gem 2
	      { 
		int gem1hittot[6] = { 0, 0, 0, 0, 0, 0 };
		for (unsigned int g2=0; g2<clusters->hits.size(); g2++)//get the # of clusters on the second GEM if the first GEM have selected the maximum cluster
		  { 	     	
		    if ((clusters->hits[g2].GEMid!=gemt)) continue; //Get the total # of clusters on GEM 2
		    gem1hittot[gemt]++;
		  };
		if (gem1hittot[gemt]==0) {// No hits on the second GEM
		  //	  continue;
		  //  printf("No hits on the second GEM %d \n",gemt);
		  if ((gem==3)&&(gemt==4))
		    {
		      zgemd=5;
		      gemdet="DS";
		    };
		  if ((gem==3)&&(gemt==5))
		    {
		      zgemd=4;
		      gemdet="MS";
		    };
		  if ((gem==4)&&(gemt==5))
		    {
		      zgemd=3;
		      gemdet="US";
		    };
		  //	      double X3 = x[0] + slopex*(zgempos[zgemd]-zgempos[gem]);
		  //	      double Y3 = y[0] + slopey*(zgempos[zgemd]-zgempos[gem]); 
		  int gem2hittot[6] = { 0, 0, 0, 0, 0, 0 };
		  for (unsigned int g3=0; g3<clusters->hits.size(); g3++)
		    {
		      if (clusters->hits[g3].GEMid!=zgemd) continue;
		      gem2hittot[zgemd]++;// search clusters on third gem (for g3)
		     
		      x[2] = clusters->hits[g3].xl*0.4-50.; y[2] = clusters->hits[g3].yl*0.4-50.;
					  //  dx=x[2]-X3;
					  //  dy=y[2]-Y3;
					  //	  printf("gem 3 after %d %d %5.2lf %5.2lf \n",g3,clusters->hits[g3].GEMid,x[2],y[2]);
					  //printf("tracks %d %d %d %d %d %d %f %f %f %f %f %f \n",gem,g1,gemt,g2,zgemd,g3,x[0],y[0],x[1],y[1],x[2],y[2]);
		      printf("see hits on the third GEM %d %d %5.2lf %5.2lf \n",zgemd,gem2hittot[zgemd],x[2],y[2]);
		    }
		  if (fabs(x[2]<40)&& fabs(y[2]<40))
		    {
		      //tracks[zgemd]++;
		      singleGEMclust->Fill(x[2],y[2]);
		    }
		};//No hits on the second GEM
	      };//start looping on gem 2
	  };// No hits on the first GEM

	}; // loop over the first GEM
/////////////////////////////////get rid of cosmic clusters
*/

	return 0;
      //  return Plugin::ok;
}

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
