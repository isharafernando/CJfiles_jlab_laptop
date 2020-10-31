#include <MUSEteleTracker.h>

//#include "cTrack.h"

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


MUSEteleTracker::MUSEteleTracker(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
}

MUSEteleTracker::~MUSEteleTracker()
{
}

Long_t MUSEteleTracker::defineHistograms()
{

	h2resmapsx=dH2(Form("MUSEtele/residualmapsx"), Form("X and Y hitmap weighted by X Residuals on MS GEM"), 50,-50,50,50,-50,50);
	h2resmapsy=dH2(Form("MUSEtele/residualmapsy"), Form("X and Y hitmap weighted by Y Residuals on MS GEM"), 50,-50,50,50,-50,50);

	std::string GEMnames0[] = {"US","4th","MI","DS"};
	std::string realornoise[] = {"Real Hits", "Noise Hits"};
	std::string particlenames[] = {"Pion", "Electron", "Muon"};
	std::string whichBHplanes[] = {" 2", " 3", "s 2 and 3"};

	for (unsigned int i=0; i<4; i++ ){
				
		totresx[i]=dH1(Form("Totals/%s Individual Hit X Residuals Total", GEMnames0[i].c_str()), "Individual Hit X Residuals;mm",200,-10,10);
		totresy[i]=dH1(Form("Totals/%s Individual Hit Y Residuals Total", GEMnames0[i].c_str()), "Individual Hit Y Residuals;mm",200,-10,10);
		
		threehitresx[i]=dH1(Form("3 GEM hits/%s Individual Hit X Residuals", GEMnames0[i].c_str()), "Individual Hit X Residuals;mm",200,-10,10);
		threehitresy[i]=dH1(Form("3 GEM hits/%s Individual Hit Y Residuals", GEMnames0[i].c_str()), "Individual Hit Y Residuals;mm",200,-10,10);
		fourhitresx[i]=dH1(Form("4 GEM hits/%s Individual Hit X Residuals", GEMnames0[i].c_str()), "Individual Hit X Residuals;mm",200,-10,10);
		fourhitresy[i]=dH1(Form("4 GEM hits/%s Individual Hit Y Residuals", GEMnames0[i].c_str()), "Individual Hit Y Residuals;mm",200,-10,10);

		firsthitresx[i]=dH1(Form("First 4 GEM hits/%s Individual Hit X Residuals", GEMnames0[i].c_str()), "Individual Hit X Residuals;mm",100,-10,10);
		firsthitresy[i]=dH1(Form("First 4 GEM hits/%s Individual Hit Y Residuals", GEMnames0[i].c_str()), "Individual Hit Y Residuals;mm",100,-10,10);
		firsthitresxcut[i]=dH1(Form("First 4 GEM hits/%s Cut Individual Hit X Residuals", GEMnames0[i].c_str()), "Individual Hit X Residuals Cut by Chi2 (1+/-0.5);mm",500,-1,1);
		firsthitresycut[i]=dH1(Form("First 4 GEM hits/%s Cut Individual Hit Y Residuals", GEMnames0[i].c_str()), "Individual Hit Y Residuals Cut by Chi2 (1+/-0.5);mm",500,-1,1);

		threechi2cutxyx[i]=dH1(Form("Chi Cut Residuals/%s X Residuals with X and Y chi2 < 3 for 3 hits",GEMnames0[i].c_str()), "X Residuals;mm",100,-1,1);
		threechi2cutxyy[i]=dH1(Form("Chi Cut Residuals/%s Y Residuals with X and Y chi2 < 3 for 3 hits",GEMnames0[i].c_str()), "Y Residuals;mm",100,-1,1);
		fourchi2cutxyx[i]=dH1(Form("Chi Cut Residuals/%s X Residuals with X and Y chi2 < 3 for 4 hits",GEMnames0[i].c_str()), "X Residuals;mm",100,-1,1);
		fourchi2cutxyy[i]=dH1(Form("Chi Cut Residuals/%s Y Residuals with X and Y chi2 < 3 for 4 hits",GEMnames0[i].c_str()), "Y Residuals;mm",100,-1,1);
		
		totxdist[i]=dH1(Form("Totals/%s x distribution",GEMnames0[i].c_str()), "x distribution;mm", 500, -50, 50);
		totydist[i]=dH1(Form("Totals/%s y distribution",GEMnames0[i].c_str()), "y distribution;mm", 500, -50, 50);
		
		tot2dmap[i]=dH2(Form("Totals/2D Track Map %s GEM", GEMnames0[i].c_str()),Form("2D Track Map %s GEM;mm;mm", GEMnames0[i].c_str()),200,-50,50,200,-50,50);

		for (unsigned int j=0; j<4; j++){

			largestremovedx[i][j]=dH1(Form("ResRemoved4hits/X Residual of %s GEM with largest residual removed %s GEM", GEMnames0[i].c_str(), GEMnames0[j].c_str()),"X Residuals;mm", 100,-10,10);
			largestremovedy[i][j]=dH1(Form("ResRemoved4hits/Y Residual of %s GEM with largest residual removed %s GEM", GEMnames0[i].c_str(), GEMnames0[j].c_str()),"Y Residuals;mm", 100,-10,10);
		
		}


		for (unsigned int j=0; j<2; j++){
			totBHresx[i][j]=dH1(Form("Totals/%s/%s Total X Residuals",realornoise[j].c_str(),GEMnames0[i].c_str()), "X Residuals;mm", 200, -10, 10);
			totBHresy[i][j]=dH1(Form("Totals/%s/%s Total Y Residuals",realornoise[j].c_str(),GEMnames0[i].c_str()), "Y Residuals;mm", 200, -10, 10);
			fourBHresx[i][j]=dH1(Form("4 GEM hits/%s/%s Total X Residuals",realornoise[j].c_str(),GEMnames0[i].c_str()), "X Residuals;mm", 200, -10, 10);
			fourBHresy[i][j]=dH1(Form("4 GEM hits/%s/%s Total Y Residuals",realornoise[j].c_str(),GEMnames0[i].c_str()), "Y Residuals;mm", 200, -10, 10);
			threeBHresx[i][j]=dH1(Form("3 GEM hits/%s/%s Total X Residuals",realornoise[j].c_str(),GEMnames0[i].c_str()), "X Residuals;mm", 200, -10, 10);
			threeBHresy[i][j]=dH1(Form("3 GEM hits/%s/%s Total Y Residuals",realornoise[j].c_str(),GEMnames0[i].c_str()), "Y Residuals;mm", 200, -10, 10);
			
			for(unsigned int N=0; N<4; N++){
				xres4mmBHcut[i][j][N]=dH1(Form("BHcuts/%i GEM hits/Plane %i/X residuals for %s GEM for 4 mm BH paddles",j+3,i,GEMnames0[N].c_str()),"X residuals;mm", 200, -9.5,10.5);
				yres4mmBHcut[i][j][N]=dH1(Form("BHcuts/%i GEM hits/Plane %i/Y residuals for %s GEM for 4 mm BH paddles",j+3,i,GEMnames0[N].c_str()),"Y residuals;mm", 200, -9.5,10.5);
				xres8mmBHcut[i][j][N]=dH1(Form("BHcuts/%i GEM hits/Plane %i/X residuals for %s GEM for 8 mm BH paddles",j+3,i,GEMnames0[N].c_str()),"X residuals;mm", 200, -9.5,10.5);
				yres8mmBHcut[i][j][N]=dH1(Form("BHcuts/%i GEM hits/Plane %i/Y residuals for %s GEM for 8 mm BH paddles",j+3,i,GEMnames0[N].c_str()),"Y residuals;mm", 200, -9.5,10.5);
										
			}

			for(unsigned int k=0; k<3; k++){
				totBHplaneresx[i][j][k]=dH1(Form("Totals/%s/Plane%s/%s Total X Residuals",realornoise[j].c_str(),whichBHplanes[k].c_str(),GEMnames0[i].c_str()), "X Residuals;mm", 200, -10, 10);
				totBHplaneresy[i][j][k]=dH1(Form("Totals/%s/Plane%s/%s Total Y Residuals",realornoise[j].c_str(),whichBHplanes[k].c_str(),GEMnames0[i].c_str()), "Y Residuals;mm", 200, -10, 10);
				fourBHplaneresx[i][j][k]=dH1(Form("4 GEM hits/%s/Plane%s/%s Total X Residuals",realornoise[j].c_str(),whichBHplanes[k].c_str(),GEMnames0[i].c_str()), "X Residuals;mm", 200, -10, 10);
				fourBHplaneresy[i][j][k]=dH1(Form("4 GEM hits/%s/Plane%s/%s Total Y Residuals",realornoise[j].c_str(),whichBHplanes[k].c_str(),GEMnames0[i].c_str()), "Y Residuals;mm", 200, -10, 10);
				threeBHplaneresx[i][j][k]=dH1(Form("3 GEM hits/%s/Plane%s/%s Total X Residuals",realornoise[j].c_str(),whichBHplanes[k].c_str(),GEMnames0[i].c_str()), "X Residuals;mm", 200, -10, 10);
				threeBHplaneresy[i][j][k]=dH1(Form("3 GEM hits/%s/Plane%s/%s Total Y Residuals",realornoise[j].c_str(),whichBHplanes[k].c_str(),GEMnames0[i].c_str()), "Y Residuals;mm", 200, -10, 10);
				

				for(unsigned int h=0; h<3; h++){
					totBHparticleresx[i][j][k][h]=dH1(Form("Totals/%s/Plane%s/%s Total X Residuals for %s hits",realornoise[j].c_str(),whichBHplanes[k].c_str(),GEMnames0[i].c_str(),particlenames[h].c_str()), "X Residuals;mm", 200, -10, 10);
					totBHparticleresy[i][j][k][h]=dH1(Form("Totals/%s/Plane%s/%s Total Y Residuals for %s hits",realornoise[j].c_str(),whichBHplanes[k].c_str(),GEMnames0[i].c_str(),particlenames[h].c_str()), "Y Residuals;mm", 200, -10, 10);
					fourBHparticleresx[i][j][k][h]=dH1(Form("4 GEM hits/%s/Plane%s/%s Total X Residuals for %s hits",realornoise[j].c_str(),whichBHplanes[k].c_str(),GEMnames0[i].c_str(),particlenames[h].c_str()), "X Residuals;mm", 200, -10, 10);
					fourBHparticleresy[i][j][k][h]=dH1(Form("4 GEM hits/%s/Plane%s/%s Total Y Residuals for %s hits",realornoise[j].c_str(),whichBHplanes[k].c_str(),GEMnames0[i].c_str(),particlenames[h].c_str()), "Y Residuals;mm", 200, -10, 10);
					threeBHparticleresx[i][j][k][h]=dH1(Form("3 GEM hits/%s/Plane%s/%s Total X Residuals for %s hits",realornoise[j].c_str(),whichBHplanes[k].c_str(),GEMnames0[i].c_str(),particlenames[h].c_str()), "X Residuals;mm", 200, -10, 10);
					threeBHparticleresy[i][j][k][h]=dH1(Form("3 GEM hits/%s/Plane%s/%s Total Y Residuals for %s hits",realornoise[j].c_str(),whichBHplanes[k].c_str(),GEMnames0[i].c_str(),particlenames[h].c_str()), "Y Residuals;mm", 200, -10, 10);
					

				}
			}
		}

		GEMXvsBHX[i]=dH2(Form("BHcuts/Plane %i/Most Upstream Gem's X coor vs. BH's hit bar's X coor",i),"GEM X vs. BH bar X;mm;mm", 200, -49.5, 50.5, 16, -39.5, 40.5);
		GEMYvsBHY[i]=dH2(Form("BHcuts/Plane %i/Most Upstream Gem's Y coor vs. BH's hit bar's Y coor",i),"GEM Y vs. BH bar Y;mm;mm", 200, -49.5, 50.5, 16, -39.5, 40.5);
		GEMXvsTotBHbar[i]=dH2(Form("BHcuts/Plane %i/Most Upstream Gem's X coor vs. BH bar hit",i),"GEM X vs. BH bar;mm", 16, -49.5, 50.5, 16, -0.5, 15.5);
		GEMYvsTotBHbar[i]=dH2(Form("BHcuts/Plane %i/Most Upstream Gem's Y coor vs. BH bar hit",i),"GEM Y vs. BH bar;mm", 16, -49.5, 50.5, 16, -0.5, 15.5);
		GEMXvsBHY[i]=dH2(Form("BHcuts/Plane %i/Most Upstream Gem's X coor vs. BH's hit bar's Y coor",i),"GEM X vs. BH bar Y;mm;mm", 200, -49.5, 50.5, 16, -39.5, 40.5);
		GEMYvsBHX[i]=dH2(Form("BHcuts/Plane %i/Most Upstream Gem's Y coor vs. BH's hit bar's X coor",i),"GEM Y vs. BH bar X;mm;mm", 200, -49.5, 50.5, 16, -39.5, 40.5);

		for(unsigned int bbar=0; bbar<16; bbar++){
			GEMXvsBHbar[bbar][i]=dH1(Form("BHcuts/Plane %i/Individual BH bars/BHbar %i Most Upstream Gem's X coor vs. BH bar hit",i,bbar),"GEM X vs. BH bar;mm", 200, -49.5, 50.5);
			GEMYvsBHbar[bbar][i]=dH1(Form("BHcuts/Plane %i/Individual BH bars/BHbar %i Most Upstream Gem's Y coor vs. BH bar hit",i,bbar),"GEM Y vs. BH bar;mm", 200, -49.5, 50.5);
		
		}

		projBHXvsBHX[i]=dH2(Form("BHcuts/Plane %i/Projected BH Bar X coor vs. BH's hit bar's X coor",i),"Projected BH Bar X vs. BH bar X;mm;mm", 16, -49.5, 50.5, 16, -49.5, 50.5);
		projBHYvsBHY[i]=dH2(Form("BHcuts/Plane %i/Projected BH Bar Y coor vs. BH's hit bar's Y coor",i),"Projected BH Bar Y vs. BH bar Y;mm;mm", 16, -49.5, 50.5, 16, -49.5, 50.5);
		
		largestremovedmx[i]=dH1(Form("ResRemoved4hits/Track slope in x if %s GEM has largest residual and removed", GEMnames0[i].c_str()),"Track slope in x", 100, -1, 1);
		largestremovedmy[i]=dH1(Form("ResRemoved4hits/Track slope in y if %s GEM has largest residual and removed", GEMnames0[i].c_str()),"Track slope in y", 100, -1, 1);
		largestremovedZX[i]=dH1(Form("ResRemoved4hits/ZX intercept if %s GEM has largest residual and removed", GEMnames0[i].c_str()),"ZX intercept", 400,-200,200);
		largestremovedZY[i]=dH1(Form("ResRemoved4hits/ZY intercept if %s GEM has largest residual and removed", GEMnames0[i].c_str()),"ZY intercept", 400,-200,200);
		
		resremovedx[i]=dH1(Form("MUSEtele/%s X Residuals with %s GEM removed",GEMnames0[i].c_str(),GEMnames0[i].c_str()), "Individual Hit X Residuals;mm",100,-10,10);
		resremovedy[i]=dH1(Form("MUSEtele/%s Y Residuals with %s GEM removed",GEMnames0[i].c_str(),GEMnames0[i].c_str()), "Individual Hit Y Residuals;mm",100,-10,10);
		resweightedx[i]=dH1(Form("MUSEtele/Weighted X distribution with %s GEM removed", GEMnames0[i].c_str()), Form("X distribution weighted by X residuals on %s GEM;mm", GEMnames0[i].c_str()), 500, -50, 50);
		resweightedy[i]=dH1(Form("MUSEtele/Weighted Y distribution with %s GEM removed", GEMnames0[i].c_str()), Form("Y distribution weighted by Y residuals on %s GEM;mm", GEMnames0[i].c_str()), 500, -50, 50);


		h2resmapx[i]=dH2(Form("MUSEtele/%s residualmapsx",GEMnames0[i].c_str()), Form("X and Y hitmap weighted by X Residuals on %s GEM;mm;mm",GEMnames0[i].c_str()), 200,-50,50,200,-50,50);
		h2resmapy[i]=dH2(Form("MUSEtele/%s residualmapsy",GEMnames0[i].c_str()), Form("X and Y hitmap weighted by Y Residuals on %s GEM;mm;mm",GEMnames0[i].c_str()), 200,-50,50,200,-50,50);

		totbest[i]=dH2(Form("Totals/Hitmap %s Best Track", GEMnames0[i].c_str()), Form("Hitmap %s GEM For the Best Track;mm;mm", GEMnames0[i].c_str()), 100, -50, 50, 100, -50, 50);


	}

	for (unsigned int i=0; i<3; i++){

		multimx[i]=dH1(Form("%i GEM hits/Track slope in x", i+2),"mx distribution", 100, -1,1);
		multimy[i]=dH1(Form("%i GEM hits/Track slope in y", i+2),"my distribution", 100, -1,1);
		multiZX[i]=dH1(Form("%i GEM hits/ZX intercept", i+2),"ZX intercept", 100, -200,200);
		multiZY[i]=dH1(Form("%i GEM hits/ZY intercept", i+2),"ZY intercept", 100, -200,200);

		for (unsigned int j=0; j<4; j++){

			multixdist[i][j]=dH1(Form("%i GEM hits/%s x distribution", i+2, GEMnames0[j].c_str()), "x distribution;mm", 500, -50, 50);
			multiydist[i][j]=dH1(Form("%i GEM hits/%s y distribution", i+2, GEMnames0[j].c_str()), "y distribution;mm", 500, -50, 50);

			multi2dmap[i][j]=dH2(Form("%i GEM hits/2D Track Map %s GEM", i+2, GEMnames0[j].c_str()),Form("2D Track Map %s GEM;mm;mm", GEMnames0[j].c_str()),200,-50,50,200,-50,50);
		
			multibest[i][j]=dH2(Form("%i GEM hits/Hitmap %s Best Track", i+2, GEMnames0[j].c_str()), Form("Hitmap %s GEM For the Best Track;mm;mm", GEMnames0[j].c_str()), 100, -50, 50, 100, -50, 50);

		}


	}

	totChi2x=dH1(Form("Totals/X Hit Reduced chi2 Total"), Form("X hit reduced chi2"), 100, 0, 10);
	totChi2y=dH1(Form("Totals/Y Hit Reduced chi2 Total"), Form("Y hit reduced chi2"), 100, 0, 10);
	threehitChi2x=dH1(Form("3 GEM hits/X Hit Reduced chi2"), Form("X hit reduced chi2"), 100, 0, 10);
	threehitChi2y=dH1(Form("3 GEM hits/Y Hit Reduced chi2"), Form("Y hit reduced chi2"), 100, 0, 10);
	fourhitChi2x=dH1(Form("4 GEM hits/X Hit Reduced chi2"), Form("X hit reduced chi2"), 100, 0, 10);
	fourhitChi2y=dH1(Form("4 GEM hits/Y Hit Reduced chi2"), Form("Y hit reduced chi2"), 100, 0, 10);

	firsthitchi2x=dH1(Form("First 4 GEM hits/X Hit Reduced chi2"), Form("X hit reduced chi2"), 100, 0, 10);
	firsthitchi2y=dH1(Form("First 4 GEM hits/Y Hit Reduced chi2"), Form("Y hit reduced chi2"), 100, 0, 10);
	

	totmx=dH1(Form("Totals/mx distribution total"),"mx distribution",100,-1,1);
	totmy=dH1(Form("Totals/my distribution total"),"my distribution",100,-1,1);
	totZX=dH1(Form("Totals/ZX intercept total"),"ZX intercept",100,-200,200);
	totZY=dH1(Form("Totals/ZY intercept total"),"ZY intercept",100,-200,200);
	
	trackmultiplicityUS=dH1("TeleTracks/US/trackmultiplicityUS","Track Multiplicity US GEM",10,0.5,10.5);
	trackmultiplicityMS=dH1("TeleTracks/MS/trackmultiplicityMS","Track Multiplicity MS GEM",10,0.5,10.5);
	trackmultiplicityDS=dH1("TeleTracks/DS/trackmultiplicityDS","Track Multiplicity DS GEM",10,0.5,10.5);

	return Plugin::ok;
}


void getLeastSquaresLine(std::vector<double> x, std::vector<double> y, std::vector<double> z, double &slopezy, double &slopezx,
 double &bzy, double &bzx){


					double avgx = std::accumulate(x.begin(),x.end(),0.0)/x.size();
					double avgy = std::accumulate(y.begin(),y.end(),0.0)/y.size();
					double avgz = std::accumulate(z.begin(),z.end(),0.0)/z.size();

					//For z-y and z-x line in target GEMs
					//form of y=mz+b
					//form of x=mz+b
					double numy = 0;
					double denom = 0;
					double numx = 0;
					//Calculating the slope
					for(unsigned int i=0; i<y.size(); i++)
					{
						numy += (y[i]-avgy)*(z[i]-avgz);
						numx += (x[i]-avgx)*(z[i]-avgz);
						denom += pow((z[i]-avgz),2);
					}
					//Now slope and intercept and you've got your line
					slopezy = numy/denom;
					slopezx = numx/denom;
					bzy = avgy-slopezy*avgz;
					bzx = avgx-slopezx*avgz;

} 


Long_t MUSEteleTracker::startup()
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


  // create output branch with tracks:
	teletracks = new TeleTracks();
	makeBranch("teletracks", (TObject**)&teletracks);
	printf(" teletracks %p\n", teletracks);

	

	return Plugin::ok;
}

int event=0;

//int trks=0;
int trk1=0;
int trk2=0;
//ofstream outf("chk_gem_rotatn.dat");

Long_t MUSEteleTracker::process()
{
	event=event+1;
  //  printf("start a new event %d \n",event);

  // char leftright[2][18] = {"downstream", "upstream"};

      // vector to store all track candidates for this event:
	std::vector <StraightTrack> TrackCands;
	
	std::vector <int>           whichclusters[4];
	std::vector <double>        chi2;
	std::string GEMnames0[] = {"US","4th","MI","DS"};

      // initialize an "empty" straight track:
	StraightTrack aTrack;
      //Target GEMs
	aTrack.x0 = -10000.0;
	aTrack.y0 = -10000.0;
	aTrack.z0 = -10000.0;
	aTrack.x1 = -10000.0;
	aTrack.y1 = -10000.0;
	aTrack.z1 = -10000.0;
	aTrack.x2 = -10000.0;
	aTrack.y2 = -10000.0;
	aTrack.z2 = -10000.0;
	aTrack.x3 = -10000.0;
	aTrack.y3 = -10000.0;
	aTrack.z3 = -10000.0;
      //IFP GEMs
	aTrack.x4 = -10000.0;
	aTrack.y4 = -10000.0;
	aTrack.z4 = -10000.0;
	aTrack.x5 = -10000.0;
	aTrack.y5 = -10000.0;
	aTrack.z5 = -10000.0;

	aTrack.mx = -10000.0;
	aTrack.my = -10000.0;
	aTrack.mxifp = -10000.0;
	aTrack.myifp = -10000.0;
	aTrack.xresidua.push_back(-10000);
	aTrack.yresidua.push_back(-10000);
	aTrack.z.push_back(-10000);
	aTrack.xchi2 = -10000.0;
	aTrack.ychi2 = -10000.0;


	TrackCands.clear();
	chi2.clear();
	teletracks->tracks.clear();
	aTrack.xresidua.clear();
	aTrack.yresidua.clear();
	aTrack.z.clear();

      // loop over all possible combinations of clusters:
	int combmulti=0, combmulti_cut=0;
	int gemmulti[4] = { 0, 0, 0, 0 };
	for (unsigned int g1=0; g1<clusters->hits.size(); g1++)
	{
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
				for (unsigned int g4=0; g4<clusters->hits.size(); g4++)
				{
					if (clusters->hits[g4].GEMid!=3) continue;
					gemmulti[3]++;
				};
			};
		};
	};

      //printf("tele %d: #clusters US: %d MI: %d DS: %d  4th: %d\n", t, gemmulti[0], gemmulti[1], gemmulti[2], gemmulti[3]);
	//I think this is what stops tracking form happening on events with multiple hits
	// if (gemmulti[0]*gemmulti[1]*gemmulti[2]*gemmulti[3]>200) {
	// //	printf("for the new event  %d, %5.3lf, %5.3lf\n",event,aTrack.mx,aTrack.my);
	// 	teletracks->tracks.push_back(aTrack);
	// };

	// BH plots replicated
	double alignment[16] = {3.0,1.25,3.8,0.75,3.5,1.75,1.6,2.45,1.2,.25,.1,0,4.5,3.0,2.5,2.5};
	int particletype = 0;
	int plane2hits = 0;
	int plane3hits = 0;
	std::string realornoise[] = {"Real Hits", "Noise Hits"};
	std::string particlenames[] = {"Pion", "Electron", "Muon"};

	std::vector< std::vector<int> > BHhits(Hits->hits.size(), std::vector<int>(2));
	std::vector< std::vector<double> > BHhitxyz(Hits->hits.size(), std::vector<double>(3));


	for(size_t i = 0; i < Hits->hits.size(); i++)
        {
			double RF = Hits->hits[i].rf;

        	//new way to extract plane and bar information after removing internal to logic function
        	int plane = Hits->hits[i].wall_id;
        	int bar = Hits->hits[i].bar_id;

			double loc[3] = {0,0,0};
			double mas[3] = {0,0,0};
			handle.BH[plane][bar]->LocalToMaster(loc,mas);

			BHhits[i][0] = plane;
			BHhits[i][1] = bar;

			BHhitxyz[i][0] = mas[0]*10;
			BHhitxyz[i][1] = mas[1]*10;
			BHhitxyz[i][2] = mas[2]*10;


			double RFmod = fmod(RF+alignment[bar],19.75);

            H1(fmod(RF+alignment[bar],19.75),TString::Format("RF/Hit in GEMs/Individual Planes/RF plane %i bar %i",plane,bar),TString::Format("RF plane %i bar %i;rf time(ns)",plane,bar),857,0,19.75);
        	H2(bar,fmod(RF+alignment[bar],19.75),Form("RF/RF plane %i",plane),Form("RF plane %i",plane),16,-0.5,15.5,857,0,19.75);
            H1(fmod(RF+alignment[bar],19.75),TString::Format("RF/Hit in GEMs/RF plane %i",plane),TString::Format("RF plane %i;rf time(ns)",plane),857,0,19.75);

			if(plane == 2)
			{
				if(RFmod > 5.0 && RFmod < 7.3)
				{
					particletype = 1;
					plane2hits++;
				}
				if(RFmod > 9.7 && RFmod < 11.5)
				{
					particletype = 2;
					plane2hits++;
				}
				if((RFmod > 18.9 && RFmod < 19.75) || (RFmod >= 0.0 && RFmod < 0.75))
				{
					particletype = 3;
					plane2hits++;
				}

			}
			if(plane == 3)
			{
				if(RFmod > 6.0 && RFmod < 8.0)
				{
					particletype = 1;
					plane3hits++;
				}
				if(RFmod > 10.5 && RFmod < 12.2)
				{
					particletype = 2;
					plane3hits++;
				}
				if((RFmod > 18.7 && RFmod < 19.75) || (RFmod >= 0.0 && RFmod < 1.5))
				{
					particletype = 3;
					plane3hits++;
				}

			}

		}
		
		int realnoise = 1;
		int planehit = 0;

		if(plane2hits == 1)
		{
			realnoise = 0;
			planehit = 1;
		}
		if(plane3hits == 1)
		{
			realnoise = 0;
			planehit = 2;
		}
		if(plane2hits == 1 && plane3hits == 1) planehit = 3;

		
		







       //if (gemmulti[0]*gemmulti[1]*gemmulti[2]*gemmulti[3]>200) return Plugin::ok; //return Plugin::ok;

       H1(gemmulti[0], Form("MUSEteleTracker/Number of possible clusters US"), Form("Number of possible clusters - US GEM"),
       	11,-0.5,10.5);
       H1(gemmulti[1], Form("MUSEteleTracker/Number of possible clusters 4TH"), Form("Number of possible clusters - 4TH GEM"),
       	11,-0.5,10.5);
       H1(gemmulti[2], Form("MUSEteleTracker/Number of possible clusters MI"), Form("Number of possible clusters - MI GEM"),
       	11,-0.5,10.5);
       H1(gemmulti[3], Form("MUSEteleTracker/Number of possible clusters DS"), Form("Number of possible clusters - DS GEM"),
       	11,-0.5,10.5);


      // The only need is the relative  distances in z here.
      //const std::vector<double> zgem0 = 
     	//{ //1846.2, 1946.2, 2046.2, 2146.2,   // OLYMPUS, left sector GEMs
     	//-536.0, -474.0 , -412.0, -350.0}; //Ethan and Ron measured for August 2018 Beamtime DS position
     	//	-836.0,-747.0,-712.0,-650};//Rough approximation of middle dowel position
		// 
	  
	  std::vector<double> zgem0 = {0, 0, 0, 0};

	  for(unsigned int k=0; k<4; k++){
		  double loc[3]={0,0,0};
		  double mas[3]={0,0,0};
		  
		  handle.GEM[k]->LocalToMaster(loc,mas);

		  zgem0[k] = mas[2]*10;
		  
	  }
	  
	  double dx, dy;
	  double xsloperec =0.0;
	  double ysloperec =0.0;


	  std::vector<double> x;
	  std::vector<double> y;
	  std::vector<double> xifp;
	  std::vector<double> yifp;
	  std::vector<double> zifp = {0,80};

	  x.resize(4);
	  y.resize(4);	

	  //Checks what gems registered a hit in the event
	  std::vector<int> GEMids;
	  GEMids.reserve(4);
	  GEMids.clear();
	  for (unsigned int i=0; i<clusters->hits.size(); i++){
	    
		if ( std::find(GEMids.begin(), GEMids.end(), clusters->hits[i].GEMid) != GEMids.end() ) continue;

		GEMids.push_back(clusters->hits[i].GEMid);

	  }
	  std::sort(GEMids.begin(), GEMids.end());


	  x.resize(GEMids.size());
	  y.resize(GEMids.size());

	  std::vector<double> zgem;
	  std::string GEMnames[GEMids.size()];

	  zgem.resize(GEMids.size());

	  for (unsigned int i=0; i<GEMids.size(); i++) {
		zgem[i]=zgem0[GEMids[i]];
		GEMnames[i]=GEMnames0[GEMids[i]];

	  }
	  


	  bool goodtrack=false;
	  int trks=0;
	  int firsthit = 0;
	  int hitsperGEM[4]={0,0,0,0};
	  std::vector<double> x1perGEM={0,0,0,0};
	  std::vector<double> y1perGEM={0,0,0,0};

      for (unsigned int g1=0; g1<clusters->hits.size(); g1++) // Check for GEMids[3] clusters
      {
	  // printf("check gem 1  %d %d \n",g1,clusters->hits[g1].GEMid);
	  	if (GEMids.size() == 4) {
      	  if (clusters->hits[g1].GEMid!=GEMids[3]) continue;
      	  x[3]=(clusters->hits[g1].xl*0.4-50.); 
      	  y[3]=(clusters->hits[g1].yl*0.4-50.);

		  hitsperGEM[3]+= 1;
		  x1perGEM[3]=x[3];
		  y1perGEM[3]=y[3];

	   //printf("gem 1 before %d %d %5.2lf %5.2lf \n",g1,clusters->hits[g1].GEMid,x[0],y[0]);
      	  // if (fabs(x[3])>48 || fabs(y[3])>48) continue;
	  // printf("gem 1 after %d %d  %5.2lf %5.2lf \n",g1,clusters->hits[g1].GEMid,x[0],y[0]);
		  }
	  for (unsigned int g2=0; g2<clusters->hits.size(); g2++) // Check for GEMids[2] clusters
	  {
	      // printf("check gem 2  %d %d \n",g2,clusters->hits[g2].GEMid);
		if (GEMids.size() >= 3) {
	  	  if (clusters->hits[g2].GEMid!=GEMids[2]) continue;
	  	  x[2]=(clusters->hits[g2].xl*0.4-50.); 
	  	  y[2]=(clusters->hits[g2].yl*0.4-50.);

		  hitsperGEM[2] += 1;
		  x1perGEM[2]=x[2];
		  y1perGEM[2]=y[2];
	      // printf("gem 2 before %d %d %5.2lf %5.2lf \n",g2,clusters->hits[g2].GEMid,x[1],y[1]);
	  	  // if (fabs(x[2])>48 || fabs(y[2])>48) continue;
	      // printf("gem 2 after %d %d %5.2lf %5.2lf \n",g2,clusters->hits[g2].GEMid,x[1],y[1]);

		}
	      for (unsigned int g3=0; g3<clusters->hits.size(); g3++) // Check for GEMids[1] clusters
		 //   for (unsigned int g3=gem3clust; g3<gem3clust+1; g3++)
	      {
		      // printf("bool min chi 2 again: %d %f\n",g3,minchi2fit);
		  // printf("check gem 3  %d %d \n",g3,clusters->hits[g3].GEMid);
		  	if(GEMids.size() < 2) break;

	      	  if (clusters->hits[g3].GEMid!=GEMids[1]) continue;
	      	  x[1]=(clusters->hits[g3].xl*0.4-50.); 
	      	  y[1]=(clusters->hits[g3].yl*0.4-50.);

			  hitsperGEM[1] += 1;
			  x1perGEM[1]=x[1];
		      y1perGEM[1]=y[1];

		  //  printf("gem 3 before %d %d %5.2lf %5.2lf \n",g3,clusters->hits[g3].GEMid,x[2],y[2]);
	      	  // if (fabs(x[1])>48 || fabs(y[1])>48) continue;


		      for (unsigned int g4=0; g4<clusters->hits.size(); g4++) // Check for GEMids[0] clusters
		      {
		      	if (clusters->hits[g4].GEMid!=GEMids[0]) continue;
		      	x[0]=(clusters->hits[g4].xl*0.4-50.); 
		      	y[0]=(clusters->hits[g4].yl*0.4-50.);

				hitsperGEM[0] += 1;
				x1perGEM[0]=x[0];
		        y1perGEM[0]=y[0];

 ///////////////////////////Need to deactivate DS GEM for tracking because some APVs did not work for some runs. Use the two clusters on US and MS GEMs to define the track and project it to the 4TH GEM. Then calculate the track residuals on the 4TH GEM./////////////////////////	     

		      	// if (fabs(x[0])>48 || fabs(y[0])>48) continue; 

				if (GEMids.size() == 4 && firsthit != -1) firsthit = 1;

		      	double slopezx = -10000;
		      	double slopezy = -10000;
				if ((xsloperec==slopezx)&&(ysloperec==slopezy))continue; 
				goodtrack=true;
				trks=trks+1;
				trackmultiplicityUS->Fill(trks);
				trackmultiplicityMS->Fill(trks);
				trackmultiplicityDS->Fill(trks);

				// printf("gem 3 after sloperec cut %d %5.2lf %5.2lf %5.2lf %5.2lf %f %f %f %f %d %f \n",g3,x[2],y[2],x[0],y[0],dx, dy,slopex,slopey,trks,minchi2fit);
				//This has to be done because of the way they initialize their track at the beginning
				

				whichclusters[0].push_back(g1);
				whichclusters[0].push_back(g2);
				whichclusters[0].push_back(g3);
				whichclusters[0].push_back(g4);
				//This should be removed from class definition
				aTrack.telescope = 1;		 

		   
		   //      //Quick crap function for tracks between IFP GEMs
	    //   //first GEM
				// for(unsigned int ifp1=0; ifp1<clusters->hits.size(); ifp1++)
				// {
				// 	if (clusters->hits[ifp1].GEMid!=4) continue;
				// 	xifp.push_back(clusters->hits[ifp1].xl*0.4-50.); 
				// 	yifp.push_back(clusters->hits[ifp1].yl*0.4-50.);
				// 	//if (fabs(xifp[0]>40)&& fabs(yifp[0]>40)) continue;
				// 	if(xifp[0]<-40 || xifp[0]>40 || yifp[0]>15 || yifp[0]<-15) continue; //cut on physical region with offset collimator
				// 	//second GEM
				// 	for(unsigned int ifp2=0; ifp2<clusters->hits.size(); ifp2++)
				// 	{
				// 		if (clusters->hits[ifp2].GEMid!=5) continue;
						
				// 		xifp.push_back(clusters->hits[ifp2].xl*0.4-50.); 
				// 		yifp.push_back(clusters->hits[ifp2].yl*0.4-50.);
				// 	  // if (fabs(xifp[1]>)&& fabs(yifp[1]>40)) continue;	
				// 	if(xifp[1]<-40 || xifp[1]>40 || yifp[1]>15 ||yifp[1]<-15) continue; //again cut on physical region
				// 	  //Straight ripped from tracking of normal GEMs below
				// 	  // Calculate track angles using the cluster coordinates on US and MS GEM data:
				// 	if(y.size()!=zgem.size() || x.size()!=zgem.size() )
				// 	{
				// 		debug(1,"Multiplicity > 1 per GEM. Skipping\n");
				// 		//std::cout << "y size x size z size "<< y.size() << " " << x.size() << " " << zgem.size() <<std::endl;
				// 		break;
				// 	}
				// 	double slopexifp = -10000;
				// 	double slopeyifp = -10000;
				// 	slopexifp = (xifp[1]-xifp[0])/(zifp[1]-zifp[0]); // x Slope between US and MI GEMS
				// 	slopeyifp = (yifp[1]-yifp[0])/(zifp[1]-zifp[0]); // y Slope between US and MI GEMS

				// 		//TA DA tracked as good as below
						
				// 	H1(xifp[0],"US X IFP distribution","X IFP Distribution",500,-50,50);	 
				// 	H1(yifp[0],"US Y IFP distribution","Y IFP Distribution",500,-50,50);
				// 	H2(xifp[0],yifp[0],"2D Track Map US - IFP","2D Track Map US - IFP",500,-50,50,500,-50,50);		 
				// 	H1(xifp[1],"DS X IFP distribution","X IFP Distribution",500,-50,50);	 
				// 	H1(yifp[1],"DS Y IFP distribution","Y IFP Distribution",500,-50,50);		 
				// 	H2(xifp[1],yifp[1],"2D Track Map DS - IFP","2D Track Map DS - IFP",500,-50,50,500,-50,50);		 
				// 	H1(slopexifp,"IFP x slope","IFP x slope",100,-1,1);
				// 	H1(slopeyifp,"IFP y slope","IFP y slope",100,-1,1);
				// 	H2(slopexifp,slopeyifp,"2D Slope Map IFP","2D Track Map IFP",100,-1,1,100,-1,1);		 

					//Only do fitting if one hit per GEM
					//At current stage GEMs can have ten clusters per event
					//Too much to deal with now, so only look if total clusters is 4 for Target GEMs and at least 1 in each IFP GEM
					//

					//Begin Least Squares line calculations for target GEMs
					//Averages
					//////////////////////////////////////////////////////////////////////////
					//NOTE: This code only finds the first valid hit for each GEM if there is ONLY ONE HIT IN THE GEM
					// and then attempts to find a least squares fit to those data points.
					//This means it is very sensitive to the accuracy of the cluster finder
					/////////////////////////////////////////////////////////////////////////
					for(unsigned int i=0; i<BHhits.size(); i++){
						GEMXvsBHX[BHhits[i][0]]->Fill(x[0],BHhitxyz[i][0]);
						GEMYvsBHY[BHhits[i][0]]->Fill(y[0],BHhitxyz[i][1]);
						
						GEMXvsTotBHbar[BHhits[i][0]]->Fill(x[0],BHhits[i][1]);
						GEMYvsTotBHbar[BHhits[i][0]]->Fill(y[0],BHhits[i][1]);
						
						GEMXvsBHY[BHhits[i][0]]->Fill(x[0],BHhitxyz[i][1]);
						GEMYvsBHX[BHhits[i][0]]->Fill(y[0],BHhitxyz[i][0]);

						GEMXvsBHbar[BHhits[i][1]][BHhits[i][0]]->Fill(x[0]);
						GEMYvsBHbar[BHhits[i][1]][BHhits[i][0]]->Fill(y[0]);
							
					}

					double bzy = 0;
					double bzx = 0;

					getLeastSquaresLine(x,y,zgem,slopezy,slopezx,bzy,bzx);
				
					
					double xchi2 = 0;
					double ychi2 = 0;
					

					//residuals calculation
					std::vector<double> xres;
					std::vector<double> yres;
					for(unsigned int i=0; i<y.size(); i++)
					{
						double residualx = x[i]-(slopezx*zgem[i]+bzx);
						double residualy = y[i]-(slopezy*zgem[i]+bzy);
						xres.push_back(residualx);

						xchi2 += pow(residualx/0.2662,2);//Non reduced chi-sq, chose 0.2 for error as estimate
						ychi2 += pow(residualy/0.2844,2);//Non reduced chi-sq, chose 0.2 for error as estimate
						
						// xchi2 += pow(residualx/clusters->hits[i].xlerr,2);
						// ychi2 += pow(residualy/clusters->hits[i].ylerr,2);
						


						yres.push_back(residualy);

						//Filling individual hit residual histograms

						if (GEMids.size() != 2){
							totresx[GEMids[i]]->Fill(residualx);
							totresy[GEMids[i]]->Fill(residualy);

							totBHresx[GEMids[i]][realnoise]->Fill(residualx);
							totBHresy[GEMids[i]][realnoise]->Fill(residualy);
							if(planehit != 0){
								totBHplaneresx[GEMids[i]][realnoise][planehit - 1]->Fill(residualx);
								totBHplaneresy[GEMids[i]][realnoise][planehit - 1]->Fill(residualy);

								if(particletype != 0){
									totBHparticleresx[GEMids[i]][realnoise][planehit - 1][particletype - 1]->Fill(residualx);
									totBHparticleresy[GEMids[i]][realnoise][planehit - 1][particletype - 1]->Fill(residualy);
									
								}
							}							
						}

						for (unsigned int b=0; b<BHhits.size(); b++){

							if (fmod(BHhits[b][0],2) == 0 && GEMids.size() != 2){

								double projectedBHhitx = 0;
								double projectedBHhity = 0;

								projectedBHhity = -(BHhitxyz[b][2]*slopezy+bzy);
								projectedBHhitx = -(BHhitxyz[b][2]*slopezx+bzx);

								projBHXvsBHX[BHhits[b][0]]->Fill(projectedBHhitx,BHhitxyz[b][0]);
								projBHYvsBHY[BHhits[b][0]]->Fill(projectedBHhity,BHhitxyz[b][1]);

								if (BHhits[b][1] > 5 && BHhits[b][1] < 12){

									if(fabs(projectedBHhity-BHhitxyz[b][1]) < 2 && fabs(projectedBHhitx-BHhitxyz[b][0]) < 5){

										xres4mmBHcut[BHhits[b][0]][GEMids.size()-3][GEMids[i]]->Fill(residualx);
										yres4mmBHcut[BHhits[b][0]][GEMids.size()-3][GEMids[i]]->Fill(residualy);
										break;
									}

								}
								if (BHhits[b][1] > 11 && BHhits[b][1] < 6){

									if(fabs(projectedBHhity-BHhitxyz[b][1]) < 4 && fabs(projectedBHhitx-BHhitxyz[b][0]) < 5){

										xres8mmBHcut[BHhits[b][0]][GEMids.size()-3][GEMids[i]]->Fill(residualx);
										yres8mmBHcut[BHhits[b][0]][GEMids.size()-3][GEMids[i]]->Fill(residualy);
										break;
									}

								}

							}
							if (fmod(BHhits[b][0],2) == 1 && GEMids.size() != 2){

								double projectedBHhitx = 0;
								double projectedBHhity = 0;

								projectedBHhitx = -(BHhitxyz[b][2]*slopezx+bzx);
								projectedBHhity = -(BHhitxyz[b][2]*slopezy+bzy);

								projBHXvsBHX[BHhits[b][0]]->Fill(projectedBHhitx,BHhitxyz[b][0]);
								projBHYvsBHY[BHhits[b][0]]->Fill(projectedBHhity,BHhitxyz[b][1]);

								if (BHhits[b][1] > 5 && BHhits[b][1] < 12){

									if(fabs(projectedBHhitx-BHhitxyz[b][0]) < 2 && fabs(projectedBHhity-BHhitxyz[b][1]) < 5){

										xres4mmBHcut[BHhits[b][0]][GEMids.size()-3][GEMids[i]]->Fill(residualx);
										yres4mmBHcut[BHhits[b][0]][GEMids.size()-3][GEMids[i]]->Fill(residualy);
										break;
									}

								}
								if (BHhits[b][1] > 11 && BHhits[b][1] < 6){

									if(fabs(projectedBHhitx-BHhitxyz[b][0]) < 2 && fabs(projectedBHhity-BHhitxyz[b][1]) < 5){

										xres8mmBHcut[BHhits[b][0]][GEMids.size()-3][GEMids[i]]->Fill(residualx);
										yres8mmBHcut[BHhits[b][0]][GEMids.size()-3][GEMids[i]]->Fill(residualy);
										break;
									}

								}

							}


						}
						

						if (firsthit == 1){
							firsthitresx[i]->Fill(residualx);
							firsthitresy[i]->Fill(residualy);

						}

						switch(GEMids.size()){

							case 3:
								threehitresx[GEMids[i]]->Fill(residualx);
								threehitresy[GEMids[i]]->Fill(residualy);
								threeBHresx[GEMids[i]][realnoise]->Fill(residualx);
								threeBHresy[GEMids[i]][realnoise]->Fill(residualy);
								if(planehit != 0){
									threeBHplaneresx[GEMids[i]][realnoise][planehit - 1]->Fill(residualx);
									threeBHplaneresy[GEMids[i]][realnoise][planehit - 1]->Fill(residualy);

									if(particletype != 0){
										threeBHparticleresx[GEMids[i]][realnoise][planehit - 1][particletype - 1]->Fill(residualx);
										threeBHparticleresy[GEMids[i]][realnoise][planehit - 1][particletype - 1]->Fill(residualy);
									
									}
								}
								break;

							case 4:
								fourhitresx[GEMids[i]]->Fill(residualx);
								fourhitresy[GEMids[i]]->Fill(residualy);
								fourBHresx[GEMids[i]][realnoise]->Fill(residualx);
								fourBHresy[GEMids[i]][realnoise]->Fill(residualy);
								if(planehit != 0){
									fourBHplaneresx[GEMids[i]][realnoise][planehit - 1]->Fill(residualx);
									fourBHplaneresy[GEMids[i]][realnoise][planehit - 1]->Fill(residualy);

									if(particletype != 0){
										fourBHparticleresx[GEMids[i]][realnoise][planehit - 1][particletype - 1]->Fill(residualx);
										fourBHparticleresy[GEMids[i]][realnoise][planehit - 1][particletype - 1]->Fill(residualy);
									
									}
								}
								break;

						}

					}
					

					if (GEMids.size() != 2) {
					ychi2 = ychi2/(GEMids.size()-2.0);//
					xchi2 = xchi2/(GEMids.size()-2.0);//reduced chi-sq
					}

					//Filling chi2 histograms

					if (GEMids.size() != 2){
						totChi2x->Fill(xchi2);
						totChi2y->Fill(ychi2);
					}

					switch(GEMids.size()){

						case 3:
							threehitChi2x->Fill(xchi2);
							threehitChi2y->Fill(ychi2);
							break;

						case 4:
							fourhitChi2x->Fill(xchi2);
							fourhitChi2y->Fill(ychi2);
							break;

					}

					// Filling first hit and chi cut histograms
					if (firsthit == 1){
						firsthitchi2x->Fill(xchi2);
						firsthitchi2y->Fill(ychi2);

					}

					if (xchi2 < 3.0 && ychi2 < 3.0){
						for(unsigned int i=0; i<y.size(); i++)
						{												
							double residualx = x[i]-(slopezx*zgem[i]+bzx);
							double residualy = y[i]-(slopezy*zgem[i]+bzy);
						
							switch (GEMids.size()){
								case 3:
									threechi2cutxyx[GEMids[i]]->Fill(residualx);
									threechi2cutxyy[GEMids[i]]->Fill(residualy);
									break;
								case 4:
									fourchi2cutxyx[GEMids[i]]->Fill(residualx);
									fourchi2cutxyy[GEMids[i]]->Fill(residualy);
									break;
							}
						}
					}

					

					//Ethan's sanity check plots
				    //Turns out I'm insane
					for(unsigned int i=0;i<y.size();i++)
					{
						if(x[i]!=-10000){
							totxdist[GEMids[i]]->Fill(x[i]);
							multixdist[GEMids.size()-2][GEMids[i]]->Fill(x[i]);
						}
						if(y[i]!=-10000){
							totydist[GEMids[i]]->Fill(y[i]);
							multiydist[GEMids.size()-2][GEMids[i]]->Fill(y[i]);						
						}
						if(x[i]!=-1e4 && y[i]!=-1e4 ){
							tot2dmap[GEMids[i]]->Fill(x[i],y[i]);
							multi2dmap[GEMids.size()-2][GEMids[i]]->Fill(x[i],y[i]);
						
						}	
					}
					
					//Filling slope and intercept distribution histograms
					totmx->Fill(slopezx);
					totmy->Fill(slopezy);
					totZX->Fill(bzx);
					totZY->Fill(bzy);
					
					multimx[GEMids.size()-2]->Fill(slopezx);
					multimy[GEMids.size()-2]->Fill(slopezy);
					multiZX[GEMids.size()-2]->Fill(bzx);
					multiZY[GEMids.size()-2]->Fill(bzy);

					// If 4 GEMs are hit the GEM hit with the largest residual is removed and the linear approximation is done with the other 3 GEMs
					
					if (GEMids.size() == 4){

						int largestresx=0;
						int largestresy=0;


						for (int i=0; i<4; i++){

							if (fabs(xres[largestresx]) < fabs(xres[i])) largestresx = i;
							if (fabs(yres[largestresy]) < fabs(yres[i])) largestresy = i;



						}

						if (fabs(xres[largestresx]) <= fabs(yres[largestresy])) largestresx = largestresy;
						if (fabs(xres[largestresx]) >= fabs(yres[largestresy])) largestresy = largestresx;

						std::vector<double> x2(4);
						std::vector<double> y2(4);
						std::vector<double> zgem2(4);
						std::vector<int> GEMids2(4);
						std::string GEMnames2[] = {"0", "0", "0"};

						std::copy (x.begin(), x.end(), x2.begin());
						std::copy (y.begin(), y.end(), y2.begin());
						std::copy (zgem.begin(), zgem.end(), zgem2.begin());
						std::copy (GEMids.begin(), GEMids.end(), GEMids2.begin());

						x2.erase(x2.begin()+largestresx);
						y2.erase(y2.begin()+largestresy);
						zgem2.erase(zgem2.begin()+largestresx);
						GEMids2.erase(GEMids2.begin()+largestresx);
						GEMnames2[0] = GEMnames0[GEMids2[0]];
						GEMnames2[1] = GEMnames0[GEMids2[1]];
						GEMnames2[2] = GEMnames0[GEMids2[2]];

						double slopezx2 = 0;
						double slopezy2 = 0;
						double bzx2=0;
						double bzy2=0;
						double residualx2=0;
						double residualy2=0;

						getLeastSquaresLine(x2, y2, zgem2, slopezy2, slopezx2, bzy2, bzx2);

						for (int i=0; i<4; i++){

							residualx2 = x[i]-(slopezx2*zgem[i]+bzx2);
							residualy2 = y[i]-(slopezy2*zgem[i]+bzy2);

							largestremovedx[i][largestresx]->Fill(residualx2);
							largestremovedy[i][largestresy]->Fill(residualy2);

						}

						largestremovedmx[largestresx]->Fill(slopezx2);
						largestremovedmy[largestresy]->Fill(slopezy2);
						largestremovedZX[largestresx]->Fill(bzx2);
						largestremovedZY[largestresy]->Fill(bzy2);

						// Checking for inefficient areas of the GEMs using the linear approximation of three GEMs and calculating the
						// residual for the excluded GEM.

						std::vector<double> x1(3);
                        std::vector<double> y1(3);
                        std::vector<double> zgem1(3);

                        double slopezy1;
                        double slopezx1;
                        double bzy1;
                        double bzx1;

                        for(unsigned int i=0; i<4; i++){

                            x1[0] = x[i];
                            x1[1] = x[(i+1)%4];
                            x1[2] = x[(i+2)%4];

                            y1[0] = y[i];
                            y1[1] = y[(i+1)%4];
                            y1[2] = y[(i+2)%4];

                            zgem1[0] = zgem[i];
                            zgem1[1] = zgem[(i+1)%4];
                            zgem1[2] = zgem[(i+2)%4];

                            getLeastSquaresLine(x1,y1,zgem1,slopezy1,slopezx1,bzy1,bzx1);
                            
                            double residualx1 = x[(i+3)%4]-(slopezx1*zgem[(i+3)%4]+bzx1);
                            double residualy1 = y[(i+3)%4]-(slopezy1*zgem[(i+3)%4]+bzy1);

							//Filling histograms of the residual for the GEM not in the linear approximation calculations

							resremovedx[(i+3)%4]->Fill(residualx1);
							resremovedy[(i+3)%4]->Fill(residualy1);

							
							residualx1 = fabs(residualx1);
							residualy1 = fabs(residualy1);

							resweightedx[(i+3)%4]->Fill(x[(i+3)%4],residualx1);
							resweightedy[(i+3)%4]->Fill(y[(i+3)%4],residualy1);
							
							h2resmapx[(i+3)%4]->Fill(x[(i+3)%4],y[(i+3)%4], residualx1);
							h2resmapy[(i+3)%4]->Fill(x[(i+3)%4],y[(i+3)%4], residualy1);

							



                         }



					}

					// For 3 GEM hits, largest residual is removed and the linear approximation is redone (Only if US and DS are hit)
					if (GEMids.size() == 3 && GEMids[0]==0 && GEMids[2]==3 ){

						int largestresx=0;
						int largestresy=0;


						for (int i=0; i<3; i++){

							if (fabs(xres[largestresx]) < fabs(xres[i])) largestresx = i;
							if (fabs(yres[largestresy]) < fabs(yres[i])) largestresy = i;

						}

						if (fabs(xres[largestresx]) <= fabs(yres[largestresy])) largestresx = largestresy;
						if (fabs(xres[largestresx]) >= fabs(yres[largestresy])) largestresy = largestresx;

						std::vector<double> x2(3);
						std::vector<double> y2(3);
						std::vector<double> zgem2(3);
						std::vector<int> GEMids2(3);
						std::string GEMnames2[] = {"0", "0"};

						std::copy (x.begin(), x.end(), x2.begin());
						std::copy (y.begin(), y.end(), y2.begin());
						std::copy (zgem.begin(), zgem.end(), zgem2.begin());
						std::copy (GEMids.begin(), GEMids.end(), GEMids2.begin());

						x2.erase(x2.begin()+largestresx);
						y2.erase(y2.begin()+largestresy);
						zgem2.erase(zgem2.begin()+largestresx);
						GEMids2.erase(GEMids2.begin()+largestresx);
						GEMnames2[0] = GEMnames0[GEMids2[0]];
						GEMnames2[1] = GEMnames0[GEMids2[1]];

						double slopezx2 = 0;
						double slopezy2 = 0;
						double bzx2=0;
						double bzy2=0;
						double residualx2=0;
						double residualy2=0;

						getLeastSquaresLine(x2, y2, zgem2, slopezy2, slopezx2, bzy2, bzx2);

						for (int i=0; i<3; i++){

							residualx2 = x[i]-(slopezx2*zgem[i]+bzx2);
							residualy2 = y[i]-(slopezy2*zgem[i]+bzy2);

							H1(residualx2, Form("ResRemoved3hits/X Residual of %s GEM with largest residual removed %s GEM", GEMnames0[GEMids[i]].c_str(), GEMnames0[GEMids[largestresx]].c_str()),"X Residuals;mm", 100,-10,10);
							H1(residualy2, Form("ResRemoved3hits/Y Residual of %s GEM with largest residual removed %s GEM", GEMnames0[GEMids[i]].c_str(), GEMnames0[GEMids[largestresy]].c_str()),"Y Residuals;mm", 100,-10,10);

						}

						H1(slopezx2, Form("ResRemoved3hits/Track Slope in X with largest residual removed %s GEM", GEMnames0[GEMids[largestresx]].c_str()),"Track Slope in X", 100,-10,10);
						H1(slopezy2, Form("ResRemoved3hits/Track Slope in Y with largest residual removed %s GEM", GEMnames0[GEMids[largestresy]].c_str()),"Track Slope in Y", 100,-10,10);
						H1(bzx2, Form("ResRemoved3hits/ZX intercept with largest residual removed %s GEM", GEMnames0[GEMids[largestresx]].c_str()),"ZX intercept;mm", 100,-200,200);
						H1(bzy2, Form("ResRemoved3hits/ZY intercept with largest residual removed %s GEM", GEMnames0[GEMids[largestresx]].c_str()),"ZY intercept;mm", 100,-200,200);


					}

					
						//Since it is a track we should put in the tracks predicted hit positions
						aTrack.x0 = slopezx*zgem0[0]+bzx;
						aTrack.y0 = slopezy*zgem0[0]+bzy;
						aTrack.z0 = zgem0[0];
						aTrack.x1 = slopezx*zgem0[1]+bzx;
						aTrack.y1 = slopezy*zgem0[1]+bzy;
						aTrack.z1 = zgem0[1];
						aTrack.x2 = slopezx*zgem0[2]+bzx;
						aTrack.y2 = slopezy*zgem0[2]+bzy;
						aTrack.z2 = zgem0[2];
						aTrack.x3 = slopezx*zgem0[3]+bzx;
						aTrack.y3 = slopezy*zgem0[3]+bzy;
						aTrack.z3 = zgem0[3];
						aTrack.bx = bzx;
						aTrack.by = bzy;

						// aTrack.x4 = xifp[0];
						// aTrack.y4 = yifp[0];
						// aTrack.z4 = zifp[0];
						// aTrack.x5 = xifp[1];
						// aTrack.y5 = yifp[1];
						// aTrack.z5 = zifp[1];
						// aTrack.bxifp = xifp[0];
						// aTrack.byifp = yifp[0];

						aTrack.mx = slopezx;
						aTrack.my = slopezy;
						
						//aTrack.mxifp = slopexifp;
						//aTrack.myifp = slopeyifp;
					

						TrackCands.push_back(aTrack);
			  /// creat 2D residual maps ////////////
			  //	  double vec2d=sqrt((dx*dx)+(dy*dy));
			  //   outf<< x[1]<<"\t"<<y[1]<<"\t"<<dx<<"\t"<<dy<<std::endl; 		 

				    h2resmapsx->Fill(x[1],y[1],dx);
				    h2resmapsy->Fill(x[1],y[1],dy);   

				// double thischi2 = 0.0;

			  // printf("before thischi2:  %d,  %d, %d,  %d\n",g1,g2,g3,aTrack.xresidua.size());

				//   for (unsigned int j=0; j<aTrack.xresidua.size(); j++){
			    //		thischi2 = pow(aTrack.xresidua[j], 2.) + pow(aTrack.yresidua[j], 2.);// THIS IS NOT CHISQ!
			    //  printf("after thischi2:  %d,  %d, %d,  %d, %d, %5.3lf, %5.3lf,%5.3lf\n",g1,g2,g3,j, aTrack.xresidua.size(),aTrack.xresidua[j],aTrack.yresidua[j],thischi2);
				//	};

			  // printf("before thischi2:  %d, %d,  %d, %d,  %5.3lf,%5.3lf, %5.3lf\n",trks,g1,g2,g3,dx,dy,thischi2);
			  // aTrack.chi_sq=minchi2fit;
			  // if (aTrack.chi_sq <1.0) teletracks->tracks.push_back(aTrack);  // put all the track to the MUSEteletracker tree
			  		teletracks->tracks.push_back(aTrack);  // put all the track to the MUSEteletracker tree
			  // if (aTrack.chi_sq <0.4) printf("use thischi2 cut : %d, %d,  %d, %d,  %5.3lf,%5.3lf, %5.3lf\n",trks,g1,g2,g3,dx,dy,aTrack.chi_sq);

					
					//chi2.push_back(thischi2);
					  
					xsloperec=slopezx;
					ysloperec=slopezy;

					firsthit = -1;

		  //TrackCands.clear(); //Ethan Commented this out because why do you clear it here if it is needed right below this to make other plots?
					//   }//IFP GEM 1

	         	//}//IFP GEM 2
			}; //DS GEM
		};//MS GEM
		if (GEMids.size() < 3) break;
	    };//4TH GEM
		if (GEMids.size() < 4) break;
	};//US GEM



	if (hitsperGEM[0] == 1 && hitsperGEM[1] == 1 && hitsperGEM[2] == 1 && hitsperGEM[3] == 1)
	{
		double slopezy1perGEM = 0;
		double slopezx1perGEM = 0;
		double bzy1perGEM = 0;
		double bzx1perGEM = 0;

		getLeastSquaresLine(x1perGEM, y1perGEM, zgem0, slopezy1perGEM, slopezx1perGEM, bzy1perGEM, bzx1perGEM);

		for (int i=0; i<4; i++)
		{
			double res1hitx = x1perGEM[i]-(slopezx1perGEM*zgem0[i]+bzx1perGEM);
			double res1hity = y1perGEM[i]-(slopezy1perGEM*zgem0[i]+bzy1perGEM);

			H1(res1hitx, Form("1 hit per GEM per Event/%s X residuals", GEMnames0[i].c_str()),Form("%s X residuals for 1 hit per GEM per Event;mm", GEMnames0[i].c_str()), 100, -10, 10);
			H1(res1hity, Form("1 hit per GEM per Event/%s Y residuals", GEMnames0[i].c_str()),Form("%s Y residuals for 1 hit per GEM per Event;mm", GEMnames0[i].c_str()), 100, -10, 10);


		}

	}

	/*
	if (!goodtrack) {
		//teletracks->tracks.push_back(aTrack);
	//	printf("for everythng after  %d, %5.3lf, %5.3lf\n",event,aTrack.mx,aTrack.my);
	};

	H1(combmulti, Form("MUSEteleTracker/CombinationMulti"),
		Form("combination multiplicity"), 1000, -0.5, 999);
	H1(combmulti_cut, 
		Form("MUSEteleTracker/CombinationMulti_cut"),
		Form("combination multiplicity cut"), 1000, -0.5, 999);
      if (TrackCands.size()==0) return Plugin::ok; //return Plugin::ok;

      // loop over all track candidates and select the best one(s):

      int best       = 0;
      int secondbest = 0;
      double best_ch2       = 10000.000;
      double secondbest_ch2 = 10000.000;

      for (unsigned int i=1; i<TrackCands.size(); i++)
      {
      	if (chi2[i]<chi2[best])
      	{
      		secondbest = best;
      		best = i;	
      		secondbest_ch2=chi2[best];
      	}
      	else {
      		if (((chi2[i]<secondbest_ch2))) 
      		{
      			secondbest_ch2=chi2[i];
      			secondbest =i;
      		};	       
      	};	  
      };	
		*/
       ///////////////////////////////////////////////
      //  teletracks->tracks.clear();  // put the best and second best track to the MUSEteletracker tree

   //    aTrack = TrackCands[best];
   //    for (unsigned int i=0; i<TrackCands[best].xresidua.size(); i++)
   //    {
   //    	trk1=trk1+1;
   //    	aTrack.xresidua.push_back(TrackCands[best].xresidua[i]);
   //    	aTrack.yresidua.push_back(TrackCands[best].yresidua[i]);
   //    	aTrack.z.push_back(TrackCands[best].z[i]);	 
	  //    //printf("fill the best tracks  %d %d  %5.3lf, %5.3lf\n",i,trk1,TrackCands[best].xresidua[i],TrackCands[best].yresidua[i]);	  
   //    };
   //     teletracks->tracks.push_back(aTrack);  // put the best track to the MUSEteletracker tree

   //    if (secondbest!=best)  //Hits with only one track gives the same track as the second best track and counts twice. THis get rid of that and fill only the second best track among all the tracks/each event.
   //    {
   //    	trk2=trk2+1;
   //    	aTrack = TrackCands[secondbest];
   //    	aTrack.xresidua = TrackCands[secondbest].xresidua;
   //    	aTrack.yresidua = TrackCands[secondbest].yresidua;
	  //  //printf("fill the second best tracks %d \n",trk2);
	  // // teletracks->tracks.push_back(aTrack);  // put the second best track to the MUSEteletracker tree
   //    };

	///////////////////////////////////////////////////////////////////

      unsigned int thissize = teletracks->tracks.size();
      //printf(" teletracks->tracks.xresidua.size() = %d\n", teletracks->tracks[thissize-1].xresidua.size());

      x.resize(4);
	  y.resize(4);
	  
	  x[0] = aTrack.x0;
      y[0] = aTrack.y0;
      x[1] = x[0] + aTrack.mx * (zgem0[1]-zgem0[0]);
      y[1] = y[0] + aTrack.my * (zgem0[1]-zgem0[0]);
      x[2] = x[0] + aTrack.mx * (zgem0[2]-zgem0[0]);
      y[2] = y[0] + aTrack.my * (zgem0[2]-zgem0[0]);
      x[3] = x[0] + aTrack.mx * (zgem0[3]-zgem0[0]);
      y[3] = y[0] + aTrack.my * (zgem0[3]-zgem0[0]);

	  //Filling in 2D track maps for the best track
	
	  for (unsigned int i=0; i<4; i++){

		  totbest[i]->Fill(x[i],y[i]);

		  multibest[GEMids.size()-2][i]->Fill(x[i],y[i]);


	  }

      //H1(chi2[best], Form("chi2"), Form("chi2 For the Best Track"), 1000., 0., 100.);
      //if (best!=secondbest)
      	//H1(chi2[secondbest], Form("2ndchi2"), Form("2ndchi2 (Chi2 For the Second Best Track"), 1000., 0., 100.);


	   x.clear();
	   y.clear();
	   xifp.clear();
	   yifp.clear();
//If these clears happen, then the output root tree is empty
//duh
//ETHAN
      // aTrack.xresidua.clear();
      // aTrack.yresidua.clear();
      // aTrack.z.clear();
	 
      // TrackCands.clear();
      // chi2.clear();
      // teletracks->tracks.clear();

      return Plugin::ok;
  }


  Long_t MUSEteleTracker::finalize()
  {

  	return Plugin::ok;

  }


  Long_t MUSEteleTracker::cmdline(char *cmd)
  {
  //add cmdline handling here

  return 0; // 0 = all ok
}


extern "C"{
	Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
	{
		return (Plugin *) new MUSEteleTracker(in,out,inf_,outf_,p);
	}
}


ClassImp(MUSEteleTracker);

