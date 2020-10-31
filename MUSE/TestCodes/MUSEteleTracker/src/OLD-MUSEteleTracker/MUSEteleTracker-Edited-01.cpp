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


//Let's define the Histograms on the root-tree
Long_t MUSEteleTracker::defineHistograms()
{
    
    h2resmapsx=dH2(Form("MUSEtele/residualmapsx"), Form("X and Y hitmap weighted by X Residuals on MS GEM"), 50,-50,50,50,-50,50);
    h2resmapsy=dH2(Form("MUSEtele/residualmapsy"), Form("X and Y hitmap weighted by Y Residuals on MS GEM"), 50,-50,50,50,-50,50);
    
    std::string GEMnames0[] = {"US","4th","MI","DS"};
    
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
            
            largestremovedx[i][j]=dH1(Form("ResRemoved/X Residual of %s GEM with largest residual removed %s GEM", GEMnames0[i].c_str(), GEMnames0[j].c_str()),"X Residuals;mm", 100,-10,10);
            largestremovedy[i][j]=dH1(Form("ResRemoved/Y Residual of %s GEM with largest residual removed %s GEM", GEMnames0[i].c_str(), GEMnames0[j].c_str()),"Y Residuals;mm", 100,-10,10);
            
        }
        
        largestremovedmx[i]=dH1(Form("ResRemoved/Track slope in x if %s GEM has largest residual and removed", GEMnames0[i].c_str()),"Track slope in x", 100, -1, 1);
        largestremovedmy[i]=dH1(Form("ResRemoved/Track slope in y if %s GEM has largest residual and removed", GEMnames0[i].c_str()),"Track slope in y", 100, -1, 1);
        largestremovedZX[i]=dH1(Form("ResRemoved/ZX intercept if %s GEM has largest residual and removed", GEMnames0[i].c_str()),"ZX intercept", 400,-200,200);
        largestremovedZY[i]=dH1(Form("ResRemoved/ZY intercept if %s GEM has largest residual and removed", GEMnames0[i].c_str()),"ZY intercept", 400,-200,200);
        
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
    double numy = 0;
    double denom = 0;
    double numx = 0;
    //Calculating the slope
    for(int i=0; i<y.size(); i++)
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
    printf(" LumiGEMhits (clusters) @%p\n", clusters);
    
    // create output branch with tracks:
    teletracks = new TeleTracks();
    makeBranch("teletracks", (TObject**)&teletracks);
    printf(" teletracks %p\n", teletracks);
    
    
    
    return Plugin::ok;
}

int event=0;
int multibesttrack=0;

int trk1=0;
int trk2=0;

Long_t MUSEteleTracker::process()
{
    event=event+1;
    
    // vector to store all track candidates for this event:
    std::vector <StraightTrack> TrackCands;
    std::vector < std::vector <StraightTrack> > BestTracks(4, std::vector <StraightTrack>(1));
    
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
    // Here, note that the GEM_ID's for US, 4th, MS, DS are assizned 0,1,2,3
    int combmulti=0, combmulti_cut=0;
    int gemmulti[4] = { 0, 0, 0, 0 };
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
                for (unsigned int g4=0; g4<clusters->hits.size(); g4++)
                {
                    if (clusters->hits[g4].GEMid!=3) continue;
                    gemmulti[3]++;
                };
            };
        };
    };
    
    printf("tele: #clusters US: %d MI: %d DS: %d  4th: %d\n", gemmulti[0], gemmulti[1], gemmulti[2], gemmulti[3]);
    //I think this is what stops tracking form happening on events with multiple hits
    // if (gemmulti[0]*gemmulti[1]*gemmulti[2]*gemmulti[3]>200) {
    // // printf("for the new event  %d, %5.3lf, %5.3lf\n",event,aTrack.mx,aTrack.my);
    // teletracks->tracks.push_back(aTrack);
    // };
    

    //if (gemmulti[0]*gemmulti[1]*gemmulti[2]*gemmulti[3]>200) return Plugin::ok; //return Plugin::ok;
    
    H1(gemmulti[0], Form("MUSEteleTracker/Number of possible clusters US"), Form("Number of possible clusters - US GEM"),
       11,-0.5,10.5);
    H1(gemmulti[1], Form("MUSEteleTracker/Number of possible clusters 4TH"), Form("Number of possible clusters - 4TH GEM"),
       11,-0.5,10.5);
    H1(gemmulti[2], Form("MUSEteleTracker/Number of possible clusters MI"), Form("Number of possible clusters - MI GEM"),
       11,-0.5,10.5);
    //H1(gemmulti[3], Form("MUSEteleTracker/Number of possible clusters DS"), Form("Number of possible clusters - DS GEM"),
    //   11,-0.5,10.5);
    
    
    // The only need is the relative  distances in z here.
    const std::vector<double> zgem0 =
    { //1846.2, 1946.2, 2046.2, 2146.2,   // OLYMPUS, left sector GEMs
        -536.0, -474.0 , -412.0, -350.0}; //Ethan and Ron measured for August 2018 Beamtime DS position
    // -836.0,-747.0,-712.0,-650};//Rough approximation of middle dowel position
    
    double dx, dy;
    double xsloperec =0.0;
    double ysloperec =0.0;
    
    
    std::vector<double> x;
    std::vector<double> y;
    //std::vector<double> xifp;
    //std::vector<double> yifp;
    //std::vector<double> zifp = {0,80};
    
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
    
    for (unsigned int g1=0; g1<clusters->hits.size(); g1++) // Check for GEMids[3] clusters
    {
        // printf("check gem 1  %d %d \n",g1,clusters->hits[g1].GEMid);
        if (GEMids.size() == 4) {
            if (clusters->hits[g1].GEMid!=GEMids[3]) continue;
            x[3]=(clusters->hits[g1].xl*0.4-50.);
            y[3]=(clusters->hits[g1].yl*0.4-50.);
            
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
                
                //  printf("gem 3 before %d %d %5.2lf %5.2lf \n",g3,clusters->hits[g3].GEMid,x[2],y[2]);
                // if (fabs(x[1])>48 || fabs(y[1])>48) continue;
                
                
                for (unsigned int g4=0; g4<clusters->hits.size(); g4++) // Check for GEMids[0] clusters
                {
                    if (clusters->hits[g4].GEMid!=GEMids[0]) continue;
                    x[0]=(clusters->hits[g4].xl*0.4-50.);
                    y[0]=(clusters->hits[g4].yl*0.4-50.);
                    
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
                    
                    
                    //Begin Least Squares line calculations for target GEMs
                    //Averages
                    //////////////////////////////////////////////////////////////////////////
                    //NOTE: This code only finds the first valid hit for each GEM if there is ONLY ONE HIT IN THE GEM
                    // and then attempts to find a least squares fit to those data points.
                    //This means it is very sensitive to the accuracy of the cluster finder
                    /////////////////////////////////////////////////////////////////////////
                    double bzy = 0;
                    double bzx = 0;
                    
                    getLeastSquaresLine(x,y,zgem,slopezy,slopezx,bzy,bzx);
                    
                    
                    double xchi2 = 0;
                    double ychi2 = 0;
                    
                    
                    //residuals calculation
                    std::vector<double> xres;
                    std::vector<double> yres;
                    for(int i=0; i<y.size(); i++)
                    {
                        double residualx = x[i]-(slopezx*zgem[i]+bzx);
                        double residualy = y[i]-(slopezy*zgem[i]+bzy);
                        xres.push_back(residualx);
                        
                        if (firsthit == 1){
                            xchi2 += pow(residualx/0.2662,2);//Non reduced chi-sq, chose 0.2 for error as estimate
                            ychi2 += pow(residualy/0.2844,2);//Non reduced chi-sq, chose 0.2 for error as estimate
                            
                        }
                        
                        else {
                            xchi2 += pow(residualx/0.2662,2);//Non reduced chi-sq, chose 0.2 for error as estimate
                            ychi2 += pow(residualy/0.2844,2);//Non reduced chi-sq, chose 0.2 for error as estimate
                            // xchi2 += pow(residualx/clusters->hits[i].xlerr,2);
                            // ychi2 += pow(residualy/clusters->hits[i].ylerr,2);
                        }
                        
                        
                        yres.push_back(residualy);
                        
                        //Filling individual hit residual histograms
                        
                        if (GEMids.size() != 2){
                            totresx[GEMids[i]]->Fill(residualx);
                            totresy[GEMids[i]]->Fill(residualy);
                        }
                        if (firsthit == 1){
                            firsthitresx[i]->Fill(residualx);
                            firsthitresy[i]->Fill(residualy);
                        }
                        
                        switch(GEMids.size()){
                                
                            case 3:
                                threehitresx[GEMids[i]]->Fill(residualx);
                                threehitresy[GEMids[i]]->Fill(residualy);
                                break;
                                
                            case 4:
                                fourhitresx[GEMids[i]]->Fill(residualx);
                                fourhitresy[GEMids[i]]->Fill(residualy);
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
                        for(int i=0; i<y.size(); i++)
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
                    for(int i=0;i<y.size();i++)
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
                    //  double vec2d=sqrt((dx*dx)+(dy*dy));
                    //   outf<< x[1]<<"\t"<<y[1]<<"\t"<<dx<<"\t"<<dy<<std::endl;
                    
                    h2resmapsx->Fill(x[1],y[1],dx);
                    h2resmapsy->Fill(x[1],y[1],dy);  
                    
                    
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
    
    
    
    if (!goodtrack) {
        //teletracks->tracks.push_back(aTrack);
        // printf("for everythng after  %d, %5.3lf, %5.3lf\n",event,aTrack.mx,aTrack.my);
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

