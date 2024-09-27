#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include <TH2.h>
#include <TCutG.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
#include <iostream>
#include <fstream>
using namespace std;

void plotElasticSingles(){

// README:
// select which matrix elements were used in replay; this will append to filename
// comment out the kinematic settings not using
// double check that histo 
// double check cuts for x


// TO DO:
// 
// convert to rad?
// pass the kinC setting into the filename automatically
// keep the hardcoded run numbers since the kinematic setting is hardcoded too
// could create a list_of_elastic_runs.dat to take away the need for hardcoding

int debug_mode = 1;
Long64_t nentries= 2000000;
std::string matrixElem="_angular6_67";


//KinC-x50-2
const int nSettings = 3;
int nRun[nSettings]={1534,1535,1536};//run numbers for LH2 Optics scans, KinC-x50-2 at -6.667 GeV/c
double centAng[nSettings]={19.145,17.995,16.850};//deg
std::string kinSetting="KinC-x50-2";

/*
//KinC-x60-3
const int nSettings = 3;
int nRun[nSettings]={1714,1715,1716};//run numbers for LH2 Optics scans, KinC-x50-2 at -6.667 GeV/c
double centAng[nSettings]={22.835,21.685,20.545};//deg
std::string kinSetting="KinC-x60-3";

*/
/*
// KinC-x36-3
const int nSettings = 3;
int nRun[nSettings]={1250,1252,1259};//run numbers for LH2 Optics scans, KinC-x50-2 at -6.667 GeV/c
double centAng[nSettings]={19.26,20.700,22.12};//deg
std::string kinSetting="KinC-x36-3";
*/
/*
// KinC_x25_3' 
const int nSettings = 3;
int nRun[nSettings]={1267,1268,1269};//run numbers for LH2 Optics scans, KinC-x25_3' at -5.639 GeV/c
double centAng[nSettings]={21.27,22.70,24.30};//deg; read from gui not camera.  camera not avail
std::string kinSetting="KinC-x25-3'";
*/

std::string out_file_name=Form("Output/elastic/hist/Singles_");
std::string file_format=".root";
TString outputhist=(out_file_name + kinSetting + matrixElem + file_format+"").c_str();


TString inputroot;
  
//set up plots
TObjArray HList(0);
TH1F *h_kinW[nSettings];
TH1F *h_delta[nSettings];
TH1F *h_yptar[nSettings];
TH1F *h_kinW2[nSettings];
TH1F *h_xfp[nSettings];
TH1F *h_xpfp[nSettings];
TH1F *h_xbcalc[nSettings];
TH1F *h_Q2calc[nSettings];
TH1F *h_nucalc[nSettings];



TH2F *h2_kinWdp[nSettings];
TH2F *h2_kinWyptar[nSettings];
TH2F *h2_kinW2yptar[nSettings];
TH2F *h2_kinW2theta[nSettings];
TH2F *h2_dpyptar[nSettings];
TH2F *h2_dptheta[nSettings];

// if needed, this is where you'd read in matrix elements
//Loop the input files

for (int ii=0; ii<nSettings; ii++){
	inputroot=Form("../ROOTfiles/OPTICS/nps_hms_optics_%d_1_-1.root",nRun[ii]);
	cout << "input root = " << inputroot <<endl;

// Define branches
TFile *fsimc = new TFile(inputroot);
TTree *tsimc = (TTree*) fsimc->Get("T");
 Double_t kinW;
   tsimc->SetBranchAddress("H.kin.W",&kinW);
 Double_t xbcalc;
   tsimc->SetBranchAddress("H.kin.x_bj",&xbcalc);
 Double_t Q2calc;
   tsimc->SetBranchAddress("H.kin.Q2",&Q2calc);
 Double_t nucalc;
   tsimc->SetBranchAddress("H.kin.nu",&nucalc);
 Double_t etracknorm;
   tsimc->SetBranchAddress("H.cal.etracknorm",&etracknorm);
 Double_t  npeSum;
   tsimc->SetBranchAddress("H.cer.npeSum",&npeSum);
// Double_t  sumhgnpe;
//   tsimc->SetBranchAddress("H.hgcer.npeSum",&sumhgnpe);
 Double_t  etottracknorm;
   tsimc->SetBranchAddress("H.cal.etottracknorm",&etottracknorm);
 Double_t  ytar;
   tsimc->SetBranchAddress("H.gtr.y",&ytar);
 Double_t  xtar;
   tsimc->SetBranchAddress("H.gtr.x",&xtar);
 Double_t  reactx;
   tsimc->SetBranchAddress("H.react.x",&reactx);
 Double_t  reacty;
   tsimc->SetBranchAddress("H.react.y",&reacty);
 Double_t  reactz;
   tsimc->SetBranchAddress("H.react.z",&reactz);
 Double_t  delta;
   tsimc->SetBranchAddress("H.gtr.dp",&delta);
 Double_t  yptar;
   tsimc->SetBranchAddress("H.gtr.ph",&yptar);
 Double_t  xptar;
   tsimc->SetBranchAddress("H.gtr.th",&xptar);
 Double_t  yfp;
   tsimc->SetBranchAddress("H.dc.y_fp",&yfp);
 Double_t  ypfp;
   tsimc->SetBranchAddress("H.dc.yp_fp",&ypfp);
 Double_t  xfp;
   tsimc->SetBranchAddress("H.dc.x_fp",&xfp);
 Double_t  xpfp;
   tsimc->SetBranchAddress("H.dc.xp_fp",&xpfp);
 Double_t  ysieve;
   tsimc->SetBranchAddress("H.extcor.ysieve",&ysieve);
 Double_t  xsieve;
   tsimc->SetBranchAddress("H.extcor.xsieve",&xsieve);
 Double_t  xbpm_tar;
   tsimc->SetBranchAddress("H.rb.raster.fr_xbpm_tar",&xbpm_tar);
 Double_t  ybpm_tar;
   tsimc->SetBranchAddress("H.rb.raster.fr_ybpm_tar",&ybpm_tar);

// Make the plot
 h_kinW[ii] = new TH1F(Form("h_kinW_%d",ii),Form("Run %d CentAngle %f; KinW; Counts",nRun[ii],centAng[ii]),100,0.5,1.5);
   HList.Add(h_kinW[ii]);
 h_delta[ii] = new TH1F(Form("h_delta_%d",ii),Form("Run %d CentAngle %f; delta; Counts",nRun[ii],centAng[ii]),100,-15,15);
   HList.Add(h_delta[ii]);
 h_yptar[ii] = new TH1F(Form("h_yptar_%d",ii),Form("Run %d CentAngle %f; yptar (deg); Counts",nRun[ii],centAng[ii]),100,-0.15,0.15);
   HList.Add(h_yptar[ii]);
 h_kinW2[ii] = new TH1F(Form("h_kinW2_%d",ii),Form("Run %d CentAngle %f; KinW2; Counts",nRun[ii],centAng[ii]),100,0.4,1.2);
   HList.Add(h_kinW2[ii]);

 h_xfp[ii] = new TH1F(Form("h_xfp_%d",ii),Form("Run %d CentAngle %f; xfp; Counts", nRun[ii],centAng[ii]),100,-50.,50.);
   HList.Add(h_xfp[ii]);
 h_xpfp[ii] = new TH1F(Form("h_xpfp_%d",ii),Form("Run %d CentAngle %f; xpfp; Counts",nRun[ii],centAng[ii]),100,-0.1,0.1);
   HList.Add(h_xpfp[ii]);
// h_xB[ii] = new TH1F(Form("h_xB_%d",ii),Form("Run %d CentAngle %f; xB; Counts",nRun[ii],centAng[ii]),100,0.9,1.15);
//   HList.Add(h_xB[ii]);
// h_xB_precut[ii] = new TH1F(Form("h_xB_precut_%d",ii),Form("Run %d; xB-precut; Counts",nRun[ii]),100,0.9,1.15);
//   HList.Add(h_xB_precut[ii]);
 h_xbcalc[ii] = new TH1F(Form("h_xbcalc_%d",ii),Form("Run %d; xB-replay; Counts",nRun[ii]),100,09.,1.1);
   HList.Add(h_xbcalc[ii]);
 h_Q2calc[ii] = new TH1F(Form("h_Q2calc_%d",ii),Form("Run %d; Q2-replay; Counts",nRun[ii]),100,2,5.5);
   HList.Add(h_Q2calc[ii]);
 h_nucalc[ii] = new TH1F(Form("h_nucalc_%d",ii),Form("Run %d; nu-replay; Counts",nRun[ii]),100,0.,10);
   HList.Add(h_nucalc[ii]);



 h2_kinWdp[ii] = new TH2F(Form("h2_kinWdp_%d",ii),Form("Run %d CentAngle %f;  KinW; Delta",nRun[ii],centAng[ii]),100,1,3,100,-15,15);
   HList.Add(h2_kinWdp[ii]);
 h2_kinWyptar[ii] = new TH2F(Form("h2_kinWyptar_%d",ii),Form("Run %d CentAngle %f; yptar (deg); kinW",nRun[ii],centAng[ii]),100,-0.15,0.15,100,0.8,1.2);
   HList.Add(h2_kinWyptar[ii]);
 h2_kinW2yptar[ii] = new TH2F(Form("h2_kinW2yptar_%d",ii),Form("Run %d CentAngle %f; yptar (deg); kinW2",nRun[ii],centAng[ii]),100,-0.15,0.15,100,1,1.25);
   HList.Add(h2_kinW2yptar[ii]);
 h2_dpyptar[ii] = new TH2F(Form("h2_dpyptar_%d",ii),Form("Run %d CentAngle %f; yptar (deg); Delta",nRun[ii],centAng[ii]),100, -0.15, 0.15, 100,-15,15);
   HList.Add(h2_dpyptar[ii]);

// ***** CHOOSE *****//
//KinC_x50_2:
  h2_kinW2theta[ii] = new TH2F(Form("h2_kinW2theta_%d",ii),Form("Run %d CentAngle %f; theta_0+yptar (deg); kinW2",nRun[ii],centAng[ii]),100,16.75,19.5,100,0.5,1.25);

//KinC_x60_3: 
//  h2_kinW2theta[ii] = new TH2F(Form("h2_kinW2theta_%d",ii),Form("Run %d CentAngle %f; theta_0+yptar (deg); kinW2",nRun[ii],centAng[ii]),100,20.25,23.,100,0.5,1.25);
//
//KinC_x36_3:
//  h2_kinW2theta[ii] = new TH2F(Form("h2_kinW2theta_%d",ii),Form("Run %d CentAngle %f; theta_0+yptar (deg); kinW2",nRun[ii],centAng[ii]),100,19.,22.5,100,0.5,1.25);
//
//KinC_x25_3':
//  h2_kinW2theta[ii] = new TH2F(Form("h2_kinW2theta_%d",ii),Form("Run %d CentAngle %f; theta_0+yptar (deg); kinW2",nRun[ii],centAng[ii]),100,21.,25,100,0.5,1.25);
 
  HList.Add(h2_kinW2theta[ii]);
// ***** CHOOSE *****//
//KinC_x50_2:
  h2_dptheta[ii] = new TH2F(Form("h2_dptheta_%d",ii),Form("Run %d CentAngle %f; theta_0+yptar (deg); Delta",nRun[ii],centAng[ii]),100, 16.75, 19.5, 100,-15,15);

//KinC_x60_3: 
// h2_dptheta[ii] = new TH2F(Form("h2_dptheta_%d",ii),Form("Run %d CentAngle %f; theta_0+yptar (deg); Delta",nRun[ii],centAng[ii]),100, 20.25, 23., 100,-15,15);

//KinC_x36_3:
//  h2_dptheta[ii] = new TH2F(Form("h2_dptheta_%d",ii),Form("Run %d CentAngle %f; theta_0+yptar (deg); Delta",nRun[ii],centAng[ii]),100, 19., 22.5, 100,-15,15);

//KinC_x25_3':
//   h2_dptheta[ii] = new TH2F(Form("h2_dptheta_%d",ii),Form("Run %d CentAngle %f; theta_0+yptar (deg); Delta",nRun[ii],centAng[ii]),100, 21.,25., 100,-15,15);

  HList.Add(h2_dptheta[ii]);



// choose events (cuts), make initial plots of events
if (debug_mode !=0){nentries = tsimc->GetEntries();}
for (int i = 0; i < nentries; i++) {
	tsimc->GetEntry(i);
	//if (i%50000==0) cout << " Entry = " << i << endl;
 	if (npeSum>6 && etracknorm>0.65 && xbcalc>0.99 && xbcalc<1.015 ){
	// && xbcalc>0.985 && xbcalc<1.015	
		double kinW2 = pow(kinW,2.);
		double theta = centAng[ii]+yptar;
		
		h_kinW[ii]->Fill(kinW);
		h_delta[ii]->Fill(delta);
		h_yptar[ii]->Fill(yptar);
		h_kinW2[ii]->Fill(kinW2);
		h_xbcalc[ii]->Fill(xbcalc);
		h_Q2calc[ii]->Fill(Q2calc);
		h_nucalc[ii]->Fill(nucalc);

		h2_kinWdp[ii]->Fill(kinW,delta);
		h2_kinWyptar[ii]->Fill(yptar,kinW);
		h2_kinW2yptar[ii]->Fill(yptar,kinW2);
		h2_kinW2theta[ii]->Fill(theta,kinW2);
		h2_dpyptar[ii]->Fill(yptar,delta);
		h2_dptheta[ii]->Fill(theta,delta);

	}//cut on npeSum, etracknorm, and kinW
   } //loop over entries
} //loop over runs


// output root file
cout << "making output root file" << endl;
TFile hsimc(outputhist,"RECREATE");
HList.Write();

}
