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

void make_fit_ntuple(Int_t nrun=1544,Int_t FileID=-2){
 Bool_t CutYtarFlag=kTRUE;
 Bool_t CutYpFpYFpFlag=kTRUE;
 Bool_t CutXpFpXFpFlag=kTRUE;

 gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);

 //  Get info for that optics run
 TString OpticsFile = Form("DATfiles/list_of_optics_run.dat");
 ifstream file_optics(OpticsFile.Data());
 TString opticsline;
 TString OpticsID="";
 Int_t RunNum=0.;
 Double_t CentAngle=0.;
 Int_t SieveFlag=1;
 Double_t ymis =0.0;
 Int_t NumFoil=0;
 TString temp;

 vector <Double_t> ztar_foil;
 Int_t ndelcut;
 vector<Double_t> delcut;

 if (file_optics.is_open()) {
    cout << " Open file = " << OpticsFile << endl;
    while (RunNum!=nrun ) {
        temp.ReadToDelim(file_optics,',');
        cout << temp << endl;
        if (temp.Atoi() == nrun) {
            RunNum = temp.Atoi();
        } else {
            temp.ReadLine(file_optics);
        }
    }
    if (RunNum==nrun) {
        temp.ReadToDelim(file_optics,',');
        OpticsID = temp;
        temp.ReadToDelim(file_optics,',');
        CentAngle = temp.Atof();
        temp.ReadToDelim(file_optics,',');
        NumFoil = temp.Atoi();
        temp.ReadToDelim(file_optics,',');
        SieveFlag = temp.Atoi();
        temp.ReadToDelim(file_optics,',');
        ndelcut = temp.Atoi();
        temp.ReadToDelim(file_optics);
        ymis = temp.Atof();
        for (Int_t nf=0;nf<NumFoil-1;nf++) {
            temp.ReadToDelim(file_optics,',');
            ztar_foil.push_back(temp.Atof());
        }
        temp.ReadToDelim(file_optics);
        ztar_foil.push_back(temp.Atof());
        for (Int_t nd=0;nd<ndelcut;nd++) {
            temp.ReadToDelim(file_optics,',');
            delcut.push_back(temp.Atof());
        }
        temp.ReadToDelim(file_optics);
        delcut.push_back(temp.Atof());
    }
 } else {
    cout << " No file = " << OpticsFile << endl;    
 }

 Double_t y_mis;
 Double_t x_mis;
 
 if (TMath::Abs(CentAngle)<40) {
    y_mis = 0.1*(0.52-0.012*TMath::Abs(CentAngle)+0.002*TMath::Abs(CentAngle)*TMath::Abs(CentAngle)); // cm
 } else {
    y_mis = 0.1*(0.52-0.012*40. + 0.002*40.*40.); // cm
 }
 
 if (TMath::Abs(CentAngle)<50) {
    x_mis = 0.1*(2.37-0.086*TMath::Abs(CentAngle)+0.0012*TMath::Abs(CentAngle)*TMath::Abs(CentAngle));
 } else {
    x_mis = 0.1*(2.37-0.086*50.+0.0012*50.*50.);
 }
 
 cout << RunNum << " " << OpticsID << " " << CentAngle << " " << NumFoil << " " << SieveFlag << " y_mis = " << y_mis << " and x_mis = " << x_mis << endl;
 if (NumFoil==0) return;
 for (Int_t nf=0;nf<NumFoil;nf++) {
    cout << nf << " foil = " << ztar_foil[nf] << endl;
 }
  
 vector <Double_t> ys_cent;
 vector <Double_t> xs_cent;
 for (Int_t nys=0;nys<9;nys++) {
    Double_t ypos=(nys-4)*0.6*2.54;
    ys_cent.push_back(ypos);
    Double_t xpos=(nys-4)*2.54;
    xs_cent.push_back(xpos);
 }
 
 TString inputroot;
 TString outputroot;
 inputroot = Form("/w/hallc-scshelf2102/nps/cploen/ROOTfiles/OPTICS/nps_hms_optics_%s_1_-1.root", OpticsID.Data(),FileID); 
 outputroot = Form("hist/Optics_%s_%d_fit_tree.root",OpticsID.Data(),FileID);
  
 TH1F *hxbpm_tar = new TH1F("hxbpm_tar",Form("Run %d ; Xbpm_tar ; Counts",nrun),100,-2.,2.);
 TH1F *hybpm_tar = new TH1F("hybpm_tar",Form("Run %d ; Ybpm_tar ; Counts",nrun),100,-2.,2.);
 

 TString YtarDeltaCutFile;
 TFile *fYtarDeltaCut;
 vector <TCutG*> ytar_delta_cut;
 if (CutYtarFlag) {
    YtarDeltaCutFile = Form("cuts/ytar_delta_%s_%d_cut.root",OpticsID.Data(),FileID);
    fYtarDeltaCut = new TFile(YtarDeltaCutFile);
    if (!fYtarDeltaCut->IsOpen()) {
        cout << "Error: Could not open file " << YtarDeltaCutFile << endl;
        return;
    }

    cout << " Cut file = " << YtarDeltaCutFile << endl;
    for (Int_t nc=0;nc<NumFoil;nc++) {
        fYtarDeltaCut->cd();
        TCutG* tempcut = (TCutG*)gROOT->FindObject(Form("delta_vs_ytar_cut_foil%d",nc));
        if (tempcut) {
            Int_t npt = tempcut->GetN();
            cout << "hYtarDelta_cut = " << nc << " npts = " << npt << endl;
            ytar_delta_cut.push_back(tempcut);
        } else {
            cout << " No hYtarDelta_cut = " << nc << endl;
            ytar_delta_cut.push_back(nullptr);
        }
    }
 }
 
 TString outCutFile;
 TFile *fcut;
 vector<vector<vector<TCutG*> > > ypfp_yfp_cut;
 vector<vector<vector<Int_t> > > ypfp_yfp_cut_flag;
 ypfp_yfp_cut.resize(NumFoil);
 ypfp_yfp_cut_flag.resize(NumFoil);
 for  (Int_t nf=0;nf<NumFoil;nf++) {
    ypfp_yfp_cut[nf].resize(ndelcut);
    ypfp_yfp_cut_flag[nf].resize(ndelcut);
 }
 if (CutYpFpYFpFlag) {
    outCutFile = Form("cuts/YpFpYFp_%s_%d_cut.root",OpticsID.Data(),FileID);
    fcut = new TFile(outCutFile);
    if (!fcut->IsOpen()) {
        cout << "Error: Could not open file " << outCutFile << endl;
        return;
    }

    cout << " Cut file = " << outCutFile << endl;
    fcut->cd();
    for  (Int_t nf=0;nf<NumFoil;nf++) {
        for  (Int_t nd=0;nd<ndelcut;nd++) {
            for (Int_t nc=0;nc<9;nc++) {
                TCutG* tempg  = (TCutG*)gROOT->FindObject(Form("hYpFpYFp_cut_yscol_%d_nfoil_%d_ndel_%d",nc,nf,nd));
                if (tempg)  {
                    ypfp_yfp_cut[nf][nd].push_back(tempg);
                    ypfp_yfp_cut_flag[nf][nd].push_back(1);
                    Int_t npt = tempg->GetN();
                    cout << "hYpFpYFp_cut_yscol = " << nc << " nfoil = " << nf << " ndel = " << nd << " npts = " << npt << endl;
                } else {
                    ypfp_yfp_cut[nf][nd].push_back(nullptr);
                    ypfp_yfp_cut_flag[nf][nd].push_back(0);
                    cout << " No hYpFpYFp_cut_yscol = " << nc << " nfoil = " << nf << " ndel = " << nd << endl;
                }
            }
        }
    }
 }
 
 TString outCutFilex;
 TFile *fcutx;
 vector<vector<vector<TCutG*> > > xpfp_xfp_cut;
 vector<vector<vector<Int_t> > > xpfp_xfp_cut_flag;
 xpfp_xfp_cut.resize(NumFoil);
 xpfp_xfp_cut_flag.resize(NumFoil);
 for  (Int_t nf=0;nf<NumFoil;nf++) {
    xpfp_xfp_cut[nf].resize(ndelcut);
    xpfp_xfp_cut_flag[nf].resize(ndelcut);
 }
 if (CutXpFpXFpFlag) {
    outCutFilex = Form("cuts/XpFpXFp_%s_%d_cut.root",OpticsID.Data(),FileID);
    fcutx = new TFile(outCutFilex);
    if (!fcutx->IsOpen()) {
        cout << "Error: Could not open file " << outCutFilex << endl;
        return;
    }

    cout << " Cut file = " << outCutFilex << endl;
    fcutx->cd();
    for  (Int_t nf=0;nf<NumFoil;nf++) {
        for  (Int_t nd=0;nd<ndelcut;nd++) {
            for (Int_t nc=0;nc<9;nc++) {
                TCutG* tempgx  = (TCutG*)gROOT->FindObject(Form("hXpFpXFp_cut_yscol_%d_nfoil_%d_ndel_%d",nc,nf,nd));
                if (tempgx)  {
                    xpfp_xfp_cut[nf][nd].push_back(tempgx);
                    xpfp_xfp_cut_flag[nf][nd].push_back(1);
                    Int_t npt = tempgx->GetN();
                    cout << "hXpFpXFp_cut_yscol = " << nc << " nfoil = " << nf << " ndel = " << nd << " npts = " << npt << endl;
                } else {
                    xpfp_xfp_cut[nf][nd].push_back(nullptr);
                    xpfp_xfp_cut_flag[nf][nd].push_back(0);
                    cout << " No hXpFpXFp_cut_yscol = " << nc << " nfoil = " << nf << " ndel = " << nd << endl;
                }
            }
        }
    }
 }

 TFile *infile = new TFile(inputroot.Data());
 if (!infile->IsOpen()) {
    cout << "Error: Could not open input file " << inputroot << endl;
    return;
 }

 cout << " Input file = " << inputroot << endl;
 infile->cd();
 TTree *T = (TTree*)infile->Get("T");
 if (!T) {
    cout << "Error: Could not find tree 'T' in input file" << endl;
    infile->Close();
    return;
 }

 ofstream file_out(Form("info/summary_ntuple_%d.dat",nrun));
 
 T->SetBranchStatus("*",0);
 T->SetBranchStatus("P.gtr.x",1);
 T->SetBranchStatus("P.gtr.y",1);
 T->SetBranchStatus("P.gtr.th",1);
 T->SetBranchStatus("P.gtr.ph",1);
 T->SetBranchStatus("P.gtr.dp",1);
 T->SetBranchStatus("P.bpm.tar.x",1);
 T->SetBranchStatus("P.bpm.tar.y",1);
 T->SetBranchStatus("P.dc.x_fp",1);
 T->SetBranchStatus("P.dc.y_fp",1);
 T->SetBranchStatus("P.dc.xp_fp",1);
 T->SetBranchStatus("P.dc.yp_fp",1);

 Double_t ytart;
 Double_t delta;
 Double_t yfp;
 Double_t xfp;
 Double_t ypfp;
 Double_t xpfp;
 Double_t xbpm;
 Double_t ybpm;

 T->SetBranchAddress("P.gtr.y",&ytart);
 T->SetBranchAddress("P.gtr.dp",&delta);
 T->SetBranchAddress("P.gtr.th",&xpfp);
 T->SetBranchAddress("P.gtr.ph",&ypfp);
 T->SetBranchAddress("P.bpm.tar.x",&xbpm);
 T->SetBranchAddress("P.bpm.tar.y",&ybpm);
 T->SetBranchAddress("P.dc.y_fp",&yfp);
 T->SetBranchAddress("P.dc.x_fp",&xfp);

 TFile *outfile = new TFile(outputroot.Data(),"RECREATE");
 if (!outfile->IsOpen()) {
    cout << "Error: Could not create output file " << outputroot << endl;
    infile->Close();
    return;
 }

 TNtuple *ntuple = new TNtuple("ntuple_optics","HMS optics ntuple","RunNum:FoilID:delta:ytar:xpfp:ypfp:xfp:yfp:xbpm:ybpm");

 Long64_t nentries = T->GetEntries();
 cout << "Number of entries: " << nentries << endl;

 for (Long64_t i=0;i<nentries;i++) {
    T->GetEntry(i);
    hxbpm_tar->Fill(xbpm);
    hybpm_tar->Fill(ybpm);
    Int_t foilid = -1;

    if (CutYtarFlag) {
        for (Int_t nf=0;nf<NumFoil;nf++) {
            if (ytar_delta_cut[nf] && ytar_delta_cut[nf]->IsInside(ytart,delta)) {
                foilid = nf;
                break;
            }
        }
    }

    if (CutYpFpYFpFlag && foilid != -1) {
        Bool_t inside = kFALSE;
        for (Int_t nd=0;nd<ndelcut;nd++) {
            for (Int_t nc=0;nc<9;nc++) {
                if (ypfp_yfp_cut_flag[foiid][nd][nc] && ypfp_yfp_cut[foiid][nd][nc]->IsInside(yfp,ypfp)) {
                    inside = kTRUE;
                    break;
                }
            }
            if (inside) break;
        }
        if (!inside) foilid = -1;
    }

    if (CutXpFpXFpFlag && foilid != -1) {
        Bool_t inside = kFALSE;
        for (Int_t nd=0;nd<ndelcut;nd++) {
            for (Int_t nc=0;nc<9;nc++) {
                if (xpfp_xfp_cut_flag[foiid][nd][nc] && xpfp_xfp_cut[foiid][nd][nc]->IsInside(xfp,xpfp)) {
                    inside = kTRUE;
                    break;
                }
            }
            if (inside) break;
        }
        if (!inside) foilid = -1;
    }

    if (foiid != -1) {
        ntuple->Fill(nrun,foiid,delta,ytart,xpfp,ypfp,xfp,yfp,xbpm,ybpm);
    }
 }

 ntuple->Write();
 hxbpm_tar->Write();
 hybpm_tar->Write();
 outfile->Close();
 infile->Close();
 file_out.close();

 cout << "Finished processing run " << nrun << endl;
}
