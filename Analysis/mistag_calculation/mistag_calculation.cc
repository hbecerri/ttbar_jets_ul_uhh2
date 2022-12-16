#include <algorithm>
#include <iterator>
#include <TROOT.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1F.h"
#include "THStack.h"
#include "TRandom.h"
#include "TUnfoldDensity.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFrame.h"
#include "TPaveLabel.h"
#include "TPad.h"
#include "TLegend.h"
#include "TRandom3.h"

void mistag_calculation(string var_name = "", string var_gen = "", string region = "", string year = "2016")
{

    gStyle->SetOptStat(0);

//----obetener_toda_la_informacion_de_entrada--------??
/*
    TChain *chreco_ttbar_semi = new TChain("AnalysisTree","");
    chreco_ttbar_semi->Add("/nfs/dust/cms/user/hugobg/output_UL/muon/uhh2.AnalysisModuleRunner.MC.TTToS*.root/AnalysisTree");
    TTree *treereco_ttbar_semi = (TTree*) chreco_ttbar_semi;

    TChain *chreco_wjets = new TChain("AnalysisTree","");
    chreco_wjets->Add("/nfs/dust/cms/user/hugobg/output_UL/muon/uhh2.AnalysisModuleRunner.MC.WJ*.root/AnalysisTree");
    TTree *treereco_wjets = (TTree*) chreco_wjets;

    TChain *chreco_diboson = new TChain("AnalysisTree","");
    chreco_diboson ->Add("/nfs/dust/cms/user/hugobg/output_UL/muon/uhh2.AnalysisModuleRunner.MC.Dib*.root/AnalysisTree");
    TTree *treereco_diboson  = (TTree*) chreco_diboson ;

    TChain *chreco_st = new TChain("AnalysisTree","");
    chreco_st->Add("/nfs/dust/cms/user/hugobg/output_UL/muon/uhh2.AnalysisModuleRunner.MC.ST*.root/AnalysisTree");
    TTree *treereco_st = (TTree*) chreco_st;

    TChain *chreco_dy = new TChain("AnalysisTree","");
    chreco_dy->Add("/nfs/dust/cms/user/hugobg/output_UL/muon/uhh2.AnalysisModuleRunner.MC.DY*.root/AnalysisTree");
    TTree *treereco_dy = (TTree*) chreco_dy;

    TChain *chreco_DATA = new TChain("AnalysisTree","");
    chreco_DATA->Add("/nfs/dust/cms/user/hugobg/output_UL/muon/uhh2.AnalysisModuleRunner.DATA.DATA_SingleMuon_Run2016.root/AnalysisTree");
    TTree *treereco_DATA = (TTree*) chreco_DATA;

    TH1D *Ttbar   = new TH1D("Ttbar","",25,400,1200);
    TH1D *WJets   = new TH1D("WJets","",25,400,1200);
    TH1D *Diboson = new TH1D("Diboson","",25,400,1200);
    TH1D *ST      = new TH1D("ST","",25,400,1200);
    TH1D *DY      = new TH1D("DY","",25,400,1200);
    TH1D *DATA      = new TH1D("DATA","",25,400,1200);

    TH1D *Ttbar_wtag     = new TH1D("Ttbar_wtag","",25,400,1200);
    TH1D *WJets_wtag     = new TH1D("WJets_wtag","",25,400,1200);
    TH1D *Diboson_wtag   = new TH1D("Diboson_wtag","",25,400,1200);
    TH1D *ST_wtag        = new TH1D("ST_wtag","",25,400,1200);
    TH1D *DY_wtag        = new TH1D("DY_wtag","",25,400,1200);
    TH1D *DATA_wtag      = new TH1D("DATA_wtag","",25,400,1200);

    TH1D *Ttbar_ttag     = new TH1D("Ttbar_ttag","",25,400,1200);
    TH1D *WJets_ttag     = new TH1D("WJets_ttag","",25,400,1200);
    TH1D *Diboson_ttag   = new TH1D("Diboson_ttag","",25,400,1200);
    TH1D *ST_ttag        = new TH1D("ST_ttag","",25,400,1200);
    TH1D *DY_ttag        = new TH1D("DY_ttag","",25,400,1200);
    TH1D *DATA_ttag      = new TH1D("DATA_ttag","",25,400,1200);


    string  weight = "(weight)*(weight_sfmu_id)*(weight_pu)*(weight_sfmu_trigger)*(weight_btag)*(muonrecSF_nominal)";
    string selcuts_boosted_central = "(btagN == 0 && chi2 > 30 && Ak8_ji_eta < 2.4 && Ak8_ji_eta > -2.4)";
    string selcuts_boosted_wtag = "(ttagN == 1 && btagN == 0 && chi2 > 30 && Ak8_ji_eta < 2.4 && Ak8_ji_eta > -2.4)";
    string selcuts_boosted_tag = "(wtagN == 1 && btagN == 0 && chi2 > 30 && Ak8_ji_eta < 2.4 && Ak8_ji_eta > -2.4)";


 
    treereco_DATA->Project("DATA","Ak8_ji_pt",Form("%s",selcuts_boosted_central.c_str()));
    treereco_ttbar_semi->Project("Ttbar","Ak8_ji_pt",Form("%s*%s",weight.c_str(),selcuts_boosted_central.c_str()));
    treereco_wjets->Project("WJets","Ak8_ji_pt",Form("%s*%s",weight.c_str(),selcuts_boosted_central.c_str()));
    treereco_diboson->Project("Diboson","Ak8_ji_pt",Form("%s*%s",weight.c_str(),selcuts_boosted_central.c_str()));
    treereco_st->Project("ST","Ak8_ji_pt",Form("%s*%s",weight.c_str(),selcuts_boosted_central.c_str()));
    treereco_dy->Project("DY","Ak8_ji_pt",Form("%s*%s",weight.c_str(),selcuts_boosted_central.c_str()));

    treereco_DATA->Project("DATA_wtag","Ak8_ji_pt",Form("%s",selcuts_boosted_wtag.c_str()));
    treereco_ttbar_semi->Project("Ttbar_wtag","Ak8_ji_pt",Form("%s*%s",weight.c_str(),selcuts_boosted_wtag.c_str()));
    treereco_wjets->Project("WJets_wtag","Ak8_ji_pt",Form("%s*%s",weight.c_str(),selcuts_boosted_wtag.c_str()));
    treereco_diboson->Project("Diboson_wtag","Ak8_ji_pt",Form("%s*%s",weight.c_str(),selcuts_boosted_wtag.c_str()));
    treereco_st->Project("ST_wtag","Ak8_ji_pt",Form("%s*%s",weight.c_str(),selcuts_boosted_wtag.c_str()));
    treereco_dy->Project("DY_wtag","Ak8_ji_pt",Form("%s*%s",weight.c_str(),selcuts_boosted_wtag.c_str()));

    treereco_DATA->Project("DATA_ttag","Ak8_ji_pt",Form("%s",selcuts_boosted_ttag.c_str()));
    treereco_ttbar_semi->Project("Ttbar_ttag","Ak8_ji_pt",Form("%s*%s",weight.c_str(),selcuts_boosted_ttag.c_str()));
    treereco_wjets->Project("WJets_ttag","Ak8_ji_pt",Form("%s*%s",weight.c_str(),selcuts_boosted_ttag.c_str()));
    treereco_diboson->Project("Diboson_ttag","Ak8_ji_pt",Form("%s*%s",weight.c_str(),selcuts_boosted_ttag.c_str()));
    treereco_st->Project("ST_ttag","Ak8_ji_pt",Form("%s*%s",weight.c_str(),selcuts_boosted_ttag.c_str()));
    treereco_dy->Project("DY_ttag","Ak8_ji_pt",Form("%s*%s",weight.c_str(),selcuts_boosted_ttag.c_str()));


    Float_t eff_data_ttag = (DATA_ttag->Integral() - Ttbar_ttag->Integral() - ST_ttag->Integral())/(DATA->Integral() - Ttbar->Integral() - ST->Integral());
    Float_t eff_data_wtag = (DATA_wtag->Integral() - Ttbar_wtag->Integral() - ST_wtag->Integral())/(DATA->Integral() - Ttbar->Integral() - ST->Integral());

    Float_t eff_MC_ttag   = (WJets_ttag->Integral())/(WJets->Inegral());
    Float_t eff_MC_wtag   = (WJets_wtag->Integral())/(WJets->Inegral());
*/
    ofstream fw("CPlusPlusSampleFile.txt", std::ofstream::out);
    fw << "Hola Mundo" << "\n";
        fw << "Me llamo Hugo" << "\n";

//    TFile out("Input_undfolding_data_.root","recreate");

}

