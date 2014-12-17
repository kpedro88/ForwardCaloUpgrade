#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TPolyMarker.h"
#include "TLine.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TFrame.h"
#include "TMinuit.h"
#include "TGaxis.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <vector>
#include <map>

#define maxHDe 16 //energy points for hadrons
//#define maxJTe 4 //energy points for jets
//#define maxJTe 8 //energy points for jets
#define maxJTe 6 //energy points for jets
#define nCballPar 7 //#pars for cball fn
#define maxZeroWts 3 //different zero layer weight variations

using namespace TMath;

//energy values
Double_t energies[] = {1., 2., 3., 5., 9., 11., 15., 20., 30., 50., 100., 150., 225., 300., 1000., 3000.};
//Double_t jet_energies[] = {20., 50., 100., 500.};
//Double_t jet_energies[] = {5., 10., 20., 50., 100., 250., 500., 750.};
Double_t jet_energies[] = {20., 50., 100., 250., 500., 750.};

class Sample {
	public:
		Sample() : dir(""), fpre(""), gname(""), gname_rat(""), hom_ecal(true), sam_ecal(1), sam(1), sam_init(false), 
		           zeroWt(0.5), zcalWt(0.0), supp_info(""), sampling_term(0.0), sampling_term_err(0.0), constant_term(0.0), constant_term_err(0.0), err_ellipse(0),
				   ecalDepth_mm(0.0), ecalDepth_lam(0.0), zcalDepth_mm(0.0), zcalDepth_lam(0.0), hcalDepth_mm(0.0), hcalDepth_lam(0.0) {}
		std::string dir;
		std::string fpre;
		std::string gname;
		std::string gname_rat;
		bool hom_ecal;
		Double_t sam_ecal;
		Double_t sam;
		bool sam_init;
		Double_t zeroWt;
		Double_t zcalWt;
		std::string supp_info;
		Double_t sampling_term;
		Double_t sampling_term_err;
		Double_t constant_term;
		Double_t constant_term_err;
		TGraph* err_ellipse;
		Double_t ecalDepth_mm, ecalDepth_lam, zcalDepth_mm, zcalDepth_lam, hcalDepth_mm, hcalDepth_lam;
};

map<int,Sample*> sample_map;

void add_sample(int snum, std::string dir_, std::string fpre_, std::string gname_, std::string gname_rat_, bool hom_ecal_, Double_t zeroWt_, Double_t zcalWt_, std::string supp_info_=""){
	Sample* sp = new Sample();
	sp->dir = dir_;
	sp->fpre = fpre_;
	sp->gname = gname_;
	sp->gname_rat = gname_rat_;
	sp->hom_ecal = hom_ecal_;
	sp->zeroWt = zeroWt_;
	sp->zcalWt = zcalWt_;
	if(!supp_info_.empty()) sp->supp_info = supp_info_;
	sample_map[snum] = sp;
}

void set_depths(int snum, Double_t e_mm, Double_t e_lam, Double_t z_mm, Double_t z_lam, Double_t h_mm, Double_t h_lam){
	Sample* sp = sample_map[snum];
	sp->ecalDepth_mm = e_mm;
	sp->ecalDepth_lam = e_lam;
	sp->zcalDepth_mm = z_mm;
	sp->zcalDepth_lam = z_lam;
	sp->hcalDepth_mm = h_mm;
	sp->hcalDepth_lam = h_lam;
}

void init_samples(int v0=0){
	if(v0<0 || v0>2) {
		cout << "No zero weight for version " << v0 << ". Samples will not be loaded." << endl;
		return;
	}
	//default ecal+hcal
	//add_sample(0,"/home/pedrok/CMSSW_4_2_8/src/ForwardCaloUpgrade/GEANT/runs","hcal","PbWO_{4} ECAL","PbWO_{4}",true);
	double w0[] = {0.5,0.49,0.68};
	add_sample(0,"/data/users/pedrok/FullSim/pbwo4","hcal","PbWO_{4} ECAL","PbWO_{4}",true,w0[v0],0.0);
	set_depths(0,220.,220./202.7,0.0,0.0,88.*17,79./153.2*17);
	//wlyso ecal w/ preshower & dead material
	add_sample(1,"/data/users/pedrok/FullSim/wlyso","wlyso","W-LYSO ECAL","W-LYSO",false,0.5,0.0);
	//wlyso ecal w/out preshower
	//add_sample(2,"/data/users/pedrok/FullSim/wlyso","wlyso2","W-LYSO ECAL v2","W-LYSO v2",false);
	add_sample(2,"/data/users/pedrok/FullSim/wlyso","wlyso2","W-LYSO ECAL","W-LYSO",false,0.5,0.0);
	//wlyso ecal w/out preshower & dead material
	//add_sample(3,"/data/users/pedrok/FullSim/wlyso","wlyso3","W-LYSO ECAL v3","W-LYSO v3",false);
	double w3[] = {0.5,0.13,0.18};
	add_sample(3,"/data/users/pedrok/FullSim/wlyso","wlyso3","W-LYSO ECAL","W-LYSO",false,w3[v0],0.0,"(no dead material)");
	set_depths(3,4*28+1.5,(4*28+1.5)/124.1,0.0,0.0,88.*17,79./153.2*17);
	//wlyso ecal, 33 layers, w/out preshower
	add_sample(4,"/data/users/pedrok/FullSim/wlyso","wlyso4","W-LYSO ECAL","W-LYSO",false,0.5,0.0,"(33 layers)");
	//wlyso ecal w/out preshower & dead material, + new dead material: 6 cm Al
	double w5[] = {0.5,0.40,0.33};
	add_sample(5,"/data/users/pedrok/FullSim/wlyso","wlyso5","W-LYSO ECAL","W-LYSO",false,w5[v0],0.0,"(6 cm Al dead mat.)");
	set_depths(5,4*28+1.5,(4*28+1.5)/124.1,0.0,0.0,88.*17,79./153.2*17);
	//wlyso ecal w/out preshower & dead material, + new dead material: 6 cm Al, HE shifted forward
	double w6[] = {0.5,0.33,0.31};
	add_sample(6,"/data/users/pedrok/FullSim/wlyso","wlyso6","W-LYSO ECAL","W-LYSO",false,w6[v0],0.0,"(6 cm Al, HE fwd)");
	set_depths(6,4*28+1.5,(4*28+1.5)/124.1,0.0,0.0,88.*17,79./153.2*17);
	//wlyso ecal w/out preshower & dead material, + new dead material: 6 cm Al, HE shifted forward + 3 extra layers
	double w7[] = {0.5,0.42,0.33};
	add_sample(7,"/data/users/pedrok/FullSim/wlyso","wlyso7","W-LYSO ECAL","W-LYSO",false,w7[v0],0.0,"(6 cm Al, HE + 3 lay)"); //previously "(Al, HE fwd + 3 lay)"
	set_depths(7,4*28+1.5,(4*28+1.5)/124.1,0.0,0.0,88.*20,79./153.2*20);
	//wlyso ecal w/out preshower & dead material, + new dead material: 6 cm Al, HE - 3 layers
	double w8[] = {0.5,0.27,0.27};
	add_sample(8,"/data/users/pedrok/FullSim/wlyso","wlyso_zcal0_hcal14","W-LYSO ECAL","W-LYSO",false,w8[v0],0.0,"(Al, HE -3 lay)");
	set_depths(8,4*28+1.5,(4*28+1.5)/124.1,0.0,0.0,88.*14,79./153.2*14);
	//wlyso ecal w/out preshower & dead material, + new dead material: 6 cm Al, HE extension w/ 4 layers of 4.9mm W & 3.7 mm scint
	//double w71[] = {0.5,1.4,1.6};
	//add_sample(71,"/data/users/pedrok/FullSim/wlyso","wlyso_z3cal","W-LYSO ECAL","W-LYSO",false,w71[v0],1.0,"(Al, ext. HE 2#lambda_{0}(W))");
	//wlyso ecal w/out preshower & dead material, plus extended HE
	double w31[] = {0.25,0.19,0.25,
	                0.25,0.14,0.25,
					0.25,0.23,0.25,
					0.25,0.25,0.25,
					0.25,0.24,0.25,
					0.25,0.20,0.25,
					0.25,0.11,0.25,
					0.25,0.19,0.25,
					0.25,0.18,0.25,
					0.25,0.16,0.25};
	double w32[] = {0.50,0.70,0.50,
	                0.50,0.58,0.50,
					0.50,0.53,0.50,
					0.50,0.51,0.50,
					0.50,0.50,0.50,
					0.50,0.50,0.50,
					0.50,0.49,0.50,
					0.50,0.48,0.50,
					0.50,0.47,0.50,
					0.50,0.46,0.50};
	for(int i = 1; i <= 10; i++){
		stringstream sn;
		sn << "wlyso_zcal" << i;
		stringstream ss;
		ss.precision(2);
		ss << "(ext. HE, " << i*1.8/9.95 << "#lambda_{0}(W))";
		add_sample(i+30,"/data/users/pedrok/FullSim/wlyso",sn.str(),"W-LYSO ECAL","W-LYSO",false,w31[(i-1)*3+v0],w32[(i-1)*3+v0],ss.str());
		set_depths(i+30,4*28+1.5,(4*28+1.5)/124.1,(double)i*18,(double)i*18/99.5,88.*17,79./153.2*17);
	}
	//wlyso ecal w/out preshower & dead material, plus extended HE, with HCAL only 14 layers (back 3 removed)
	double w41[] = {0.25,0.11,0.25,
	                0.25,0.19,0.25,
					0.25,0.2,0.25,
					0.25,0.19,0.25,
					0.25,0.23,0.25,
					0.25,0.24,0.25,
					0.25,0.18,0.25,
					0.25,0.16,0.25,
					0.25,0.19,0.25,
					0.25,0.22,0.25};
	double w42[] = {0.50,0.67,0.50,
	                0.50,0.55,0.50,
					0.50,0.52,0.50,
					0.50,0.50,0.50,
					0.50,0.48,0.50,
					0.50,0.47,0.50,
					0.50,0.47,0.50,
					0.50,0.46,0.50,
					0.50,0.45,0.50,
					0.50,0.45,0.50};
	for(int i = 1; i <= 10; i++){
		stringstream sn;
		sn << "wlyso_zcal" << i << "_hcal14";
		stringstream ss;
		ss.precision(2);
		ss << "(ext. HE, " << i*1.8/9.95 << "#lambda_{0}, H14)";
		add_sample(i+40,"/data/users/pedrok/FullSim/wlyso",sn.str(),"W-LYSO ECAL","W-LYSO",false,w41[(i-1)*3+v0],w42[(i-1)*3+v0],ss.str());
		set_depths(i+40,4*28+1.5,(4*28+1.5)/124.1,(double)i*18,(double)i*18/99.5,88.*14,79./153.2*14);
	}	
//	for(int i = 1; i <= 5; i++){
//		stringstream sn;
//		sn << "wlyso_z2cal" << i;
//		stringstream ss;
//		ss.precision(2);
//		ss << "(extended HE, " << i*4.8/9.95 << "#lambda_{0}(W))";
//		add_sample(i+40,"/data/users/pedrok/FullSim/wlyso",sn.str(),"W-LYSO ECAL","W-LYSO",false,1.0,ss.str());		
//	}
	//wlyso ecal w/out preshower & dead material, + new dead material: 6 cm Al, Tungsten HE 18 layers (physics list QGSP_BERT_EMV)
	double w9[] = {0.5,0.37,0.37};
	add_sample(9,"/data/users/pedrok/FullSim/wlyso/qgsp_bert_emv","wlyso_whcal","W-LYSO ECAL","W-LYSO",false,w9[v0],0.0,"(WHCAL, Q_B_E)");
	set_depths(9,4*28+1.5,(4*28+1.5)/124.1,0.0,0.0,58.7*18,49.7/99.46*18);
	//wlyso ecal w/out preshower & dead material, + new dead material: 6 cm Al, Tungsten HE 18 layers (physics list QGSP_BERT_HP)
	double w10[] = {0.5,0.35,0.35};
	add_sample(10,"/data/users/pedrok/FullSim/wlyso/qgsp_bert_hp","wlyso_whcal","W-LYSO ECAL","W-LYSO",false,w10[v0],0.0,"(WHCAL, Q_B_H)");
	set_depths(10,4*28+1.5,(4*28+1.5)/124.1,0.0,0.0,58.7*18,49.7/99.46*18);
	//default ecal+hcal, no dead material
	double w100[] = {0.5,0.3,0.2};
	add_sample(100,"/data/users/pedrok/FullSim/pbwo4","hcal_nodead","PbWO_{4} ECAL","PbWO_{4}",true,w100[v0],0.0,"(no dead material)");
	set_depths(100,220.,220./202.7,0.0,0.0,88.*17,79./153.2*17);
	//default ecal+hcal, liquid tile nonuniformity
	add_sample(1000,"/data/users/pedrok/FullSim/pbwo4","hcal_nonuni","PbWO_{4} ECAL","nonuni",true,w0[v0],0.0,"(w/ nonuniformity)");
	set_depths(1000,220.,220./202.7,0.0,0.0,88.*17,79./153.2*17);
	//default ecal+hcal, step fn nonuniformity
	add_sample(1001,"/data/users/pedrok/FullSim/pbwo4","hcal_nonuni2","PbWO_{4} ECAL","nonuni",true,w0[v0],0.0,"(w/ step fn nonuni)");
	set_depths(1001,220.,220./202.7,0.0,0.0,88.*17,79./153.2*17);
	//default ecal+hcal, 20% light yield
	add_sample(1002,"/data/users/pedrok/FullSim/pbwo4","hcal_nonuni3","PbWO_{4} ECAL","HCAL 20%",true,w0[v0],0.0,"(HCAL 20% LY)");
	set_depths(1002,220.,220./202.7,0.0,0.0,88.*17,79./153.2*17);		
}

std::string pdir = "plots";
std::string pformat = "png";

void set_pinfo(std::string pd, std::string pf) { pdir = pd; pformat = pf; }

void set_zeroWt(int snum, Double_t z) { Sample* sp = sample_map[snum]; sp->zeroWt = z; }
void set_zcalWt(int snum, Double_t z) { Sample* sp = sample_map[snum]; sp->zcalWt = z; }

//-------------------------------------
//class to store energy resolution data
class energyRes {
	//vars
	private:
		Double_t energy;
		Int_t imip;
		std::vector<Double_t> stats; //stats are N, mean, RMS
		std::vector<Double_t> stat_errs;
		TF1 *fit;
		TH1F *hist;

	public:
		//constructors
		energyRes(Double_t en, Int_t im) : energy(en), imip(im), fit(NULL), hist(NULL) {}
		
		//set members
		//0: Nevents, 1: mean, 2: stdev
		void setStats(std::vector<Double_t>& _stats, std::vector<Double_t>& _stat_errs){
			stats = _stats;
			stat_errs = _stat_errs;
		}
		void setFit(TF1* _fit) { fit = _fit; }
		void setEnergy(Double_t en) { energy = en; }
		void setMip(Int_t im) { imip = im; }
		void setHist(TH1F* _hist) { hist = _hist; }
	
		//access members
		Double_t getEnergy() { return energy; }
		Int_t getMip() { return imip; }
		Double_t getStat(int istat) { return stats[istat]; }
		Double_t getStatErr(int istat) { return stat_errs[istat]; }
		TF1* getFit() { return fit; }
		TH1F* getHist() { return hist; }

};

//------------------------------------------
//Double-sided Crystal Ball function
//parameters:
//N, mu, sigma, aL, nL, aR, nR
//0,  1,     2,  3,  4,  5,  6
Double_t cball(Double_t *x, Double_t *par){
	//ensure sigma > 0 and a > 0
	//Double_t N = 1/(sigma*(n/a*1/(n-1)*Exp(-a*a/2) + Sqrt(Pi()/2)*(1+Erf(a/Sqrt(2))))); //normalized N
	Double_t N = par[0]; //let N float
	Double_t mu = par[1];
	par[2] = Abs(par[2]);
	Double_t sigma = par[2];
	par[3] = Abs(par[3]);
	Double_t aL = par[3];
	par[4] = (par[4]>1) ? par[4] : 1.01; //n>1 required
	Double_t nL = par[4];
	par[5] = Abs(par[5]);
	Double_t aR = par[5];
	par[6] = (par[6]>1) ? par[6] : 1.01; //n>1 required
	Double_t nR = par[6];	
	Double_t arg = (x[0]-mu)/sigma;
	
	//left tail
	if(arg <= -aL){
		return N*Power(nL/aL,nL)*Exp(-Power(aL,2)/2)*Power(nL/aL-aL-arg,-nL);
	}
	//right tail
	else if(arg >= aR){
		return N*Power(nR/aR,nR)*Exp(-Power(aR,2)/2)*Power(nR/aR-aR+arg,-nR);
	}
	//core
	else{
		return N*Exp(-Power(arg,2)/2);
	}

}

//------------------------------------------
//function to calculate ecal sampling factor
Double_t g4_ecal_sample(int snum, Double_t energy, bool do_show, bool do_print=false){

	Sample* sp = sample_map[snum];
	if(!sp) { std::cout << "Sample " << snum << " is not loaded." << std::endl; return 0; }

	//make filenames
	std::stringstream drawname, fname, oname, piname;
	fname << sp->dir << "/" << sp->fpre << "_elec_" << energy << "gev_10k.root";
	oname << pdir << "/" << sp->fpre << "_ecal_sample_" << energy << "gev_peak.png";
	piname << "e- " << energy << " GeV";

	//open file and tree
	TFile* _file;
	_file = TFile::Open((fname.str()).c_str());
	TTree* totalTree = (TTree*)_file->Get("Total");

	//get histo from tree (no display)
	std::string cutname;
	//drawname << "(ecal)/1000>>hsam(200,0," << energy/5 << ")";
	drawname << "(ecal)/1000>>hsam";
	cutname = "";
	totalTree->Draw((drawname.str()).c_str(),cutname.c_str(),"hist goff");
	TH1F* hsam = (TH1F*)gDirectory->Get("hsam");
	
	//use parameters from histo to start fit
	TSpectrum* spec = new TSpectrum(5);
	spec->Search(hsam,5,"nodraw goff");
	Float_t* xpos = spec->GetPositionX();
	Float_t* ypos = spec->GetPositionY();

	Double_t m = xpos[0];
	Double_t me = hsam->GetMeanError();
	Double_t N = hsam->GetEntries();
	std::stringstream s_mean;
	Double_t f = energy/m;
	s_mean << f;

	TPolyMarker* pm = new TPolyMarker(1, xpos, ypos);
	hsam->GetListOfFunctions()->Add(pm);
	pm->SetMarkerStyle(23);
	pm->SetMarkerColor(kRed);
	pm->SetMarkerSize(1.3);

	//plotting and printing
	if (do_show){
		TCanvas* can = new TCanvas("sample","sample",700,500);
		can->cd();
		
		hsam->SetTitle("");
		hsam->GetXaxis()->SetTitle("Energy [GeV]");
		hsam->Draw();
		
		std::stringstream Nname;
		Nname << "N = " << N;
		
		//determine placing of pave
		Double_t xmin;
		if (m/((hsam->GetXaxis()->GetXmax() + hsam->GetXaxis()->GetXmin())/2) < 1) xmin = 0.7;
		else xmin = 0.2;
		
		//legend
		TPaveText *pave = new TPaveText(xmin,0.65,xmin+0.2,0.9,"NDC");
		pave->AddText((piname.str()).c_str());
		pave->AddText((Nname.str()).c_str());
		pave->AddText("Peak sampling factor:");
		pave->AddText((s_mean.str()).c_str());
		pave->SetFillColor(0);
		pave->SetBorderSize(0);
		pave->SetTextFont(42);
		pave->SetTextSize(0.04);
		pave->Draw("same");

		std::cout.precision(6);
		std::cout << "f = " << f << " +/- " << energy*me/(m*m) << std::endl;

		if(do_print) can->Print((oname.str()).c_str(),"png");
	}
	else _file->Close();

	return f;
}

//--------------------------------------
//function to calculate sampling factors
Double_t g4_sample(int snum, int num, Double_t sam_ecal, bool do_show, bool do_print=false, bool var_range=false){
	if (num>=maxHDe || num<0) {
		std::cout << "num must be between 0 and " << maxHDe - 1 << std::endl;
		return 0;
	}
	
	Sample* sp = sample_map[snum];
	if(!sp) { std::cout << "Sample " << snum << " is not loaded." << std::endl; return 0; }

	//make filenames
	std::stringstream drawname, cutname, fname, oname, piname;
	fname << sp->dir << "/" << sp->fpre << "_pion_" << energies[num] << "gev_10k.root";
	oname << pdir << "/" << sp->fpre << "_sample_" << energies[num] << "gev_peak.png";
	piname << "#pi^{-} " << energies[num] << " GeV";

	//open file and tree
	TFile* _file;
	_file = TFile::Open((fname.str()).c_str());
	TTree* totalTree = (TTree*)_file->Get("Total");

	//get histo from tree (no display)
	//define mip as sam_ecal*ecal < 1 gev = 1000 mev
	if(var_range) {
		drawname << "(hcal+";
		if(sp->zcalWt>0) drawname << sp->zcalWt << "*zcal+";
		drawname << sp->zeroWt << "*zero)/1000>>hsam(200)";
	}
	else {
		drawname << "(hcal+";
		if(sp->zcalWt>0) drawname << sp->zcalWt << "*zcal+";
		drawname << sp->zeroWt << "*zero)/1000>>hsam(200,0," << energies[num]/100 << ")";
	}
	cutname << "(" << sam_ecal << "*ecal)<1000";
	totalTree->Draw((drawname.str()).c_str(),(cutname.str()).c_str(),"hist goff");
	TH1F* hsam = (TH1F*)gDirectory->Get("hsam");
	
	//use parameters from histo to start fit
	TSpectrum* spec = new TSpectrum(5);
	spec->Search(hsam,5,"nodraw goff");
	Float_t* xpos = spec->GetPositionX();
	Float_t* ypos = spec->GetPositionY();

	Double_t m = xpos[0];
	Double_t me = hsam->GetMeanError();
	Double_t N = hsam->GetEntries();
	std::stringstream s_mean;
	s_mean.precision(3);
	Double_t f = energies[num]/m;
	s_mean << f << " #pm " << energies[num]*(me/(m*m));

	TPolyMarker* pm = new TPolyMarker(1, xpos, ypos);
	hsam->GetListOfFunctions()->Add(pm);
	pm->SetMarkerStyle(23);
	pm->SetMarkerColor(kRed);
	pm->SetMarkerSize(1.3);

	//plotting and printing
	if (do_show){
		TCanvas* can = new TCanvas("sample","sample",700,500);
		can->cd();
		TPad* pad = new TPad("graph","",0,0,1,1);
		pad->SetMargin(0.12,0.05,0.15,0.05);
		pad->Draw();
		pad->cd();
		
		//formatting
		hsam->SetTitle("");
		hsam->GetXaxis()->SetTitle("E_{hcal} [GeV]");
		//hsam->SetStats(kTRUE);
		//gStyle->SetOptStat("mr");
		hsam->SetLineWidth(2);
		hsam->SetLineColor(kBlack);
		hsam->GetYaxis()->SetTitleSize(32/(pad->GetWh()*pad->GetAbsHNDC()));
		hsam->GetYaxis()->SetLabelSize(28/(pad->GetWh()*pad->GetAbsHNDC()));
		hsam->GetXaxis()->SetTitleSize(32/(pad->GetWh()*pad->GetAbsHNDC()));
		hsam->GetXaxis()->SetLabelSize(28/(pad->GetWh()*pad->GetAbsHNDC()));
		hsam->GetYaxis()->SetTickLength(12/(pad->GetWh()*pad->GetAbsHNDC()));
		hsam->GetXaxis()->SetTickLength(12/(pad->GetWh()*pad->GetAbsHNDC()));
		
		hsam->Draw();
		
		std::stringstream Nname;
		Nname << "N = " << N;
		
		//determine placing of pave
		Double_t xmin;
		if (m/((hsam->GetXaxis()->GetXmax() + hsam->GetXaxis()->GetXmin())/2) < 1) xmin = 0.65;
		else xmin = 0.2;
		
		//legend
		TPaveText *pave = new TPaveText(xmin,0.65,xmin+0.2,0.85,"NDC");
		pave->AddText((piname.str()).c_str());
		pave->AddText((Nname.str()).c_str());
		pave->AddText("Peak sampling factor:");
		pave->AddText((s_mean.str()).c_str());
		pave->SetFillColor(0);
		pave->SetBorderSize(0);
		pave->SetTextFont(42);
		pave->SetTextSize(0.05);
		pave->Draw("same");

		std::cout.precision(6);
		std::cout << "f = " << f << " +/- " << energies[num]*(me/(m*m)) << std::endl;

		if(do_print) can->Print((oname.str()).c_str(),"png");
	}
	else _file->Close();

	return f;
}

//---------------------------------------------
//function to get sampling factors for a sample
void get_sampling_factors(int snum, bool do_show = false){
	Sample* sp = sample_map[snum];
	if(!sp) { std::cout << "Sample " << snum << " is not loaded." << std::endl; return; }
	if(!sp->sam_init){
		//ecal sampling factor
		if(sp->hom_ecal) sp->sam_ecal = 1; //homogeneous ecal, no sampling factor needed
		else sp->sam_ecal = g4_ecal_sample(snum,50,do_show); //get sampling factor from 50 gev electrons

		//get sampling factor from 50 gev pions
		//zeroWt = 0.5 is standard
		sp->sam = g4_sample(snum,9,sp->sam_ecal,do_show);
		//std::cout << sam << std::endl;
		sp->sam_init = true;
	}
}

//----------------------------------
//function to fit energy resolutions
energyRes* get_res(int snum, int num, int imip, bool do_fit, bool do_jet, bool do_show, bool do_print=false, bool do_batch=false){
	int maxE = do_jet ? maxJTe : maxHDe;

	if (num>=maxE || num<0) { std::cout << "num must be between 0 and " << maxE - 1 << std::endl; energyRes* theRes = new energyRes(0,0); return theRes; }

	double energy = 0;
	if(do_jet) energy = jet_energies[num];
	else energy = energies[num];
	
	if (imip>2 || imip<0) { std::cout << "imip must be between 0 and 2" << std::endl; energyRes* theRes = new energyRes(0,0); return theRes; }
	
	if(do_jet) imip = 2; //no mip options for jets

	Sample* sp = sample_map[snum];
	if(!sp) { std::cout << "Sample " << snum << " is not loaded." << std::endl; energyRes* theRes = new energyRes(0,0); return theRes; }
	
	//make filenames
	std::string rtype = "pion";
	if(do_jet) rtype = "monojet";
	std::stringstream fname, piname;
	fname << sp->dir << "/" << sp->fpre << "_" << rtype << "_" << energy << "gev_10k.root";
	if(do_jet) piname << "d jet " << energy << " GeV";
	else piname << "#pi^{-} " << energy << " GeV";

	//open file and tree
	TFile* _file;
	_file = TFile::Open((fname.str()).c_str());
	TTree* totalTree = (TTree*)_file->Get("Total");

	//make tree drawing expressions
	//define mip as ecal < 1 gev = 1000 mev
	double mipcut = 0.8;
	std::stringstream drawname;
	std::stringstream cutname;
	std::stringstream etacut;
	//default histo settings
	//double Emin = 0.1*energies[num]; //lower cut to remove weird peaks near E=zero
	double Emin = 0;
	double Emax = 3*energy;
	if(energy>1000. || do_jet) Emax = 2*energy;
	int nbins = 100;
	
	//ecal & hcal energies are already "calibrated" at sim level (sampling factors for HBHEHO, poisson smearing & pe conversion for HF)
	get_sampling_factors(snum);
	drawname << "(" << sp->sam_ecal << "*ecal+" << sp->sam << "*(hcal+";
	if(sp->zcalWt>0) drawname << sp->zcalWt << "*zcal+";
	drawname << sp->zeroWt << "*zero))/1000>>htemp(" << nbins << "," << Emin << "," << Emax << ")";
	//0 is mip, 1 is nomip, 2 is total
	if(imip==0) cutname << sp->sam_ecal << "*ecal/1000 < " << mipcut;
	else if(imip==1) cutname << sp->sam_ecal << "*ecal/1000 > " << mipcut;
	else cutname << "";

	TH1F* h_res; //to store histos drawn from tree
	TF1* gfit;
	TF1* gsnL;
	TF1* gsnR;

	//plotting variables
	TCanvas* can;
	TPad* pad;
	TLegend* leg;
	TPaveText* pave;
	TPaveText* pave_par;
	TLine *aLline;
	TLine *aRline;

	//create instance of energyRes object
	energyRes* theRes = new energyRes(energy,imip);

	//draw w/ appropriate cut
	totalTree->Draw((drawname.str()).c_str(),(cutname.str()).c_str(),"hist goff");
	h_res = (TH1F*)gDirectory->Get("htemp");
	h_res->SetTitle("");
	h_res->GetXaxis()->SetTitle("Energy [GeV]");
	h_res->SetLineWidth(2);
	h_res->SetLineColor(kBlack);

	//names
	std::string omip, ofit;
	if (imip==0) omip = "mip";
	else if (imip==1) omip = "nomip";
	else omip = "tot";
	//if(do_fit) ofit = "cballD";
	if(do_fit) ofit = "fit";
	else ofit = "nofit";
	std::stringstream oname;
	oname << pdir << "/" << rtype << "_response_" << sp->fpre << "_" << omip << "_" << ofit << "_" << energy << "gev";
	
	//get values from histo
	Double_t m = h_res->GetMean();
	Double_t me = h_res->GetMeanError();
	//Double_t m = h_res->GetBinCenter(h_res->GetMaximumBin()); //peak
	Double_t s = h_res->GetRMS();
	Double_t se = h_res->GetRMSError();
	Int_t N = h_res->GetEntries();
	
	std::vector<Double_t> stats(3,0);
	std::vector<Double_t> stat_err(3,0);
	stats[0] = N;
	stat_err[0] = 0;
	stats[1] = m;
	stat_err[1] = me;
	stats[2] = s;
	stat_err[2] = se;

	//find peak
	TSpectrum *spec = new TSpectrum(5);
	if(nbins < 100) spec->Search(h_res,6,"nobackground nodraw goff"); //turn off background removal when nbins too small
	else spec->Search(h_res,6,"nodraw goff");
	Float_t* xpos = spec->GetPositionX();
	Float_t* ypos = spec->GetPositionY();
	Double_t p = xpos[0];
	Double_t ph = ypos[0];
	if(do_show) std::cout << "peak: " << p << std::endl;
	
	//setup fitting function & do fit
	if (do_fit){
		/*gfit = new TF1((oname.str()).c_str(),cball,Emin,Emax,nCballPar);
		gfit->SetParameters(ph,p,s,1,1.1,1,1.1);
		
		//limits on parameters: 0 < a < 10, 1 < n < 200
		gfit->SetParLimits(3,0,10);
		gfit->SetParLimits(5,0,10);
		gfit->SetParLimits(4,1.01,200);
		gfit->SetParLimits(6,1.01,200);*/
		
		gfit = new TF1(omip.c_str(),"gaus",0,h_res->GetXaxis()->GetXmax());
		if(do_jet){
			gfit->SetParameters(ph,p,s);
			if(m > p) gfit->SetRange(p-1.5*s,p+1.0*s); //high tail
			else gfit->SetRange(p-1.0*s,p+1.5*s); //low tail
		}
		else{
			gfit->SetParameters((Double_t)N,m,s);
			gfit->SetRange(m-2*s,m+1*s); //fit within 2 std devs
		}
		
		//formatting
		gfit->SetLineColor(kRed);
		gfit->SetMarkerColor(kRed);
		gfit->SetLineWidth(2);
		//fit
		h_res->Fit(gfit,"LNQR");
	}
	
	//store parameters
	theRes->setStats(stats,stat_err);
	if(do_fit) theRes->setFit(gfit);
	//store histo
	h_res->SetDirectory(0);
	theRes->setHist(h_res);
	
	std::stringstream muname, signame, musigname, aLname, nLname, aRname, nRname, Nname, chiname;
	muname.precision(2);
	signame.precision(2);
	musigname.precision(3);
	aLname.precision(2);
	nLname.precision(2);
	aRname.precision(2);
	nRname.precision(2);
	chiname.precision(5);
	if (do_fit) {
		muname << fixed << "#mu = " << gfit->GetParameter(1) << " #pm " << gfit->GetParError(1);
		signame << fixed << "#sigma = " << gfit->GetParameter(2) << " #pm " << gfit->GetParError(2);
		musigname << fixed << "#sigma/#mu = " << gfit->GetParameter(2)/gfit->GetParameter(1) << " #pm " << 
		gfit->GetParameter(2)/gfit->GetParameter(1) * sqrt( Power(gfit->GetParError(1),2)/Power(gfit->GetParameter(1),2) +  Power(gfit->GetParError(2),2)/Power(gfit->GetParameter(2),2) );
		//aLname << fixed << "a_{L} = " << gfit->GetParameter(3) << " #pm " << gfit->GetParError(3);
		//nLname << fixed << "n_{L} = " << gfit->GetParameter(4) << " #pm " << gfit->GetParError(4);
		//aRname << fixed << "a_{R} = " << gfit->GetParameter(5) << " #pm " << gfit->GetParError(5);
		//nRname << fixed << "n_{R} = " << gfit->GetParameter(6) << " #pm " << gfit->GetParError(6);
		chiname << fixed << "#chi^{2}/ndf = " << gfit->GetChisquare()/gfit->GetNDF();
	}
	else {
		muname << fixed << "Mean = " << m << " #pm " << me;
		signame << fixed << "RMS = " << s << " #pm " << se;
		musigname << fixed << "RMS/Mean = " << s/m << " #pm " << s/m*sqrt((me*me)/(m*m)+(se*se)/(s*s));
	}
	Nname << "N = " << N; 

	//plotting
	if (do_show){
		can = new TCanvas(omip.c_str(),omip.c_str(),700,500);
		can->cd();
		pad = new TPad("graph","",0,0,1,1);
		pad->SetMargin(0.12,0.05,0.15,0.05);
		pad->Draw();
		pad->cd();
		
		//formatting
		h_res->SetStats(kTRUE);
		gStyle->SetOptStat("mr");
		h_res->GetYaxis()->SetTitleSize(32/(pad->GetWh()*pad->GetAbsHNDC()));
		h_res->GetYaxis()->SetLabelSize(28/(pad->GetWh()*pad->GetAbsHNDC()));
		h_res->GetXaxis()->SetTitleSize(32/(pad->GetWh()*pad->GetAbsHNDC()));
		h_res->GetXaxis()->SetLabelSize(28/(pad->GetWh()*pad->GetAbsHNDC()));
		h_res->GetYaxis()->SetTickLength(12/(pad->GetWh()*pad->GetAbsHNDC()));
		h_res->GetXaxis()->SetTickLength(12/(pad->GetWh()*pad->GetAbsHNDC()));
		
		//plot histo and fit
		h_res->Draw("hist");
		if(do_fit) gfit->Draw("same");	
	
		//determine placing of legend and paves - par pave goes on side with more space
		Double_t xmin;
		if (m/((h_res->GetXaxis()->GetXmax() + h_res->GetXaxis()->GetXmin())/2) < 1) xmin = 0.65;
		else xmin = 0.2;
	
		if(do_fit) { //legend
			leg = new TLegend(xmin,0.78,xmin+0.2,0.88);
			leg->AddEntry(h_res,"Standalone");
			leg->AddEntry(gfit,"Fit");
			leg->SetFillColor(0);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.05);
			leg->SetTextFont(42);
			leg->Draw("same");
			
			can->Update();
			/*
			//left line
			Double_t bndL = gfit->GetParameter(1) - gfit->GetParameter(2)*gfit->GetParameter(3);
			aLline = new TLine(bndL,pad->GetUymin(),bndL,pad->GetUymax());
			aLline->SetLineStyle(2);
			aLline->SetLineWidth(3);
			aLline->SetLineColor(kBlue);
			aLline->Draw("same");
			
			//left gaussian
			gsnL = new TF1("gsn","gaus",Emin,bndL);
			gsnL->SetParameters(gfit->GetParameter(0),gfit->GetParameter(1),gfit->GetParameter(2));
			gsnL->SetLineColor(kRed);
			gsnL->SetMarkerColor(kRed);
			gsnL->SetLineWidth(2);
			gsnL->SetLineStyle(2);
			gsnL->Draw("same");

			//line
			Double_t bndR = gfit->GetParameter(1) + gfit->GetParameter(2)*gfit->GetParameter(5);
			aRline = new TLine(bndR,pad->GetUymin(),bndR,pad->GetUymax());
			aRline->SetLineStyle(2);
			aRline->SetLineWidth(3);
			aRline->SetLineColor(kBlue);
			aRline->Draw("same");
			
			//right gaussian
			gsnR = new TF1("gsn","gaus",bndR,Emax);
			gsnR->SetParameters(gfit->GetParameter(0),gfit->GetParameter(1),gfit->GetParameter(2));
			gsnR->SetLineColor(kRed);
			gsnR->SetMarkerColor(kRed);
			gsnR->SetLineWidth(2);
			gsnR->SetLineStyle(2);
			gsnR->Draw("same");			
			*/
		}
		
		//name
		std::stringstream mipname;
		if(imip==0) mipname << "mip, E_{ecal} < " << mipcut << " GeV";
		else if(imip==1) mipname << "nomip, E_{ecal} > " << mipcut << " GeV";
		else mipname << "total";
		
		//pave
		pave = new TPaveText(xmin,0.68,xmin+0.2,0.78,"NDC");
		pave->AddText(("("+sp->gname+")").c_str());
		pave->AddText((piname.str()).c_str());
		if(imip!=2) pave->AddText((mipname.str()).c_str());
		pave->SetFillColor(0);
		pave->SetBorderSize(0);
		pave->SetTextFont(42);
		pave->SetTextSize(0.05);
		pave->Draw("same");

		//par pave
		Double_t ymin1;
		//if(do_fit) ymin1 = 0.26;
		//else ymin1 = 0.51;
		ymin1 = 0.47;
		pave_par = new TPaveText(xmin,ymin1,xmin+0.2,ymin1+0.05*4,"NDC");
		pave_par->AddText((Nname.str()).c_str());
		pave_par->AddText((muname.str()).c_str());
		pave_par->AddText((signame.str()).c_str());
		pave_par->AddText((musigname.str()).c_str());
		//if(do_fit){
		//	pave_par->AddText((aLname.str()).c_str());
		//	pave_par->AddText((nLname.str()).c_str());
		//	pave_par->AddText((aRname.str()).c_str());
		//	pave_par->AddText((nRname.str()).c_str());
		//	pave_par->AddText((chiname.str()).c_str());
		//}
		pave_par->SetFillColor(0);
		pave_par->SetBorderSize(0);
		pave_par->SetTextFont(42);
		pave_par->SetTextSize(0.05);
		pave_par->Draw("same");
		
		if (imip==0) std::cout << "mip:" << std::endl;
		else if (imip==1) std::cout << "nomip:" << std::endl;
		else std::cout << "total:" << std::endl;
		
		std::cout << Nname.str() << std::endl;
		std::cout << muname.str() << std::endl;
		std::cout << signame.str() << std::endl;
		std::cout << musigname.str() << std::endl;
		if(do_fit){
		//	std::cout << "aL = " << gfit->GetParameter(3) << " +/- " << gfit->GetParError(3) << std::endl;
		//	std::cout << "nL = " << gfit->GetParameter(4) << " +/- " << gfit->GetParError(4) << std::endl;
		//	std::cout << "aR = " << gfit->GetParameter(5) << " +/- " << gfit->GetParError(5) << std::endl;
		//	std::cout << "nR = " << gfit->GetParameter(6) << " +/- " << gfit->GetParError(6) << std::endl;
			std::cout << "chi^2/ndf = " << gfit->GetChisquare()/gfit->GetNDF() << std::endl;
		}
		
		if(do_print) can->Print((oname.str()+"."+pformat).c_str(),pformat.c_str());
		if(do_batch) _file->Close();
	}
	else { _file->Close(); }
	
	//return data structure with relevant info
	return theRes;
}

//---------------------------------------------------------------------------
//function to plot res vs. energy for pions (uses current dir, etc. settings)
//imip: 0 = mip, 1 = nomip, 2 = total
//qty: 0 = response, 1 = resolution, 2 = mip percentage
TGraphErrors* g4_plot_res(int snum, int imip, int qty, bool do_fit, bool do_jet, bool do_show, bool do_print=false){
	int maxE = do_jet ? maxJTe : maxHDe;

	Sample* sp = sample_map[snum];
	if(!sp) { std::cout << "Sample " << snum << " is not loaded." << std::endl; return 0; }

	//store values from get_res
	Double_t* vals = new Double_t[maxE]; //sigma or mean
	Double_t* xvals = new Double_t[maxE]; //sigma or mean
	Double_t* y_errors = new Double_t[maxE]; //errors on pars
	Double_t* x_errors = new Double_t[maxE]; //errors on energy (usually 0)

	//for storage of output info
	energyRes* res_temp;

	for (int i = 0; i < maxE; i++){
		double energy;
		if(do_jet) energy = jet_energies[i];
		else energy = energies[i];
	
		Double_t v, ve;
	
		if(qty==2){
			res_temp = get_res(snum,i,2,do_fit,do_jet,0,0,1); //tot
			Double_t N = res_temp->getStat(0);
			res_temp = get_res(snum,i,0,do_fit,do_jet,0,0,1); //mip
			Double_t Nmip = res_temp->getStat(0);
			v = Nmip/N;
			ve = 0;
		}
		else {
			res_temp = get_res(snum,i,imip,do_fit,do_jet,0,0,1);

			Double_t m, me, s, se;
			if(do_fit){
				TF1* fit = res_temp->getFit();
				m = fit->GetParameter(1);
				me = fit->GetParError(1);
				s = fit->GetParameter(2);
				se = fit->GetParError(2);
			}
			else{
				m = res_temp->getStat(1);
				me = res_temp->getStatErr(1);
				s = res_temp->getStat(2);
				se = res_temp->getStatErr(2);
			}
			
			if(qty==1){
				v = s/m;
				ve = v*sqrt(se*se/(s*s)+me*me/(m*m));
			}
			else if(qty==0){
				v = m/energy;
				ve = me/energy;
			}
		}
		
		xvals[i] = energy;
		//xvals[i] = 1./sqrt(energy);
		vals[i] = v;
		y_errors[i] = ve;
		x_errors[i] = 0;
	}

	TCanvas* can;
	TPaveText* pave;
	TGraphErrors* val_graph;

	Int_t col, mrk;
	col = kBlack; mrk = 20;

	//graph values
	std::string qtyaxes[] = {"Response (#mu/E_{true})","Resolution (#sigma/#mu)","MIP fraction"};
	if(do_jet) { qtyaxes[0] = "Jet " + qtyaxes[0]; qtyaxes[1] = "Jet " + qtyaxes[1]; }
	else { qtyaxes[0] = "Pion " + qtyaxes[0]; qtyaxes[1] = "Pion " + qtyaxes[1]; }
	val_graph = new TGraphErrors(maxE,xvals,vals,x_errors,y_errors);
	val_graph->GetXaxis()->SetTitle("Energy [GeV]");
	val_graph->GetYaxis()->SetTitle(qtyaxes[qty].c_str());
	val_graph->SetTitle("");
	val_graph->SetMarkerStyle(mrk);
	val_graph->SetMarkerColor(col);
	val_graph->SetMarkerSize(1.5);
	val_graph->SetLineColor(col);
	val_graph->SetFillColor(0);

	//fit jet resolution
	TF1* gfit = 0;
	std::stringstream sname, cname, chiname;
	sname.precision(2);
	cname.precision(2);
	chiname.precision(2);	
	if(do_jet && qty==1){
		gfit = new TF1("jetres","sqrt([0]*[0]/x+[1]*[1])",val_graph->GetXaxis()->GetXmin(),val_graph->GetXaxis()->GetXmax());
		//val_graph->Fit(gfit,"NQR");
		val_graph->Fit(gfit,"NR");
		std::cout.precision(2);
		std::cout << "chi^2/ndf = " << gfit->GetChisquare() << " / " << gfit->GetNDF() << " = " << gfit->GetChisquare()/gfit->GetNDF() << std::endl;
		double a, ae, b, be;
		a = gfit->GetParameter(0);
		ae = gfit->GetParError(0);
		b = gfit->GetParameter(1);
		be = gfit->GetParError(1);
		std::cout << sp->gname << ((!sp->supp_info.empty()) ? " " + sp->supp_info : "") << ": " << a*100 << "%/sqrt(E) (+) " << b*100 << "%" << std::endl;
		//store in sample
		sp->sampling_term = a;
		sp->sampling_term_err = ae;
		sp->constant_term = b;
		sp->constant_term_err = be;
		//std::cout << "(+/- " << sp->sampling_term_err << ", +/- " << sp->constant_term_err << ", chi2/ndf = " << gfit->GetChisquare() << " / " << gfit->GetNDF() << ")" << std::endl;
		
		//get error ellipse from minuit
		gMinuit->SetErrorDef(1); //1-sigma
		TGraph* g_ellipse = (TGraph*)gMinuit->Contour(500,0,1); //100 pts, correlation of params 0 and 1
		sp->err_ellipse = g_ellipse;
		
		//for pave
		sname << "s = " << a*100 << "% #pm " << ae*100 << "%";
		cname << "c = " << b*100 << "% #pm " << be*100 << "%";
		chiname << "#chi^{2}/ndf = " << gfit->GetChisquare() << " / " << gfit->GetNDF() << " = " << gfit->GetChisquare()/gfit->GetNDF();
	}
	
	if(do_show){
		can = new TCanvas("res","res",700,500);
		can->cd();
		can->SetLogx();

		if(qty) val_graph->GetYaxis()->SetRangeUser(0,0.4);
		else val_graph->GetYaxis()->SetRangeUser(0,1.1);
		val_graph->Draw("APZ");

		//legend, pave coords
		double y1;
		if(qty) y1 = 0.5;
		else y1 = 0.2;	
		
		if(gfit){
			gfit->SetLineColor(kRed);
			gfit->SetMarkerColor(kRed);
			gfit->SetLineWidth(2);
			gfit->Draw("same");
		}
		
		double y2 = y1 + 0.1;
		if(gfit) y2 = y1 + 6*0.06;
		pave = new TPaveText(0.5,y1,0.9,y2,"NDC");
		pave->AddText(sp->gname.c_str());
		if(!sp->supp_info.empty()) pave->AddText(sp->supp_info.c_str());
		else pave->AddText(" ");
		if(gfit){
			pave->AddText("f(E) = s/#surdE #oplus c = #sqrt{s^{2}/E + c^{2}}");
			pave->AddText((sname.str()).c_str());
			pave->AddText((cname.str()).c_str());
			pave->AddText((chiname.str()).c_str());
		}
		if(qty==2 || do_jet) {} //no mip pave for nmip
		else if(imip==0) pave->AddText("mip, E_{ecal} < 1 GeV");
		else if(imip==1) pave->AddText("nomip, E_{ecal} > 1 GeV");
		else pave->AddText("total, N = 10000");
		pave->SetFillColor(0);
		pave->SetBorderSize(0);
		pave->SetTextFont(42);
		pave->SetTextSize(0.05);
		pave->Draw("same");
		
		if(do_print){
			//names
			std::string omip;
			if (imip==0) omip = "mip";
			else if (imip==1) omip = "nomip";
			else omip = "tot";
			std::string stype = sp->fpre;
			std::string rtype = "pion";
			if(do_jet) rtype = "monojet";
			std::string qtyname[] = {"mu","sigma","nmip"};

			std::stringstream oname;
			oname << pdir << "/" << stype << "_" << rtype << "_" << qtyname[qty] << "_" << omip << "." << pformat;
			can->Print((oname.str()).c_str(),pformat.c_str());
		}
	}

	return val_graph;
}

//--------------------------------------------------------------------------
//function to compare plots of response or resolution from previous function
void g4_comp(int s1, int s2, int imip, int qty, bool do_fit, bool do_jet, bool do_print=false){
	int maxE = do_jet ? maxJTe : maxHDe;

	TGraphErrors* the_graphs[2];
	Double_t* vals[2];
	Double_t* errs[2];

	Sample* sp1 = sample_map[s1];
	Sample* sp2 = sample_map[s2];
	if(!sp1) { std::cout << "Sample " << s1 << " is not loaded." << std::endl; return; }
	else if(!sp2) { std::cout << "Sample " << s2 << " is not loaded." << std::endl; return; }
	int legentries = 0;
	if(!sp1->supp_info.empty()) legentries += 2;
	else legentries += 1;
	if(!sp2->supp_info.empty()) legentries += 2;
	else legentries += 1;
	
	//setup canvas with histo and ratio areas
	std::string qtycans[] = {"mu-comp","sigma-comp","nmip-comp"};
	TCanvas* can = new TCanvas(qtycans[qty].c_str(),qtycans[qty].c_str(),700,750);
	TPad* pad1;
	can->cd();
	pad1 = new TPad("graph","",0,0.27,1.0,1.0);
	pad1->SetMargin(0.15,0.05,0.02,0.05);
	pad1->Draw();
	pad1->SetTicks(1,1);
	pad1->SetLogx();
	TPad* pad2;
	can->cd();
	pad2 = new TPad("dmc","",0,0,1.0,0.25);
	pad2->SetMargin(0.15,0.05,0.35,0.05);
	pad2->Draw();
	pad2->SetTicks(1,1);
	pad2->SetLogx();
	pad1->cd();

	//legend, pave coords
	double y1, y2;
	if(qty==0) { y1 = 0.05; y2 = y1 + 0.05*legentries; }
	else { y2 = 0.9; y1 = y2 - 0.05*legentries; }

	TLegend* leg = new TLegend(0.5,y1,0.9,y2);
	//std::string gname[] = {"PbWO_{4} ECAL","W-LYSO ECAL"};
	//std::string gname_rat[] = {"PbWO_{4}","W-LYSO"};
	std::string gname[] = {sp1->gname,sp2->gname};
	std::string gname_supp[] = {sp1->supp_info,sp2->supp_info};
	std::string gname_rat[] = {sp1->gname_rat,sp2->gname_rat};
	
	//get and plot graphs
	//set_dir("/home/pedrok/CMSSW_4_2_8/src/ForwardCaloUpgrade/GEANT/runs","hcal","plots");
	//hom_ecal = true;
	//sam_init = false;
	the_graphs[0] = g4_plot_res(s1,imip,qty,do_fit,do_jet,0);//default
	//set_dir("/data/users/pedrok/FullSim/wlyso",wname,"plots");
	//hom_ecal = false;
	//sam_init = false;	
	the_graphs[1] = g4_plot_res(s2,imip,qty,do_fit,do_jet,0);//shashlik
	
	Color_t color[2] = {kBlue, kRed};
	Int_t marker[2] = {21, 33};
	
	for(int i = 0; i < 2; i++){
		vals[i] = the_graphs[i]->GetY();
		errs[i] = the_graphs[i]->GetEY();

		leg->AddEntry(the_graphs[i],(gname[i]).c_str(),"pl");
		if(!gname_supp[i].empty()) leg->AddEntry((TObject*)0,gname_supp[i].c_str(),"");
		//if(qty) the_graphs[i]->GetYaxis()->SetRangeUser(0,1);
		//else the_graphs[i]->GetYaxis()->SetRangeUser(0,1.3);

		//formatting
		if(qty==1) the_graphs[i]->GetYaxis()->SetRangeUser(0,0.5);
		else if(qty==2) the_graphs[i]->GetYaxis()->SetRangeUser(0,1.0);
		else if(qty==0) the_graphs[i]->GetYaxis()->SetRangeUser(0,1.1);
		the_graphs[i]->GetXaxis()->SetLabelOffset(999);
		the_graphs[i]->GetYaxis()->SetTitleOffset(1.1);
		the_graphs[i]->GetYaxis()->SetTitleSize(32/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[i]->GetYaxis()->SetLabelSize(28/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[i]->GetYaxis()->SetTickLength(12/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[i]->GetXaxis()->SetTickLength(12/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[i]->SetMarkerStyle(marker[i]);
		the_graphs[i]->SetMarkerColor(color[i]);
		the_graphs[i]->SetMarkerSize(1.5);
		the_graphs[i]->SetLineColor(color[i]);
		the_graphs[i]->SetFillColor(0);
	}
	the_graphs[0]->Draw("APZ");
	the_graphs[1]->Draw("PZ same");

	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.05);
	leg->SetTextFont(42);
	leg->SetTextSize(0.05);
	leg->Draw("same");
	
	//setup prelim text
	double umax3 = .95;
	TPaveText* pave2 = new TPaveText(umax3-0.3,0.95,umax3,1.0,"NDC");
	pave2->SetFillColor(0);
	pave2->SetBorderSize(0);
	pave2->SetTextFont(42);
	pave2->SetTextSize(0.04);
	//pave2->AddText("CMS preliminary");
	//pave2->AddText("(no radiation damage)");
	//pave2->Draw("same");

	//plot ratios
	pad2->cd();

	Double_t* rat = new Double_t[maxE]; //sigma or mean
	Double_t* y_errors = new Double_t[maxE]; //errors on pars
	Double_t* x_errors = new Double_t[maxE]; //errors on energy (usually 0)

	//take ratio point-by-point w/ error progression
	for(int i=0; i < maxE; i++){
		rat[i] = vals[1][i]/vals[0][i];
		y_errors[i] = rat[i]*rat[i]*((errs[1][i]*errs[1][i])/(vals[1][i]*vals[1][i]) + (errs[0][i]*errs[0][i])/(vals[0][i]*vals[0][i]));
		x_errors[i] = 0;
	}

	TGraphErrors* rat_graph = new TGraphErrors(maxE,the_graphs[0]->GetX(),rat,x_errors,y_errors);
	rat_graph->GetXaxis()->SetTitle("Energy [GeV]");
	rat_graph->GetYaxis()->SetTitle(("#frac{" + gname_rat[1] + "}{" + gname_rat[0] + "}").c_str());
	rat_graph->SetTitle("");
	rat_graph->SetMarkerStyle(20);
	rat_graph->SetMarkerColor(kBlack);
	rat_graph->SetMarkerSize(1.25);
	rat_graph->SetLineColor(kBlack);
	rat_graph->SetFillColor(0);
	rat_graph->GetYaxis()->SetTitleOffset(0.35);
	rat_graph->GetXaxis()->SetLabelColor(1);
	rat_graph->GetYaxis()->SetLabelColor(1);
	rat_graph->GetXaxis()->SetTitleSize(32/(pad2->GetWh()*pad2->GetAbsHNDC()));
	rat_graph->GetXaxis()->SetLabelSize(28/(pad2->GetWh()*pad2->GetAbsHNDC()));
	rat_graph->GetYaxis()->SetTitleSize(32/(pad2->GetWh()*pad2->GetAbsHNDC()));
	rat_graph->GetYaxis()->SetLabelSize(28/(pad2->GetWh()*pad2->GetAbsHNDC()));
	rat_graph->GetYaxis()->SetTickLength(6/(pad2->GetWh()*pad2->GetAbsHNDC()));
	rat_graph->GetXaxis()->SetTickLength(12/(pad2->GetWh()*pad2->GetAbsHNDC()));
	rat_graph->GetYaxis()->SetNdivisions(503);
	if(do_jet) rat_graph->GetYaxis()->SetRangeUser(0.75,1.25);
	else rat_graph->GetYaxis()->SetRangeUser(0.5,1.5);
	rat_graph->Draw("APZ");

	//line
	TLine *line = new TLine(0,1,rat_graph->GetXaxis()->GetXmax(),1);
	line->SetLineStyle(2);
	line->SetLineWidth(1);
	line->SetLineColor(kBlack);
	line->Draw("same");

	if(do_print){
		//names
		std::string omip;
		if (imip==0) omip = "mip";
		else if (imip==1) omip = "nomip";
		else omip = "tot";
		std::string rtype = "pion";
		if(do_jet) rtype = "monojet";
		std::string qtyname[] = {"mu","sigma","nmip"};

		std::stringstream oname;
		oname << pdir << "/" << rtype << "_" << qtyname[qty] << "_comp_" << sp1->fpre << "_" << sp2->fpre << "_";
		if(do_fit) oname << "fit_";
		oname << omip << "." << pformat;
		can->Print((oname.str()).c_str(),pformat.c_str());
	}

}

//----------------------------------------------------------------------------------------
//function to compare multiple (>2) plots of response or resolution from previous function
void g4_multicomp(int s1, std::vector<int> s2, int imip, int qty, bool do_fit, bool do_jet, bool do_print=false){
	int maxE = do_jet ? maxJTe : maxHDe;

	if(s2.size()+1 > 6) std::cout << "Note: more markers and colors need to be added (only 6 currently defined)" << std::endl;
	
	TGraphErrors** the_graphs = new TGraphErrors*[s2.size()+1];
	Double_t** vals = new Double_t*[s2.size()+1];
	Double_t** errs = new Double_t*[s2.size()+1];
	Sample** sp = new Sample*[s2.size()+1];
	int legentries = 0;
	
	//get s1
	sp[0] = sample_map[s1];
	if(!sp[0]) { std::cout << "Sample " << s1 << " is not loaded." << std::endl; return; }
	if(!sp[0]->supp_info.empty()) legentries += 2;
	else legentries += 1;
	the_graphs[0] = g4_plot_res(s1,imip,qty,do_fit,do_jet,0); //default
	
	//get s2's	
	for(int i = 1; i <= s2.size(); i++){
		sp[i] = sample_map[s2[i-1]];
		if(!sp[i]) { std::cout << "Sample " << s2[i-1] << " is not loaded." << std::endl; return; }
		if(!sp[i]->supp_info.empty()) legentries += 2;
		else legentries += 1;
		the_graphs[i] = g4_plot_res(s2[i-1],imip,qty,do_fit,do_jet,0); //shashlik
	}
	
	//setup canvas with histo and ratio areas
	std::string qtycans[] = {"mu-comp","sigma-comp","nmip-comp"};
	TCanvas* can = new TCanvas(qtycans[qty].c_str(),qtycans[qty].c_str(),700,750);
	TPad* pad1;
	can->cd();
	pad1 = new TPad("graph","",0,0.27,1.0,1.0);
	pad1->SetMargin(0.15,0.05,0.02,0.05);
	pad1->Draw();
	pad1->SetTicks(1,1);
	pad1->SetLogx();
	TPad* pad2;
	can->cd();
	pad2 = new TPad("dmc","",0,0,1.0,0.25);
	pad2->SetMargin(0.15,0.05,0.35,0.05);
	pad2->Draw();
	pad2->SetTicks(1,1);
	pad2->SetLogx();
	pad1->cd();

	//legend, pave coords
	double y2, y1;
	if(qty==0) { y1 = 0.05; y2 = y1 + 0.05*legentries; }
	else { y2 = 0.9; y1 = y2 - 0.05*legentries; }
	
	TLegend* leg = new TLegend(0.5,y1,0.9,y2);
	
	Color_t color[] = {kBlack, kBlue, kMagenta, kRed, kYellow+2, kOrange+7};
	Int_t marker[] = {20, 21, 22, 23, 33, 29};
	
	for(int i = 0; i < s2.size()+1; i++){
		vals[i] = the_graphs[i]->GetY();
		errs[i] = the_graphs[i]->GetEY();

		leg->AddEntry(the_graphs[i],sp[i]->gname.c_str(),"p");
		if(!sp[i]->supp_info.empty()) leg->AddEntry((TObject*)0,sp[i]->supp_info.c_str(),"");
		//if(qty) the_graphs[i]->GetYaxis()->SetRangeUser(0,1);
		//else the_graphs[i]->GetYaxis()->SetRangeUser(0,1.3);

		//formatting
		if(qty==1) the_graphs[i]->GetYaxis()->SetRangeUser(0,0.5);
		else if(qty==2) the_graphs[i]->GetYaxis()->SetRangeUser(0,1.0);
		else if(qty==0) the_graphs[i]->GetYaxis()->SetRangeUser(0,1.1);
		//if(qty==1) the_graphs[i]->GetYaxis()->SetRangeUser(0,0.2);
		//else if(qty==2) the_graphs[i]->GetYaxis()->SetRangeUser(0,1.0);
		//else if(qty==0) the_graphs[i]->GetYaxis()->SetRangeUser(0.7,1.05);

		the_graphs[i]->GetXaxis()->SetLabelOffset(999);
		the_graphs[i]->GetYaxis()->SetTitleOffset(1.1);
		the_graphs[i]->GetYaxis()->SetTitleSize(32/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[i]->GetYaxis()->SetLabelSize(28/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[i]->GetYaxis()->SetTickLength(12/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[i]->GetXaxis()->SetTickLength(12/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[i]->SetMarkerStyle(marker[i]);
		the_graphs[i]->SetMarkerColor(color[i]);
		the_graphs[i]->SetMarkerSize(1.5);
		the_graphs[i]->SetLineColor(color[i]);
		the_graphs[i]->SetFillColor(0);
		if(i==0) the_graphs[i]->Draw("APZ");
		else the_graphs[i]->Draw("PZ same");
	}

	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.05);
	leg->SetTextFont(42);
	leg->SetTextAlign(11); //testing...
	leg->Draw("same");
	
	//setup prelim text
	double umax3 = .95;
	TPaveText* pave2 = new TPaveText(umax3-0.3,0.95,umax3,1.0,"NDC");
	pave2->SetFillColor(0);
	pave2->SetBorderSize(0);
	pave2->SetTextFont(42);
	pave2->SetTextSize(0.04);
	//pave2->AddText("CMS preliminary");
	//pave2->AddText("(no radiation damage)");
	//pave2->Draw("same");

	//plot ratios
	pad2->cd();
	
	TGraphErrors** rat_graph = new TGraphErrors*[s2.size()];
	
	for(int i = 1; i <= s2.size(); i++){
		Double_t* rat = new Double_t[maxE]; //sigma or mean
		Double_t* y_errors = new Double_t[maxE]; //errors on pars
		Double_t* x_errors = new Double_t[maxE]; //errors on energy (usually 0)

		//take ratio point-by-point w/ error progression
		for(int j = 0; j < maxE; j++){
			rat[j] = vals[i][j]/vals[0][j];
			y_errors[j] = rat[j]*rat[j]*((errs[i][j]*errs[i][j])/(vals[i][j]*vals[i][j]) + (errs[0][j]*errs[0][j])/(vals[0][j]*vals[0][j]));
			x_errors[j] = 0;
		}

		rat_graph[i-1] = new TGraphErrors(maxE,the_graphs[0]->GetX(),rat,x_errors,y_errors);
		rat_graph[i-1]->GetXaxis()->SetTitle("Energy [GeV]");
		rat_graph[i-1]->GetYaxis()->SetTitle(("#frac{various}{" + sp[0]->gname_rat + "}").c_str());
		rat_graph[i-1]->SetTitle("");
		rat_graph[i-1]->SetMarkerStyle(marker[i]);
		rat_graph[i-1]->SetMarkerColor(color[i]);
		rat_graph[i-1]->SetMarkerSize(1.25);
		rat_graph[i-1]->SetLineColor(color[i]);
		rat_graph[i-1]->SetFillColor(0);
		rat_graph[i-1]->GetYaxis()->SetTitleOffset(0.35);
		rat_graph[i-1]->GetXaxis()->SetLabelColor(1);
		rat_graph[i-1]->GetYaxis()->SetLabelColor(1);
		rat_graph[i-1]->GetXaxis()->SetTitleSize(32/(pad2->GetWh()*pad2->GetAbsHNDC()));
		rat_graph[i-1]->GetXaxis()->SetLabelSize(28/(pad2->GetWh()*pad2->GetAbsHNDC()));
		rat_graph[i-1]->GetYaxis()->SetTitleSize(32/(pad2->GetWh()*pad2->GetAbsHNDC()));
		rat_graph[i-1]->GetYaxis()->SetLabelSize(28/(pad2->GetWh()*pad2->GetAbsHNDC()));
		rat_graph[i-1]->GetYaxis()->SetTickLength(6/(pad2->GetWh()*pad2->GetAbsHNDC()));
		rat_graph[i-1]->GetXaxis()->SetTickLength(12/(pad2->GetWh()*pad2->GetAbsHNDC()));
		rat_graph[i-1]->GetYaxis()->SetNdivisions(503);
		if(do_jet) {
			rat_graph[i-1]->GetYaxis()->SetRangeUser(0.75,1.25);
			//if(qty==0) rat_graph[i-1]->GetYaxis()->SetRangeUser(0.95,1.25);
			//else if(qty==1) rat_graph[i-1]->GetYaxis()->SetRangeUser(0.75,1.05);
		}
		else rat_graph[i-1]->GetYaxis()->SetRangeUser(0.5,1.5);
		
		if(i==1) rat_graph[i-1]->Draw("APZ");
		else rat_graph[i-1]->Draw("PZ same");
	}

	//line
	TLine *line = new TLine(0,1,rat_graph[0]->GetXaxis()->GetXmax(),1);
	line->SetLineStyle(2);
	line->SetLineWidth(1);
	line->SetLineColor(kBlack);
	line->Draw("same");

	if(do_print){
		//names
		std::string omip;
		if (imip==0) omip = "mip";
		else if (imip==1) omip = "nomip";
		else omip = "tot";
		std::string rtype = "pion";
		if(do_jet) rtype = "monojet";
		std::string qtyname[] = {"mu","sigma","nmip"};

		std::stringstream oname;
		oname << pdir << "/" << rtype << "_" << qtyname[qty] << "_comp_";
		for(int i = 0; i <= s2.size(); i++){
			oname << sp[i]->fpre << "_";
		}
		if(do_fit) oname << "fit_";
		oname << omip << "." << pformat;
		can->Print((oname.str()).c_str(),pformat.c_str());
	}

}

//--------------------------------------------------------
//macro to run g4_multicomp with predefined set of samples
void run_g4_multicomp(int imip, int qty, bool do_fit, bool do_jet, bool do_print=false){
	int s1 = 0;
	std::vector<int> s2;
	s2.push_back(100);
	s2.push_back(3);
	s2.push_back(5);
	s2.push_back(7);
	
	g4_multicomp(s1,s2,imip,qty,do_fit,do_jet,do_print);
}

//----------------------------------------------------------------------------------------
//function to compare multiple (>2) plots of response or resolution from previous function
//performance plot style (Phase 2 TP)
void g4_multicomp_tp(int s1, std::vector<int> s2, int imip, int qty, bool do_fit, bool do_jet, bool do_print=false){
	int maxE = do_jet ? maxJTe : maxHDe;

	if(s2.size()+1 > 6) std::cout << "Note: more markers and colors need to be added (only 6 currently defined)" << std::endl;
	
	TGraphErrors** the_graphs = new TGraphErrors*[s2.size()+1];
	Double_t** vals = new Double_t*[s2.size()+1];
	Double_t** errs = new Double_t*[s2.size()+1];
	Sample** sp = new Sample*[s2.size()+1];
	
	//get s1
	sp[0] = sample_map[s1];
	if(!sp[0]) { std::cout << "Sample " << s1 << " is not loaded." << std::endl; return; }
	the_graphs[0] = g4_plot_res(s1,imip,qty,do_fit,do_jet,0); //default
	
	//get s2's	
	for(int i = 1; i <= s2.size(); i++){
		sp[i] = sample_map[s2[i-1]];
		if(!sp[i]) { std::cout << "Sample " << s2[i-1] << " is not loaded." << std::endl; return; }
		the_graphs[i] = g4_plot_res(s2[i-1],imip,qty,do_fit,do_jet,0); //shashlik
	}
	
	//setup canvas with histo and ratio areas
	std::string qtycans[] = {"mu-comp","sigma-comp","nmip-comp"};
	TCanvas* can = new TCanvas(qtycans[qty].c_str(),qtycans[qty].c_str(),700,750);
	TPad* pad1;
	can->cd();
	pad1 = new TPad("graph","",0,0.27,1.0,1.0);
	pad1->SetMargin(0.15,0.05,0.02,0.05);
	pad1->Draw();
	pad1->SetTicks(1,1);
	pad1->SetLogx();
	TPad* pad2;
	can->cd();
	pad2 = new TPad("dmc","",0,0,1.0,0.25);
	pad2->SetMargin(0.15,0.05,0.35,0.05);
	pad2->Draw();
	pad2->SetTicks(1,1);
	pad2->SetLogx();
	pad1->cd();

	//legend, pave coords
	double y2;
	if(qty==0) y2 = 0.4;
	else y2 = 0.9;
	
	double y1 = y2 - 0.05*(s2.size()+1);

	//TLegend* leg = new TLegend(0.4,y1,0.9,y2);
	TLegend* leg = new TLegend(0.5,y1,0.9,y2);
	
	Color_t color[] = {kBlack, kBlue, kMagenta, kRed, kYellow+2, kOrange+7};
	Int_t marker[] = {20, 21, 22, 23, 33, 29};
	
	for(int i = 0; i < s2.size()+1; i++){
		vals[i] = the_graphs[i]->GetY();
		errs[i] = the_graphs[i]->GetEY();

		string legname;
		legname = sp[i]->gname;
		leg->AddEntry(the_graphs[i],(legname).c_str(),"p");
		//if(qty) the_graphs[i]->GetYaxis()->SetRangeUser(0,1);
		//else the_graphs[i]->GetYaxis()->SetRangeUser(0,1.3);

		//formatting
		if(qty==1) the_graphs[i]->GetYaxis()->SetRangeUser(0,0.2);
		else if(qty==2) the_graphs[i]->GetYaxis()->SetRangeUser(0,1.0);
		else if(qty==0) the_graphs[i]->GetYaxis()->SetRangeUser(0.7,1.05);
		the_graphs[i]->GetXaxis()->SetLabelOffset(999);
		the_graphs[i]->GetYaxis()->SetTitleOffset(1.1);
		the_graphs[i]->GetYaxis()->SetTitleSize(32/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[i]->GetYaxis()->SetLabelSize(28/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[i]->GetYaxis()->SetTickLength(12/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[i]->GetXaxis()->SetTickLength(12/(pad1->GetWh()*pad1->GetAbsHNDC()));
		the_graphs[i]->SetMarkerStyle(marker[i]);
		the_graphs[i]->SetMarkerColor(color[i]);
		the_graphs[i]->SetMarkerSize(1.5);
		the_graphs[i]->SetLineColor(color[i]);
		the_graphs[i]->SetFillColor(0);
		if(i==0) the_graphs[i]->Draw("APZ");
		else the_graphs[i]->Draw("PZ same");
	}

	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.05);
	leg->SetTextFont(42);
	leg->SetTextAlign(11); //testing...
	leg->Draw("same");
	
	//setup prelim text
	double umax3 = .95;
	TPaveText* pave2 = new TPaveText(umax3-0.3,0.95,umax3,1.0,"NDC");
	pave2->SetFillColor(0);
	pave2->SetBorderSize(0);
	pave2->SetTextFont(42);
	pave2->SetTextSize(0.04);
	//pave2->AddText("CMS preliminary");
	//pave2->AddText("(no radiation damage)");
	//pave2->Draw("same");

	//plot ratios
	pad2->cd();
	
	TGraphErrors** rat_graph = new TGraphErrors*[s2.size()];
	
	for(int i = 1; i <= s2.size(); i++){
		Double_t* rat = new Double_t[maxE]; //sigma or mean
		Double_t* y_errors = new Double_t[maxE]; //errors on pars
		Double_t* x_errors = new Double_t[maxE]; //errors on energy (usually 0)

		//take ratio point-by-point w/ error progression
		for(int j = 0; j < maxE; j++){
			rat[j] = vals[i][j]/vals[0][j];
			y_errors[j] = rat[j]*rat[j]*((errs[i][j]*errs[i][j])/(vals[i][j]*vals[i][j]) + (errs[0][j]*errs[0][j])/(vals[0][j]*vals[0][j]));
			x_errors[j] = 0;
		}

		rat_graph[i-1] = new TGraphErrors(maxE,the_graphs[0]->GetX(),rat,x_errors,y_errors);
		rat_graph[i-1]->GetXaxis()->SetTitle("Energy [GeV]");
		if(s2.size()>1) rat_graph[i-1]->GetYaxis()->SetTitle(("#frac{various}{" + sp[0]->gname_rat + "}").c_str());
		else rat_graph[i-1]->GetYaxis()->SetTitle(("#frac{" + sp[1]->gname_rat + "}{" + sp[0]->gname_rat + "}").c_str());
		rat_graph[i-1]->SetTitle("");
		rat_graph[i-1]->SetMarkerStyle(marker[i]);
		rat_graph[i-1]->SetMarkerColor(color[i]);
		rat_graph[i-1]->SetMarkerSize(1.25);
		rat_graph[i-1]->SetLineColor(color[i]);
		rat_graph[i-1]->SetFillColor(0);
		rat_graph[i-1]->GetYaxis()->SetTitleOffset(0.35);
		rat_graph[i-1]->GetXaxis()->SetLabelColor(1);
		rat_graph[i-1]->GetYaxis()->SetLabelColor(1);
		rat_graph[i-1]->GetXaxis()->SetTitleSize(32/(pad2->GetWh()*pad2->GetAbsHNDC()));
		rat_graph[i-1]->GetXaxis()->SetLabelSize(28/(pad2->GetWh()*pad2->GetAbsHNDC()));
		rat_graph[i-1]->GetYaxis()->SetTitleSize(32/(pad2->GetWh()*pad2->GetAbsHNDC()));
		rat_graph[i-1]->GetYaxis()->SetLabelSize(28/(pad2->GetWh()*pad2->GetAbsHNDC()));
		rat_graph[i-1]->GetYaxis()->SetTickLength(6/(pad2->GetWh()*pad2->GetAbsHNDC()));
		rat_graph[i-1]->GetXaxis()->SetTickLength(12/(pad2->GetWh()*pad2->GetAbsHNDC()));
		rat_graph[i-1]->GetYaxis()->SetNdivisions(503);
		if(do_jet) {
			rat_graph[i-1]->GetYaxis()->SetRangeUser(0.75,1.25);
			if(qty==0) rat_graph[i-1]->GetYaxis()->SetRangeUser(0.95,1.25);
			else if(qty==1) rat_graph[i-1]->GetYaxis()->SetRangeUser(0.75,1.05);
		}
		else rat_graph[i-1]->GetYaxis()->SetRangeUser(0.5,1.5);
		
		if(i==1) rat_graph[i-1]->Draw("APZ");
		else rat_graph[i-1]->Draw("PZ same");
	}

	//line
	TLine *line = new TLine(0,1,rat_graph[0]->GetXaxis()->GetXmax(),1);
	line->SetLineStyle(2);
	line->SetLineWidth(1);
	line->SetLineColor(kBlack);
	line->Draw("same");

	if(do_print){
		//names
		std::string omip;
		if (imip==0) omip = "mip";
		else if (imip==1) omip = "nomip";
		else omip = "tot";
		std::string rtype = "pion";
		if(do_jet) rtype = "monojet";
		std::string qtyname[] = {"mu","sigma","nmip"};

		std::stringstream oname;
		oname << pdir << "/tp_" << rtype << "_" << qtyname[qty] << "_comp_";
		for(int i = 0; i <= s2.size(); i++){
			oname << sp[i]->fpre << "_";
		}
		if(do_fit) oname << "fit_";
		oname << omip << "." << pformat;
		can->Print((oname.str()).c_str(),pformat.c_str());
	}

}

//--------------------------------------------------------
//macro to run g4_multicomp with predefined set of samples
void run_g4_multicomp_tp(int imip, int qty, bool do_fit, bool do_jet, bool do_print=false){
	//modify sample display names
	Sample* sp = sample_map[0];
	if(!sp) { std::cout << "Sample " << 0 << " is not loaded." << std::endl; return; }
	sp->gname = "PbWO_{4} ECAL + HE";
	sp = sample_map[5];
	if(!sp) { std::cout << "Sample " << 5 << " is not loaded." << std::endl; return; }
	sp->gname = "W-LYSO ECAL + HE";
	//sp = sample_map[7];
	//if(!sp) { std::cout << "Sample " << 7 << " is not loaded." << std::endl; return; }
	//sp->gname = "W-LYSO ECAL + ext. HE";
	
	int s1 = 0;
	std::vector<int> s2;
	s2.push_back(5);
	//s2.push_back(7);
	
	g4_multicomp_tp(s1,s2,imip,qty,do_fit,do_jet,do_print);
}

//----------------------------------------------------------------------
//function to compare multiple response distributions for a given energy
void plot_overlay(std::vector<int> snum, int num, int imip, bool do_fit, bool do_jet, bool do_print=false){
	int maxE = do_jet ? maxJTe : maxHDe;

	if(snum.size() > 6) std::cout << "Note: more markers and colors need to be added (only 6 currently defined)" << std::endl;

	if (num>=maxE || num<0) { std::cout << "num must be between 0 and " << maxE - 1 << std::endl; energyRes* theRes = new energyRes(0,0); return; }

	double energy = 0;
	if(do_jet) energy = jet_energies[num];
	else energy = energies[num];
	
	if (imip>2 || imip<0) { std::cout << "imip must be between 0 and 2" << std::endl; energyRes* theRes = new energyRes(0,0); return; }
	
	if(do_jet) imip = 2; //no mip options for jets
	
	TH1F** hist = new TH1F*[snum.size()];
	Sample** sp = new Sample*[snum.size()];
	
	//get samples
	for(int i = 0; i < snum.size(); i++){
		sp[i] = sample_map[snum[i]];
		if(!sp[i]) { std::cout << "Sample " << snum[i-1] << " is not loaded." << std::endl; return; }
		energyRes* res_temp = get_res(snum[i],num,imip,do_fit,do_jet,0);
		hist[i] = res_temp->getHist();
	}
	
	//setup canvas with histo and ratio areas
	std::stringstream cname;
	cname << (do_jet ? "jet" : "pion") << "_res_" << energy;
	TCanvas* can = new TCanvas((cname.str()).c_str(),(cname.str()).c_str(),900,500);
	can->cd();
	TPad* pad = new TPad("graph","",0,0,1,1);
	pad->SetMargin(0.075,0.35,0.15,0.05);
	pad->SetTicks(1,1);
	pad->SetLogy();
	pad->Draw();
	pad->cd();

	//determine placing of legend and paves - legend goes on side with more space
	//judge by checking how many means are left vs. right
	//int m_left, m_right;
	//m_left = m_right = 0;
	//for(int i = 0; i < snum.size(); i++){
	//	//check if mean is to the left of x-axis middle
	//	if(hist[i]->GetMean()/((hist[i]->GetXaxis()->GetXmax() + hist[i]->GetXaxis()->GetXmin())/2) < 1) ++m_left;
	//	else ++m_right;
	//}
	
	//legend, pave coords
	double x1 = 0.675;
	double y2 = 0.9;
	double y1 = y2 - 0.05 - 0.1*(snum.size());
	//double y2 = y1 + 0.1*(snum.size());

	TLegend* leg = new TLegend(x1,y1,x1+0.325,y2-0.05);
	
	Color_t color[] = {kBlack, kBlue, kMagenta, kRed, kYellow+2, kOrange+7};
	Int_t marker[] = {20, 21, 22, 23, 33, 29};
	
	for(int i = 0; i < snum.size(); i++){
		string legname;
		if(!sp[i]->supp_info.empty()) legname = "#splitline{" + sp[i]->gname + "}{" + sp[i]->supp_info + "}";
		else legname = sp[i]->gname;
		leg->AddEntry(hist[i],(legname).c_str(),"l");
		//if(qty) the_graphs[i]->GetYaxis()->SetRangeUser(0,1);
		//else the_graphs[i]->GetYaxis()->SetRangeUser(0,1.3);

		//formatting
		hist[i]->GetYaxis()->SetTitleSize(32/(pad->GetWh()*pad->GetAbsHNDC()));
		hist[i]->GetYaxis()->SetLabelSize(28/(pad->GetWh()*pad->GetAbsHNDC()));
		hist[i]->GetXaxis()->SetTitleSize(32/(pad->GetWh()*pad->GetAbsHNDC()));
		hist[i]->GetXaxis()->SetLabelSize(28/(pad->GetWh()*pad->GetAbsHNDC()));
		hist[i]->GetYaxis()->SetTickLength(12/(pad->GetWh()*pad->GetAbsHNDC()));
		hist[i]->GetXaxis()->SetTickLength(12/(pad->GetWh()*pad->GetAbsHNDC()));
		hist[i]->SetLineColor(color[i]);

		if(i==0) hist[i]->Draw("hist");
		else hist[i]->Draw("hist same");
	}

	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.05);
	leg->SetTextFont(42);
	leg->Draw("same");
	
	//setup prelim text
	std::stringstream piname;
	if(do_jet) piname << "jet " << energy << " GeV";
	else piname << "#pi^{-} " << energy << " GeV";	
	
	TPaveText* pave2 = new TPaveText(x1,y2-0.05,x1+0.325,y2,"NDC");
	pave2->SetFillColor(0);
	pave2->SetBorderSize(0);
	pave2->SetTextFont(42);
	pave2->SetTextSize(0.05);
	//pave2->AddText("CMS preliminary");
	//pave2->AddText("(no radiation damage)");
	pave2->AddText((piname.str()).c_str());
	pave2->Draw("same");

	if(do_print){
		//names
		std::string omip;
		if (imip==0) omip = "mip";
		else if (imip==1) omip = "nomip";
		else omip = "tot";
		std::string rtype = "pion";
		if(do_jet) rtype = "monojet";

		std::stringstream oname;
		oname << pdir << "/" << rtype << "_" << energy << "gev_overlay_";
		for(int i = 0; i < snum.size(); i++){
			oname << sp[i]->fpre << "_";
		}
		if(do_fit) oname << "fit_";
		oname << omip << "." << pformat;
		can->Print((oname.str()).c_str(),pformat.c_str());
	}	
}

//--------------------------------------------------------
//macro to run plot_overlay with predefined set of samples
void run_plot_overlay(int num, int imip, bool do_fit, bool do_jet, bool do_print=false){
	std::vector<int> snum;
	snum.push_back(0);
	//snum.push_back(2);
	//snum.push_back(3);
	snum.push_back(3);
	snum.push_back(5);
	//snum.push_back(6);
	snum.push_back(7);
	//snum.push_back(71);
	
	plot_overlay(snum,num,imip,do_fit,do_jet,do_print);
}

//------------------------------------------------------------------------------------
//function to compare jet quantities for different geom configs and zero layer weights
//qty: 0 = zero layer weight, 1 = sampling term, 2 = constant term
void comp_jet_qty(std::vector<int> snum, int qty, bool do_fit, bool do_print=false){
	int maxE = maxJTe;
	
	TGraphErrors** graph = new TGraphErrors*[maxZeroWts];
	Sample** sp = new Sample*[snum.size()];
	
	double ymin = 1e10;
	double ymax = 0;
	
	for(int z = 0; z < maxZeroWts; z++){
		init_samples(z); //initialize samples with specified type of zero weight

		//arrays for graph
		Double_t* vals = new Double_t[snum.size()]; //qty
		Double_t* xvals = new Double_t[snum.size()]; //geom config number
		Double_t* y_errors = new Double_t[snum.size()]; //errors on qty
		
		//get qty values from samples
		for(int i = 0; i < snum.size(); i++){
			sp[i] = sample_map[snum[i]];
			if(!sp[i]) { std::cout << "Sample " << snum[i-1] << " is not loaded." << std::endl; return; }
			xvals[i] = i;
			if(qty==0){
				vals[i] = sp[i]->zeroWt;
				y_errors[i] = 0;
			}
			else {
				g4_plot_res(snum[i],2,1,0,1,0); //fit jet resolution
				if(qty==1){
					vals[i] = sp[i]->sampling_term*100;
					y_errors[i] = sp[i]->sampling_term_err*100;
				}
				else if(qty==2){
					vals[i] = sp[i]->constant_term*100;
					y_errors[i] = sp[i]->constant_term_err*100;
				}
			}
			
			//check min/max
			if(vals[i]<ymin) {ymin = vals[i];}
			if(vals[i]>ymax) {ymax = vals[i];}
		}

		//create graph
		graph[z] = new TGraphErrors(snum.size(),xvals,vals,0,y_errors);
	}

	std::string ynames[] = {"zero layer weight","sampling term [%]","constant term [%]"}; //vs. qty
	std::string legnames[] = {"w_{0} = 0.5","w_{0}* v1 (pions)","w_{0}* v2 (jets)"}; //vs. z/maxZeroWts

	Color_t color[] = {kBlack, kBlue, kMagenta, kRed, kYellow+2, kOrange+7};
	Int_t marker[] = {20, 21, 22, 23, 33, 29};
	
	//setup canvas
	std::stringstream cname;
	cname << "comp_jet_" << qty;
	TCanvas* can = new TCanvas((cname.str()).c_str(),(cname.str()).c_str(),900,500);
	can->cd();
	TPad* pad = new TPad("graph","",0,0,1,1);
	pad->SetMargin(0.125,0.3,0.15,0.05);
	pad->SetTicks(1,1);
	pad->Draw();
	pad->cd();

	//leg, pave coords
	double x1 = 0.725;
	double x2 = x1 + 0.275;
	double y2 = 0.95;
	double y1 = y2 - 0.05*maxZeroWts;
	double y0 = y1 - 0.1*(snum.size());

	TLegend* leg = new TLegend(x1,y1,x2,y2);
	TPaveText* pave = new TPaveText(x1,y0,x2,y1,"NDC");
	
	//second loop for formatting, etc.
	for(int z = 0; z < maxZeroWts; z++){
		//formatting
		graph[z]->SetTitle("");
		graph[z]->GetYaxis()->SetTitle(ynames[qty].c_str());
		graph[z]->GetXaxis()->SetTitle("Geometry configuration");
		graph[z]->GetXaxis()->SetNdivisions(snum.size());
		graph[z]->GetYaxis()->SetTitleOffset(0.9);
		graph[z]->GetYaxis()->SetTitleSize(32/(pad->GetWh()*pad->GetAbsHNDC()));
		graph[z]->GetYaxis()->SetLabelSize(28/(pad->GetWh()*pad->GetAbsHNDC()));
		graph[z]->GetXaxis()->SetTitleSize(32/(pad->GetWh()*pad->GetAbsHNDC()));
		graph[z]->GetXaxis()->SetLabelSize(28/(pad->GetWh()*pad->GetAbsHNDC()));
		graph[z]->GetYaxis()->SetTickLength(12/(pad->GetWh()*pad->GetAbsHNDC()));
		graph[z]->GetXaxis()->SetTickLength(12/(pad->GetWh()*pad->GetAbsHNDC()));
		graph[z]->GetYaxis()->SetRangeUser(ymin*0.9,ymax*1.1);
		graph[z]->SetLineColor(color[z]);
		graph[z]->SetMarkerColor(color[z]);
		graph[z]->SetMarkerStyle(marker[z]);
		graph[z]->SetMarkerSize(1.5);
		
		//add to legend
		leg->AddEntry(graph[z],legnames[z].c_str(),"pl");
		
		if(z==0) graph[z]->Draw("APZ");
		else graph[z]->Draw("PZ same");
	}

	//make list of samples
	for(int i = 0; i < snum.size(); i++){
		//sp[i] = sample_map[snum[i]];
		//if(!sp[i]) { std::cout << "Sample " << snum[i-1] << " is not loaded." << std::endl; return; }

		stringstream sname;
		if(!sp[i]->supp_info.empty()) sname << "#splitline{" << i << ". " << sp[i]->gname << "}{   " << sp[i]->supp_info << "}";
		else sname << i << ". " << sp[i]->gname;
		pave->AddText((sname.str()).c_str());
	}
	
	//drawing
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.05);
	leg->SetTextFont(42);
	leg->Draw("same");

	pave->SetFillColor(0);
	pave->SetBorderSize(0);
	pave->SetTextFont(42);
	pave->SetTextSize(0.05);
	pave->SetTextAlign(11);
	pave->Draw("same");	

	if(do_print){
		//names
		std::string qtyname[] = {"zeroWt","sampling","constant"};
		std::stringstream oname;
		oname << pdir << "/comp_jet_" << qtyname[qty] << "_";
		for(int i = 0; i < snum.size(); i++){
			oname << sp[i]->fpre << "_";
		}
		if(do_fit) oname << "fit";
		oname << "." << pformat;
		can->Print((oname.str()).c_str(),pformat.c_str());
	}	
}

//--------------------------------------------------------
//macro to run comp_jet_qty with predefined set of samples
void run_comp_jet_qty(int qty, bool do_fit, bool do_print=false){
	std::vector<int> snum;
	snum.push_back(0);
	//snum.push_back(2);
	//snum.push_back(3);
	snum.push_back(3);
	snum.push_back(5);
	//snum.push_back(6);
	snum.push_back(7);
	//snum.push_back(71);
	
	comp_jet_qty(snum,qty,do_fit,do_print);
}


//----------------------------------------------------------------
//macro to compare resolution terms for different zcal thicknesses
//qty=0: sampling term, qty=1: constant term, qty=2: response for jet with energy=num, qty=3: resolution for jet with energy=num
//qty=4: zero layer weight, qty=5: zcal weight
TGraphErrors* zcal_comp(int qty, int num, bool do_fit, int version=1, bool do_show=true, bool do_print=false){
	double energy = jet_energies[num];
	int nlay = 11;
	//int nlay = 1;
	//if(version==1) nlay += 10;
	//else if(version==2) nlay += 5;
	
	//store values from g4_plot_res
	Double_t* vals = new Double_t[nlay]; //sigma or mean
	Double_t* xvals = new Double_t[nlay]; //sigma or mean
	Double_t* y_errors = new Double_t[nlay]; //errors on pars
	Double_t* x_errors = new Double_t[nlay]; //errors on energy (usually 0)

	//get values for specified qty
	for (int i = 0; i < nlay; i++){
		int snum;
		//if(i==0) snum = 3;
		//use version with 6cm Al dead material
		if(i==0) {
			if(version==1) snum = 5;
			else if(version==2) snum = 8;
		}
		else {
			if(version==1) snum = 30+i; //zcal samples numbered 3x, x = # layers
			else if(version==2) snum = 40+i; //z2cal samples numbered 4x
		}
		Sample* sp = sample_map[snum];
		
		if(qty==0 || qty==1){
			g4_plot_res(snum,2,1,do_fit,1,0); //plots and fits jet resolution
			
			//graph in %
			if(qty==0){
				vals[i] = sp->sampling_term*100;
				y_errors[i] = sp->sampling_term_err*100;
			}
			else if(qty==1){
				vals[i] = sp->constant_term*100;
				y_errors[i] = sp->constant_term_err*100;
			}
		}
		else if(qty==2 || qty==3){
			energyRes* res_temp = get_res(snum,num,2,do_fit,1,0,0,1);
	
			Double_t m, me, s, se;
			if(do_fit){
				TF1* fit = res_temp->getFit();
				m = fit->GetParameter(1);
				me = fit->GetParError(1);
				s = fit->GetParameter(2);
				se = fit->GetParError(2);
			}
			else{
				m = res_temp->getStat(1);
				me = res_temp->getStatErr(1);
				s = res_temp->getStat(2);
				se = res_temp->getStatErr(2);
			}
			
			Double_t v, ve;
			if(qty==3){
				v = s/m;
				ve = v*sqrt(se*se/(s*s)+me*me/(m*m));
			}
			else if(qty==2){
				v = m/energy;
				ve = me/energy;
			}
			
			vals[i] = v;
			y_errors[i] = ve;
		}
		else if(qty==4){
			vals[i] = sp->zeroWt;
			y_errors[i] = 0;
		}
		else if(qty==5){
			vals[i] = sp->zcalWt;
			y_errors[i] = 0;
		}
		
		double abs_thick = 1.8;
		xvals[i] = i*abs_thick/9.95; //units of interaction lengths of tungsten
		x_errors[i] = 0;
	}

	std::stringstream jetname;
	jetname << energy << " GeV jet";

	//graph values
	std::string qtyaxes[] = {"sampling term [%]","constant term [%]","Response (#mu/E_{true})","Resolution (#sigma/#mu)","zero layer weight","HE extension weight"};
	TGraphErrors* val_graph = new TGraphErrors(nlay,xvals,vals,x_errors,y_errors);
	val_graph->GetXaxis()->SetTitle("HE extension depth [#lambda_{0}(W)]");
	val_graph->GetYaxis()->SetTitle(qtyaxes[qty].c_str());
	val_graph->SetTitle("");
	val_graph->SetMarkerStyle(20);
	val_graph->SetMarkerColor(kBlack);
	val_graph->SetMarkerSize(1.5);
	val_graph->SetLineColor(kBlack);
	val_graph->SetFillColor(0);
	
	stringstream cname;
	cname << "zcal_comp_qty" << qty;
	TCanvas* can = new TCanvas((cname.str()).c_str(),(cname.str()).c_str(),700,500);
	can->cd();
	TPad* pad = new TPad("graph","",0,0,1,1);
	pad->SetMargin(0.15,0.05,0.15,0.05);
	pad->Draw();
	pad->cd();

	//formatting
	val_graph->GetYaxis()->SetTitleOffset(1.1);
	val_graph->GetYaxis()->SetTitleSize(32/(pad->GetWh()*pad->GetAbsHNDC()));
	val_graph->GetXaxis()->SetTitleSize(32/(pad->GetWh()*pad->GetAbsHNDC()));
	val_graph->GetYaxis()->SetLabelSize(28/(pad->GetWh()*pad->GetAbsHNDC()));
	val_graph->GetXaxis()->SetLabelSize(28/(pad->GetWh()*pad->GetAbsHNDC()));
	val_graph->GetYaxis()->SetTickLength(12/(pad->GetWh()*pad->GetAbsHNDC()));
	val_graph->GetXaxis()->SetTickLength(12/(pad->GetWh()*pad->GetAbsHNDC()));

	//draw
	if(do_show){
		val_graph->Draw("APZ");
		
		TPaveText* pave;
		if(qty==2 || qty==3){
			pave = new TPaveText(0.3,0.85,0.5,0.9,"NDC");
			pave->SetFillColor(0);
			pave->SetBorderSize(0);
			pave->SetTextFont(42);
			pave->SetTextSize(0.05);
			pave->AddText((jetname.str()).c_str());
			pave->Draw("same");
		}
		
		if(do_print){
			//names
			std::string rtype = "monojet";
			std::string qtyname[] = {"sterm","cterm","mu","sigma","zeroWt","zcalWt"};

			std::stringstream oname;
			oname << pdir << "/" << rtype << "_" << qtyname[qty];
			if(qty==2 || qty==3) oname << "_" << energy << "gev";
			oname << "_";
			if(version==1) oname << "zcal";
			else if(version==2) oname << "zcal_hcal14";
			oname << "_comp." << pformat;
			can->Print((oname.str()).c_str(),pformat.c_str());
		}
	}
	else can->Close();
	
	return val_graph;
}

//--------------------------------------------------------------------------
//function to compare plots of zcal qtys from previous function
void zcal_multicomp(int qty, int num, bool do_fit, bool do_print=false){
	double energy = jet_energies[num];
	std::stringstream jetname;
	jetname << energy << " GeV jet";

	TGraphErrors* the_graphs[2];
	Double_t* vals[2];

	//get graphs
	the_graphs[0] = zcal_comp(qty,num,do_fit,1,0);
	the_graphs[1] = zcal_comp(qty,num,do_fit,2,0);

	//setup canvas with histo and ratio areas
	stringstream cname;
	cname << "zcal_multicomp_qty" << qty;
	TCanvas* can = new TCanvas((cname.str()).c_str(),(cname.str()).c_str(),700,500);
	can->cd();
	TPad* pad = new TPad("graph","",0,0,1,1);
	pad->SetMargin(0.15,0.05,0.15,0.05);
	pad->Draw();
	pad->cd();

	//legend, pave coords
	double y1 = 0.75;
	double y2 = y1 + 0.1;

	TLegend* leg = new TLegend(0.6,y1,0.9,y2);
	std::string gname[] = {"HCAL 17 layers","HCAL 14 layers"};
	
	Color_t color[2] = {kBlue, kRed};
	Int_t marker[2] = {21, 33};
	
	Double_t ymin = 1e10;
	Double_t ymax = 0;
	
	//find max and min for y range
	for(int i = 0; i < 2; i++){
		vals[i] = the_graphs[i]->GetY();

		for(int j = 0; j < the_graphs[i]->GetN(); j++){
			if(vals[i][j]<ymin) ymin = vals[i][j];
			if(vals[i][j]>ymax) ymax = vals[i][j];
		}
	}
	
	//formatting
	for(int i = 0; i < 2; i++){
		leg->AddEntry(the_graphs[i],(gname[i]).c_str(),"p");
		//if(qty) the_graphs[i]->GetYaxis()->SetRangeUser(0,1);
		//else the_graphs[i]->GetYaxis()->SetRangeUser(0,1.3);

		//formatting
		the_graphs[i]->GetYaxis()->SetRangeUser(ymin*0.9,ymax*1.1);
		the_graphs[i]->SetMarkerStyle(marker[i]);
		the_graphs[i]->SetMarkerColor(color[i]);
		the_graphs[i]->SetMarkerSize(1.5);
		the_graphs[i]->SetLineColor(color[i]);
		the_graphs[i]->SetFillColor(0);
	}
	the_graphs[0]->Draw("APZ");
	the_graphs[1]->Draw("PZ same");

	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.05);
	leg->SetTextFont(42);
	leg->Draw("same");

	TPaveText* pave;
	if(qty==2 || qty==3){
		pave = new TPaveText(0.6,y1-0.06,0.9,y1-0.01,"NDC");
		pave->SetFillColor(0);
		pave->SetBorderSize(0);
		pave->SetTextFont(42);
		pave->SetTextSize(0.05);
		pave->AddText((jetname.str()).c_str());
		pave->Draw("same");
	}

	if(do_print){
		//names
		std::string rtype = "monojet";
		std::string qtyname[] = {"sterm","cterm","mu","sigma","zeroWt","zcalWt"};

		std::stringstream oname;
		oname << pdir << "/" << rtype << "_" << qtyname[qty];
		if(qty==2 || qty==3) oname << "_" << energy << "gev";
		oname << "_zcal_multicomp." << pformat;
		can->Print((oname.str()).c_str(),pformat.c_str());
	}

}

//--------------------------------------------
//macro to make ECAL-HCAL energy "banana plot"
void g4_banana_plot(int snum, int num, int imip, bool do_jet, bool do_print=false){
	gStyle->SetPalette(1);

	int maxE = do_jet ? maxJTe : maxHDe;

	if (num>=maxE || num<0) { std::cout << "num must be between 0 and " << maxE - 1 << std::endl; return; }

	double energy = 0;
	if(do_jet) energy = jet_energies[num];
	else energy = energies[num];
	
	if (imip>2 || imip<0) { std::cout << "imip must be between 0 and 2" << std::endl; return; }
	
	if(do_jet) imip = 2; //no mip options for jets

	Sample* sp = sample_map[snum];
	if(!sp) { std::cout << "Sample " << snum << " is not loaded." << std::endl; return; }
	
	//make filenames
	std::string rtype = "pion";
	if(do_jet) rtype = "monojet";
	std::stringstream fname, piname;
	fname << sp->dir << "/" << sp->fpre << "_" << rtype << "_" << energy << "gev_10k.root";
	if(do_jet) piname << "d jet " << energy << " GeV";
	else piname << "#pi^{-} " << energy << " GeV";

	//open file and tree
	TFile* _file;
	_file = TFile::Open((fname.str()).c_str());
	TTree* totalTree = (TTree*)_file->Get("Total");

	//make tree drawing expressions
	//define mip as ecal < 1 gev = 1000 mev
	double mipcut = 0.8;
	std::stringstream drawname;
	std::stringstream cutname;
	std::stringstream etacut;
	//default histo settings
	//double Emin = 0.1*energies[num]; //lower cut to remove weird peaks near E=zero
	double Emin = 0;
	double Emax = 3*energy;
	if(energy>1000. || do_jet) Emax = 2*energy;
	double Emax_ecal = Emax;
	if(imip==0) Emax_ecal = mipcut;
	int nbins = 100;
	
	//ecal & hcal energies need to be calibrated
	get_sampling_factors(snum);
	drawname << sp->sam << "*(hcal+";
	if(sp->zcalWt>0) drawname << sp->zcalWt << "*zcal+";
	drawname << sp->zeroWt << "*zero)/1000:" << "(" << sp->sam_ecal << "*ecal)/1000>>htemp(" 
			<< nbins << "," << Emin << "," << Emax_ecal << "," << nbins << "," << Emin << "," << Emax << ")";

	//0 is mip, 1 is nomip, 2 is total
	if(imip==0) cutname << sp->sam_ecal << "*ecal/1000 < " << mipcut;
	else if(imip==1) cutname << sp->sam_ecal << "*ecal/1000 > " << mipcut;
	else cutname << "";

	TH2F* h_res; //to store histos drawn from tree
	
	//draw w/ appropriate cut
	totalTree->Draw((drawname.str()).c_str(),(cutname.str()).c_str(),"hist goff");
	h_res = (TH2F*)gDirectory->Get("htemp");

	h_res->SetTitle("");
	h_res->GetXaxis()->SetTitle("ECAL energy [GeV]");
	h_res->GetYaxis()->SetTitle("HCAL energy [GeV]");
	h_res->GetZaxis()->SetTitle("number of events");

	//names
	std::string omip, ofit;
	if (imip==0) omip = "mip";
	else if (imip==1) omip = "nomip";
	else omip = "tot";
	std::stringstream oname;
	oname << pdir << "/" << rtype << "_banana_" << sp->fpre << "_" << omip << "_" << energy << "gev";
	
	std::stringstream Nname, mipname;
	Nname << "N = " << h_res->GetEntries();
	if(imip==0) mipname << "mip, E_{ecal} < " << mipcut << " GeV";
	else if(imip==1) mipname << "nomip, E_{ecal} > " << mipcut << " GeV";
	
	TCanvas* can = new TCanvas((oname.str()).c_str(),(oname.str()).c_str(),800,500);
	TPad* pad1 = new TPad("graph","",0,0,1,1);
	pad1->SetMargin(0.125,0.225,0.15,0.05);
	pad1->SetTicks(1,1);
	pad1->Draw();
	pad1->cd();
	
	h_res->GetYaxis()->SetTitleOffset(0.75);
	h_res->GetZaxis()->SetTitleOffset(1);
	h_res->GetXaxis()->SetTitleOffset(0.85);
	h_res->GetZaxis()->SetTitleSize(38/(pad1->GetWh()*pad1->GetAbsHNDC()));
	h_res->GetZaxis()->SetLabelSize(32/(pad1->GetWh()*pad1->GetAbsHNDC()));
	h_res->GetYaxis()->SetTitleSize(38/(pad1->GetWh()*pad1->GetAbsHNDC()));
	h_res->GetYaxis()->SetLabelSize(32/(pad1->GetWh()*pad1->GetAbsHNDC()));
	h_res->GetXaxis()->SetTitleSize(38/(pad1->GetWh()*pad1->GetAbsHNDC()));
	h_res->GetXaxis()->SetLabelSize(32/(pad1->GetWh()*pad1->GetAbsHNDC()));
	h_res->GetYaxis()->SetTickLength(12/(pad1->GetWh()*pad1->GetAbsHNDC()));
	h_res->GetXaxis()->SetTickLength(12/(pad1->GetWh()*pad1->GetAbsHNDC()));
	
	h_res->Draw("COLZ");
	
	int pavenum = 4;
	if(snum==3||snum==4) pavenum++;
	TPaveText* pave = new TPaveText(0.4,0.85-pavenum*0.05,0.6,0.85,"NDC");
	pave->AddText((piname.str()).c_str());
	if(imip==0 || imip==1) pave->AddText((mipname.str()).c_str());
	pave->AddText((Nname.str()).c_str());
	pave->AddText(sp->gname.c_str());
	if(!sp->supp_info.empty()) pave->AddText(sp->supp_info.c_str());
	//if(snum==3) pave->AddText("(no dead material)");
	//else if(snum==4) pave->AddText("(33 layers)");
	//else if(snum==5) pave->AddText("(extended HE)");
	pave->SetFillColor(0);
	pave->SetBorderSize(0);
	pave->SetTextFont(42);
	pave->SetTextSize(0.05);
	pave->Draw("same");

	TLine* line = new TLine(Emin,energy,energy,Emin);
	//line->SetLineStyle(2);
	line->SetLineWidth(3);
	line->SetLineColor(kBlack);
	line->Draw("same");
	
	if(do_print) can->Print((oname.str()+"."+pformat).c_str(),pformat.c_str());

}

//--------------------------------------------------------------------------------------------------
//macro to draw longitudinal shower distributions (ECAL + ZCAL + HCAL)
//dtype: 0 = nlay, 1 = mm, 2 = lambda
void long_dist(int snum, int num, bool do_jet, int dtype=0, bool do_print=false){
	int maxE = do_jet ? maxJTe : maxHDe;

	if (num>=maxE || num<0) { std::cout << "num must be between 0 and " << maxE - 1 << std::endl; return; }

	double energy = 0;
	if(do_jet) energy = jet_energies[num];
	else energy = energies[num];
	
	Sample* sp = sample_map[snum];
	if(!sp) { std::cout << "Sample " << snum << " is not loaded." << std::endl; return; }
	
	//make filenames
	std::string rtype = "pion";
	if(do_jet) rtype = "monojet";
	std::stringstream fname, piname;
	fname << sp->dir << "/" << sp->fpre << "_" << rtype << "_" << energy << "gev_10k.root";
	if(do_jet) piname << "d jet " << energy << " GeV";
	else piname << "#pi^{-} " << energy << " GeV";

	//open file
	TFile* _file;
	_file = TFile::Open((fname.str()).c_str());
	std::string hname[3] = {"EcalLoShape","ZcalLoShape","HcalLoShape"};
	std::string htitle[3] = {"ECAL layers","HCAL ext. layers","HCAL layers"};
	std::string legname[3] = {sp->gname,"HCAL ext.","HCAL"};
	Color_t colors[3] = {kBlue,kMagenta,kRed};
	TH1F* hist[3];
	double depth_mm[3] = {sp->ecalDepth_mm,sp->zcalDepth_mm,sp->hcalDepth_mm};
	double depth_lam[3] = {sp->ecalDepth_lam,sp->zcalDepth_lam,sp->hcalDepth_lam};
	double depth[3] = {};
	double depth_tot = 0;
	double marginL[3] = {}; //in px
	double marginR[3] = {}; //in px
	double margin_tot = 0; //in px
	
	//open histos
	for(int h = 0; h < 3; h++){
		hist[h] = 0;
		if(h!=1 || (h==1 && sp->zcalWt!=0)) hist[h] = (TH1F*)_file->Get(hname[h].c_str());
		if(!hist[h]) continue;
		
		if(dtype==0) depth[h] = hist[h]->GetNbinsX();
		else if(dtype==1) depth[h] = depth_mm[h];
		else if(dtype==2) depth[h] = depth_lam[h];
		depth_tot += depth[h];
		//marginL[h] = 70;
		marginL[h] = 5;
		marginR[h] = 5;
		margin_tot += marginL[h] + marginR[h];
		//std::cout << "nlay " << h << ": " << nlay[h] << std::endl;
	}
	//std::cout << "nlay_tot: " << nlay_tot << std::endl;
	
	//make canvas with pads relative in size to # of layers per subdet
	double c_width = 900;
	TCanvas* can = new TCanvas("long_dist","long_dist",c_width,500);
	TPad* pad[3];
	//pad location
	double p_width = 0;
	double p_y1 = 0;
	
	//legend
	TLegend* leg = new TLegend(0.5,0.7,0.9,0.85);
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.05);
	leg->SetTextFont(42);
	
	//setup pads and draw histos
	for(int h = 0; h < 3; h++){
		//pad location
		p_width = (marginL[h] + marginR[h] + depth[h]/depth_tot*(c_width-margin_tot))/c_width;

		can->cd();
		pad[h] = 0;
		if(!hist[h]) continue;

		std::stringstream padname;
		padname << "p" << h+1;
		pad[h] = new TPad((padname.str()).c_str(),"",p_y1,0,min(p_y1+p_width,1.),1);
		pad[h]->SetMargin(marginL[h]/(pad[h]->GetWw()*pad[h]->GetAbsWNDC()),marginR[h]/(pad[h]->GetWw()*pad[h]->GetAbsWNDC()),0.05,0.05);
		pad[h]->SetTicks(1,1);
		pad[h]->Draw();
		pad[h]->cd();
		//std::cout << "p_y1 = " << p_y1 << ", p_width = " << p_width << ", Ww = " << pad[h]->GetWw() << ", AbsWNDC = " << pad[h]->GetAbsWNDC() << std::endl;
		
		//histo formatting
		hist[h]->SetTitle("");
		//hist[h]->GetXaxis()->SetTitle(htitle[h].c_str());
		hist[h]->GetXaxis()->SetTitle("");
		hist[h]->GetYaxis()->SetTitle("");
		hist[h]->GetYaxis()->SetRangeUser(0,hist[h]->GetMaximum()*1.05);
		hist[h]->GetYaxis()->SetTitleSize(32/(pad[h]->GetWw()*pad[h]->GetAbsWNDC()));
		hist[h]->GetYaxis()->SetLabelSize(28/(pad[h]->GetWw()*pad[h]->GetAbsWNDC()));
		hist[h]->GetXaxis()->SetTitleSize(32/(pad[h]->GetWw()*pad[h]->GetAbsWNDC()));
		hist[h]->GetXaxis()->SetLabelSize(28/(pad[h]->GetWw()*pad[h]->GetAbsWNDC()));
		hist[h]->GetYaxis()->SetTickLength(6/(pad[h]->GetWw()*pad[h]->GetAbsWNDC()));
		hist[h]->GetXaxis()->SetTickLength(6/(pad[h]->GetWh()*pad[h]->GetAbsHNDC()));
		hist[h]->GetXaxis()->SetLabelOffset(999);
		hist[h]->GetYaxis()->SetLabelOffset(999);
		hist[h]->SetFillColor(colors[h]);
		//hist[h]->GetXaxis()->SetTitleOffset(5/(pad[h]->GetWw()*pad[h]->GetAbsWNDC()));
		//hist[h]->GetYaxis()->SetTitleOffset(5/(pad[h]->GetWw()*pad[h]->GetAbsWNDC()));
		
		//draw on current pad
		leg->AddEntry(hist[h],legname[h].c_str(),"f");
		hist[h]->Draw("hist");
		
		//hack to remove zcal unfilled bin
		if(h==1) {
			hist[h]->GetXaxis()->SetRange(1,hist[h]->GetNbinsX()-1);
			pad[h]->Update();
		}
		
		//increment y1 position for next pad
		p_y1 += p_width;
	}
	
	//draw legend on last pad
	leg->Draw("same");
	
	//pave with scaling info
	std::string dname[3] = {"scaled by number of layers","scaled by depth in mm","scaled by depth in #lambda_{0}"};
	TPaveText* pave = new TPaveText(0.5,0.49,0.9,0.69,"NDC");
	pave->AddText("1/E dE/d[layer]");
	pave->AddText("(energy from sensitive layers)");
	pave->AddText(dname[dtype].c_str());
	pave->AddText((piname.str()).c_str());
	pave->SetFillColor(0);
	pave->SetBorderSize(0);
	pave->SetTextFont(42);
	pave->SetTextSize(0.05);
	pave->Draw("same");

	if(do_print){
		//names
		std::string rtype = "pion";
		if(do_jet) rtype = "monojet";

		std::stringstream oname;
		oname << pdir << "/" << rtype << "_" << energy << "gev_longdist_" << sp->fpre << "." << pformat;
		can->Print((oname.str()).c_str(),pformat.c_str());
	}
	
}

//-----------------------------------------------------------------------------
//function to compare multiple plots of error ellipses for jet resolution terms
void err_multicomp(std::vector<int> snum, bool do_print=false){
	int maxE = maxJTe;

	if(snum.size() > 6) std::cout << "Note: more markers and colors need to be added (only 6 currently defined)" << std::endl;
	
	TGraph** the_graphs = new TGraph*[snum.size()];
	TGraph** the_centers = new TGraph*[snum.size()];
	Double_t** errs = new Double_t*[snum.size()];
	Sample** sp = new Sample*[snum.size()];
	
	//get samples
	for(int i = 0; i < snum.size(); i++){
		sp[i] = sample_map[snum[i]];
		if(!sp[i]) { std::cout << "Sample " << snum[i] << " is not loaded." << std::endl; return; }
		g4_plot_res(snum[i],2,1,0,1,0); //fit jet resolution
		the_graphs[i] = sp[i]->err_ellipse;
		
		//make separate graph for center point
		Double_t* xt = new Double_t[1];
		xt[0] = sp[i]->sampling_term;
		Double_t* yt = new Double_t[1];
		yt[0] = sp[i]->constant_term;
		the_centers[i] = new TGraph(1,xt,yt);
	}
	
	//setup canvas
	TCanvas* can = new TCanvas("err_comp","err_comp",900,500);
	can->cd();
	TPad* pad = new TPad("graph","",0,0,1,1);
	pad->SetMargin(0.125,0.35,0.15,0.075);
	pad->SetTicks(1,1);
	pad->Draw();
	pad->cd();

	//legend, pave coords
	double x1 = 0.675;
	double y2 = 0.9;
	double y1 = y2 - 0.05 - 0.1*(snum.size());

	TLegend* leg = new TLegend(x1,y1,x1+0.325,y2-0.05);
	
	Color_t color[] = {kBlack, kBlue, kMagenta, kRed, kYellow+2, kOrange+7};
	Int_t marker[] = {20, 21, 22, 23, 33, 29};
	
	//axis ranges
	double ymin = 1e10;
	double ymax = -1e10;
	double xmin = 1e10;
	double xmax = -1e10;
	
	for(int i = 0; i < snum.size(); i++){
		Double_t* xvals = the_graphs[i]->GetX();
		Double_t* yvals = the_graphs[i]->GetY();
		Int_t npts = the_graphs[i]->GetN();
		
		for(int j = 0; j < npts; j++){
			if(xvals[i]<xmin) xmin = xvals[i];
			if(xvals[i]>xmax) xmax = xvals[i];
			if(yvals[i]<ymin) ymin = yvals[i];
			if(yvals[i]>ymax) ymax = yvals[i];
		}

		string legname;
		if(!sp[i]->supp_info.empty()) legname = "#splitline{" + sp[i]->gname + "}{" + sp[i]->supp_info + "}";
		else legname = sp[i]->gname;
		leg->AddEntry(the_centers[i],(legname).c_str(),"p");
		//if(qty) the_graphs[i]->GetYaxis()->SetRangeUser(0,1);
		//else the_graphs[i]->GetYaxis()->SetRangeUser(0,1.3);

		//formatting
		the_centers[i]->SetMarkerStyle(marker[i]);
		the_centers[i]->SetMarkerColor(color[i]);
		the_centers[i]->SetMarkerSize(1.0);
		the_graphs[i]->SetLineColor(color[i]);
		the_graphs[i]->SetLineWidth(2);
	}
	
	//histo for axis ranges
	//std::cout << "axis range: (" << xmin << "," << ymin << "), (" << xmax << "," << ymax << ")" << std::endl;
	TH1F* h_axis = new TH1F("axis","",1,xmin*0.95,xmax*1.05);
	h_axis->GetYaxis()->SetRangeUser(ymin*0.95,ymax*1.05);
	//TGaxis::SetMaxDigits(2); //force scientific notation
	h_axis->GetYaxis()->SetTitleOffset(0.9);
	h_axis->GetXaxis()->SetTitleSize(32/(pad->GetWh()*pad->GetAbsHNDC()));
	h_axis->GetXaxis()->SetLabelSize(28/(pad->GetWh()*pad->GetAbsHNDC()));
	h_axis->GetYaxis()->SetTitleSize(32/(pad->GetWh()*pad->GetAbsHNDC()));
	h_axis->GetYaxis()->SetLabelSize(28/(pad->GetWh()*pad->GetAbsHNDC()));
	h_axis->GetYaxis()->SetTickLength(12/(pad->GetWh()*pad->GetAbsHNDC()));
	h_axis->GetXaxis()->SetTickLength(12/(pad->GetWh()*pad->GetAbsHNDC()));
	h_axis->GetXaxis()->SetTitle("sampling term");
	h_axis->GetYaxis()->SetTitle("constant term");
	h_axis->Draw();
	
	//draw ellipses and centers
	for(int i = 0; i < snum.size(); i++){
		the_graphs[i]->Draw("l same");
		the_centers[i]->Draw("p same");
	}

	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.05);
	leg->SetTextFont(42);
	leg->Draw("same");

	if(do_print){
		//names
		std::stringstream oname;
		oname << pdir << "/jet_res_err_ellipse_comp_";
		for(int i = 0; i < snum.size(); i++){
			oname << sp[i]->fpre;
			if(i < snum.size()-1) oname << "_";
		}
		oname << "." << pformat;
		can->Print((oname.str()).c_str(),pformat.c_str());
	}

}

//-------------------------------------------------------
//overload of above function for case of only 2 samples
void err_multicomp(int s1, int s2, bool do_print=false){
	std::vector<int> snum;
	snum.push_back(s1);
	snum.push_back(s2);
	err_multicomp(snum,do_print);
}

//--------------------------------------------------------
//macro to run err_multicomp with predefined set of samples
void run_err_multicomp(bool do_print=false){
	std::vector<int> snum;
	snum.push_back(0);
	snum.push_back(3);
	snum.push_back(5);
	//snum.push_back(6);
	snum.push_back(7);
	
	err_multicomp(snum,do_print);
}