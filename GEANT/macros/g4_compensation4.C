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

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <vector>
#include <map>

#define maxHDe 7 //energy points for hadrons (and electrons)
#define nCballPar 7 //#pars for cball fn
#define maxZeroWts 3 //different zero layer weight variations

using namespace TMath;

//energy values
Double_t energies[] = {20., 30., 50., 100., 150., 225., 300.};

enum Detector { Ecal=1, Hcal=2 };

class Sample {
	public:
		Sample() : dir(""), fpre(""), name(""), name_rat(""), sam_elec(1), sam_pion(1), sam_init(false), zeroWt(0.5), det(Ecal), supp_info(""),
				   eh(-1), eh_err(-1), k(-1), k_err(-1) {}
		std::string dir;
		std::string fpre;
		std::string name, name_rat;
		Double_t sam_elec, sam_pion;
		bool sam_init;
		Double_t zeroWt;
		Detector det;
		std::string supp_info;
		Double_t eh, eh_err, k, k_err;
};

std::map<int,Sample*> sample_map;

void add_sample(int snum, std::string dir_, std::string fpre_, std::string name_, std::string name_rat_, Double_t zeroWt_, Detector det_, std::string supp_info_=""){
	Sample* sp = new Sample();
	sp->dir = dir_;
	sp->fpre = fpre_;
	sp->name = name_;
	sp->name_rat = name_rat_;
	sp->zeroWt = zeroWt_;
	sp->det = det_;
	if(!supp_info_.empty()) sp->supp_info = supp_info_;
	sample_map[snum] = sp;
}

void init_samples(int v0=0){
	if(v0<0 || v0>2) {
		cout << "No zero weight for version " << v0 << ". Samples will not be loaded." << endl;
		return;
	}
	//default hcal
	double w0[] = {0.5,0.49,0.68};
	add_sample(0,"/data/users/pedrok/FullSim/compensation","hcal_only","HCAL","HCAL",w0[v0],Hcal);
	//pbwo4 ecal
	add_sample(1,"/data/users/pedrok/FullSim/compensation","pbwo4_only","PbWO_{4} ECAL","PbWO_{4}",0.5,Ecal);
	//w-lyso ecal
	add_sample(2,"/data/users/pedrok/FullSim/compensation","wlyso_only","W-LYSO ECAL","W-LYSO",0.5,Ecal);
}

std::string pdir = "plots";
std::string pformat = "png";

void set_pinfo(std::string pd, std::string pf) { pdir = pd; pformat = pf; }

void set_zeroWt(int snum, Double_t z) { Sample* sp = sample_map[snum]; sp->zeroWt = z; }

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

//--------------------------------------
//function to calculate sampling factors
std::pair<Double_t,Double_t> g4_sample(int snum, Double_t energy, bool do_pion, bool do_show, bool do_print=false, bool set_val=true){
	Sample* sp = sample_map[snum];
	if(!sp) { std::cout << "Sample " << snum << " is not loaded." << std::endl; return std::pair<Double_t,Double_t>(0.,0.); }

	//select correct file
	std::string fpre = sp->fpre;
	if(do_pion) fpre += "_pion";
	else fpre += "_elec";

	//make filenames
	std::stringstream drawname, fname, piname;
	fname << sp->dir << "/" << fpre << "_" << energy << "gev_10k.root";
	if(do_pion) piname << "#pi^{-} " << energy << " GeV";
	else piname << "e^{-} " << energy << " GeV";

	//open file and tree
	TFile* _file;
	_file = TFile::Open((fname.str()).c_str());
	TTree* totalTree = (TTree*)_file->Get("Total");

	//get histo from tree (no display)
	//define mip as sam_ecal*ecal < 1 gev = 1000 mev (for pions in HCAL)
	if(sp->det==Hcal) drawname << "(hcal+" << sp->zeroWt << "*zero)/1000>>hsam(200)";
	else drawname << "(ecal)/1000>>hsam(200)";
	
	totalTree->Draw((drawname.str()).c_str(),"","hist goff");
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
	Double_t f = energy/m;
	Double_t f_err = energy*(me/(m*m));
	s_mean << f << " #pm " << f_err;

	TPolyMarker* pm = new TPolyMarker(1, xpos, ypos);
	hsam->GetListOfFunctions()->Add(pm);
	pm->SetMarkerStyle(23);
	pm->SetMarkerColor(kRed);
	pm->SetMarkerSize(1.3);

	std::cout.precision(6);
	std::cout << "f_" << (do_pion ? "pion" : "elec") << " = " << f << " +/- " << f_err << std::endl;
	
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
		hsam->GetXaxis()->SetTitle("Energy [GeV]");
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

		if(do_print) {
			std::stringstream oname;
			oname << pdir << "/" << fpre << "_sample_" << energy << "gev_peak.png";
			can->Print((oname.str()).c_str(),"png");
		}
	}
	else _file->Close();

	//store value in sample
	if(set_val){
		if(do_pion) sp->sam_pion = f;
		else sp->sam_elec = f;
	}

	return std::pair<Double_t,Double_t>(f,f_err);
}

//---------------------------------------------
//function to get sampling factors for a sample
void get_sampling_factors(int snum, bool do_show = false){
	Sample* sp = sample_map[snum];
	if(!sp) { std::cout << "Sample " << snum << " is not loaded." << std::endl; return; }
	if(!sp->sam_init){
		Double_t energy = 50;
		//Double_t energy = 300;
	
		//sampling factor from elec
		g4_sample(snum,energy,0,0);
		//sampling factor from pion
		g4_sample(snum,energy,1,0);

		sp->sam_init = true;
	}
}

//------------------------------------
//function to fit energy distributions
energyRes* get_res(int snum, Double_t energy, bool do_pion, bool use_f_pion, bool do_fit, bool do_show, bool do_print=false, bool do_batch=false){
	Sample* sp = sample_map[snum];
	if(!sp) { std::cout << "Sample " << snum << " is not loaded." << std::endl; energyRes* theRes = new energyRes(0,0); return theRes; }
	
	//select correct file
	std::string fpre = sp->fpre;
	if(do_pion) fpre += "_pion";
	else fpre += "_elec";

	//make filenames
	std::stringstream drawname, fname, piname;
	fname << sp->dir << "/" << fpre << "_" << energy << "gev_10k.root";
	if(do_pion) piname << "#pi^{-} " << energy << " GeV";
	else piname << "e^{-} " << energy << " GeV";

	//open file and tree
	TFile* _file;
	_file = TFile::Open((fname.str()).c_str());
	TTree* totalTree = (TTree*)_file->Get("Total");

	//default histo settings
	//double Emin = 0.1*energies[num]; //lower cut to remove weird peaks near E=zero
	double Emin = 0;
	double Emax = 2*energy;
	int nbins = 100;
	
	//ecal & hcal energies need to be calibrated
	get_sampling_factors(snum);

	//make tree drawing expressions
	//define mip as ecal < 1 gev = 1000 mev
	if(use_f_pion) drawname << sp->sam_pion;
	else drawname << sp->sam_elec;
	
	if(sp->det==Hcal) drawname << "*(hcal+" << sp->zeroWt << "*zero)/1000";
	else drawname << "*ecal/1000";

	drawname << ">>htemp(" << nbins << "," << Emin << "," << Emax << ")";
	//std::cout << drawname.str() << std::endl;

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
	energyRes* theRes = new energyRes(energy,2);

	//draw w/ appropriate cut
	totalTree->Draw((drawname.str()).c_str(),"","hist goff");
	h_res = (TH1F*)gDirectory->Get("htemp");
	h_res->SetTitle("");
	h_res->GetXaxis()->SetTitle("Energy [GeV]");
	h_res->SetLineWidth(2);
	h_res->SetLineColor(kBlack);

	//names
	std::string ofit;
	if(do_fit) ofit = "fit";
	else ofit = "nofit";
	std::stringstream oname;
	oname << pdir << "/" << fpre;
	if(use_f_pion) oname << "_fpion";
	oname << "_response_" << ofit << "_" << energy << "gev";
	
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
		gfit = new TF1("resp","gaus",0,h_res->GetXaxis()->GetXmax());
		//if(do_jet){
		//	gfit->SetParameters(ph,p,s);
		//	if(m > p) gfit->SetRange(p-1.5*s,p+1.0*s); //high tail
		//	else gfit->SetRange(p-1.0*s,p+1.5*s); //low tail
		//}
		//else{
			gfit->SetParameters((Double_t)N,m,s);
			gfit->SetRange(m-2*s,m+1*s); //fit within 2 std devs
			//if(m > p) gfit->SetRange(p-2*s,p+1*s); //high tail
			//else gfit->SetRange(p-1*s,p+2*s); //low tail
		//}
		
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
		can = new TCanvas((oname.str()).c_str(),(oname.str()).c_str(),700,500);
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
		
		//pave
		pave = new TPaveText(xmin,0.68,xmin+0.2,0.78,"NDC");
		pave->AddText(sp->name.c_str());
		pave->AddText((piname.str()).c_str());
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
		
		std::cout << "response:" << std::endl;
		
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
//qty: 0 = response, 1 = resolution, 2 = sampling factor
TGraphErrors* g4_plot_res(int snum, int qty, bool do_pion, bool use_f_pion, bool do_fit, bool do_show, bool do_print=false){
	Sample* sp = sample_map[snum];
	if(!sp) { std::cout << "Sample " << snum << " is not loaded." << std::endl; return 0; }

	//store values from get_res
	Double_t* vals = new Double_t[maxHDe]; //sigma or mean
	Double_t* xvals = new Double_t[maxHDe]; //sigma or mean
	Double_t* y_errors = new Double_t[maxHDe]; //errors on pars
	Double_t* logxvals = new Double_t[maxHDe]; //sigma or mean

	//for storage of output info
	energyRes* res_temp;

	for (int i = 0; i < maxHDe; i++){
		double energy = energies[i];
	
		Double_t v, ve;
	
		if(qty==2){
			//get sampling factor for energy without setting value
			std::pair<Double_t,Double_t> f_temp = g4_sample(snum,energy,do_pion,0);
			v = f_temp.first;
			ve = f_temp.second;
		}
		else {
			res_temp = get_res(snum,energy,do_pion,use_f_pion,do_fit,0,0,1);

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
		logxvals[i] = log(energy);
		vals[i] = v;
		y_errors[i] = ve;
	}

	TCanvas* can;
	TPaveText* pave;
	TGraphErrors* val_graph;
	TGraphErrors* fit_graph;

	Int_t col, mrk;
	col = kBlue; mrk = 21;

	//graph values
	std::string qtyaxes[] = {"Response (#mu/E_{true})","Resolution (#sigma/#mu)","sampling factor"};
	if(do_pion) qtyaxes[0] = sp->name_rat + " #pi^{-} Response (R_{cal}/E_{gen})";
	else qtyaxes[0] = sp->name_rat + " e^{-} Response (R_{cal}/E_{gen})";
	fit_graph = new TGraphErrors(maxHDe,logxvals,vals,0,y_errors);
	val_graph = new TGraphErrors(maxHDe,xvals,vals,0,y_errors);
	val_graph->GetXaxis()->SetTitle("Energy [GeV]");
	val_graph->GetYaxis()->SetTitle(qtyaxes[qty].c_str());
	val_graph->SetTitle("");
	val_graph->SetMarkerStyle(mrk);
	val_graph->SetMarkerColor(col);
	val_graph->SetMarkerSize(1.5);
	val_graph->SetLineColor(col);
	val_graph->SetFillColor(0);

	//fit response for e/pi from pions
	TF1* gfit = 0;
	TF1* gline = 0;
	double c, ce, k, ke;
	std::stringstream epiname, ehname, kname;
	if(qty==0 && do_pion){
		gfit = new TF1("epi_fit","pol1",fit_graph->GetXaxis()->GetXmin(),fit_graph->GetXaxis()->GetXmax());
		//gfit = new TF1("epi_fit","(1+([0]-1)*([1]*x))/[0]",val_graph->GetXaxis()->GetXmin(),val_graph->GetXaxis()->GetXmax());
		//gfit = new TF1("epi","(1+([0]-1)*([1]*x^[2]))/[0]",val_graph->GetXaxis()->GetXmin(),val_graph->GetXaxis()->GetXmax());
		gfit->SetParameter(0,0.9);
		gfit->SetParameter(1,0.01);
		//gfit->SetParameter(2,-2.8);
		fit_graph->Fit(gfit,"NR");
		
		//results
		double a, ae, b, be;
		a = gfit->GetParameter(0);
		ae = gfit->GetParError(0);
		b = gfit->GetParameter(1);
		be = gfit->GetParError(1);
		
		//transform to desired params
		c = 1./a;
		ce = ae/(a*a);
		k = b/(1.-a);
		ke = sqrt(pow(be/(1.-a),2) + pow(ae*b/pow(1-a,2),2));
		std::cout.precision(2);
		std::cout << "e/h = " << c << " +/- " << ce << ", k = " << k << " +/- " << ke << std::endl;
		
		//store params in sample
		sp->eh = c;
		sp->eh_err = ce;
		sp->k = k;
		sp->k_err = ke;

		epiname.precision(2);
		epiname << "#frac{R}{E} = #frac{#pi}{e}(E) = #frac{1 + (e/h - 1) #upoint k ln(E)}{e/h}";
		ehname.precision(2);
		ehname << "e/h = " << c << " #pm " << ce;
		kname.precision(2);
		kname << "k = " << k << " #pm " << ke;
		
		//line for E instead of log(E)
		gline = new TF1("epi","(1+([0]-1)*([1]*log(x)))/[0]",val_graph->GetXaxis()->GetXmin(),val_graph->GetXaxis()->GetXmax());
		gline->SetParameter(0,c);
		gline->SetParameter(1,k);
		//formatting
		gline->SetLineColor(kRed);
		gline->SetMarkerColor(kRed);
		gline->SetLineWidth(2);
	}
	
	if(do_show){
		std::string cname;
		cname = "res";
		can = new TCanvas(cname.c_str(),cname.c_str(),700,500);
		can->cd();
		//can->SetLogx();

		//if(qty) val_graph->GetYaxis()->SetRangeUser(0,0.4);
		//else val_graph->GetYaxis()->SetRangeUser(0,1.1);
		val_graph->Draw("APZ");

		//legend, pave coords
		double y1;
		if(qty) y1 = 0.5;
		else y1 = 0.2;
		
		std::string pavename = sp->name;
		if(do_pion) pavename += " #pi^{-}";
		else pavename += " e^{-}";
		
		pave = new TPaveText(0.5,y1,0.95,y1+0.2,"NDC");
		if(qty==0 && do_pion){
			pave->AddText((epiname.str()).c_str());
			pave->AddText((ehname.str()).c_str());
			pave->AddText((kname.str()).c_str());
		}
		else{
			pave->AddText(pavename.c_str());
		}
		pave->SetFillColor(0);
		pave->SetBorderSize(0);
		pave->SetTextFont(42);
		pave->SetTextSize(0.05);
		pave->Draw("same");
		
		if(gline) gline->Draw("same");
		
		if(do_print){
			std::string fpre = sp->fpre;
			if(do_pion) fpre += "_pion";
			else fpre += "_elec";
		
			//names
			std::string ofit;
			if(do_fit) ofit = "fit";
			else ofit = "nofit";
			std::string qtyname[] = {"mu","sigma","sam"};
			std::stringstream oname;
			oname << pdir << "/" << fpre;
			if(use_f_pion) oname << "_fpion";
			oname << "_" << qtyname[qty] << "_" << ofit;
			oname << "." << pformat;
			can->Print((oname.str()).c_str(),pformat.c_str());
		}
	}

	return val_graph;
}

void calc_eh(int snum, bool do_fit=true, bool do_print=false){
	Sample* sp = sample_map[snum];
	if(!sp) { std::cout << "Sample " << snum << " is not loaded." << std::endl; return; }

	//get eh from response
	g4_plot_res(snum, 0, 1, 0, do_fit, 1, do_print);
}