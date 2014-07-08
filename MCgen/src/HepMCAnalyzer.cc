#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

// essentials !!!
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

// Global FWCore clases
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// user include files
#include "ForwardCaloUpgrade/MCgen/interface/HepMCAnalyzer.h"
#include "TROOT.h"
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

//HepMC headers
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

using namespace std;
using namespace edm;

HepMCAnalyzer::HepMCAnalyzer(const edm::ParameterSet& iConfig) : ctr(0) { 
	outname = iConfig.getParameter<string>("fileName");
	nparts = iConfig.getParameter<int>("nparts");
}

HepMCAnalyzer::~HepMCAnalyzer() { }

// ------------ method called for each event  ------------
void
HepMCAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // here's an example of accessing particles in the event record (HepMCProduct)
  //
  Handle< HepMCProduct > EvtHandle ;
  
  // find initial (unsmeared, unfiltered,...) HepMCProduct
  //
  iEvent.getByLabel( "generator", EvtHandle ) ;
  
  const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;

  //split events into nparts files
  asciiOutput[ctr++]->write_event(Evt);
  if(ctr>=nparts) ctr = 0;
}

// ------------ method called once each job just before starting event loop  ------------
void 
HepMCAnalyzer::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HepMCAnalyzer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
HepMCAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&) {
  //open HepMC IO class
  
  //no part number in filename if only 1 part
  if(nparts==1){
    asciiOutput.push_back(new HepMC::IO_GenEvent(outname, std::ios::out));
  }
  else{
	std::string outpre;
	std::size_t pos = outname.find(".dat");
	//std::cout << "outname.length() = " << outname.length() << ", pos = " << pos << ", -> " << outname.length()-pos-1 << std::endl;
	if(pos!=std::string::npos) outpre = outname.substr(0,pos);
	else outpre = outname;
	
	for(int i = 0; i < nparts; i++){
	  std::stringstream outpart;
	  outpart << outpre << "_part" << i+1 << ".dat";
	  //std::cout << outpart.str() << std::endl;
	  asciiOutput.push_back(new HepMC::IO_GenEvent(outpart.str(), std::ios::out));
	}
  }
}

// ------------ method called when ending the processing of a run  ------------
void 
HepMCAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) { 
  //close HepMC IO class
  for(int i = 0; i < nparts; i++){
    delete asciiOutput[i];
  }
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HepMCAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { }

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HepMCAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { }

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HepMCAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



