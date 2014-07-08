#ifndef HepMCAnalyzer_h
#define HepMCAnalyzer_h

// system include files
#include <memory>
#include <string>
#include <vector>

// user include files
//#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TROOT.h"

//HepMC headers
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

class HepMCAnalyzer : public edm::EDAnalyzer {
	public:
		explicit HepMCAnalyzer(const edm::ParameterSet&);
		~HepMCAnalyzer();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
		
	private:
		virtual void beginJob();
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob();
	
		virtual void beginRun(edm::Run const&, edm::EventSetup const&);
		virtual void endRun(edm::Run const&, edm::EventSetup const&);
		virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
		virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

		//member variables
		std::vector<HepMC::IO_GenEvent*> asciiOutput;
		std::string outname;
		int nparts;
		int ctr;

};

//define this as a plug-in
DEFINE_FWK_MODULE(HepMCAnalyzer);

#endif