//
// ********************************************************************
// ********************************************************************
//

#include "EventAction.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VTrajectory.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

// Constructor of EventAction-class and
// assignment of pinter runAct(run) <=> runAct=run

EventAction::EventAction(RunAction* run,DetectorConstruction* det,HistoManager* hist)
:runAct(run),detCon(det),myana(hist),
dEdL(0), dEdR(0), dEdRHcal(0), dEdLHcal(0), dEdRZcal(0), dEdLZcal(0), dEdLAbs(0), dEdRAbs(0), 
RangeEcalLay(0), RangeHcalLay(0), dECellsEcal(0), dEHitsEcal(0),
printModulo(20)
{}

EventAction::~EventAction()
{}

// Member function which run at the start of each event

void EventAction::BeginOfEventAction(const G4Event* evt)
{  

// get bin number for histograms
//------------------------------

  nLayEcal  = detCon->GetNbOfEcalLayers();
  if(nLayEcal!=1) nLayEcal += 1; //+1 for zero layer
  nLayHcal  = detCon->GetNbOfHcalLayers();
  nLayZcal  = detCon->GetNbOfZcalLayers()+1; //+1 for zero layer
  nRtot    = myana->GetnRtot();
  nLtot    = myana->GetnLtot();
  if( nLayEcal != 1 ) nLtot = nLayEcal;
  nRtoth   = myana->GetHcalnRtot();
  nRtotz   = myana->GetHcalnRtot();
  nRtotAbs = myana->GetAbsnRtot();
  nLtotAbs = myana->GetAbsnLtot();
  if( nLayEcal != 1 ) nLtotAbs = nLayEcal-1; //no abs zero layer
  nEcalCells   = detCon->GetNbOfEcalCells();

// initialize dynamic bin arrays
//-----------------------------
  dEdL     = new G4double[nLtot];
  dEdR     = new G4double[nRtot];
  dEdLAbs  = new G4double[nLtotAbs];
  dEdRAbs  = new G4double[nRtotAbs];
  dEdRHcal = new G4double[nRtoth];
  dEdLHcal = new G4double[nLayHcal];
  dEdRZcal = new G4double[nRtotz];
  dEdLZcal = new G4double[nLayZcal];
  RangeEcalLay = new G4double[nLayEcal];
  RangeHcalLay = new G4double[nLayHcal];
  dECellsEcal  = new G4double[nEcalCells];
  dEHitsEcal  = new G4double[nRtot*nRtot];

  G4int evtNb = evt->GetEventID();
  evtNbOld    = evt->GetEventID();

  if (evtNb%printModulo == 0) { 
    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
//    CLHEP::HepRandom::showEngineStatus();
  }
 
// initialisation per event

  EnergyEcal = EnergyHcal = EnergyZero = EnergyZcal = EnergyAbs = 0. ;
  EnergyEcal05 = EnergyHcal05 = EnergyZero05 = EnergyZcal05 = 0. ;
  RangeHcal = RangeEcal = 0. ;
  for (G4int i=0; i<nLtot; ++i)    { dEdL[i] = 0.; }
  for (G4int j=0; j<nRtot; ++j)    { dEdR[j] = 0.; } 
  for (G4int k=0; k<nRtoth; ++k)   { dEdRHcal[k] = 0.; }
  for (G4int k=0; k<nRtotz; ++k)   { dEdRZcal[k] = 0.; }
  for (G4int l=0; l<nLayHcal; ++l)  { dEdLHcal[l] = 0.; }     
  for (G4int l=0; l<nLayZcal; ++l)  { dEdLZcal[l] = 0.; }     
  for (G4int l=0; l<nLayEcal; ++l)  { RangeEcalLay[l] = 0.; }     
  for (G4int l=0; l<nLayHcal; ++l)  { RangeHcalLay[l] = 0.; }     
  for (G4int m=0; m<nLtotAbs; ++m) { dEdLAbs[m] = 0.; }
  for (G4int n=0; n<nRtotAbs; ++n) { dEdRAbs[n] = 0.; }     
  for (G4int i=0; i<nEcalCells; ++i) { dECellsEcal[i] = 0.; }
  for (G4int j=0; j<nRtot*nRtot; ++j) { dEHitsEcal[j] = 0.; }
}

// Member function at the end of each event

void EventAction::EndOfEventAction(const G4Event* evt)
{

//accumulates statistics
//----------------------
  runAct->fillPerEvent(EnergyEcal, EnergyHcal, EnergyZero, EnergyZcal);
  
// fill histos
//------------
  myana-> FillEnergy(EnergyEcal, EnergyHcal, EnergyZero, EnergyZcal, EnergyAbs, EnergyEcal05, EnergyHcal05, EnergyZero05, EnergyZcal05);

// fill Ecal and Hcal range of charged particles
//----------------------------------------------
  myana-> FillRange(RangeHcal, RangeEcal);

// fill Hcal longitudinal shower profile
//----------------------------------------
  myana-> FillHcalLongShape(dEdLHcal);  

// fill Hcal transverse shower profile
//-------------------------------------
  myana-> FillHcalTransShape(dEdRHcal);
  
// fill Zcal longitudinal shower profile
//----------------------------------------
  myana-> FillZcalLongShape(dEdLZcal);  

// fill Zcal transverse shower profile
//-------------------------------------
  myana-> FillZcalTransShape(dEdRZcal);

// fill Ecal transverse shower profile
//-------------------------------------
  myana-> FillTransShape(dEdR);

// fill longitudinal shower profile
//----------------------------------
  myana-> FillLongShape(dEdL);

// fill Ecal absorber transverse shower profile
//----------------------------------------------
  myana-> FillAbsTransShape(dEdRAbs);
  
// fill Ecal absorber longitudinal shower profile
//------------------------------------------------
  myana-> FillAbsLongShape(dEdLAbs);

// fill Ecal cells energy
//------------------------
  myana-> FillCells(nEcalCells,dECellsEcal);   

// fill Ecal transverse hits energy
//---------------------------------
  myana-> FillEcalTransHits(dEHitsEcal);

//print per event (modulo n)
//---------------------------
  G4int evtNb = evt->GetEventID();

  if (evtNb%printModulo == 0) {
    G4cout << "---> End of event: " << evtNb << G4endl;	

    G4cout
       << " Total deposited energy: " 
       << "  ECAL = " << std::setw(7)
                     << G4BestUnit(EnergyEcal,"Energy")
       << ", HCAL = " << std::setw(7)
                     << G4BestUnit(EnergyHcal,"Energy")
       << ", ZERO = " << std::setw(7)
                     << G4BestUnit(EnergyZero,"Energy")
       << G4endl;
	  
  }

// delete dynamic bin arrays
//-------------------------
  delete dEdL;
  delete dEdR;
  delete dEdLAbs;
  delete dEdRAbs;
  delete dEdRHcal;
  delete dEdLHcal;
  delete dEdRZcal;
  delete dEdLZcal;
  delete RangeEcalLay;  
  delete RangeHcalLay;  
  delete dECellsEcal;
  delete dEHitsEcal;

}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
