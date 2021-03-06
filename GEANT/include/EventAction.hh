//
// ********************************************************************
// ********************************************************************
//

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
//#include "HistoManager.hh"
//#include "DetectorConstruction.hh"

class RunAction;
class DetectorConstruction;
class HistoManager;

// Determination of EventAction-class

class EventAction : public G4UserEventAction
{
public:
  EventAction(RunAction*, DetectorConstruction*, HistoManager*);
  virtual ~EventAction();

  void  BeginOfEventAction(const G4Event*);
  void    EndOfEventAction(const G4Event*);
    
  void AddEcal(G4double de) { EnergyEcal += de; };
  void AddEcal05(G4double de) { EnergyEcal05 += de; };
  void AddZero(G4double de) { EnergyZero += de; };
  void AddZero05(G4double de) { EnergyZero05 += de; };
  void AddZcal(G4double de) { EnergyZcal += de; };
  void AddZcal05(G4double de) { EnergyZcal05 += de; };  
  void AddHcal(G4double de) { EnergyHcal += de; };   
  void AddHcal05(G4double de) { EnergyHcal05 += de; };
  void AddAbs(G4double de)  { EnergyAbs  += de; };

  void AddHcalRange(G4double dl, G4int lay_hcal) {
                                 RangeHcal += dl;
                                 RangeHcalLay[lay_hcal] +=dl; };

  void AddEcalRange(G4double dl, G4int lay_ecal) {
                                 RangeEcal += dl;
                                 RangeEcalLay[lay_ecal] +=dl; };

  void fillEcalStep(G4double dEstep, G4int Lbin, G4int Rbin) {
                                     dEdL[Lbin] += dEstep; 
                                     dEdR[Rbin] += dEstep; };

  void fillHcalStep(G4double dEstep, G4int Lbin, G4int Rbin) {
                                     dEdLHcal[Lbin] += dEstep;
                                     dEdRHcal[Rbin] += dEstep; };

  void fillZcalStep(G4double dEstep, G4int Lbin, G4int Rbin) {
                                     dEdLZcal[Lbin] += dEstep;
                                     dEdRZcal[Rbin] += dEstep; };
									 
  void fillAbsStep(G4double dEstep, G4int Lbin, G4int Rbin) {
                                    dEdLAbs[Lbin] += dEstep;
                                    dEdRAbs[Rbin] += dEstep; };

  void fillEcalCell(G4int Lcell, G4double dEstep) {dECellsEcal[Lcell] += dEstep;};

  void fillEcalHits(G4int Lhits, G4double dEstep) {dEHitsEcal[Lhits] += dEstep;};

  G4int   GetEventNb()              {return evtNbOld;};

private:

   RunAction*  runAct;
   DetectorConstruction* detCon;
   HistoManager* myana;

   G4double    EnergyEcal, EnergyHcal, EnergyZero, EnergyZcal, EnergyAbs;
   G4double    EnergyEcal05, EnergyHcal05, EnergyZcal05, EnergyZero05;
   G4double    RangeHcal,  RangeEcal;
   G4int       nRtot, nLtot, nRtoth, nRtotz;
   G4int       nLayEcal, nLayHcal, nLayZcal;
   G4int       nRtotAbs, nLtotAbs;
   G4int       nEcalCells;
   G4double    *dEdL, *dEdR, *dEdRHcal, *dEdLHcal, *dEdRZcal, *dEdLZcal;
   G4double    *dEdLAbs, *dEdRAbs;
   G4double    *RangeEcalLay, *RangeHcalLay;
   G4double    *dECellsEcal, *dEHitsEcal;
   G4int       printModulo;                     
   G4int       evtNbOld;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
