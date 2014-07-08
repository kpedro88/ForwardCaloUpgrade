//
// ********************************************************************
// ********************************************************************
//

#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin,
                     HistoManager* histo):Det(det),Kin(kin), myana(histo)
{}

RunAction::~RunAction()
{}


void RunAction::BeginOfRunAction(const G4Run* aRun)
{ 

// Set Job Run Number from input cards

  G4int JobRunNumber = myana->GetJobRunNumber();
  ((G4Run *)(aRun))->SetRunID(JobRunNumber);

  G4cout << "### Run " << JobRunNumber << " start." << G4endl;
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

//inform the runManager to save random number seed

  G4RunManager::GetRunManager()->SetRandomNumberStore(true);

// histos cleaning and after definition

  myana-> Clear();

  G4double Ekin = Kin->GetParticleGun()->GetParticleEnergy();

  myana-> Book(Ekin);

// initialize cumulative quantities

  sumEcal  = sumHcal  = sumZero  = sumZcal  = 0.;
  sum2Ecal = sum2Hcal = sum2Zero = sum2Zcal = 0.;

}


void RunAction::fillPerEvent(G4double eEcal, G4double eHcal, G4double eZero, G4double eZcal)
{

// accumulate statistics

  sumEcal += eEcal;  sum2Ecal += eEcal*eEcal;
  sumHcal += eHcal;  sum2Hcal += eHcal*eHcal;
  sumZero += eZero;  sum2Zero += eZero*eZero;
  sumZcal += eZcal;  sum2Zcal += eZcal*eZcal;

}

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;

// save histos

  myana-> Save();
  
// compute statistics: mean and rms

  sumEcal /= NbOfEvents; sum2Ecal /= NbOfEvents;
  G4double rmsEcal = sum2Ecal - sumEcal*sumEcal;
  if (rmsEcal >0.) rmsEcal = std::sqrt(rmsEcal); else rmsEcal = 0.;
  
  sumHcal /= NbOfEvents; sum2Hcal /= NbOfEvents;
  G4double rmsHcal = sum2Hcal - sumHcal*sumHcal;
  if (rmsHcal >0.) rmsHcal = std::sqrt(rmsHcal); else rmsHcal = 0.;

  sumZero /= NbOfEvents; sum2Zero /= NbOfEvents;
  G4double rmsZero = sum2Zero - sumZero*sumZero;
  if (rmsZero >0.) rmsZero = std::sqrt(rmsZero); else rmsZero = 0.;

  sumZcal /= NbOfEvents; sum2Zcal /= NbOfEvents;
  G4double rmsZcal = sum2Zcal - sumZcal*sumZcal;
  if (rmsZcal >0.) rmsZcal = std::sqrt(rmsZcal); else rmsZcal = 0.;
  
// print
  
  G4cout
     << "\n--------------------End of Run------------------------------\n"
     << "\n mean Energy in ECAL : " << G4BestUnit(sumEcal,"Energy")
     << " +- "                      << G4BestUnit(rmsEcal,"Energy")  
     << "\n mean Energy in HCAL : " << G4BestUnit(sumHcal,"Energy")
     << " +- "                      << G4BestUnit(rmsHcal,"Energy")
     << "\n mean Energy in Zero : " << G4BestUnit(sumZero,"Energy")
     << " +- "                      << G4BestUnit(rmsZero,"Energy");
	 
	 if(sumZcal>0){
	    G4cout << "\n mean Energy in Zcal : " << G4BestUnit(sumZcal,"Energy")
               << " +- "                      << G4BestUnit(rmsZcal,"Energy");
	 }
	 
     G4cout << "\n------------------------------------------------------------\n"
     << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
