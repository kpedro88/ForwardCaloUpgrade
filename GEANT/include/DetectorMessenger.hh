//
// ********************************************************************

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWith3Vector;
class G4UIcmdWithABool;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorMessenger: public G4UImessenger
{
  public:
    DetectorMessenger(DetectorConstruction* );
   ~DetectorMessenger();

    void SetNewValue(G4UIcommand*, G4String);

  private:
    DetectorConstruction* Detector;

    G4UIdirectory*             ecalDir;
    G4UIdirectory*             detDir;
    G4UIcmdWithAnInteger*      NbLayersCmd;
    G4UIcmdWithAnInteger*      ZcalNbLayersCmd;
    G4UIcmdWithAnInteger*      HcalNbLayersCmd;
    G4UIcmdWithAString*        HcalAbsMaterCmd;
    G4UIcmdWithADoubleAndUnit* HcalAbsThickCmd;
    G4UIcmdWithAString*        HcalSensMaterCmd;
    G4UIcmdWithAString*        EcalAbsMaterCmd;
    G4UIcmdWithAString*        EcalSensMaterCmd;
    G4UIcmdWithADoubleAndUnit* EcalAbsThickCmd;
    G4UIcmdWithADoubleAndUnit* EcalSensThickCmd;
    G4UIcmdWithAString*        ZcalAbsMaterCmd;
    G4UIcmdWithAString*        ZcalSensMaterCmd;
    G4UIcmdWithADoubleAndUnit* ZcalAbsThickCmd;
    G4UIcmdWithADoubleAndUnit* ZcalSensThickCmd;
	G4UIcmdWithADoubleAndUnit* DeadThickCmd;
	G4UIcmdWithADoubleAndUnit* NewDeadThickCmd;
	G4UIcmdWithADoubleAndUnit* HcalZeroThickCmd;
	G4UIcmdWithADoubleAndUnit* HcalOffsetCmd;
	G4UIcmdWithABool*          PreshowerCmd;
    G4UIcmdWithADoubleAndUnit* MagFieldCmd;
    G4UIcmdWithoutParameter*   UpdateCmd;
    G4UIcmdWith3Vector*        NbCellEcalCmd;
    G4UIcmdWith3Vector*        BirksConsHcalCmd;
    G4UIcmdWith3Vector*        BirksConsEcalCmd;
    G4UIcmdWith3Vector*        BirkL3ConsEcalCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

