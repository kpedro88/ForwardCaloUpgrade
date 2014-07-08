//
// ********************************************************************
// ********************************************************************
//

#include "DetectorMessenger.hh"

#include <sstream>

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:Detector(Det)
{ 
  ecalDir = new G4UIdirectory("/ecal/");
  ecalDir->SetGuidance("UI commands specific to this example");
  
  detDir = new G4UIdirectory("/ecal/det/");
  detDir->SetGuidance("detector construction commands");

  NbLayersCmd = new G4UIcmdWithAnInteger("/ecal/det/setNbOfLayers",this);
  NbLayersCmd->SetGuidance("Set number of Ecal layers.");
  NbLayersCmd->SetParameterName("NbLayers",false);
  NbLayersCmd->SetRange("NbLayers>0");
  NbLayersCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  HcalAbsMaterCmd = new G4UIcmdWithAString("/ecal/det/setHcalAbsMat",this);
  HcalAbsMaterCmd->SetGuidance("Select Material of the Hcal Absorber.");
  HcalAbsMaterCmd->SetParameterName("choice",false);
  HcalAbsMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  HcalAbsThickCmd = new G4UIcmdWithADoubleAndUnit("/ecal/det/setHcalAbsThick",this);
  HcalAbsThickCmd->SetGuidance("Set Thickness of the Hcal Absorber.");
  HcalAbsThickCmd->SetParameterName("Size",false);
  HcalAbsThickCmd->SetRange("Size>=0. && Size<=400.");
  HcalAbsThickCmd->SetUnitCategory("Length");
  HcalAbsThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  HcalSensMaterCmd = new G4UIcmdWithAString("/ecal/det/setHcalSensMat",this);
  HcalSensMaterCmd->SetGuidance("Select Sensitive Material of Hcal");
  HcalSensMaterCmd->SetParameterName("choice",false);
  HcalSensMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  HcalNbLayersCmd = new G4UIcmdWithAnInteger("/ecal/det/setHcalNbOfLayers",this);
  HcalNbLayersCmd->SetGuidance("Set number of Hcal layers.");
  HcalNbLayersCmd->SetParameterName("HcalNbLayers",false);
  HcalNbLayersCmd->SetRange("HcalNbLayers>=0");
  HcalNbLayersCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  HcalOffsetCmd = new G4UIcmdWithADoubleAndUnit("/ecal/det/setHcalOffset",this);
  HcalOffsetCmd->SetGuidance("Select offset for Hcal");
  HcalOffsetCmd->SetParameterName("Size",false);
  HcalOffsetCmd->SetRange("Size>=400.");
  HcalOffsetCmd->SetUnitCategory("Length");
  HcalOffsetCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  HcalZeroThickCmd = new G4UIcmdWithADoubleAndUnit("/ecal/det/setHcalZeroThick",this);
  HcalZeroThickCmd->SetGuidance("Select zero layer thickness for Hcal");
  HcalZeroThickCmd->SetParameterName("Size",false);
  HcalZeroThickCmd->SetRange("Size>=0. && Size<=13");
  HcalZeroThickCmd->SetUnitCategory("Length");
  HcalZeroThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  EcalAbsMaterCmd = new G4UIcmdWithAString("/ecal/det/setEcalAbsMat",this);
  EcalAbsMaterCmd->SetGuidance("Select Material of the Ecal Absorber.");
  EcalAbsMaterCmd->SetParameterName("choice",false);
  EcalAbsMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  EcalAbsThickCmd = new G4UIcmdWithADoubleAndUnit("/ecal/det/setEcalAbsThick",this);
  EcalAbsThickCmd->SetGuidance("Set Thickness of the Ecal Absorber");
  EcalAbsThickCmd->SetParameterName("Size",false);
  EcalAbsThickCmd->SetRange("Size>=0. && Size<=400.");
  EcalAbsThickCmd->SetUnitCategory("Length");
  EcalAbsThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  EcalSensMaterCmd = new G4UIcmdWithAString("/ecal/det/setEcalSensMat",this);
  EcalSensMaterCmd->SetGuidance("Select Sensitive Material of Ecal");
  EcalSensMaterCmd->SetParameterName("choice",false);
  EcalSensMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  EcalSensThickCmd = new G4UIcmdWithADoubleAndUnit("/ecal/det/setEcalSensThick",this);
  EcalSensThickCmd->SetGuidance("Set Thickness of the Ecal Sensitive");
  EcalSensThickCmd->SetParameterName("Size",false);
  EcalSensThickCmd->SetRange("Size>=0.0");
  EcalSensThickCmd->SetUnitCategory("Length");  
  EcalSensThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ZcalNbLayersCmd = new G4UIcmdWithAnInteger("/ecal/det/setZcalNbOfLayers",this);
  ZcalNbLayersCmd->SetGuidance("Set number of Zcal layers.");
  ZcalNbLayersCmd->SetParameterName("ZcalNbLayers",false);
  ZcalNbLayersCmd->SetRange("ZcalNbLayers>=0");
  ZcalNbLayersCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  ZcalAbsMaterCmd = new G4UIcmdWithAString("/ecal/det/setZcalAbsMat",this);
  ZcalAbsMaterCmd->SetGuidance("Select Material of the Zcal Absorber.");
  ZcalAbsMaterCmd->SetParameterName("choice",false);
  ZcalAbsMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ZcalAbsThickCmd = new G4UIcmdWithADoubleAndUnit("/ecal/det/setZcalAbsThick",this);
  ZcalAbsThickCmd->SetGuidance("Set Thickness of the Zcal Absorber");
  ZcalAbsThickCmd->SetParameterName("Size",false);
  ZcalAbsThickCmd->SetRange("Size>=0. && Size<=400.");
  ZcalAbsThickCmd->SetUnitCategory("Length");
  ZcalAbsThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ZcalSensMaterCmd = new G4UIcmdWithAString("/ecal/det/setZcalSensMat",this);
  ZcalSensMaterCmd->SetGuidance("Select Sensitive Material of Zcal");
  ZcalSensMaterCmd->SetParameterName("choice",false);
  ZcalSensMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ZcalSensThickCmd = new G4UIcmdWithADoubleAndUnit("/ecal/det/setZcalSensThick",this);
  ZcalSensThickCmd->SetGuidance("Set Thickness of the Zcal Sensitive");
  ZcalSensThickCmd->SetParameterName("Size",false);
  ZcalSensThickCmd->SetRange("Size>=0.0 && Size<=400.");
  ZcalSensThickCmd->SetUnitCategory("Length");  
  ZcalSensThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  UpdateCmd = new G4UIcmdWithoutParameter("/ecal/det/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);

  MagFieldCmd = new G4UIcmdWithADoubleAndUnit("/ecal/det/setMagField",this);  
  MagFieldCmd->SetGuidance("Define magnetic field.");
  MagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  MagFieldCmd->SetParameterName("Bz",false);
  MagFieldCmd->SetUnitCategory("Magnetic flux density");
  MagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  NbCellEcalCmd = new G4UIcmdWith3Vector("/ecal/det/setEcalCells",this);
  NbCellEcalCmd->SetGuidance("set number and transverse size Ecal cell");
  NbCellEcalCmd->SetGuidance("nb of cells; transverse cell size [mm]");
  NbCellEcalCmd->SetParameterName("nCells","dxCell"," ",true);
  NbCellEcalCmd->SetRange("nCells>=1 && nCells < 26 && dxCell>0");
  NbCellEcalCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  BirksConsHcalCmd = new G4UIcmdWith3Vector("/ecal/det/setHcalBirks",this);
  BirksConsHcalCmd->SetGuidance("set Hcal Birks constant");
  BirksConsHcalCmd->SetGuidance("birk1; birk2; birk3");
  BirksConsHcalCmd->SetParameterName("birk1","birk2","birk3",true);
  BirksConsHcalCmd->SetRange("birk1>=0 && birk2>=0 && birk3>=0");
  BirksConsHcalCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  BirksConsEcalCmd = new G4UIcmdWith3Vector("/ecal/det/setEcalBirks",this);
  BirksConsEcalCmd->SetGuidance("set Ecal Birks constant");
  BirksConsEcalCmd->SetGuidance("birk1; birk2; birk3");
  BirksConsEcalCmd->SetParameterName("birk1","birk2","birk3",true);
  BirksConsEcalCmd->SetRange("birk1>=0 && birk2>=0 && birk3>=0");
  BirksConsEcalCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  BirkL3ConsEcalCmd = new G4UIcmdWith3Vector("/ecal/det/setEcalBirkL3",this);
  BirkL3ConsEcalCmd->SetGuidance("set Ecal Birk L3 constants");
  BirkL3ConsEcalCmd->SetGuidance("birk1; birk2; birk3");
  BirkL3ConsEcalCmd->SetParameterName("birk1","birk2","birk3",true);
  BirkL3ConsEcalCmd->SetRange("birk1>=0 && birk2>=0 && birk3>=0");
  BirkL3ConsEcalCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DeadThickCmd = new G4UIcmdWithADoubleAndUnit("/ecal/det/setDeadThick",this);
  DeadThickCmd->SetGuidance("Set Thickness of the Dead Material");
  DeadThickCmd->SetParameterName("Size",false);
  DeadThickCmd->SetRange("Size>=0. && Size<=800.");
  DeadThickCmd->SetUnitCategory("Length");
  DeadThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NewDeadThickCmd = new G4UIcmdWithADoubleAndUnit("/ecal/det/setNewDeadThick",this);
  NewDeadThickCmd->SetGuidance("Set Thickness of the New Dead Material");
  NewDeadThickCmd->SetParameterName("Size",false);
  NewDeadThickCmd->SetRange("Size>=0. && Size<=200.");
  NewDeadThickCmd->SetUnitCategory("Length");
  NewDeadThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  PreshowerCmd = new G4UIcmdWithABool("/ecal/det/preshower",this);
  PreshowerCmd->SetGuidance("Enable or disable preshower");
  PreshowerCmd->SetParameterName("preshower",false);
  PreshowerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete NbLayersCmd;
  delete EcalAbsMaterCmd;
  delete EcalAbsThickCmd; 
  delete EcalSensMaterCmd;
  delete EcalSensThickCmd;
  delete HcalAbsMaterCmd;
  delete HcalSensMaterCmd;
  delete HcalZeroThickCmd;
  delete HcalOffsetCmd;
  delete HcalNbLayersCmd;
  delete ZcalNbLayersCmd;
  delete ZcalAbsMaterCmd;
  delete ZcalAbsThickCmd; 
  delete ZcalSensMaterCmd;
  delete ZcalSensThickCmd;
  delete DeadThickCmd;
  delete NewDeadThickCmd;
  delete PreshowerCmd;
  delete MagFieldCmd;
  delete NbCellEcalCmd;
  delete BirksConsHcalCmd;
  delete BirksConsEcalCmd;
  delete BirkL3ConsEcalCmd;
  delete UpdateCmd;
  delete detDir;
  delete ecalDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == NbLayersCmd )
   { Detector->SetNbOfEcalLayers(NbLayersCmd->GetNewIntValue(newValue));}

  if( command == HcalAbsMaterCmd )
   { Detector->SetHcalAbsMaterial(newValue);}

  if( command == HcalAbsThickCmd )
   { Detector->SetHcalAbsThickness(HcalAbsThickCmd->GetNewDoubleValue(newValue));}
   
  if( command == HcalSensMaterCmd )
   { Detector->SetHcalSensMaterial(newValue);}

  if( command == HcalNbLayersCmd )
   { Detector->SetNbOfHcalLayers(HcalNbLayersCmd->GetNewIntValue(newValue));}   

  if( command == HcalOffsetCmd )
   { Detector->SetHcalOffset(HcalOffsetCmd->GetNewDoubleValue(newValue));}
   
  if( command == HcalZeroThickCmd )
   { Detector->SetHcalZeroThickness(HcalZeroThickCmd->GetNewDoubleValue(newValue));}

  if( command == EcalSensMaterCmd )
   { Detector->SetEcalSensMaterial(newValue);}

  if( command == EcalAbsMaterCmd )
   { Detector->SetEcalAbsMaterial(newValue);}

  if( command == EcalAbsThickCmd )
   { Detector->SetEcalAbsThickness(EcalAbsThickCmd->GetNewDoubleValue(newValue));}
   
  if( command == EcalSensThickCmd )
   { Detector->SetEcalSensThickness(EcalSensThickCmd->GetNewDoubleValue(newValue));}

  if( command == ZcalNbLayersCmd )
   { Detector->SetNbOfZcalLayers(ZcalNbLayersCmd->GetNewIntValue(newValue));}   

  if( command == ZcalSensMaterCmd )
   { Detector->SetZcalSensMaterial(newValue);}

  if( command == ZcalAbsMaterCmd )
   { Detector->SetZcalAbsMaterial(newValue);}

  if( command == ZcalAbsThickCmd )
   { Detector->SetZcalAbsThickness(ZcalAbsThickCmd->GetNewDoubleValue(newValue));}
   
  if( command == ZcalSensThickCmd )
   { Detector->SetZcalSensThickness(ZcalSensThickCmd->GetNewDoubleValue(newValue));}
   
  if( command == DeadThickCmd )
   { Detector->SetDeadThickness(DeadThickCmd->GetNewDoubleValue(newValue));}

  if( command == NewDeadThickCmd )
   { Detector->SetNewDeadThickness(NewDeadThickCmd->GetNewDoubleValue(newValue));}
   
  if( command == PreshowerCmd )
   { Detector->SetPreshower(PreshowerCmd->GetNewBoolValue(newValue));}
   
  if( command == MagFieldCmd )
   { Detector->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));}
   
  if( command == NbCellEcalCmd )
   { Detector->SetEcalCells(NbCellEcalCmd->GetNew3VectorValue(newValue));}
   
  if( command == BirksConsHcalCmd )
   { Detector->SetHcalBirksConstant(BirksConsHcalCmd->GetNew3VectorValue(newValue));}
   
  if( command == BirksConsEcalCmd )
   { Detector->SetEcalBirksConstant(BirksConsEcalCmd->GetNew3VectorValue(newValue));}
   
  if( command == BirkL3ConsEcalCmd )
   { Detector->SetEcalBirkL3Constant(BirkL3ConsEcalCmd->GetNew3VectorValue(newValue));}
   
  if( command == UpdateCmd )
   { Detector->UpdateGeometry(); }
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
