//
// ====================================================================

#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"

// ====================================================================

class TH1D;
class TH2D;
class TTree;

class HistoMessenger;
class DetectorConstruction;

 const G4int  nhist = 13; 

class HistoManager 
{

  public:

    HistoManager(DetectorConstruction*);
   ~HistoManager();

    void Book(G4double);
    void Clear();
    void Save();

    void FillEnergy(G4double, G4double, G4double, G4double, G4double, G4double, G4double, G4double, G4double);
    void FillRange(G4double, G4double);

    void FillHcalTransShape(G4double*);
    void FillHcalLongShape(G4double*);

    void FillZcalTransShape(G4double*);
    void FillZcalLongShape(G4double*);

    void FillTransShape(G4double*);
    void FillLongShape(G4double*);

    void FillAbsTransShape(G4double*);
    void FillAbsLongShape(G4double*);

    void FillCells(G4int, G4double*);
    void FillEcalTransHits(G4double*);

    void SetFileName(G4String);

    void SetSensLBining(G4ThreeVector);
    void SetSensRBining(G4ThreeVector);

    void SetAbsLBining(G4ThreeVector);
    void SetAbsRBining(G4ThreeVector);

    void SetHcalRBining(G4ThreeVector);

    void SetEcalResponse(G4ThreeVector);
    void SetEcalCellNoise(G4double);

    void SetJobRunNumber(G4int);

    G4int       GetnLtot()           {return nLtot;};
    G4int       GetnRtot()           {return nRtot;};
    G4double    GetdLbin()           {return dLbin;};
    G4double    GetdRbin()           {return dRbin;};

    G4int       GetAbsnLtot()        {return nLtotAbs;};
    G4int       GetAbsnRtot()        {return nRtotAbs;};
    G4double    GetAbsdLbin()        {return dLbinAbs;};
    G4double    GetAbsdRbin()        {return dRbinAbs;};

    G4int       GetHcalnRtot()       {return nRtotHcal;};
    G4double    GetHcaldRbin()       {return dRbinHcal;};
    G4int       GetZcalnRtot()       {return nRtotZcal;};
    G4double    GetZcaldRbin()       {return dRbinZcal;};

    G4int       GetJobRunNumber()    {return RunNumber;};
 
    G4String    GetfileName()        {return fileName;};

  private:
    DetectorConstruction* detCon;
	G4int nLayEcal, nLayHcal, nLayZcal;

    TH1D*  histo[nhist];
    TH2D*  hits;

    TTree*    tree_tot;
    TTree*    tree_vec;
    TTree*    tree_ran;
    TTree*    tree_cell;

    G4double  EdepEcalRad, EdepEcalLong;
    G4double  EdepAbsRad, EdepAbsLong;
    G4double  EdepHcalRad, EdepHcalLong;
    G4double  EdepZcalRad, EdepZcalLong;
    G4double  EdepEcalHits;

    G4double  e_ecal, e_hcal, e_zero, e_zcal, e_abs;
    G4double  e_ecal05, e_hcal05, e_zcal05, e_zero05;
    G4double  r_hcal, r_ecal;
    G4double* e_vec;

    G4int     n_cells;
    G4double  e_dep[25], e_phot[25], e_unif[25], e_eff[25];
    G4double  LightYield, LightCollEff, LightCollUnif, CellNoise;

    G4int    nLtot,  nRtot,  nLtotAbs,  nRtotAbs,  nRtotHcal,  nRtotZcal;       
    G4double dLbin,  dRbin,  dLbinAbs,  dRbinAbs,  dRbinHcal,  dRbinZcal;      

    G4int    RunNumber;

    G4String  fileName ;
    HistoMessenger* histoMessenger;
};


#endif
//===============================================================
