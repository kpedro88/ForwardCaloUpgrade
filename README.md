This documentation describes the FCAL standalone simulation, which approximates the CMS endcap specifications in a simple rectangular geometry.

1) GEOMETRY
-----------------------------------------------------------------------------------------------------------------------------------------------

	The transverse size of all objects in this geometry is 540cm x 540cm.
	The x and y directions are transverse, while the z direction is longitudinal.
	The total size of the "world" is 600cm x 600cm x 2400cm.

	a) The default HCAL Endcap (HE) is made up of 17 layers of 79mm brass absorber and 3.7 mm plastic scintillator.
	Each scintillator layer is contained in an 8mm aluminum wrapper, inside a 9mm air gap behind the absorber layer.
	The total thickness of the default HE is 1496mm. The default offset from the front of the ECAL is 818mm.
	At the front of the HE is a "zero layer", 9mm of plastic scintillator contained in a 15mm aluminum wrapper.
	Birks' Law is simulated for energy deposits in the sensitive layers (including the zero layer), using the Chou-Birks parameterization:
	double rkb = birk1/density;
	double c  = birk2*rkb*rkb;
	if (abs(charge) >= 2.) rkb /= birk3; // based on alpha particle data
	weight = 1./(1.+rkb*dedx+c*dedx*dedx);
	
	Adjustable parameters:
	* number of HCAL layers: 0 or higher. If 0, the HCAL (including the zero layer) is not created.
	  /ecal/det/setHcalNbOfLayers
	* zero layer thickness: 0 to 13mm
	  /ecal/det/setHcalZeroThick
	* absorber material: brass (default) or tungsten
	  /ecal/det/setHcalAbsMat
	* absorber layer thickness: 0 to 400mm
	  /ecal/det/setHcalAbsThick
	* sensitive material: plastic (default), LAG, YAG, LSO, YSO, LYSO, PbWO, LYSOc
	  /ecal/det/setHcalSensMat
	* offset from ecal: > 400mm
	  /ecal/det/setHcalOffset
	* Chou-Birks' law constants: 3 arguments, birks1 birks2 birks3
	  /ecal/det/setHcalBirks
        * Non-uniform light yield: integer case (0, 1, 2, 3...), with function defined in SteppingAction
         /ecal/det/LTnonuniform
	  
	b) Optionally, an HE extension can be added. This extension is a sampling calorimeter placed between the HE and the zero layer.
	It is called "ZCAL" in the input commands. It can only be added if the old dead material is disabled.
	Note that the HCAL offset is not automatically changed when the ZCAL is turned on, therefore the user must be careful to set
	the total ZCAL thickness so that it will fit in the gap between ECAL and HCAL.
	Birks' Law is simulated for energy deposits in the sensitive layers, using the HCAL constants.
	
	Adjustable parameters:
	* number of ZCAL layers: 0 or higher. If 0, the ZCAL is not created.
	  /ecal/det/setZcalNbOfLayers
	* absorber material: tungsten (default), etc.
	  /ecal/det/setZcalAbsMat
	* absorber layer thickness: 0 to 400mm
	  /ecal/det/setZcalAbsThick
	* sensitive material: plastic (default), LAG, YAG, LSO, YSO, LYSO, PbWO, LYSOc
	  /ecal/det/setZcalSensMat
	* sensitive layer thickness: 0 to 400mm
	  /ecal/det/setZcalSensThick
	  
	c) Dead material is simulated between the ECAL and the HCAL.
	The old dead material is made of a material containing copper, hydrogen, oxygen, and carbon (to simulate cables).
	The default thickness is 234mm = 0.6 interaction lengths, based on studies of the material budget in CMSSW.
	The old dead material is placed to be centered between the ECAL and HCAL.
	Alternatively, a new dead material made of aluminum and placed directly behind the ECAL can be simulated.
	(This simulates the "strongback" for the Shashlik ECAL.)
	It will only be turned on if the old dead material is turned off.
	
	Adjustable parameters:
	* Old dead material thickness: 0 to 800mm (234mm default). If 0, the old dead material is not created.
	  /ecal/det/setDeadThick
	* New dead material thickness: 0 (default) to 200mm. Only created if the old dead material thickness is 0.
	  /ecal/det/setNewDeadThick
	  
	d) The default ECAL Endcap (EE) is a homogeneous volume of PbWO4 (lead tungstate) crystal, 220mm thick.
	It is offset from the interaction point (IP) by a distance of 3154mm.
	The default configuration for the Shashlik ECAL (an alternative calorimeter for the Phase 2 upgrade) is
	28 layers of 2.5 mm tungsten absorber and 29 layers of 1.5 mm LYSO scintillator.
	There is a 2.5mm G10 plate behind the ECAL and a 4.5mm aluminum support in front of the ECAL.
	Birks' Law is simulated for energy deposits in the sensitive layers (including the zero layer), using the L3 parameterization:
	weight = 1. - birkSlope*log(rkb*dedx);
	if (weight < birkCut) weight = birkCut;
	else if (weight > 1.) weight = 1.;

	Adjustable parameters:
	* number of ECAL layers: 1 or higher. 1 is a homogeneous calorimeter, 2+ is a sampling calorimeter.
	  /ecal/det/setNbOfLayers
	* absorber material: lead (default), brass, tungsten, etc.
	  /ecal/det/setEcalAbsMat
	* absorber layer thickness: 0 to 400mm
	  /ecal/det/setEcalAbsThick
	* sensitive material: PbWO (default), LAG, YAG, LSO, YSO, LYSO, LYSOc
	  /ecal/det/setEcalSensMat
	* sensitive layer thickness: 0+ mm (unlimited to allow "infinite" ECALs for e/h compensation measurements)
	  /ecal/det/setEcalSensThick
	* L3 Birks' law constants: 3 arguments, on/off birkSlope birkCut
	  /ecal/det/setEcalBirkL3
	* Chou-Birks' law constants: (used if L3 is turned off) 3 arguments, birks1 birks2 birks3
	  /ecal/det/setEcalBirks
	* Number of cells for segmented readout: two arguments, nCells dxCell
	  0 < nCells < 26, dxCell > 0
	  /ecal/det/setEcalCells
	
	e) The preshower detector (ES) is simulated in front of the ECAL. It is enabled by default but can be turned off. It consists of:
	27.2mm aluminum, 165mm in front of EE
	13mm lead, 187.5mm in front of EE
	2.5mm G10 plate, 196.25mm in front of EE
	
	Adjustable parameters:
	* enable or disable (1 or 0)
	  /ecal/det/preshower
	  
	f) A uniform magnetic field can be simulated along the z-axis. It is disabled by default.
	
	Adjustable parameters:
	* field strength: in units of kG, 0 (default) disables
	  /ecal/det/setMagField

2) OUTPUT
-----------------------------------------------------------------------------------------------------------------------------------------------

	a) Several histograms are booked and filled. These are:
	* HCAL deposited energy (sensitive layers, MeV)
	* HCAL transverse shower profile (mm)
	* HCAL longitudinal shower profile (layer)
	
	* Zero layer deposited energy (MeV)
	
	* ZCAL transverse shower profile (mm)
	* ZCAL longitudinal shower profile (layer)
	
	* ECAL deposited energy (GeV)
	* ECAL transverse shower profile (mm)
	* ECAL longitudinal shower profile (mm or layer)
	* ECAL Moliere radius (mm)
	* ECAL absorber transverse shower profile (mm)
	* ECAL absorber longitudinal shower profile (mm or layer)
	* ECAL absorber Moliere radius (mm)
	* ECAL transverse hit point distributions (mm vs. mm)
	
	b) Several trees are filled. These are:
	* Total (energy per event), with branches:
	  * ecal (ECAL sensitive layer energy, MeV)
	  * hcal (HCAL sensitive layer energy, MeV)
	  * zero (zero layer energy, MeV)
	  * zcal (ZCAL sensitive layer energy, MeV)
	  * ecal05 (ECAL sensitive layer energy in a cone phi<0.5, MeV)
	  * hcal05 (HCAL sensitive layer energy in a cone phi<0.5, MeV)
	  * zero05 (zero layer energy in a cone phi<0.5, MeV)
	  * zcal05 (ZCAL sensitive layer energy in a cone phi<0.5, MeV)
	  * abse (ECAL absorber layer energy, MeV)
	
	* Range (total range of charged particles in sensitive layers), with branches:
	  * hran (HCAL range)
	  * eran (ECAL range)
	  
	* Vector (energy deposits for each layer of HCAL), with branches:
	  * nLayHcal (layer number)
	  * e_vec (energy, MeV)
	
	* Cell (energy deposits in cells of ECAL), with branches:
	  * n_cells (number of cells)
	  * e_dep (deposited energy, MeV)
	  * e_phot (number of photons)
	  * e_unif (number of photons with uniform light collection)
	  * e_eff (effective number of photons)
	  
	c) Adjustable parameters:
	* Output ROOT file name
	  /test/histo/setRootName
	* job run number (necessary to save random number seed)
	  /test/histo/setRunNumber
	* binning for HCAL transverse shower profile: 2 arguments, nbins binwidth (mm)
	  /test/histo/setHcalRbin
	* binning for ECAL transverse shower profile: 2 arguments, nbins binwidth (mm)
	  /test/histo/setSensRbin
	* binning for ECAL longitudinal shower profile (homogeneous only): 2 arguments, nbins binwidth (mm)
	  /test/histo/setSensLbin
	* binning for ECAL absorber transverse shower profile: 2 arguments, nbins binwidth (mm)
	  /test/histo/setAbsRbin
	* binning for ECAL absorber longitudinal shower profile (homogeneous only): 2 arguments, nbins binwidth (mm)
	  /test/histo/setAbsLbin

3) EVENT GENERATION
-----------------------------------------------------------------------------------------------------------------------------------------------
	
	a) The default event generation uses the particle gun.
	See http://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides/ForApplicationDeveloper/html/AllResources/Control/UIcommands/_gun_.html
	
	b) There is also an option to simulate complex physics events in the HepMC ASCII format.
	The folder ForwardCaloUpgrade/MCgen shows an example of how to generate monojet events using Pythia in CMSSW,
	and convert the gen-level info in the output ROOT file to the HepMC format in a .dat text file.
	To use this option, remove the particle gun settings from the input configuration, and replace them with these commands:
	/generator/select hepmcAscii
	/generator/hepmcAscii/open your_file.dat
	/generator/hepmcAscii/verbose 0

4) MISCELLANEOUS OPTIONS
-----------------------------------------------------------------------------------------------------------------------------------------------

	a) Production cuts can be set for all media or for specific media:
	* all media: cut (mm)
	  /run/setCut
	* specific medium: 2 arguments, material cut (mm)
	  /run/setCutForRegion
	  
	b) The default physics list is QGSP_BERT_EMV.
	A different physics list can be specified by using a second command line argument (after the input configuration),
	or by defining the environment variable PHYSLIST.
	scripts/turn_on_hp.sh shows an example script to switch to the QSGP_BERT_HP physics list.
	Note that some physics lists may require external data files.

5) RUNNING
-----------------------------------------------------------------------------------------------------------------------------------------------

	These instructions assume that you have an installation of Geant4.
	The scripts in setup/ and the GNUmakefile will need to be modified to point to the correct directories for your Geant4 installation.
	(Both .csh and .sh versions of the scripts in setup/ are provided.)
	Some information about installing Geant4 can be found at:
	https://twiki.cern.ch/twiki/pub/CMS/FCALSimSLHCStandalone/geant4_instructions.txt

	cmsrel CMSSW_4_2_8
	cd CMSSW_4_2_8/src
	cmsenv
	git clone git@github.com:kpedro88/ForwardCaloUpgrade
	cd ForwardCaloUpgrade/MCgen
	scram b
	cd ../GEANT
	source setup/setup.csh
	source setup/setup_hepmc.csh
	source setup/setup_geant_opengl.csh
	gmake
	$G4BIN/Linux-g++/fcalor hadr01.in > hadr01.log

6) EXAMPLES
-----------------------------------------------------------------------------------------------------------------------------------------------

	In the main GEANT directory, hadr01.in is a usable example configuration (as shown above).
	If the executable is run without any command-line arguments, it will open the interactive Geant4 terminal, and (if configured correctly)
	create a visualization of the geometry specified in vis.mac.

	The folder configs/ contains a large number of template input files.
	These template input files demonstrate different geometries and types of events.
	To make usable input files, the template files are given to the scripts in scripts/, which fill in certain pieces of information
	(e.g. particle energy).
	The user should modify any directories in the files in configs/ or scripts/ as needed.
	
	The different geometry configurations are as follows:
	(default ZCAL layer configuration is 18mm tungsten absorber, 3.7mm plastic scintillator)
	* hcal: default ECAL + default HCAL
	* hcal_only: default HCAL, particle gun position set to skip ECAL
	* pbwo4_only: "infinite" (20 interaction lengths) PbWO4 ECAL, no HCAL or preshower or dead material
	* wlyso3: W-LYSO ECAL + default HCAL, no preshower, no dead material
	* wlyso5: W-LYSO ECAL + default HCAL, no preshower, new dead material (6cm Al)
	* wlyso7: W-LYSO ECAL + HCAL shifted forward 30cm w/ 20 layers, no preshower, new dead material (6cm Al)
	* wlyso_zcal: W-LYSO ECAL + default HCAL + default ZCAL with 1-10 layers, no preshower, new dead material (6cm Al)
	* wlyso_zcal_hcal14: W-LYSO ECAL + HCAL w/ 14 layers + default ZCAL with 1-10 layers, no preshower, new dead material (6cm Al)
	* wlyso_only: "infinite" (20 interaction lengths) W-LYSO ECAL, no HCAL or preshower or dead material
	* wlyso_whcal: W-LYSO ECAL + HCAL with 18 layers of 49.7mm tungsten absorber and default scintillator, no preshower, new dead material (6cm Al)
	* whcal_only: HCAL with 18 layers of 49.7mm tungsten absorber and default scintillator, particle gun position set to skip ECAL
	
	In scripts/, the scripts *temp* use these template files to create usable input files for the simulation.
	The scripts *sub* loop over the *temp* scripts.
	These scripts are currently configured for batch submission to Condor.

7) MACROS
-----------------------------------------------------------------------------------------------------------------------------------------------

	The folder macros/ shows some ROOT macros that can be used to analyze the output ROOT files from the simulation.
	Documentation of these macros is in progress. (The code is commented to some degree.)

8) OTHER VERSIONS
-----------------------------------------------------------------------------------------------------------------------------------------------

	A twiki page with info about the standalone simulation can be found here:
	https://twiki.cern.ch/twiki/bin/viewauth/CMS/FCALSimSLHCStandalone
	
	The version preceding this version can be found here:
	https://github.com/cms-fcal-upgrade/ForwardCaloUpgrade/tree/master/GEANT_PROJECTS/Shashlik
	
	A parallel version, developed by Todd Adams et al. at FSU, can be found here:
	https://github.com/cms-fcal-upgrade/ForwardCaloUpgrade/tree/master/GEANT_PROJECTS/UpdatedShashlik
	This version supports the simulation of pileup events, as well as improved segmented readout for both ECAL and HCAL.
