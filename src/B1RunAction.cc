//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B1RunAction.cc 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
// #include "B1Run.hh"
#include "G4String.hh"

#include "G4ThreeVector.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
//#include "G4ParameterManager.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "B4Analysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction(long int maxEvents, G4String codeName)
  : G4UserRunAction(),
    fGlandularEdep("Edep", 0.),
    fGlandularEdep2("Edep2", 0.),
    event_counter(0),
    fMaxEvents(maxEvents),
    fCodeName{codeName}
{
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;
  const G4double picogray  = 1.e-12*gray;

  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);

  //Creating analyis manager to analyse hit collections
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(0);

  //Setting prefix of all analysed file names
  analysisManager-> SetFileName(fCodeName);

  // Creating Ntuple for hit info
  analysisManager->CreateNtuple("EdepTable","EdepTable");
  analysisManager->CreateNtupleIColumn("EventID");
  analysisManager->CreateNtupleFColumn("x/mm");
  analysisManager->CreateNtupleFColumn("y/mm");
  analysisManager->CreateNtupleFColumn("z/mm");
  analysisManager->CreateNtupleFColumn("Energy/keV");
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run*)
{
  //Informing  runManager to not save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  //Reseting parameters to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

  //Opening analysis files to write to
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager -> OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EndOfRunAction(const G4Run* run)
{

  // Get detector constuction class for geoemtry information
  const B1DetectorConstruction* detectorConstruction
   = static_cast<const B1DetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  // Calculate the primary flux (needed for I0)
  G4int nofEvents = run->GetNumberOfEvent(); // Total primaries
  G4double front_phantom_area = detectorConstruction->GetFrontPhantomArea();
  G4double front_phantom_distance = detectorConstruction->GetFrontPhantomDistance();
  G4double primary_flux_front_phantom = (G4double)nofEvents / front_phantom_area;
  G4double source_detector_distance = 700.0*mm;
  G4double primary_flux_detector = (primary_flux_front_phantom*front_phantom_distance*front_phantom_distance)/(source_detector_distance*source_detector_distance);
  long int I0_primaries = static_cast<long int>(primary_flux_detector*400.0);

  //If no events in run, simulation ends here
  if (nofEvents == 0) return;

  /////////////// WRITE HITS INFORMATION TO FILE ////////////////

  //Merge parameters
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Run conditions
  //note: There is no primary generator action object for "master"
  //run manager for multi-threaded mode.
  const B1PrimaryGeneratorAction* generatorAction
   = static_cast<const B1PrimaryGeneratorAction*>
    (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

  //Writing and closing analysis files
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
  ////////////////////////////////////////////////////////////////


  /////////////// Density / Volume / Dose Calculations /////////////////
  
  // Get the energy deposited in glandular voxels
  G4double edep = fGlandularEdep.GetValue(); ///joule; // MeV to J
  //std::cout << "Glandular Energy = " << edep << " J" << std::endl;
  G4double edep2 = fGlandularEdep2.GetValue(); ///joule; // Mev to J
  G4double rms = edep2 - edep*edep/nofEvents; // Unsure on this but does it even matter?
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;

  // Get the voxel count information from detectorConstruction
  G4int nGlandularVoxelsTotal = detectorConstruction->GetnGlandularVoxelsTotal();
  G4int nVoxelsTotal = detectorConstruction->GetnVoxelsTotal();
  G4int nGlandularVoxelsRegion = detectorConstruction->GetnGlandularVoxelsRegion();
  G4int nVoxelsRegion = detectorConstruction->GetnVoxelsRegion();

  // Calculate the mass of glanduar voxels (all in kg)
  G4double glandularVoxelMass = detectorConstruction->GetGlandularVoxelMass();
  G4double glandularMassTotal = (G4double)nGlandularVoxelsTotal * glandularVoxelMass;


  // Calculate the dose
  G4double meanGlandularDose = edep / glandularMassTotal;
  G4double meanGlandularDoseError = rms / glandularMassTotal;

  // Density Calculations
  G4double phantomDensity = (G4double)nGlandularVoxelsTotal / (G4double)nVoxelsTotal;
  G4double regionDensity = (G4double)nGlandularVoxelsRegion / (G4double)nVoxelsRegion;

  // Volume (cm3)
  G4double voxelVolume = detectorConstruction->GetVoxelVolume();
  G4double phantomVolume = voxelVolume * (G4double)nVoxelsTotal;
  G4double regionVolume = voxelVolume * (G4double)nVoxelsRegion;

  // Information output
  std::ofstream infoFile;
  G4String infoFileID = "info_" + fCodeName + ".txt";
  infoFile.open(infoFileID);

  if (infoFile.is_open()) {
      // Print to both console and file
      G4cout << "===== Phantom Information =====" << G4endl;
      infoFile << "===== Phantom Information =====" << G4endl;

      G4cout << "Phantom Volume (cm3): " << phantomVolume / cm3 << G4endl;
      infoFile << "Phantom Volume (cm3): " << phantomVolume / cm3 << G4endl;

      G4cout << "Region Volume (cm3): "  << regionVolume / cm3 << G4endl;
      infoFile << "Region Volume (cm3): "  << regionVolume / cm3 << G4endl;

      G4cout << "===============================" << G4endl;
      infoFile << "===============================" << G4endl;

      G4cout << "Density (phantom): "  << 100*phantomDensity << G4endl;
      infoFile << "Density (phantom): "  << 100*phantomDensity << G4endl;

      G4cout << "Density (region): "  << 100*regionDensity << G4endl;
      infoFile << "Density (region): "  << 100*regionDensity << G4endl;

      G4cout << "===============================" << G4endl;
      infoFile << "===============================" << G4endl;

      G4cout << "MGD (microgray): " << meanGlandularDose / CLHEP::microgray << G4endl;
      infoFile << "MGD (microgray): " << meanGlandularDose / CLHEP::microgray << G4endl;

      G4cout << "MGD error (microgray): "  << meanGlandularDoseError / CLHEP::microgray << G4endl;
      infoFile << "MGD error (microgray): "  << meanGlandularDoseError / CLHEP::microgray << G4endl;

      G4cout << "===============================" << G4endl;
      infoFile << "===============================" << G4endl;

      G4cout << "Primary detector flux (ph/mm2): " << primary_flux_detector << G4endl;
      infoFile << "Primary detector flux (ph/mm2): " << primary_flux_detector << G4endl;

      G4cout << "I0 primaries (ph): " << I0_primaries << G4endl << G4endl;
      infoFile << "I0 primaries (ph): " << I0_primaries << G4endl;

      infoFile.close();  // Close the file
  } else {
      G4cerr << "Error opening file!" << G4endl;
  }
}

// Accumulating the energy deposited in glandular volumes in each event for dose calculation
// fEdep2 is the square, later used for rms
void B1RunAction::AccumulateGlandularEdep(G4double edep)
{
  fGlandularEdep  += edep;
  fGlandularEdep2 += edep*edep;
}

// Increment the event counter by 1 each time there is a new event detected
void B1RunAction::IncrementEventCounter(){
  event_counter ++;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
