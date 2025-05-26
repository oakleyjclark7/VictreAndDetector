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
// $Id: B1EventAction.cc 93886 2015-11-03 08:28:26Z gcosmo $
//


#include "B1EventAction.hh"
#include "B1RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "B4Analysis.hh"

#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "B2TrackerHit.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4UImanager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::B1EventAction(B1RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fTHC1ID(-1),
  fEdepGlandular(0.0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::BeginOfEventAction(const G4Event* event)
{
  //If 1st event in run define location of hit collection
  if (fTHC1ID==-1){
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    fTHC1ID = sdManager->GetCollectionID("tracker1/hitscollection1");
 }

  // Set fEdepGlandular to zero (this will track energy depositied in glandular tissue)
  fEdepGlandular = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::EndOfEventAction(const G4Event* event)
{

  // Report energy deposited in glandular volumes to accumulate in run action
  fRunAction->AccumulateGlandularEdep(fEdepGlandular);

  //Retrieving hit collection
  G4HCofThisEvent* hce = event->GetHCofThisEvent();

  //Retrieving detector specific hit collections
  B2TrackerHitsCollection* tHC1
    = static_cast<B2TrackerHitsCollection*>(hce->GetHC(fTHC1ID));

  //Processing results
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  //Processing detector 1 results
  G4int n_hit = tHC1->entries();
  for (G4int i=0; i<n_hit;i++)
   {
    // Get hits information
     B2TrackerHit* hit = (*tHC1)[i];
     G4ThreeVector localPos = hit->GetPos() / 1000.0; // Convert from mm to m
     G4double Edeposited = hit->GetEdep() * 1000.0; // change to kev
     G4int EventID = hit->GetEventID();

     // Fill the ntuple with the hits information
     analysisManager->FillNtupleIColumn(0,EventID);
     analysisManager->FillNtupleFColumn(1,localPos[0] + 0.01); // Shift x by 1cm so goes from 0cm to 2cm
     analysisManager->FillNtupleFColumn(2,localPos[1] + 0.01); // Shift y by 1cm so goes from 0cm to 2cm
     analysisManager->FillNtupleFColumn(3,localPos[2]);
     analysisManager->FillNtupleFColumn(4,Edeposited);
     analysisManager->AddNtupleRow();
   }

  // If there were any detector hits, then its a recorded event, increment count
  if (n_hit > 0){
    // Increment counter
    fRunAction->IncrementEventCounter();

    // Check whether there have been enough events
    if (fRunAction->get_event_counter() >= fRunAction->get_max_events()){
      G4UImanager* UImanager = G4UImanager::GetUIpointer();
      UImanager->ApplyCommand("/run/abort");
    }

    // Progress
    if (fRunAction->get_event_counter()%1000 == 0){
      // Calculate progress
      double percentComplete = 100.0*(double)fRunAction->get_event_counter()/(double)fRunAction->get_max_events();
      G4cout << percentComplete << "%" << G4endl; 
    }

  }
}

// Accumlate energy deposited in glandular volumes during this event
void B1EventAction::AddGlandularEdep(G4double glandularEdep){
    fEdepGlandular += glandularEdep; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
