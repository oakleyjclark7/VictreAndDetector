#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4VPhysicalVolume.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class B1DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    B1DetectorConstruction(const G4String& VICTRE_filepath);
    virtual ~B1DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
    //G4LogicalVolume* GetScoringDoseVolume() const { return fScoringDoseVolume; }

    G4int GetnGlandularVoxelsTotal() const { return nGlandularVoxelsTotal;}
    G4int GetnVoxelsTotal() const { return nVoxelsTotal;}
    G4int GetnGlandularVoxelsRegion() const {return nGlandularVoxelsRegion;}
    G4int GetnVoxelsRegion() const {return nVoxelsRegion;}

    G4double GetGlandularVoxelMass() const {return glandularVoxelMass;}
    G4double GetVoxelVolume() const {return voxelVolume;}

    // Front of phantom information for calculating primary flux
    G4double GetFrontPhantomDistance() const {return front_phantom_distance;}
    G4double GetFrontPhantomArea() const {return front_phantom_area;}

  protected:
    G4LogicalVolume*  fScoringVolume;
    //G4LogicalVolume*  fScoringDoseVolume;

  private:
    G4String fVICTREFilePath;

    G4int nGlandularVoxelsTotal = 0; // Number of glandular voxels
    G4int nVoxelsTotal = 0; // Total number of breast voxels (not air or paddle)
    G4int nGlandularVoxelsRegion = 0; // Number of glandular voxels within the HEXITEC cross section
    G4int nVoxelsRegion = 0; // Total number of breast voxels (not air or paddle) within the HEXITEC cross section

    G4double glandularVoxelMass;
    G4double voxelVolume;

    // Front of phantom info for calculating flux
    G4double front_phantom_distance;
    G4double front_phantom_area;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif