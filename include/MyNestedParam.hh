#include "G4VNestedParameterisation.hh"

#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4VTouchable.hh"

#include <vector>

class G4VPhysicalVolume;

class MyNestedParam : public G4VNestedParameterisation
{
    public:

        // Constructor
        MyNestedParam(const std::vector<G4Material*>& materials,
                    size_t* materialIndices,
                    const G4ThreeVector& voxelSize,
                    G4int nX, G4int nY, G4int nZ);

        G4Material* ComputeMaterial(G4VPhysicalVolume *currentVol,
                              const G4int repNo,
                              const G4VTouchable *parentTouch);

        void ComputeTransformation(const G4int no,
                    G4VPhysicalVolume* currentPV) const;

        G4int GetNumberOfMaterials() const;

        G4Material* GetMaterial(G4int idx) const;

    private:

        G4double fdX, fdY, fdZ; // Voxel half dimensions!!! HALF
        G4int fnX, fnY, fnZ; // Number of voxels
        std::vector<G4Material*> fMaterials; // Materials vector
        size_t* fMaterialIndices; // Index in materials corresponding to each voxel

};