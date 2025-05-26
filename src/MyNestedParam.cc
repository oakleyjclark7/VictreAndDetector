#include "MyNestedParam.hh"

#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"


// Constructor
MyNestedParam::MyNestedParam(const std::vector<G4Material*>& materials,
                    size_t* materialIndices,
                    const G4ThreeVector& voxelSize,
                    G4int nX, G4int nY, G4int nZ)
:
    fMaterials(materials),
    fMaterialIndices(materialIndices),
    fdX(voxelSize.x()),
    fdY(voxelSize.y()),
    fdZ(voxelSize.z()),
    fnX(nX),
    fnY(nY),
    fnZ(nZ)
{}


G4Material* MyNestedParam::ComputeMaterial(G4VPhysicalVolume* physVol, const G4int iz,
                const G4VTouchable* parentTouch){
    //G4cout << "Nested: ComputeMaterial called" << G4endl;
    if (parentTouch == nullptr){
        return fMaterials[0];
    }

    // Copy number of voxels.
    // Copy number of X and Y are obtained from replication number.
    // Copy nymber of Z is the copy number of current voxel.
    G4int ix = parentTouch->GetReplicaNumber(0);
    G4int iy = parentTouch->GetReplicaNumber(1);
    G4int copyID = ix + fnX*iy + fnX*fnY*iz;

    // Get the material
    //G4cout << "Nested: ComputeMaterial finished" << G4endl;
    //G4cout << "CopyID: " << copyID << G4endl;
    //G4cout << "Material indices: " << fMaterialIndices.size() << G4endl;
    //G4cout << "Material index: " << fMaterialIndices[copyID] << G4endl;
    //G4cout << fMaterials[fMaterialIndices[copyID]] << G4endl;
    
    G4Material* mat = fMaterials[fMaterialIndices[copyID]];
    
    if (G4VVisManager::GetConcreteInstance() && physVol) {
        // Create visualization attributes (initialized once)
        static G4VisAttributes* matPaddleAttr = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.1));    // Blue, alpha = 0.1
        static G4VisAttributes* matWhiteOpaqueAttr = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 1.0)); // White, fully opaque
        static G4VisAttributes* invisibleAttr = new G4VisAttributes(G4VisAttributes::GetInvisible());   // Invisible

        // Check the material name and apply the corresponding visualization attributes
        if (mat->GetName() == "matPaddle") {
            // If it's the paddle, set blue with alpha = 0.1
            physVol->GetLogicalVolume()->SetVisAttributes(matPaddleAttr);
        }
        else if (mat->GetName() == "matSkin" || mat->GetName() == "matGlandular") {
            // If it's skin or glandular, set white and fully opaque
            physVol->GetLogicalVolume()->SetVisAttributes(matWhiteOpaqueAttr);
        }
        else {
            // Anything else is invisible
            physVol->GetLogicalVolume()->SetVisAttributes(invisibleAttr);
        }
    }

    return mat;

    //return fMaterials[fMaterialIndices[copyID]];
}

void MyNestedParam::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
    // Position of voxels.
    // x and y positions are already defined in DetectorConstruction by using
    // replicated volume. Here only we need to define is z positions of voxels.
    physVol->SetTranslation(G4ThreeVector(0.,0.,(2.*static_cast<double>(copyNo)
                                                +1.)*fdZ - fdZ*fnZ));
}

G4int MyNestedParam::GetNumberOfMaterials() const{
    return G4int(fMaterials.size());
}

G4Material* MyNestedParam::GetMaterial(G4int idx) const{
    return fMaterials[idx];
}


