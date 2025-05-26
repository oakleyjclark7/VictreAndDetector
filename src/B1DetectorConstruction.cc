#include "B1DetectorConstruction.hh"
#include "MyNestedParam.hh"


#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "B2TrackerSD.hh"
#include "G4SDManager.hh"
#include "G4PVParameterised.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <random>
#include <cmath>

#include <zlib.h>

#include <iostream>


// Splits up a string by a delimeter - used to help parse header file
std::vector<std::string> tokenize(std::string const &str, const char delim) 
{ 
    // construct a stream from the string 
    std::stringstream ss(str); 
    
    std::vector<std::string> out;
    std::string sstring; 
    while (std::getline(ss, sstring, delim)) { 
        out.push_back(sstring); 
    } 
    return out;
} 


// Constructor - sets VICTRE_filepath
B1DetectorConstruction::B1DetectorConstruction(const G4String& VICTRE_filepath)
: G4VUserDetectorConstruction(),
  fScoringVolume(0),
  fVICTREFilePath(VICTRE_filepath)
{ }


// Default destructor
B1DetectorConstruction::~B1DetectorConstruction()
{ }

// MAIN Construct method - building the geometry
G4VPhysicalVolume* B1DetectorConstruction::Construct(){

  // Variable name to use for checking overlaps when placing a logical volume
  G4bool checkOverlaps = true;

  // ------------------------------------------------------------

  ///////////////////////// Elements ////////////////////////////

  // For the phantom
  G4Element* elH = 
    new G4Element("Hydrogen","H",1.0,1.008*g/mole);
  G4Element* elC = 
    new G4Element("Carbon","C",6.0,12.01*g/mole);
  G4Element* elN = 
    new G4Element("Nitrogen","N",7.0,14.007*g/mole);
  G4Element* elO = 
    new G4Element("Oxygen","O",8.0,16.0*g/mole);
  G4Element* elNa = 
    new G4Element("Sodium","Na",11.0,22.99*g/mole);
  G4Element* elP = 
    new G4Element("Phosphorus","P",15.0,30.974*g/mole);
  G4Element* elS = 
    new G4Element("Sulfur","P",16.0,32.065*g/mole);
  G4Element* elCl = 
    new G4Element("Chlorine","Cl",17.0,35.453*g/mole);
  G4Element* elK = 
    new G4Element("Potassium","P",19.0,39.098*g/mole);
  G4Element* elFe = 
    new G4Element("Iron","Fe",26.0,55.845*g/mole);

  // For the CZT detector
  G4Element* elCd =
    new G4Element("Cadmium","Cd",48.0,112.4*g/mole);
  G4Element* elZn =
    new G4Element("Zinc","Zn",30.0,65.4*g/mole);
  G4Element* elTe =
    new G4Element("Tellurium","Te",52.0,127.6*g/mole);

  ///////////////////////////////////////////////////////////////

  // ------------------------------------------------------------

  ////////////////////////// MATERIALS //////////////////////////

  // Get nist manager for G4_AIR material
  G4NistManager* nist = G4NistManager::Instance();

  // Air material - Geant4 defualt
  G4Material* matAir = nist->FindOrBuildMaterial("G4_AIR");

  // Adipose Material - from https://iopscience.iop.org/article/10.1088/0031-9155/60/16/N311/pdf
  G4double Adipose_density = 0.9301*g/cm3;
  G4Material* matAdipose =
    new G4Material("matAdipose",Adipose_density,5);
  matAdipose->AddElement(elH,11.2*0.01);
  matAdipose->AddElement(elC,61.9*0.01);
  matAdipose->AddElement(elN,1.7*0.01);
  matAdipose->AddElement(elO,25.1*0.01);
  matAdipose->AddElement(elP,0.1*0.01);

  // Skin Material (same as glandular but different name) - from https://iopscience.iop.org/article/10.1088/0031-9155/60/16/N311/pdf
  G4double Skin_density = 1.04*g/cm3;
  G4Material* matSkin =
    new G4Material("matSkin",Skin_density,5);
  matSkin->AddElement(elH,10.2*0.01);
  matSkin->AddElement(elC,18.4*0.01);
  matSkin->AddElement(elN,3.2*0.01);
  matSkin->AddElement(elO,67.7*0.01);
  matSkin->AddElement(elP,0.5*0.01);

  // Glandular Material - from https://iopscience.iop.org/article/10.1088/0031-9155/60/16/N311/pdf
  G4double Glandular_density = 1.04*g/cm3;
  G4Material* matGlandular =
    new G4Material("matGlandular",Glandular_density,5);
  matGlandular->AddElement(elH,10.2*0.01);
  matGlandular->AddElement(elC,18.4*0.01);
  matGlandular->AddElement(elN,3.2*0.01);
  matGlandular->AddElement(elO,67.7*0.01);
  matGlandular->AddElement(elP,0.5*0.01);

  // Muscle Material - from NIST Muscle, Skeletal (ICRU-44)
  G4double Muscle_density = 1.05*g/cm3;
  G4Material* matMuscle =
    new G4Material("matMuscle",Muscle_density,9);
  matMuscle->AddElement(elH,10.2*0.01);
  matMuscle->AddElement(elC,14.3*0.01);
  matMuscle->AddElement(elN,3.4*0.01);
  matMuscle->AddElement(elO,71.0*0.01);
  matMuscle->AddElement(elNa,0.1*0.01);
  matMuscle->AddElement(elP,0.2*0.01);
  matMuscle->AddElement(elS,0.3*0.01);
  matMuscle->AddElement(elCl,0.1*0.01);
  matMuscle->AddElement(elK,0.4*0.01);

  // Polycarbonate Material (Paddle)
  G4double Paddle_density = 1.20*g/cm3;
  G4Material* matPaddle =
    new G4Material("matPaddle",Paddle_density,3);
  matPaddle->AddElement(elH,5.2*0.01);
  matPaddle->AddElement(elC,71.1*0.01);
  matPaddle->AddElement(elO,23.7*0.01);

  // Blood Material - (ICRU44) https://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html
  G4double Blood_density = 1.06*g/cm3;
  G4Material* matBlood =
    new G4Material("matBlood",Blood_density,10);
  matBlood->AddElement(elH,10.2*0.01);
  matBlood->AddElement(elC,11.0*0.01);
  matBlood->AddElement(elN,3.3*0.01);
  matBlood->AddElement(elO,74.5*0.01);
  matBlood->AddElement(elNa,0.1*0.01);
  matBlood->AddElement(elP,0.1*0.01);
  matBlood->AddElement(elS,0.2*0.01);
  matBlood->AddElement(elCl,0.3*0.01);
  matBlood->AddElement(elK,0.2*0.01);
  matBlood->AddElement(elFe,0.1*0.01);

  // CZT materal  
  G4double CZT_density = 5.8*g/cm3;
  G4Material* matCZT =
    new G4Material("CZT",CZT_density,3);
  matCZT->AddElement(elCd, 43.0*perCent);
  matCZT->AddElement(elZn, 2.77*perCent);
  matCZT->AddElement(elTe, 54.23*perCent);


  // Create a materials vector for the phantom
  std::vector<G4Material*> phantom_materials_vector {matAir,matAdipose,matSkin, matGlandular,
                                            matMuscle,matPaddle,matBlood};

  ///////////////////////////////////////////////////////////////

  // ------------------------------------------------------------

  /////////////////////////// WORLD /////////////////////////////
  
  // Dimensions 60x60x200cm made of air
  G4double world_sizeX = 60*cm;
  G4double world_sizeY = 60*cm;
  G4double world_sizeZ = 200*cm;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  // World solid
  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ);     //its size
  
  // World logical 
  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
  // WORLD PHYSICAL
  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0,0,0),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  ///////////////////////////////////////////////////////////////

  // ------------------------------------------------------------

  ////////////////////////// HEXITEC ////////////////////////////

  // Dimensions
  G4double HEX_side_length = 20*mm;
  G4double HEX_thickness = 2*mm;

  // HEXITEC solid
  G4Box* HEXITEC_solid = 
    new G4Box("HEXITEC",
        0.5*HEX_side_length, 0.5*HEX_side_length, 0.5*HEX_thickness);
  
  // HEXITEC logical
  G4LogicalVolume* HEXITEC_logical = 
    new G4LogicalVolume(HEXITEC_solid,         //its solid
                        matCZT,          //its material==
                        "HEXITEC");           //its name

  // HEXITEC physical
  // Place HEXITEC (always at same position - (0,0,z=+1mm)
  // This places the very front edge of HEXITEC at z=0
  G4VPhysicalVolume* HEXITEC_physical = 
    new G4PVPlacement(0,
                G4ThreeVector(0*mm,0*mm,1*mm),
                HEXITEC_logical,
                "HEXITEC",
                logicWorld,
                false,
                0,
                false);

  ///////////////////////////////////////////////////////////////

  // ------------------------------------------------------------

  ////////////////// Read VICTRE config file ////////////////////

  // Get the VICTRE file IDs from the VICTRE_config.txt file
  std::string line;
  std::string headerfid, imagefid;

  // The VICTRE config file
  std::ifstream inputFile(fVICTREFilePath);

  // Read each line from the file
  while (std::getline(inputFile, line)) {
      // Find the position of the equal sign
      size_t equalPos = line.find('=');

      if (equalPos != std::string::npos) {
          // Extract the variable name and value
          std::string varName = line.substr(0, equalPos);
          std::string varValue = line.substr(equalPos + 1);

          // Trim leading and trailing whitespaces
          varName.erase(0, varName.find_first_not_of(" \t\n\r\f\v"));
          varName.erase(varName.find_last_not_of(" \t\n\r\f\v") + 1);

          varValue.erase(0, varValue.find_first_not_of(" \t\n\r\f\v"));
          varValue.erase(varValue.find_last_not_of(" \t\n\r\f\v") + 1);

          // Use the variable name to decide which variable to assign
          if (varName == "header_file") {
              headerfid = varValue;
          } else if (varName == "image_file") {
              imagefid = varValue;
          }
      }
  }

  ///////////////////////////////////////////////////////////////

  // ------------------------------------------------------------

  ////////////////// Read VICTRE header file ////////////////////

  // Read the VICTRE header file to get dimensions and spacing (read VICTRE geometry)
  std::ifstream headerFile(headerfid);
  G4int xdim, ydim, zdim;
  G4double xspacing, yspacing, zspacing;
  while (std::getline(headerFile,line)){
      
    // find dimension
      if ( line.find("DimSize") != std::string::npos){	
          std::stringstream ss(line);
          std::vector<std::string> dims = tokenize(line,' ');
          xdim = std::stoi(dims[2]);
          ydim = std::stoi(dims[3]);
          zdim = std::stoi(dims[4]);
      }
      // find voxel size - given in mm
      if (line.find("ElementSpacing") != std::string::npos){
          std::stringstream ss(line);
          std::vector<std::string> spacingVec = tokenize(line,' ');
          xspacing = std::stod(spacingVec[2])*mm;
          yspacing = std::stod(spacingVec[3])*mm;
          zspacing = std::stod(spacingVec[4])*mm;
      }
  }

  ///////////////////////////////////////////////////////////////

  // ------------------------------------------------------------

  ////////////////// Read VICTRE image file /////////////////////

  // Read from .raw.gz
  size_t* voxelData = new size_t[xdim * ydim * zdim];
  gzFile imageFile = gzopen(imagefid.c_str(), "rb");
  if (!imageFile) {
      std::cerr << "Error opening .gz file!" << std::endl;
  }
  unsigned char value;
  G4int lineNum = 0;
  while (gzread(imageFile, &value, sizeof(unsigned char)) == sizeof(unsigned char)) {
      voxelData[lineNum] = static_cast<size_t>(value);  // Cast to size_t
      lineNum++;
  }
  gzclose(imageFile);

  ///////////////////////////////////////////////////////////////

  // ------------------------------------------------------------

  /////////// Phantom parameterisation and placement ////////////

  //----- Define the volume that contains all the voxels
  G4VSolid* cont_solid = new G4Box("phantomContainer",xdim*0.5*xspacing,
                               ydim*0.5*yspacing,
                               zdim*0.5*zspacing);
  G4LogicalVolume* cont_logic =
    new G4LogicalVolume( cont_solid,
                         phantom_materials_vector[0], // Material doesn't matter
                         "phantomContainer",
                         0, 0, 0 );
  // make container invisible ---- MAY WISH TO UNDO THIS IF YOU WANT TO SEE CONTAINER OUTLINE      
  cont_logic->SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));

  // Place the phantom container
  G4double minX = -20*mm; // Chest wall inline with the beam at -2cm
  G4double maxZ = -25*mm; // back of the breast paddle 2.5cm from the detector
  G4ThreeVector posCentreVoxels(minX + 0.5*xspacing*xdim,
                                0*mm,
                                maxZ - 0.5*zspacing*zdim);

  G4VPhysicalVolume* cont_phys =
  new G4PVPlacement(0,  // rotation
                    posCentreVoxels,
                    cont_logic,     // The logic volume
                    "phantomContainer",  // Name
                    logicWorld,  // Mother
                    false,           // No op. bool.
                    1);              // Copy number

  // For RunAction - need to know the area of front of phantom, and the distance to it
  front_phantom_area = xdim*0.5*xspacing*ydim*0.5*yspacing;
  G4double source_z = -700*mm;
  front_phantom_distance = (maxZ - zspacing*zdim) - source_z;

  // --- Y Slice
  G4String yRepName("RepY");
  // Solid
  G4VSolid* solYRep = new G4Box(yRepName,xdim*0.5*xspacing,
                                yspacing*0.5,
                                zdim*0.5*zspacing);
  // Logical
  G4LogicalVolume* logYRep = new G4LogicalVolume(solYRep,phantom_materials_vector[0],yRepName);
  // Replicate
  new G4PVReplica(yRepName,logYRep,cont_logic,kYAxis,ydim,yspacing);
  // Set invisible
  logYRep->SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));



  // --- X Slice
  G4String xRepName("RepX");
  // Solid
  G4VSolid* solXRep = new G4Box(xRepName,0.5*xspacing,
                                yspacing*0.5,
                                zdim*0.5*zspacing);
  // Logical
  G4LogicalVolume* logXRep = new G4LogicalVolume(solXRep,phantom_materials_vector[0],xRepName);
  // Replicate
  new G4PVReplica(xRepName,logXRep,logYRep,kXAxis,xdim,xspacing);
  // Set invisible
  logXRep->SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));


  //--- Z Slice
  // Solid
  G4VSolid* solVoxel = new G4Box("phantom",0.5*xspacing,
                                0.5*yspacing,0.5*zspacing);
  // Logical
  G4LogicalVolume* logicVoxel = new G4LogicalVolume(solVoxel,phantom_materials_vector[0],"phantom");
  // Set invisible 
  //logicVoxel->SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));


  // Parameterisation
  G4ThreeVector voxelSize(0.5*xspacing,0.5*yspacing,0.5*zspacing);
  MyNestedParam* param = 
    new MyNestedParam(phantom_materials_vector,
                      voxelData,
                      voxelSize,
                      xdim, ydim, zdim);

  new G4PVParameterised("phantom",    // their name
                        logicVoxel, // their logical volume
                        logXRep,      // Mother logical volume
                        kZAxis,       // Are placed along this axis
                        //kUndefined,
                        // Are placed along this axis
                        zdim,      // Number of cells
                        param);       // Parameterisation.

  ///////////////////////////////////////////////////////////////

  // ------------------------------------------------------------

  /////////////////// Phantom information /////////////////////// ----------------------- NEED TO CHANGE THIS - AIR OR PADDLE - + REDO ALLIGNMENT

  lineNum = 0;
  for (size_t idx = 0; idx < xdim * ydim * zdim; ++idx){
    size_t value = voxelData[idx];
    //G4cout << lineNum << G4endl;
    // If the voxel is not air or paddle
    if (value != 0 && value != 5){

      // Add to voxel count
      nVoxelsTotal++;

      // Find if the voxel is in the detector region
      bool inDetectorRegion = false;
      int i = lineNum % xdim;
      int j = (lineNum / xdim) % ydim;
      int k = lineNum / (xdim * ydim);

      G4double xPosition = i * xspacing + minX; // mm
      G4double minY = -0.5 * ydim * yspacing;
      G4double yPosition = j*yspacing + minY; // mm
      G4double minZ = maxZ - zdim * zspacing;
      G4double zPosition = k*zspacing + minZ; //mm

      //G4cout << "minX: " << minX << ", maxX: " << minX+xdim*xspacing << G4endl;
      //G4cout << "minY: " << minY << ", maxY: " << minY+ydim*yspacing << G4endl;
      //G4cout << "minZ: " << minZ << ", maxZ: " << maxZ <<G4endl;

      //G4cout << "Pos: " << xPosition << ", " << yPosition << ", " << zPosition << G4endl;

      // Find the chest wall angle -> theta
      double theta_vox = std::atan((xPosition - minX)/(zPosition + 700.0));

      // Find the other angle -> phi - absolute value because this angle is symmetric along the zaxis
      double phi_vox = std::atan((std::abs(yPosition))/(zPosition + 700.0));

      // Calculate the min and max theta for region
      double theta_max = std::atan((3.0)/(70.0));
      double theta_min = std::atan((1.0)/(70.0));

      // Calculate the max phi - symmetric so can use absolute value of this angle
      double phi_max = std::atan((1.0)/(70.0));

      //G4cout << "theta voxel: " << theta_vox << ", (" << theta_min << " - " << theta_max << ")" << G4endl;
      //G4cout << "phi voxel: " << phi_vox << ", max -> " << phi_max << G4endl;

      // Check if in region -> angle to voxel less than maximum angles
      if (theta_vox <= theta_max && theta_vox >= theta_min && phi_vox <= phi_max){
        inDetectorRegion = true;
      }

      // If it is in the region
      if (inDetectorRegion) {
          nVoxelsRegion++;
      }

      // If the voxel is glandular
      if (value == 3){
          nGlandularVoxelsTotal++;
          
          // If the voxel is glandular and in region
          if (inDetectorRegion) {
              nGlandularVoxelsRegion++;
          }
      }
    }
    lineNum++;
  }

  // Volume and mass of a voxel
  voxelVolume = solVoxel->GetCubicVolume();
  glandularVoxelMass = Glandular_density * voxelVolume;

  ///////////////////////////////////////////////////////////////

  // ------------------------------------------------------------

  fScoringVolume = HEXITEC_logical;

  // Always return the physical world
  return physWorld;

}


void B1DetectorConstruction::ConstructSDandField()
{
  //Signposting generated detectors as sensitive detectors
  //SD Manager allows hit collection of each detector to be obtained individually
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  B2TrackerSD* aTrackerSD = new B2TrackerSD("tracker1","hitscollection1");
  SetSensitiveDetector("HEXITEC", aTrackerSD, true);
  SDman->AddNewDetector(aTrackerSD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......




