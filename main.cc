#include "B1DetectorConstruction.hh"
#include "B1ActionInitialization.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "PhysicList.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4String.hh"

#include "Randomize.hh"
#include <string>
#include <fstream>
#include <sstream>

int main(int argc, char** argv) {
    // Ensure at least 1 argument (VICTRE file path) is provided
    if (argc < 2) {
        G4cerr << "Error: VICTRE file must be provided as the first argument." << G4endl;
        G4cerr << "Usage: " << argv[0] << " <VICTRE_filepath> [Macro_filepath]" << G4endl;
        return 1;  // Exit with an error if no VICTRE file is provided
    }

    // Detect interactive mode (if only one argument - VICTRE file) and define UI session
    G4UIExecutive* ui = nullptr;

    // File paths from the command-line arguments
    G4String VICTRE_filepath = argv[1];  // The first argument (VICTRE file path)
    G4String Macro_filepath = "";        // The second argument (macro file) - optional

    if (argc == 2) {
        // Only the VICTRE file is provided, run in interactive mode
        ui = new G4UIExecutive(argc, argv);
    } else if (argc == 3) {
        // Both VICTRE and Macro file paths are provided, use batch mode
        Macro_filepath = argv[2];
    } else {
        // Too many arguments provided
        G4cerr << "Error: Too many arguments." << G4endl;
        G4cerr << "Usage: " << argv[0] << " <VICTRE_filepath> [Macro_filepath]" << G4endl;
        return 1;
    }

    // Choose the Random engine
    G4Random::setTheEngine(new CLHEP::RanecuEngine);

    // Construct the run manager  
    G4RunManager* runManager = new G4RunManager;
    runManager->SetVerboseLevel(0);

    // Read the max number of events from the VICTRE file
    long int maxEvents = 1000;  // Default value for max number of events
    std::ifstream file(VICTRE_filepath);
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            // Check if the line contains "n_target_detector_events"
            if (line.find("n_target_detector_events") != std::string::npos) {
                // Find the position of the equal sign
                size_t pos = line.find("=");
                if (pos != std::string::npos) {
                    // Extract the value after the equal sign
                    std::string value = line.substr(pos + 1);
                    
                    // Trim any leading/trailing spaces and convert to integer
                    std::stringstream valueStream(value);
                    valueStream >> maxEvents;
                }
            }
        }
        file.clear();
        file.seekg(0, std::ios::beg);
        //file.close();
    } else {
        G4cerr << "Unable to open VICTRE file: " << VICTRE_filepath << G4endl;
        return 1;  // Exit if the file can't be opened
    }

    // Set mandatory initialization classes
    // Pass the VICTRE file path to the Detector Construction
    runManager->SetUserInitialization(new B1DetectorConstruction(VICTRE_filepath));

    // Physics list
    G4VModularPhysicsList* physicsList = new QBBC;
    physicsList->SetVerboseLevel(0);
    runManager->SetUserInitialization(physicsList);

    // Find the codename for file saving
    //file.open();
    if (!file.is_open()) {
        G4cerr << "Error: Could not open file: " << VICTRE_filepath << G4endl;
        return 1;
    }
    G4String line;
    G4String headerFilePath;
    // Loop through the lines to find the one containing "header_file"
    while (std::getline(file, line)) {
        if (line.find("header_file") != G4String::npos) {
            // Extract the file path after "header_file ="
            size_t pos = line.find('=');
            if (pos != G4String::npos) {
                headerFilePath = line.substr(pos + 1);  // Get the part after '='
                headerFilePath.erase(0, headerFilePath.find_first_not_of(" \t"));  // Trim leading whitespace
            }
            break;
        }
    }
    file.close();
    if (headerFilePath.empty()) {
        G4cerr << "Error: 'header_file' not found in the file." << G4endl;
        return 1;
    }

    // Extract the directory name from headerFilePath
    size_t lastSlash = headerFilePath.rfind('/');
    if (lastSlash == G4String::npos) {
        G4cerr << "Error: headerFilePath contains no directory. File is in current folder." << G4endl;
        return 1;
    }
    size_t secondLastSlash = headerFilePath.rfind('/', lastSlash - 1);
    G4String codeName;
    if (secondLastSlash == G4String::npos) {
        // Only one slash — file is in a direct subdirectory
        codeName = headerFilePath.substr(0, lastSlash);
    } else {
        // Extract the directory name between the two slashes
        codeName = headerFilePath.substr(secondLastSlash + 1, lastSlash - secondLastSlash - 1);
    }



    // User action initialization
    runManager->SetUserInitialization(new B1ActionInitialization(maxEvents,codeName));

    // Initialize visualization
    G4VisManager* visManager = new G4VisExecutive("Quiet");
    visManager->Initialize();

    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    // Process macro or start UI session
    if (!ui) {
        // Batch mode (Macro file provided)
        G4String command = "/control/execute ";
        UImanager->ApplyCommand(command + Macro_filepath);  // Use the macro file path
    } else {
        // Interactive mode (only VICTRE file provided)
        UImanager->ApplyCommand("/control/execute interactive_mode.mac");
        ui->SessionStart();
        delete ui;
    }

    // Job termination
    delete visManager;
    delete runManager;

    return 0;  // Return 0 to indicate success
}
