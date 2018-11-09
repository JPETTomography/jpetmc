#include "DetectorConstruction.h"
#include "PhysicsList.h"
#include "G4RunManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "ActionInitialization.h"
#include "Randomize.hh"
#include "time.h"
#include "G4INCLRandom.hh"

int main (int argc,char** argv)
{

    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
    G4long seed = time(NULL);
    CLHEP::HepRandom::setTheSeed(seed);
    long seeds[2];
    time_t systime = time(NULL);
    seeds[0] = (long) systime;
    seeds[1] = (long) (systime * G4UniformRand());
    G4Random::setTheSeeds(seeds);



    G4UIExecutive* ui = 0;
    if ( argc == 1 ) {
        ui = new G4UIExecutive(argc, argv);
    }


    G4RunManager* runManager = new G4RunManager;

    runManager->SetUserInitialization(DetectorConstruction::GetInstance());
    runManager->SetUserInitialization(new PhysicsList);
    runManager->SetUserInitialization(new ActionInitialization);

    //runManager->Initialize();
    //int numberOfEvent = 10;
    //runManager->BeamOn(numberOfEvent);


    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();


    if ( ! ui ) {
        // batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command+fileName);
    }
    else {
        // interactive mode
        UImanager->ApplyCommand("/control/execute init_vis.mac");
        ui->SessionStart();
        delete ui;
    }


    delete visManager;
    delete runManager;
    return 0;
}
