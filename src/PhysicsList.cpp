#include "PhysicsList.h"
#include "G4EmLivermorePolarizedPhysics.hh"
//#include "G4EmStandardPhysics_option4.hh"


PhysicsList::PhysicsList() 
: G4VModularPhysicsList()
{
     SetVerboseLevel(1);
     RegisterPhysics(new G4EmLivermorePolarizedPhysics());
     //RegisterPhysics(new G4EmStandardPhysics_option4());
}

PhysicsList::~PhysicsList()
{
}

