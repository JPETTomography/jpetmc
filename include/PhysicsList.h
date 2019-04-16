#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh" 

/**
 * \class PhysicsList
 * \brief standard GEANT4 package is used for physics 
 * used G4EmLivermorePolarizedPhysics deals properly with 
 * polarized particles
 */
class PhysicsList: public G4VModularPhysicsList
{
    public:
    PhysicsList();
    virtual ~PhysicsList();
};

#endif
