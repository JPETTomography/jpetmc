#ifndef ActionInitialization_h
#define ActionInitialization_h 1
#include "G4VUserActionInitialization.hh"


/**
 * \class ActionInitialization
 * \brief function inherited from GEANT4; be careful while implementing multithread mode
 */
class ActionInitialization : public G4VUserActionInitialization 
{
    public:
        ActionInitialization(); 
        virtual ~ActionInitialization();
        virtual void BuildForMaster() const; //< put here functions called only once during multithread mode as histograms
        virtual void Build() const; //< functions called for each thread
};


#endif
