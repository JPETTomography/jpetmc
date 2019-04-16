#ifndef EventAction_h
#define EventAction_h

#include "G4UserEventAction.hh"
#include "globals.hh" 
#include "HistoManager.h"
#include "JPetGeantDecayTree.h"
#include "JPetGeantScinHits.h"
#include <memory>

/**
 * \class EventAction
 * \brief class managing events 
 *
 * Here also one can reach for information
 * after propagation of all tracks in detector.
 */
class EventAction : public G4UserEventAction
{
    public:
        EventAction(HistoManager* histo); 
        virtual ~EventAction();
        virtual void BeginOfEventAction(const G4Event*);
        virtual void EndOfEventAction(const G4Event* anEvent);

    private:
        //std::unique_ptr<SavedEvent> fEvent {new SavedEvent()};
        HistoManager* fHisto;
        G4int fScinCollID;
};

#endif
