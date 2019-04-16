#ifndef RunAction_h 
#define RunAction_h 1
#include "G4UserRunAction.hh"
#include "globals.hh"

#include "HistoManager.h"


class G4Run;

/**
 * \class RunAction
 * \brief class manages the runs = set of events
 */
class RunAction : public G4UserRunAction 
{
     public:
         RunAction(HistoManager*); 
         virtual ~RunAction();

         virtual void BeginOfRunAction(const G4Run*);
         virtual void   EndOfRunAction(const G4Run*);

     private:
         HistoManager* fHistoManager;


};

#endif
