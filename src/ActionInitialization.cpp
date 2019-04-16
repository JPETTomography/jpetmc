#include "ActionInitialization.h"
#include "PrimaryGeneratorAction.h"
#include "RunAction.h"
#include "HistoManager.h"
#include "EventAction.h"
#include "TrackingAction.h"
#include "SteppingAction.h"


ActionInitialization::ActionInitialization()
     : G4VUserActionInitialization()
{}

ActionInitialization::~ActionInitialization()
{}

void ActionInitialization::BuildForMaster() const
{
//    HistoManager*  histo = new HistoManager();
//    SetUserAction(new RunAction(histo));

}

void ActionInitialization::Build() const
{


    HistoManager*  histo = new HistoManager();

    SetUserAction(new EventAction(histo));
    SetUserAction(new RunAction(histo));
    SetUserAction(new PrimaryGeneratorAction(histo));

    SetUserAction(new TrackingAction);
    SetUserAction(new SteppingAction);

}

