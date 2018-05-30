#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"


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

