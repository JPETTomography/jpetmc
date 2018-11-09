#include "SteppingAction.h"
#include "G4TransportationManager.hh" 

SteppingAction::SteppingAction()
{
     G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->SetPushVerbosity(0);
}


SteppingAction::~SteppingAction()
{}
