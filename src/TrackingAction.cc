#include "G4Track.hh"
#include "G4TrackingManager.hh"
#include "TrackingAction.hh"
#include "Trajectory.hh"

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{

// Create trajectory only for primaries 
    fpTrackingManager->SetStoreTrajectory(true); 
    fpTrackingManager->SetTrajectory(new Trajectory(aTrack));
//    if(aTrack->GetParentID()==0)
//    { fpTrackingManager->SetStoreTrajectory(true); }
//    else
//    { fpTrackingManager->SetStoreTrajectory(false); }
 }


void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{

}
