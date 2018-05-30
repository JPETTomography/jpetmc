#include "EventAction.hh"
#include "G4PrimaryParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Event.hh"

#include "Trajectory.hh"
#include "G4TrajectoryContainer.hh"

#include "G4EventManager.hh"
#include "G4SDManager.hh"
#include "DetectorHit.hh"


EventAction::EventAction(HistoManager* histo )
     :G4UserEventAction(), 
     fScinCollID(-1)
{

    fHisto = histo;
}

EventAction::~EventAction()
{}


void EventAction::BeginOfEventAction(const G4Event*)
{

    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    if(fScinCollID<0)
    {
        G4String colNam;
        fScinCollID = SDman->GetCollectionID(colNam="detectorCollection");
    }

     fHisto->Clear();
}


void EventAction::EndOfEventAction(const G4Event* anEvent)
{


     if(anEvent->GetNumberOfPrimaryVertex()==0) return;

     G4int id =  anEvent->GetEventID();
     fHisto->SetEventNumber(id);


     /*
     G4TrajectoryContainer* trajectoryContainer = anEvent->GetTrajectoryContainer();
     G4int trackNum = 0; 
     if (trajectoryContainer) trackNum = trajectoryContainer->entries();

     for (G4int i=0; i<trackNum; i++)
     {
         Trajectory* vec = (Trajectory*)((*( anEvent->GetTrajectoryContainer()))[i]); 
         printf("in %i  track ID %i, pdg %i parentID %i, x y x %4.2f %4.2f %4.2f \n",
                 i,
                 vec->GetTrackID(),
                 vec->GetPDGEncoding(),
                 vec->GetParentID(),
                 vec->GetVertexPosition().x(),
                 vec->GetVertexPosition().y(),
                 vec->GetVertexPosition().z()
               );
     }
     */

    G4HCofThisEvent * HCE = anEvent->GetHCofThisEvent();
    DetectorHitsCollection* DHC = 0;
    if(HCE)
    {
        DHC = (DetectorHitsCollection*)(HCE->GetHC(fScinCollID));
        int n_hit = DHC->entries();
        //G4cout << "Detector hits "<< n_hit <<    
        //"-------------------------------------------------" << G4endl;
        for (int i=0; i<n_hit; i++)
        {
          fHisto->AddNewHit((DetectorHit*)DHC->GetHit(i) );
           //DetectorHit* dh =  (DetectorHit*)DHC->GetHit(i);
           //printf("dh %i trID %i  gammMuliplicity %i gammaindex %i  \n", 
           //        i,
           //        dh->GetTrackID(),
           //        dh->GetGenGammaMultiplicity(),
           //        dh->GetGenGammaIndex()
           //        );
        }

    }


    for( int i=0; i<anEvent->GetNumberOfPrimaryVertex(); i++)
    {
        VtxInformation* info =  (VtxInformation*) anEvent->GetPrimaryVertex(i)->GetUserInformation();    
        if( info != 0 )
        {
            fHisto->AddGenInfo(info);
        }
    }

    // save full information about event in final ntuples
    fHisto->SaveEvtPack();

     if ( id % 1000 == 0)
     {
         printf (" === Processed %i events \n", id);
     }

}
