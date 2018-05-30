#include "Trajectory.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4Polyline.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"
#include "G4DynamicParticle.hh"
#include "G4PrimaryParticle.hh"
#include "G4SystemOfUnits.hh"


G4ThreadLocal G4Allocator<Trajectory> * myTrajectoryAllocator = 0;

Trajectory::Trajectory()
:G4VTrajectory(),  fPositionRecord(0) 
{
}

Trajectory::Trajectory(const G4Track* aTrack) 
:G4VTrajectory(), fPositionRecord(0)
{
    fParticleDefinition = aTrack->GetDefinition();
    fParticleName = fParticleDefinition->GetParticleName();
    fPDGCharge = fParticleDefinition->GetPDGCharge();
    fPDGEncoding = fParticleDefinition->GetPDGEncoding();

    fTrackID = aTrack->GetTrackID(); 


    fPositionRecord =new TrajectoryPointContainer();
    fPositionRecord->push_back(new G4TrajectoryPoint(aTrack->GetPosition()));
    fMomentum = aTrack->GetMomentumDirection();
    fVertexPosition = aTrack->GetPosition(); 
    fGlobalTime = aTrack->GetGlobalTime();

    fParentID = aTrack->GetParentID();

}

Trajectory::~Trajectory() 
{
    size_t i; 
     for(i=0;i<fPositionRecord->size();i++){
         delete  (*fPositionRecord)[i];

     }
     fPositionRecord->clear();  
     delete fPositionRecord;
}


void Trajectory::MergeTrajectory(G4VTrajectory* secondTrajectory)  
{ 
//  if(!secondTrajectory) return; 
//  
//  Trajectory* seco = (Trajectory*)secondTrajectory;            
//  G4int ent = seco->GetPointEntries();                                 
//  //
//  // initial point of the second trajectory should not be merged       
//  for(int i=1;i<ent;i++) 
//  {                                                                    
//    fPositionRecord->push_back((*(seco->fPositionRecord))[i]);         
//  }
//  delete (*seco->fPositionRecord)[0];
//  seco->fPositionRecord->clear();
//  
} 


void Trajectory::AppendStep(const G4Step* aStep)
{
    fPositionRecord->push_back( new G4TrajectoryPoint(aStep->GetPostStepPoint()->GetPosition() ));
}
    
//void Trajectory::DrawTrajectory() const
//{
//   G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
//   G4ThreeVector pos;
//
//   G4Polyline pPolyline;
//   for (size_t i = 0; i < fPositionRecord->size() ; i++) {
//     G4TrajectoryPoint* aTrajectoryPoint = 
//       (G4TrajectoryPoint*)((*fPositionRecord)[i]);
//     pos = aTrajectoryPoint->GetPosition();
//     pPolyline.push_back( pos );
//   }
//
//   G4Colour colour(0.2,0.2,0.2);
//   if(fParticleDefinition==G4Gamma::GammaDefinition())
//      colour = G4Colour(0.,0.,1.);
//   else if(fParticleDefinition==G4Electron::ElectronDefinition()
//         ||fParticleDefinition==G4Positron::PositronDefinition())
//      colour = G4Colour(1.,1.,0.);
//   else if(fParticleDefinition==G4MuonMinus::MuonMinusDefinition()
//         ||fParticleDefinition==G4MuonPlus::MuonPlusDefinition())
//      colour = G4Colour(0.,1.,0.);
//   else if(fParticleDefinition->GetParticleType()=="meson")
//   {
//      if(fPDGCharge!=0.)
//         colour = G4Colour(1.,0.,0.);
//      else
//         colour = G4Colour(0.5,0.,0.);
//   }
//   else if(fParticleDefinition->GetParticleType()=="baryon")
//   {
//      if(fPDGCharge!=0.)
//         colour = G4Colour(0.,1.,1.);
//      else 
//         colour = G4Colour(0.,0.5,0.5);
//   }
//
//   G4VisAttributes attribs(colour);
//   pPolyline.SetVisAttributes(attribs);
//   if(pVVisManager) pVVisManager->Draw(pPolyline);
//
//
//}
