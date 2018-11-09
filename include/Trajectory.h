#ifndef Trajectory_h
#define Trajectory_h 1 

#include "G4VTrajectory.hh"
#include "globals.hh" 
#include "G4Track.hh" 
#include "G4ParticleDefinition.hh"
#include "G4TrajectoryPoint.hh"
#include "G4Allocator.hh" 


// from geant10.4 manual:
// When a customized version of G4Trajectory declares any new class variables, operator new and operator delete must be provided. 


typedef std::vector<G4VTrajectoryPoint*> TrajectoryPointContainer;

class Trajectory : public G4VTrajectory 
{
    public:
        Trajectory();
        Trajectory(const G4Track* aTrack); 
        virtual ~Trajectory();

        inline void* operator new(size_t);
        inline void  operator delete(void*);  
        inline int operator == (const Trajectory& right) const  {return (this==&right);} 

        virtual void MergeTrajectory(G4VTrajectory* secondTrajectory);
        virtual void AppendStep(const G4Step* aStep);
        virtual int GetPointEntries() const { return fPositionRecord->size();}
        virtual G4VTrajectoryPoint* GetPoint(G4int i) const { return (*fPositionRecord)[i]; }
//        virtual void DrawTrajectory() const;
        
        virtual G4int GetTrackID() const { return fTrackID; } 
        virtual G4int GetParentID() const { return fParentID; }
        virtual G4String GetParticleName() const { return fParticleName; } 
        virtual G4double GetCharge() const { return fPDGCharge; } 
        virtual G4int GetPDGEncoding() const { return fPDGEncoding; } 
        virtual G4ThreeVector GetInitialMomentum() const { return fMomentum;}
        virtual G4ThreeVector GetVertexPosition() const { return fVertexPosition;}
//
        G4double GetTime() const { return fGlobalTime; } 
//
    private:
     TrajectoryPointContainer* fPositionRecord;
     G4int                        fTrackID;                           
     G4int                        fParentID;                     
     G4int                        fTrackStatus;                   
     G4ParticleDefinition*        fParticleDefinition;              
     G4String                     fParticleName;                    
     G4double                     fPDGCharge;                        
     G4int                        fPDGEncoding;                      
     G4ThreeVector                fMomentum;                          
     G4ThreeVector                fVertexPosition;                    
     G4double                     fGlobalTime; 
};


extern G4ThreadLocal G4Allocator<Trajectory> * myTrajectoryAllocator;

inline void* Trajectory::operator new(size_t) 
{
    if(!myTrajectoryAllocator) 
        myTrajectoryAllocator = new G4Allocator<Trajectory>;  
    return (void*)myTrajectoryAllocator->MallocSingle();          
}

inline void Trajectory::operator delete(void* aTrajectory)        
{
     myTrajectoryAllocator->FreeSingle((Trajectory*)aTrajectory);
}

#endif
