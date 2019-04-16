#ifndef DetectorSD_h
#define DetectorSD_h 1

#include "G4VSensitiveDetector.hh"
#include "DetectorHit.h"

class DetectorSD : public G4VSensitiveDetector
{
    public:
        DetectorSD(G4String name, G4int scinSum);
        virtual ~DetectorSD();

        virtual void Initialize(G4HCofThisEvent* HCE);

    private:
        const G4double timeIntervals = 0.08; // 80 ps
        G4int totScinNum; 
        std::vector<G4int> previousHitHistory;
        std::vector<G4double> previousHitTimeHistory;

        DetectorHitsCollection *fDetectorCollection;

    protected:
         virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);


};


#endif
