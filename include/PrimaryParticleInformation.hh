#ifndef PRIMARY_PARTICLE_INFORMATION_H
#define PRIMARY_PARTICLE_INFORMATION_H

#include "G4VUserPrimaryParticleInformation.hh" 
#include "globals.hh"
#include "G4ThreeVector.hh"

class PrimaryParticleInformation : public G4VUserPrimaryParticleInformation {
    public:
        PrimaryParticleInformation();
        virtual ~PrimaryParticleInformation();
        void Clear();
        virtual void Print() const;

        G4int GetIndex(){return fIndex;}
        G4int GetGammaMultiplicity(){return fDecayMultiplicity;}

        void SetIndex(G4int i){fIndex = i;}
        void SetGammaMultiplicity(G4int i){fDecayMultiplicity =i;}
    private:
        G4int fIndex;
        G4int fDecayMultiplicity;
        // 1 - prompt gamma
        // 2 = 2gamma
        // 3 = 3gamma

};



#endif
