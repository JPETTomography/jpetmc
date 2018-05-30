#ifndef MaterialExtension_h
#define MaterialExtension_h


#include "G4VMaterialExtension.hh"
#include "G4Material.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

/**
 * \class MaterialExtension
 * \brief extended material keeps information about 2g/3g creation fraction
 */

const G4double fTauoPsVaccum = 142*ns;

const G4double fTauBulk = 0.2*ns;
// data for AL  https://doi.org/10.1007/s00339-015-9214-0
const G4double fTauoPsAl = 0.32*ns;
const G4double foPsProbabilityAl = 0.12;

//class MaterialExtension : public G4VMaterialExtension
class MaterialExtension : public G4Material
{
    public:
      MaterialExtension(const G4String& name, const G4Material* baseMaterial);

      ~MaterialExtension();
      G4Material* GetMaterial() {return fMaterial;};

      void SetoPsLifetime(G4double);
      void Set3gProbability(G4double);

      G4double Get3gFraction(); 
      G4double GetoPsLifetime(){return foPslifetime;};


      G4bool IsTarget() {return fTarget;};
      void AllowsAnnihilations(G4bool tf){ fTarget = tf;};

      void Print() const {;};
      G4bool IsExtended() const { return true; }

    private:
      G4Material* fMaterial;
      G4bool fTarget;

      G4double f3gFraction; ///<  3g/2g events fraction 
      G4double foPsPobability; 
      G4double foPslifetime;

};

#endif
