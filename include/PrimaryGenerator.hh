#ifndef PrimaryGenerator_h 
#define PrimaryGenerator_h 1

#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4VPrimaryGenerator.hh"
#include "MaterialExtension.hh"
#include<TGenPhaseSpace.h> 
#include<TLorentzVector.h>
#include "BeamParams.hh"
#include "SourceParams.hh"


class PrimaryGenerator : public G4VPrimaryGenerator 
{
    public:
        PrimaryGenerator();
        ~PrimaryGenerator();
        void GenerateBeam(BeamParams*, G4Event* );
        void GenerateIsotope(SourceParams*, G4Event* );
        void GenerateEvtChamberRun3(G4Event*);

    public:
         virtual void GeneratePrimaryVertex(G4Event*);

     private:

         // return: vtx position, 2/3g ratio, meanlifetime
         std::tuple<G4ThreeVector,G4double,G4double> GetVerticesDistribution(); 
        void GenerateTwoGammaVertex(G4PrimaryVertex*);
        void GenerateThreeGammaVertex(G4PrimaryVertex*);
        void GeneratePromptGammaSodium(G4PrimaryVertex*);

        G4ThreeVector VertexUniformInCylinder(G4double, G4double);

        Double_t calculate_mQED(Double_t mass_e, Double_t w1, Double_t w2, Double_t w3);
       
};


#endif 
