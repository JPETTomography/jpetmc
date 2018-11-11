#ifndef PrimaryGenerator_h 
#define PrimaryGenerator_h 1

#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4VPrimaryGenerator.hh"
#include "MaterialExtension.h"
#include<TGenPhaseSpace.h> 
#include<TLorentzVector.h>
#include "BeamParams.h"
#include "SourceParams.h"


class PrimaryGenerator : public G4VPrimaryGenerator 
{
    public:
        PrimaryGenerator();
        ~PrimaryGenerator();
        void GenerateBeam(BeamParams*, G4Event* );
        void GenerateIsotope(SourceParams*, G4Event* );
        void GenerateNema(G4int, G4Event* );
        void GenerateEvtChamberWithSodiumAndPorousMaterial(G4Event* event, G4double maxXhalf, G4double maxYhalf, G4double maxZhalf); ///< as input the maximal dimension(/2) of annihilation chamber are passed (to speed up)

    public:
         virtual void GeneratePrimaryVertex(G4Event*);

     private:

         std::tuple<G4ThreeVector,G4double,G4double> GetVerticesDistribution(G4double,G4double,G4double);  ///< return: vtx position, 2/3g ratio, meanlifetime; as input the maximal dimension(/2) of annihilation chamber are taken (to speed up simulatons)
        void GenerateTwoGammaVertex(G4PrimaryVertex*);
        void GenerateThreeGammaVertex(G4PrimaryVertex*);
        void GeneratePromptGammaSodium(G4PrimaryVertex*);

        G4ThreeVector VertexUniformInCylinder(G4double, G4double);

        Double_t calculate_mQED(Double_t mass_e, Double_t w1, Double_t w2, Double_t w3);
       
};


#endif 
