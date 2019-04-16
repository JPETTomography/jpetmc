#include "PrimaryGeneratorAction.h"
#include "DetectorConstruction.h" 
#include "G4PrimaryVertex.hh"



PrimaryGeneratorAction::PrimaryGeneratorAction(HistoManager* histo) 
    :G4VUserPrimaryGeneratorAction(),
    fPrimaryGenerator(0), fHisto(histo)
{
    fPrimaryGenerator = new PrimaryGenerator();
    fBeam = new BeamParams();
    fIsotope = new SourceParams();
    fMessenger = new PrimaryGeneratorActionMessenger(this);

}

PrimaryGeneratorAction::~PrimaryGeneratorAction() 
{
    delete fIsotope;
    delete fBeam;
    delete fPrimaryGenerator;
    delete fMessenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) 
{

    // if setup for dedicated run is set then ignore its modifications made by user 
    G4int nRun = DetectorConstruction::GetInstance()->GetRunNumber();
    if(  nRun != 0){
        if( GetSourceTypeInfo() != "run"){
          //G4Exception("PrimaryGeneratorAction","PG03",JustWarning,
          //       "User can not modify the predefined run geometry");
            SetSourceTypeInfo("run");
        }
    }


    if( GetSourceTypeInfo() == ("run")) {
        if (nRun == 3) {
            fPrimaryGenerator->GenerateEvtChamberWithSodiumAndPorousMaterial(event,2.5,2.5,7.6);
        } else if (nRun == 5) {
            fPrimaryGenerator->GenerateEvtChamberWithSodiumAndPorousMaterial(event,10.0,10.0,20.0);
        } else if (nRun == 6) {
            fPrimaryGenerator->GenerateEvtChamberWithSodiumAndPorousMaterial(event,15.0,15.0,31.0);
        } else {
          G4Exception("PrimaryGeneratorAction","PG05",FatalException,
                 "Called run with non-exisitng geometry");
        }
    } else if (GetSourceTypeInfo() == ("beam")) {
        fPrimaryGenerator->GenerateBeam(fBeam,event);
    } else if (GetSourceTypeInfo() == ("isotope")) {
        fPrimaryGenerator->GenerateIsotope(fIsotope,event);
    } else if (GetSourceTypeInfo() == ("nema")) {
        fPrimaryGenerator->GenerateNema(GetNemaPoint(),event);
    } else {
          G4Exception("PrimaryGeneratorAction","PG04",FatalException,
                 "Required source type is not allowed");
    }


       // fPrimaryGenerator->GeneratePrimaryVertex(event);


}

void PrimaryGeneratorAction::SetSourceTypeInfo(G4String newSourceType)
{
    if (std::find(std::begin(fAllowedSourceTypes), std::end(fAllowedSourceTypes), newSourceType) != std::end(fAllowedSourceTypes))
    {
        // setup found
        G4int nRun = DetectorConstruction::GetInstance()->GetRunNumber();
        if( (nRun == 0) && (newSourceType != "run")){ 
            fGenerateSourceType = newSourceType;
        } else if (nRun >0) {
            fGenerateSourceType = "run";
        } else {
            G4Exception("PrimaryGeneratorAction","PG01",JustWarning,
                    "Chosen detector geometry corresponds to run number and it can not be changed");
        }
    } else {
        G4Exception("PrimaryGeneratorAction","PG02",JustWarning,
                "Please pick from avaliable setups: beam/isotope");
    }

}



