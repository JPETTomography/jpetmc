#include "PrimaryGeneratorActionMessenger.hh"
#include "PrimaryGeneratorAction.hh"

PrimaryGeneratorActionMessenger::PrimaryGeneratorActionMessenger(PrimaryGeneratorAction* primGeneratorAction)
    :fPrimGen(primGeneratorAction)
{

    fDirectory = new G4UIdirectory("/jpetmc/source/"); 
    fDirectory->SetGuidance("Commands for controling the gamma quanta source (beam/target) and its parameters");

    fSourceType = new G4UIcmdWithAString("/jpetmc/source/setType",this);
    fSourceType->SetCandidates("beam isotope");
    fSourceType->SetDefaultValue("beam"); 

    fGammaBeamSetEnergy = new G4UIcmdWithADoubleAndUnit("/jpetmc/source/gammaBeam/setEnergy",this);
    fGammaBeamSetEnergy->SetGuidance("Set energy (value and unit) for beam of gamma quanta");
    fGammaBeamSetEnergy->SetDefaultValue(511);
    fGammaBeamSetEnergy->SetDefaultUnit("keV");
    fGammaBeamSetEnergy->SetUnitCandidates("keV");


    fGammaBeamSetPosition = new G4UIcmdWith3VectorAndUnit("/jpetmc/source/gammaBeam/setPosition",this);
    fGammaBeamSetPosition->SetGuidance("Set vertex position of the gamma quanta beam");
    fGammaBeamSetPosition->SetDefaultValue(G4ThreeVector(0,0,0));
    fGammaBeamSetPosition->SetDefaultUnit("cm");
    fGammaBeamSetPosition->SetUnitCandidates("cm");
    fGammaBeamSetPosition->SetParameterName("Xvalue","Yvalue","Zvalue",false);

    fGammaBeamSetMomentum = new G4UIcmdWith3VectorAndUnit("/jpetmc/source/gammaBeam/setMomentum",this);
    fGammaBeamSetMomentum->SetGuidance("Set momentum of the gamma quanta beam");
    fGammaBeamSetMomentum->SetDefaultValue(G4ThreeVector(1,0,0));
    fGammaBeamSetMomentum->SetDefaultUnit("keV");
    fGammaBeamSetMomentum->SetUnitCandidates("keV");
    fGammaBeamSetMomentum->SetParameterName("Xvalue","Yvalue","Zvalue",false);

    fGammaBeamSetPolarization = new G4UIcmdWith3Vector("/jpetmc/source/gammaBeam/setPolarization",this);
    fGammaBeamSetPolarization->SetGuidance("Set initial polarization of the gamma quanta beam");
    fGammaBeamSetPolarization->SetDefaultValue(G4ThreeVector(0,0,0));
    fGammaBeamSetPolarization->SetParameterName("Xvalue","Yvalue","Zvalue",false);


    fIsotopeSetShape = new G4UIcmdWithAString("/jpetmc/source/isotope/setShape",this);
    fIsotopeSetShape->SetCandidates("cylinder");

    fIsotopeSetGenGammas = new G4UIcmdWithAnInteger("/jpetmc/source/isotope/setNGamma",this);
    fIsotopeSetGenGammas->SetGuidance("Give number of gamma quanta to generate 1 / 2 / 3 ");

    fIsotopeSetShapeDimCylinderRadius = new G4UIcmdWithADoubleAndUnit("/jpetmc/source/isotope/setShape/cylinderRadius",this);
    fIsotopeSetShapeDimCylinderRadius->SetGuidance("For cylindrical shape - set radius"); 
    fIsotopeSetShapeDimCylinderRadius->SetDefaultValue(10);
    fIsotopeSetShapeDimCylinderRadius->SetDefaultUnit("cm");
    fIsotopeSetShapeDimCylinderRadius->SetUnitCandidates("cm"); 

    fIsotopeSetShapeDimCylinderZ = new G4UIcmdWithADoubleAndUnit("/jpetmc/source/isotope/setShape/cylinderZ",this);
    fIsotopeSetShapeDimCylinderZ->SetGuidance("For cylindrical shape - set z (half length)"); 
    fIsotopeSetShapeDimCylinderZ->SetDefaultValue(10);
    fIsotopeSetShapeDimCylinderZ->SetDefaultUnit("cm");
    fIsotopeSetShapeDimCylinderZ->SetUnitCandidates("cm"); 


}


PrimaryGeneratorActionMessenger::~PrimaryGeneratorActionMessenger()
{
    delete fIsotopeSetShape;
    delete fIsotopeSetShapeDimCylinderRadius;
    delete fIsotopeSetShapeDimCylinderZ;
    delete fSourceType;
    delete fGammaBeamSetEnergy;
    delete fGammaBeamSetPosition;
    delete fGammaBeamSetMomentum;
    delete fGammaBeamSetPolarization;
}


void PrimaryGeneratorActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if(command==fSourceType){
        fPrimGen->SetSourceTypeInfo(newValue);
    }

    if(command==fGammaBeamSetEnergy){
        CheckIfBeam();
        fPrimGen->GetBeamParams()->SetEnergy(fGammaBeamSetEnergy->GetNewDoubleRawValue(newValue));
    }

    if(command==fGammaBeamSetPosition){
        CheckIfBeam();
        fPrimGen->GetBeamParams()->SetVtxPosition(fGammaBeamSetPosition->GetNew3VectorValue(newValue));
    }

    if(command==fGammaBeamSetMomentum){
        CheckIfBeam();
        fPrimGen->GetBeamParams()->SetMomentum(fGammaBeamSetMomentum->GetNew3VectorValue(newValue));
    }


    if(command==fIsotopeSetShape){
        fPrimGen->GetIsotopeParams()->SetShape(newValue);
    }

    if(command==fIsotopeSetGenGammas){
        fPrimGen->GetIsotopeParams()->SetGammasNumber(fIsotopeSetGenGammas->GetNewIntValue(newValue));
    }


    if(command==fIsotopeSetShapeDimCylinderRadius){
        CheckIfIsotope();
        fPrimGen->GetIsotopeParams()->SetShapeDim(0,fIsotopeSetShapeDimCylinderRadius->GetNewDoubleRawValue(newValue));
    }


    if(command==fIsotopeSetShapeDimCylinderZ){
        CheckIfIsotope();
        fPrimGen->GetIsotopeParams()->SetShapeDim(1,fIsotopeSetShapeDimCylinderRadius->GetNewDoubleRawValue(newValue));
    }

}


void PrimaryGeneratorActionMessenger::CheckIfIsotope()
{
        if(fPrimGen->GetSourceTypeInfo() != "isotope")
        {
            G4Exception("PrimaryGeneratorActionMessenger","PGM02",JustWarning,
                    "Changed sourceType to isotope");
            fPrimGen->SetSourceTypeInfo("isotope");
        }

}

void PrimaryGeneratorActionMessenger::CheckIfBeam()
{
        // check if we really changing parameters corresponding to the beam
        if(fPrimGen->GetSourceTypeInfo() != "beam")
        {
            G4Exception("PrimaryGeneratorActionMessenger","PGM01",JustWarning,
                    "Changed sourceType to beam");
            fPrimGen->SetSourceTypeInfo("beam");
        }
}
