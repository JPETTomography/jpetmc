#include "VtxInformation.h"

VtxInformation::VtxInformation()
    :fVtxPosition(0), fTwoGammaGen(false), fThreeGammaGen(false), fPromptGammaGen(false), fnRun(0), fLifetime(0)
{}


VtxInformation::~VtxInformation()
{}


void VtxInformation::SetVtxPosition(G4double x, G4double y, G4double z)
{ 
    fVtxPosition.setX(x);
    fVtxPosition.setY(y);
    fVtxPosition.setZ(z);
}


void VtxInformation::Clear()
{
    fTwoGammaGen = false;
    fThreeGammaGen = false;
    fPromptGammaGen = false;
    fnRun = 0;
    fLifetime =0.;
    fVtxPosition.setX(0);
    fVtxPosition.setY(0);
    fVtxPosition.setZ(0);

}

void VtxInformation::Print() const
{
    
}

