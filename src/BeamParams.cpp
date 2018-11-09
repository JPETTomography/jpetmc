#include "BeamParams.h"

BeamParams::BeamParams()
{
    energy=511*keV;
    this->SetVtxPosition(0,0,0);
    this->SetPolarization(0,0,0);
    this->SetMomentum(1,0,0);
}

void BeamParams::SetVtxPosition(G4double x, G4double y, G4double z)
{
    vtxPosition.setX(x);
    vtxPosition.setY(y);
    vtxPosition.setZ(z);
}

void BeamParams::SetPolarization(G4double x, G4double y, G4double z)
{
    polarization.setX(x);
    polarization.setY(y);
    polarization.setZ(z);
    polarization = polarization.unit();
}

void BeamParams::SetMomentum(G4double x, G4double y, G4double z)
{
    momentum.setX(x);
    momentum.setY(y);
    momentum.setZ(z);
    momentum = momentum.unit();
}

void BeamParams::SetVtxPosition(G4ThreeVector v)
{
    vtxPosition =v;
}

void BeamParams::SetPolarization(G4ThreeVector v)
{
    polarization=v;
    polarization = polarization.unit();
}

void BeamParams::SetMomentum(G4ThreeVector v)
{
    momentum=v;
    momentum = momentum.unit();
}


