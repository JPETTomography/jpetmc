#include "PrimaryParticleInformation.hh"

PrimaryParticleInformation::PrimaryParticleInformation()
    :fIndex(0), fDecayMultiplicity(0)
{}


PrimaryParticleInformation::~PrimaryParticleInformation()
{}



void PrimaryParticleInformation::Clear()
{
    fIndex = 0;
    fDecayMultiplicity = 0;
}

void PrimaryParticleInformation::Print() const
{
    
}

