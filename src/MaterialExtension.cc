#include "MaterialExtension.hh"


//MaterialExtension::MaterialExtension(G4Material* mat,const G4String& name):
//    G4VMaterialExtension(name),
//    fMaterial(mat)
MaterialExtension::MaterialExtension(const G4String& name, const G4Material* baseMaterial)
: G4Material(name,baseMaterial->GetDensity(),baseMaterial,    
   baseMaterial->GetState(),baseMaterial->GetTemperature(),
   baseMaterial->GetPressure())
{
    foPsPobability = 0.;
    foPslifetime = 0.;
    fTarget = false ;
}



MaterialExtension::~MaterialExtension()
{;}


G4double MaterialExtension::Get3gFraction()
{
    /// f3g = (1-Pi)/372 + 3/4 tau_i/tau_vac \times Pi
    //  Pi - propability of Ps formation
    //  f^{oPs}_{3g} = tau_i/tau_vac 
    
    // foPsProbability = 4/3 I_3 measured in PALS
    // 3g = direct + oPs
    f3gFraction = (1.-foPsPobability)/372 + (3/4)*(foPslifetime/fTauoPsVaccum)*foPsPobability; 

    return f3gFraction;
}

void MaterialExtension::Set3gProbability(G4double x)
{
    if ( x <0 || x > 1)
    {
        G4Exception ("MaterialExtension", "ME01", FatalException,
                "3gamma fraction can not extend 1");

    } else {
        foPsPobability = x;
    };
}


void MaterialExtension::SetoPsLifetime(G4double x)
{
    if ( x <0 || x > fTauoPsVaccum )
    {
        G4Exception ("MaterialExtension", "ME02", FatalException,
                "oPs lifetime can not excided 142*ns");

    } else {
        foPslifetime = x;
    };

}


