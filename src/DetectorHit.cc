#include "DetectorHit.hh"

DetectorHit::DetectorHit()
    : G4VHit(), fScinID(0), fTrackID(-1), fEdep(0.0), fTime(0), fPos(0),
   fGenGammaMultiplicity(0), fGenGammaIndex(0) 
{}

DetectorHit::~DetectorHit()
{}


