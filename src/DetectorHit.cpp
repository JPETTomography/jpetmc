#include "DetectorHit.h"

DetectorHit::DetectorHit()
    : G4VHit(), fScinID(0), fTrackID(-1), fTrackPDG(0), fEdep(0.0), fTime(0), fPos(0), fNumInteractions(0), fPolarizationIn(0,0,0), fPolarizationOut(0,0,0), 
    fMomentumIn(0,0,0), fMomentumOut(0,0,0),
   fGenGammaMultiplicity(0), fGenGammaIndex(0) 
{}

DetectorHit::~DetectorHit()
{}


