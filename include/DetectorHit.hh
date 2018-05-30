#ifndef DetectorHit_h
#define DetectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Types.hh"


/**
 * \class DetectorHit
 * \brief Class containing hits in the sensitive part of the detector
 *
 * class is filled in DetectorSD class and rewritten into ROOT tree in 
 * HistoManager::FillScin  
 */
class DetectorHit : public G4VHit
{
    public:
    DetectorHit(); //< standard constructor
    virtual ~DetectorHit();
    void SetEdep(G4double de)  { fEdep = de; }
    void SetTime(G4double val) {fTime = val;}
    void SetTrackID(G4int i) { fTrackID =i;}
    void SetTrackPDG(G4int i) {fTrackPDG = i;}
    void SetScinID(G4int i) {fScinID = i;}
    void SetPosition(G4ThreeVector xyz) {fPos =xyz;}
    void SetInteractionNumber() {fNumInteractions = 1;}
    void AddInteraction() { fNumInteractions += 1;}
    void AddEdep(G4double de) { fEdep += de;}

    void SetPolarizationIn(G4ThreeVector xyz) {fPolarizationIn =xyz;}
    void SetPolarizationOut(G4ThreeVector xyz) {fPolarizationOut =xyz;}
    void SetMomentumIn(G4ThreeVector xyz)  {fMomentumIn =xyz;}
    void SetMomentumOut(G4ThreeVector xyz) {fMomentumOut =xyz;}
    void SetProcessName(G4String str) {fName = str;}


    G4int GetScinID() {return fScinID;}
    G4int GetTrackID() {return fTrackID;}
    G4double GetEdep() {return fEdep;}
    G4double GetTime() {return fTime;}
    G4ThreeVector GetPosition() {return fPos;}
    G4int GetTrackPDG() {return fTrackPDG;}
    G4ThreeVector GetPolarizationIn()  {return fPolarizationIn;}
    G4ThreeVector GetPolarizationOut() {return fPolarizationOut;}
    G4ThreeVector GetMomentumIn()  {return fPolarizationIn;}
    G4ThreeVector GetMomentumOut() {return fPolarizationOut;}
    G4int GetNumInteractions() {return fNumInteractions;}
    G4String GetProcessName() {return fName;}

    void SetGenGammaMultiplicity(G4int i){fGenGammaMultiplicity=i;}
    void SetGenGammaIndex(G4int i){fGenGammaIndex=i;}

    G4int GetGenGammaMultiplicity(){return fGenGammaMultiplicity;}
    G4int GetGenGammaIndex(){return fGenGammaIndex;}


    private:
    G4int fScinID; //< scintillator number (arbitrary!; not following convension used in laboratory )
    G4int fTrackID; //< track identificator 
    G4int fTrackPDG; //< Particle Data Group numbering for particles
    G4double fEdep; //< total energy deposited in the strip
    G4double fTime; //< interaction time
    G4ThreeVector fPos; //< 3D interaction position
    G4int fNumInteractions; //< counters for interaction classified as single hit
    G4ThreeVector fPolarizationIn;  //< polarization vector before interaction
    G4ThreeVector fPolarizationOut; //< polarization vector after interaction
    G4ThreeVector fMomentumIn;  //< momentum vector before interaction
    G4ThreeVector fMomentumOut; //< momentum vector after interaction
    G4String fName; //< process name

    // from PrimaryParticleInformation
    G4int fGenGammaMultiplicity;
    G4int fGenGammaIndex;

};

typedef G4THitsCollection<DetectorHit> DetectorHitsCollection;

#endif
