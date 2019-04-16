#ifndef EvtInfo_h
#define EvtInfo_h 1

#include "TObject.h"
#include "TVector3.h"
#include "globals.hh"


class EvtInfo : public TObject
{
    public:
        EvtInfo();
        ~EvtInfo();
        void FillEvtInfo(G4double id, G4double x, G4double y, G4double z, G4double t);


    private:
        G4int fEvtId;
        G4double fX;     
        G4double fY;     
        G4double fZ;
        G4double fTime;
            

        ClassDef(EvtInfo,1)
};


#endif
