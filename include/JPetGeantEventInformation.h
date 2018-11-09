#ifndef JPET_GEANT_EVENT_INFORMATION_H
#define JPET_GEANT_EVENT_INFORMATION_H 1

#include "TObject.h"
#include <vector>
#include "TVector3.h" 
#include "TBits.h"

/**
 * \class JPetGeantEventInformation
 * \brief keeps information about initial simulation parameters
 * e.g. vertices and times distributions for annihilation 
 * and prompt gamma photons
 */ 
class JPetGeantEventInformation : public TObject {
    public:
        JPetGeantEventInformation();
        ~JPetGeantEventInformation();
        void Clear();

        void SetThreeGammaGen(bool tf){fgenGammaNum.SetBitNumber(2,tf);};
        void SetTwoGammaGen(bool tf){fgenGammaNum.SetBitNumber(1,tf);};
        void SetPromptGammaGen(bool tf){fgenGammaNum.SetBitNumber(0,tf);};
        void SetRunNr(int x){fnRun =x;};
        void SetVtxPosition(double x, double y, double z){fVtxPosition.SetXYZ(x,y,z);};
        void SetVtxPromptPosition(double x, double y, double z){fVtxPromptPosition.SetXYZ(x,y,z);};
        void SetLifetime(double x){fLifetime=x;};
        void SetPromptLifetime(double x){fPromptLifetime=x;};

        bool GetThreeGammaGen(){return fgenGammaNum.TestBitNumber(2);};
        bool GetTwoGammaGen(){return fgenGammaNum.TestBitNumber(1);};
        bool GetPromptGammaGen(){return fgenGammaNum.TestBitNumber(0);};
        int GetRunNr(){return fnRun;};
        double GetVtxPositionX(){return fVtxPosition.X();};
        double GetVtxPositionY(){return fVtxPosition.Y();};
        double GetVtxPositionZ(){return fVtxPosition.Z();};
        double GetVtxPromptPositionX(){return fVtxPromptPosition.X();};
        double GetVtxPromptPositionY(){return fVtxPromptPosition.Y();};
        double GetVtxPromptPositionZ(){return fVtxPromptPosition.Z();};

        double GetLifetime(){return fLifetime;};
        double GetPromptLifetime(){return fPromptLifetime;};


    private:
        TVector3 fVtxPosition; ///< xyz annihilation coordinated
        TVector3 fVtxPromptPosition; ///< xyz of prompt photon emmision
        TBits fgenGammaNum; ///< bitNR 0-prompt; 1-back-to-back; 2- oPs 
        int fnRun; ///< number should follow the JPet run numbering scheme
        double fLifetime; ///< lifetime of generated bound state or direct annihilation; see specific simulation details
        double fPromptLifetime; ///< generated lifetime of emmited prompt photon; filled only if prompt gamma is generated

    private:
     ClassDef(JPetGeantEventInformation,3)

};


#endif
