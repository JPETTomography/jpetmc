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

        void SetThreeGammaGen(bool tf){fGenGammaNum.SetBitNumber(2,tf);};
        void SetTwoGammaGen(bool tf){fGenGammaNum.SetBitNumber(1,tf);};
        void SetPromptGammaGen(bool tf){fGenGammaNum.SetBitNumber(0,tf);};
        void SetRunNr(int x){fnRun =x;};
        void SetVtxPosition(double x, double y, double z){fVtxPosition.SetXYZ(x,y,z);};
        void SetVtxPromptPosition(double x, double y, double z){fVtxPromptPosition.SetXYZ(x,y,z);};
        void SetLifetime(double x){fLifetime=x;};
        void SetPromptLifetime(double x){fPromptLifetime=x;};

        bool GetThreeGammaGen() const {return fGenGammaNum.TestBitNumber(2);};
        bool GetTwoGammaGen() const {return fGenGammaNum.TestBitNumber(1);};
        bool GetPromptGammaGen() const {return fGenGammaNum.TestBitNumber(0);};
        int GetRunNr() const {return fnRun;};
        double GetVtxPositionX() const {return fVtxPosition.X();};
        double GetVtxPositionY() const {return fVtxPosition.Y();};
        double GetVtxPositionZ() const {return fVtxPosition.Z();};
        double GetVtxPromptPositionX() const {return fVtxPromptPosition.X();};
        double GetVtxPromptPositionY() const {return fVtxPromptPosition.Y();};
        double GetVtxPromptPositionZ() const {return fVtxPromptPosition.Z();};

        double GetLifetime() const {return fLifetime;};
        double GetPromptLifetime() const {return fPromptLifetime;};


    private:
        const unsigned int fMaxGammaNumberIndex = 3;
        TVector3 fVtxPosition; ///< xyz annihilation coordinated
        TVector3 fVtxPromptPosition; ///< xyz of prompt photon emmision
        TBits fGenGammaNum; ///< bitNR 0-prompt; 1-back-to-back; 2- oPs 
        int fnRun = -1; ///< number should follow the JPet run numbering scheme
        double fLifetime = -1.0; ///< lifetime of generated bound state or direct annihilation; see specific simulation details
        double fPromptLifetime = -1.0; ///< generated lifetime of emmited prompt photon; filled only if prompt gamma is generated

    private:
     ClassDef(JPetGeantEventInformation,3)

};


#endif
