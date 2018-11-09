#ifndef JPETGEANTEVENTPACK_H
#define JPETGEANTEVENTPACK_H

#include <TObject.h>
#include <TVector3.h> 
#include <TClonesArray.h>
#include <TBits.h>
#include "JPetGeantScinHits.h"
#include "JPetGeantDecayTree.h"
#include "JPetGeantEventInformation.h"

/**
 * \class JPetGeantEventPack
 * \brief container that keeps information about single event: 
 *  initial parameters (in JPetGeantEventInformation) 
 *  detector hits (in JPetGeantScinHits - true information, without detector resolution)
 *  and decay tree (in JPetGeantDecayTree)
 *  Class is directly processed in JPetGeantParser
 */
class JPetGeantEventPack : public TObject
{
    public:
        JPetGeantEventPack();
        ~JPetGeantEventPack();
        void Clear();


        JPetGeantScinHits* ConstructNextHit();
        JPetGeantDecayTree* ConstructNextDecayTree();

        JPetGeantScinHits* GetHit(int i){return dynamic_cast<JPetGeantScinHits*>(fMCHits[i]);};
        JPetGeantDecayTree* GetDecayTree(int i){return dynamic_cast<JPetGeantDecayTree*>(fMCDecayTrees[i]);};

        JPetGeantEventInformation* GetEventInformation(){return fGenInfo;};

        unsigned int GetNumberOfHits(){ return fHitIndex;};
        unsigned int GetNumberOfDecayTrees(){ return fMCDecayTreesIndex;};
        unsigned int GetEventNumber(){ return fEvtIndex;};
        void SetEventNumber(int x){fEvtIndex = x;};
          


    private:
        TClonesArray fMCHits;
        TClonesArray fMCDecayTrees;
        JPetGeantEventInformation* fGenInfo;

        unsigned int fEvtIndex;
        unsigned int fHitIndex;
        unsigned int fMCDecayTreesIndex;


        ClassDef(JPetGeantEventPack,2)
};
#endif
