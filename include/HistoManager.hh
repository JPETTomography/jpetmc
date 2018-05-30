#ifndef HistoManager_h 
#define HistoManager_h 1 

#include "globals.hh"
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include "JPetGeantDecayTree.h"
#include "JPetGeantScinHits.h"
#include "JPetGeantEventPack.h"
#include "JPetGeantEventInformation.h"
#include "DetectorHit.hh"
#include "VtxInformation.hh"


class TFile;
class TTree;
const G4int MaxHisto = 4;

/**
 * \class HistoManager
 * \brief class reach for informations stored during sumulations
 * and saves them into ROOT tree structures
 */
class HistoManager
{
    public:
        HistoManager();
        ~HistoManager();

        void Book(); //< call once; book all trees and histograms
        void Save(); //< call once; save all trees and histograms
        void SaveEvtPack(){fTree->Fill();}; 
        void Clear(){fEventPack->Clear();}; 

        void AddGenInfo(VtxInformation* info);
        void AddNewHit(DetectorHit*);
        void SetEventNumber(int x){fEventPack->SetEventNumber(x);};


    private:
        TFile*   fRootFile;
        TH1D*    fHisto[MaxHisto];
        TTree*   fTree; 
        TBranch* fBranchTrk;
        TBranch* fBranchScin;
        TBranch* fBranchEventPack;

        JPetGeantEventPack* fEventPack;
        JPetGeantEventInformation* fGeantInfo; 

};

#endif
