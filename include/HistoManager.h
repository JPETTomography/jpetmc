#ifndef HistoManager_h 
#define HistoManager_h 1 

#include "globals.hh"
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include "JPetGeantDecayTree.h"
#include "JPetGeantScinHits.h"
#include "JPetGeantEventPack.h"
#include "JPetGeantEventInformation.h"
#include "DetectorHit.h"
#include "VtxInformation.h"


class TFile;
class TTree;
const int MaxHisto = 6;
const int MaxHisto2D = 7;

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

        void SetHistogramCreation(bool tf){fMakeControlHisto=tf;};
        bool MakeControlHisto(){return fMakeControlHisto;};


    private:
        bool     fMakeControlHisto; 
        TH1F*    fHisto[MaxHisto];
            // 0 - generated gamma
            // 1 - time
            // 2 - energy
            // 3 - z position
            // 4 - lifetime
            // 5 - prompt lifetime

        TH2F*    fHisto2D[MaxHisto2D];
            // 0 - XY hits
            // 1 - XY annihilation
            // 2 - XZ annihilation
            // 3 - YZ annihilation
            // 4 - XY prompt 
            // 5 - XZ prompt 
            // 6 - YZ prompt 

        TFile*   fRootFile;
        TTree*   fTree; 
        TBranch* fBranchTrk;
        TBranch* fBranchScin;
        TBranch* fBranchEventPack;

        JPetGeantEventPack* fEventPack;
        JPetGeantEventInformation* fGeantInfo; 

        void BookHistograms();
};

#endif
