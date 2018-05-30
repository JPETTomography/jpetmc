#include "HistoManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include <vector>



HistoManager::HistoManager()
{
    fEventPack = new JPetGeantEventPack();
    fGeantInfo =  fEventPack->GetEventInformation();
}

HistoManager::~HistoManager()
{
}

void HistoManager::Book()
{
    G4String fileName = "mcGeant.root";
    fRootFile = new TFile(fileName,"RECREATE");  
    if (! fRootFile) { 
        G4cout << " HistoManager::Book :"                    
            << " problem creating the ROOT TFile "  
            << G4endl;   
        return; 
    }

    Int_t bufsize=32000;
    Int_t splitlevel=2;


    fTree = new TTree("T", "Tree keeps output from Geant simulation",splitlevel);
    fTree->SetAutoSave(1000000000); // autosave when 1 Gbyte written
    fBranchEventPack = fTree->Branch("eventPack", &fEventPack, bufsize, splitlevel);

}

void HistoManager::AddGenInfo(VtxInformation* info)
{
    bool is3g = info->GetThreeGammaGen();
    bool is2g = info->GetTwoGammaGen();
    bool isprompt = info->GetPromptGammaGen();

    if ( is2g || is3g )
    {
        fGeantInfo->SetThreeGammaGen(is3g);
        fGeantInfo->SetTwoGammaGen(is2g);
        fGeantInfo->SetVtxPosition(info->GetVtxPositionX(),info->GetVtxPositionY(),info->GetVtxPositionZ());
        fGeantInfo->SetLifetime(info->GetLifetime());
        fGeantInfo->SetRunNr(info->GetRunNr());
    }

    if (isprompt)
    {
        fGeantInfo->SetPromptGammaGen(isprompt);
        fGeantInfo->SetPromptLifetime(info->GetLifetime());
        fGeantInfo->SetVtxPromptPosition(info->GetVtxPositionX(),info->GetVtxPositionY(),info->GetVtxPositionZ());
    }


}

void HistoManager::AddNewHit(DetectorHit* hit)
{

  JPetGeantScinHits* geantHit =  fEventPack->ConstructNextHit();

   geantHit->Fill(
           fEventPack->GetEventNumber(),    //int evID, 
           hit->GetScinID(),    //int scinID, 
           hit->GetTrackID(),    //int trkID, 
           hit->GetTrackPDG(),    //int trkPDG, 
           hit->GetNumInteractions(),    //int nInter,
           hit->GetEdep()/keV,    //float ene, 
           hit->GetTime()/ps    //float time, 
               );

         // ugly way but there is no easy way to cast g4vector into root vector
    geantHit->SetHitPosition(
             hit->GetPosition().getX()/cm,
             hit->GetPosition().getY()/cm,
             hit->GetPosition().getZ()/cm);

    geantHit->SetPolarizationIn(
             hit->GetPolarizationIn().getX(),
             hit->GetPolarizationIn().getY(),
             hit->GetPolarizationIn().getZ());

    geantHit->SetPolarizationOut(
             hit->GetPolarizationOut().getX(),
             hit->GetPolarizationOut().getY(),
             hit->GetPolarizationOut().getZ());

    geantHit->SetMomentumIn(
             hit->GetMomentumIn().getX()/keV,
             hit->GetMomentumIn().getY()/keV,
             hit->GetMomentumIn().getZ()/keV);

    geantHit->SetMomentumOut(
             hit->GetMomentumOut().getX()/keV,
             hit->GetMomentumOut().getY()/keV,
             hit->GetMomentumOut().getZ()/keV);

    geantHit->SetGenGammaMultiplicity(hit->GetGenGammaMultiplicity());
    geantHit->SetGenGammaIndex(hit->GetGenGammaIndex());
}




void HistoManager::Save()
{
    if (! fRootFile) return;
     //fRootFile->Write(); 
     fTree->Write();
     fRootFile->Close(); 

    G4cout << "\n----> Histograms and ntuples are saved\n" << G4endl;
}

