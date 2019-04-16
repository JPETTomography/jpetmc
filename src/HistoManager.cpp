#include "HistoManager.h"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include <vector>



HistoManager::HistoManager()
: fMakeControlHisto(true)
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

      if ( MakeControlHisto()){
            BookHistograms(); 
      }
}

void HistoManager::BookHistograms()
{

      /// histograms - registered hits parameters
      fHisto[0] = new TH1F("gen_gamma_multiplicity","Generated gammas multiplicity",10,0,10);
    fHisto[1] = new TH1F("gen_hit_time","Gen hit time",100, 0.0, 15000.0);
    fHisto[2] = new TH1F("gen_hit_eneDepos","Gen hit ene deposition",750, 0.0, 1500.0);
    fHisto[3] = new TH1F("gen_hits_z_pos", "gen hits Z position",100, -60.0, 60.0);
    fHisto2D[0] = new TH2F("gen_hits_xy_pos","GEN hits XY pos",
      121, -60.5, 60.5,121, -60.5, 60.5);


    /// histograms - source position
    fHisto[4] = new TH1F("gen_lifetime","Gen lifetime",100, 0.0, 1500.0);
    fHisto[5] = new TH1F("gen_prompt_lifetime","Gen prompt lifetime",100, 0.0, 1500.0);


    fHisto2D[1] = new TH2F("gen_XY","GEN XY coordinates of annihilation point",
    121, -21.5, 21.5,121, -21.5, 21.5);
    fHisto2D[2] = new TH2F("gen_XZ","GEN XZ coordinates of annihilation point",
          121, -21.5, 21.5,121, -60.5, 60.5);
    fHisto2D[3] = new TH2F("gen_YZ","GEN YZ coordinates of  annihilation point",
          121, -21.5, 21.5,121, -60.5, 60.5);
    fHisto2D[4] = new TH2F("gen_prompt_XY","GEN prompt XY coordinates of annihilation point",
    121, -21.5, 21.5,121, -21.5, 21.5);
    fHisto2D[5] = new TH2F("gen_prompt_XZ","GEN prompt XZ coordinates of annihilation point",
          121, -21.5, 21.5,121, -60.5, 60.5);
    fHisto2D[6] = new TH2F("gen_prompt_YZ","GEN prompt YZ coordinates of  annihilation point",
            121, -21.5, 21.5,121, -60.5, 60.5);

            

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

            if ( MakeControlHisto()){
                  if(is2g) fHisto[0]->Fill(2);
                  if(is3g) fHisto[0]->Fill(3);
            
                  fHisto[4]->Fill(info->GetLifetime());
                  fHisto2D[1]->Fill(info->GetVtxPositionX(),info->GetVtxPositionY());
                  fHisto2D[2]->Fill(info->GetVtxPositionX(),info->GetVtxPositionZ());
                  fHisto2D[3]->Fill(info->GetVtxPositionY(),info->GetVtxPositionZ());
            }


    }

    if (isprompt)
    {
        fGeantInfo->SetPromptGammaGen(isprompt);
        fGeantInfo->SetPromptLifetime(info->GetLifetime());
        fGeantInfo->SetVtxPromptPosition(info->GetVtxPositionX(),info->GetVtxPositionY(),info->GetVtxPositionZ());
        fGeantInfo->SetRunNr(info->GetRunNr());

            if ( MakeControlHisto()){
                  fHisto[0]->Fill(1);
                  fHisto[5]->Fill(info->GetLifetime());
                  fHisto2D[4]->Fill(info->GetVtxPositionX(),info->GetVtxPositionY());
                  fHisto2D[5]->Fill(info->GetVtxPositionX(),info->GetVtxPositionZ());
                  fHisto2D[6]->Fill(info->GetVtxPositionY(),info->GetVtxPositionZ());
            }

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


      if ( MakeControlHisto()){
            fHisto[1]->Fill(hit->GetTime()/ps);
            fHisto[2]->Fill(hit->GetEdep()/keV);
            fHisto[3]->Fill(hit->GetPosition().getZ()/cm);
            fHisto2D[0]->Fill(hit->GetPosition().getX()/cm,
            hit->GetPosition().getY()/cm);

      }
}




void HistoManager::Save()
{
    if (! fRootFile) return;
     //fRootFile->Write(); 
     fTree->Write();

      if ( MakeControlHisto()){
            for(int i=0; i<MaxHisto; i++) fHisto[i]->Write();
            for(int i=0; i<MaxHisto2D; i++) fHisto2D[i]->Write();
      }
     fRootFile->Close(); 

    G4cout << "\n----> Histograms and ntuples are saved\n" << G4endl;
}

