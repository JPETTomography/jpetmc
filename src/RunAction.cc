#include "RunAction.hh"
#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


RunAction::RunAction(HistoManager* histo)
    : G4UserRunAction(),  fHistoManager(histo)
{

}

RunAction::~RunAction()
{
}

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
    fHistoManager->Book();
}

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
     fHistoManager->Save();
}

