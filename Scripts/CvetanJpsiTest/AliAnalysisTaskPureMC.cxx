/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                       *
* Author: Baldo Sahlmueller, Friederike Bock                     *
* Version 1.0                                 *
*                                       *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its    *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee, provided that the above copyright notice appears in all    *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims    *
* about the suitability of this software for any purpose. It is      *
* provided "as is" without express or implied warranty.               *
**************************************************************************/

//////////////////////////////////////////////////////////////////
//----------------------------------------------------------------
// Class used to do analysis on conversion photons + calo photons
//----------------------------------------------------------------
//////////////////////////////////////////////////////////////////
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliAnalysisTaskPureMC.h"
#include "AliVParticle.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliEventplane.h"
#include "AliInputEventHandler.h"
#include <vector>
#include <map>

ClassImp(AliAnalysisTaskPureMC)

//________________________________________________________________________
AliAnalysisTaskPureMC::AliAnalysisTaskPureMC(): AliAnalysisTaskSE(),
  fOutputContainer(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fMCStack(NULL),
  fHistNEvents(NULL),
  fHistXSection(NULL),
  fHistPtJpsi(NULL),
  fHistPtJpsiSel(NULL),
  fHistYJpsi(NULL),
  fHistPhiVsPt(NULL),
  fHistPhiVsPt05(NULL),
  fIsMC(1),
  fPrompt(kTRUE)
{
  
}

//________________________________________________________________________
AliAnalysisTaskPureMC::AliAnalysisTaskPureMC(const char *name):
  AliAnalysisTaskSE(name),
  fOutputContainer(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fMCStack(NULL),
  fHistNEvents(NULL),
  fHistXSection(NULL),
  fHistPtJpsi(NULL),
  fHistPtJpsiSel(NULL),
  fHistYJpsi(NULL),
  fHistPhiVsPt(NULL),
  fHistPhiVsPt05(NULL),
  fIsMC(1),
  fPrompt(kTRUE)
{
  // Define output slots here
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskPureMC::~AliAnalysisTaskPureMC()
{

}

//________________________________________________________________________
void AliAnalysisTaskPureMC::UserCreateOutputObjects(){

  printf("%s\n",fPrompt?"Prompt J/psi":"B -> J/psi");

  // Create histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
    fOutputContainer          = NULL;
  }
  if(fOutputContainer == NULL){
    fOutputContainer          = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }
  
  fHistNEvents                = new TH1F("NEvents", "NEvents", 3, -0.5, 2.5);
  fHistNEvents->Sumw2();
  fOutputContainer->Add(fHistNEvents);

  fHistXSection               = new TH1D("XSection", "XSection", 1000000, 0, 1e4);
  fHistXSection->Sumw2();
  fOutputContainer->Add(fHistXSection);
  
  fHistPtJpsi                 = new TH1F("PtJpsi", "PtJpsi", 400, 0, 200);
  fHistPtJpsi->Sumw2();
  fOutputContainer->Add(fHistPtJpsi);

  fHistPtJpsiSel              = new TH1F("PtJpsiSel", "PtJpsi", 120, 0, 12);
  fHistPtJpsiSel->Sumw2();
  fOutputContainer->Add(fHistPtJpsiSel);
  
  fHistYJpsi                 = new TH1F("YJpsi", "YJpsi", 400, -10, 10);
  fHistYJpsi->Sumw2();
  fOutputContainer->Add(fHistYJpsi);

  fHistPhiVsPt = new TH2F("PhiVsPt","",120,0,12,12,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fHistPhiVsPt->Sumw2();
  fOutputContainer->Add(fHistPhiVsPt);

  fHistPhiVsPt05 = new TH2F("PhiVsPt05","",120,0,12,12,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fHistPhiVsPt05->Sumw2();
  fOutputContainer->Add(fHistPhiVsPt05);
    
  PostData(1, fOutputContainer);
}

//_____________________________________________________________________________
void AliAnalysisTaskPureMC::UserExec(Option_t *)
{

  fInputEvent = InputEvent();
//   cout << "I found an Event" << endl;
  
  fMCEvent = MCEvent();
  if(fMCEvent == NULL) fIsMC = 0;
  
  if (fIsMC==0) return;
//   cout << "I found an MC header" << endl;
    
  fMCStack = fMCEvent->Stack();
  if(fMCStack == NULL) fIsMC = 0;
  if (fIsMC==0) return;
  
//   cout << "the stack is intact" << endl;
  
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();  

  if (TMath::Abs(mcProdVtxZ) < 10 ){
    fHistNEvents->Fill(0);
  } else {
    fHistNEvents->Fill(1);
  }
  
  AliGenEventHeader* mcEH = fMCEvent->GenEventHeader();
  AliGenPythiaEventHeader *pyH  = dynamic_cast<AliGenPythiaEventHeader*>(mcEH);
  AliGenHijingEventHeader *hiH  = 0;
  AliGenDPMjetEventHeader *dpmH = 0;
  
  // it can be only one save some casts
  // assuming PYTHIA and HIJING are the most likely ones...
  if(!pyH){
    hiH = dynamic_cast<AliGenHijingEventHeader*>(mcEH);
    if(!hiH){
      dpmH = dynamic_cast<AliGenDPMjetEventHeader*>(mcEH);
    }
  }
  
  // fetch the trials on a event by event basis, not from pyxsec.root otherwise 
  // we will get a problem when running on proof since Notify may be called 
  // more than once per file
  // consider storing this information in the AOD output via AliAODHandler
  Float_t ntrials = 0;
  if (!pyH || !hiH || dpmH) {
    AliGenCocktailEventHeader *ccEH = dynamic_cast<AliGenCocktailEventHeader *>(mcEH);
    if (ccEH) {
      TList *genHeaders = ccEH->GetHeaders();
      for (int imch=0; imch<genHeaders->GetEntries(); imch++) {
        if(!pyH)pyH = dynamic_cast<AliGenPythiaEventHeader*>(genHeaders->At(imch));
        if(!hiH)hiH = dynamic_cast<AliGenHijingEventHeader*>(genHeaders->At(imch));
        if(!dpmH)dpmH = dynamic_cast<AliGenDPMjetEventHeader*>(genHeaders->At(imch));
      }
    }
  }

  // take the trials from the p+p event
  if(hiH)ntrials = hiH->Trials();
  if(dpmH)ntrials = dpmH->Trials();
  if(pyH)ntrials = pyH->Trials();
  if(ntrials)fHistNEvents->Fill(2,ntrials); 
  
  Double_t xSection = 0;
  Double_t ptHard = 0;
  if (pyH) xSection = pyH->GetXsection();
  if (pyH) ptHard = pyH->GetPtHard();
  if (xSection) fHistXSection->Fill(xSection);
  
  ProcessMCParticles();

  
  PostData(1, fOutputContainer);
}


//________________________________________________________________________
void AliAnalysisTaskPureMC::ProcessMCParticles()
{

  // Loop over all primary MC particle 
  for(Long_t i = 0; i < fMCStack->GetNtrack(); i++) {
    // fill primary histograms
    TParticle* particle         = (TParticle *)fMCStack->Particle(i);
    if (!particle) continue;

    // discard unphysical particles from some generators
    if (particle->Pt() == 0 || particle->Energy() <= 0)
      continue;

    // J/psi
    if (TMath::Abs(particle->GetPdgCode()) != 443 ) continue;

    Bool_t hasMother            = kFALSE;
    if (particle->GetMother(0)>-1) 
      hasMother                 = kTRUE;
    TParticle* motherParticle   = NULL;
    if( hasMother ) 
      motherParticle            = (TParticle *)fMCStack->Particle(particle->GetMother(0));
    if (motherParticle) 
      hasMother                 = kTRUE;
    else 
      hasMother                 = kFALSE;

    TParticle *daughter=fMCStack->Particle(particle->GetDaughter(0));
    //    if (TMath::Abs(daughter->GetPdgCode()) == 443) continue;

    if (0) {
      printf("%ld %f %f %d     ",i,particle->Pt(),particle->Y(),fMCStack->IsPhysicalPrimary(i));
      if (hasMother) printf("%d %f %f    ",particle->GetMother(0),motherParticle->Pt(),motherParticle->Y());
      if (daughter) printf("%d",fMCStack->IsPhysicalPrimary(particle->GetDaughter(0)));
      printf("\n");
    }

    if (!fMCStack->IsPhysicalPrimary(particle->GetDaughter(0))) continue;

    // B0 511  -511 -> B*_2c 545 - 545
    //    if (motherParticle) motherParticle->Print();
    // if (!fPrompt) {
    //   if (!hasMother) continue;
    //   Int_t motherPdg = TMath::Abs(motherParticle->GetPdgCode());
    //   if (motherPdg < 511 || motherPdg > 545) continue;
    //   if (1) printf("From B: %d %f %f    %f %f    %f %f\n",
    // 		    motherPdg,
    // 		    motherParticle->Pt(),particle->Pt(),
    // 		    motherParticle->Y(),particle->Y(),
    // 		    motherParticle->Phi(),particle->Phi());
    // }
    
    fHistPtJpsi->Fill(particle->Pt());
    fHistYJpsi->Fill(particle->Y());

    if (TMath::Abs(particle->Y())<2.5 || TMath::Abs(particle->Y())>4.0) continue;
    
    Double_t pT = particle->Pt();
    Double_t y = particle->Y();
    Double_t phi = particle->Phi();

    CorrelateWithHadrons(i,pT,y,phi);
  }
  if (0) printf("end of event\n");
}

void AliAnalysisTaskPureMC::CorrelateWithHadrons(Int_t iJpsi, Double_t pT, Double_t y, Double_t phi)
{
  fHistPtJpsiSel->Fill(pT);
  for(Long_t i = 0; i < fMCStack->GetNtrack(); i++) {
    TParticle* particle         = (TParticle *)fMCStack->Particle(i);
    if (!particle) continue;
    if (particle->Pt() == 0 || particle->Energy() <= 0) continue;

    if (0) printf("DEBUG: %d  %d  %s  %d  %d\n",i,particle->GetPdgCode(),particle->GetName(),
		  particle->GetUniqueID(),fMCStack->IsPhysicalPrimary(i));

    if (!fMCStack->IsPhysicalPrimary(i)) continue;

    AliMCParticle *track = static_cast<AliMCParticle*>(fMCEvent->GetTrack(i));
    if (!track) {
      printf("ERROR: Could not receive track %ld", i);
      continue;
    }
    if (TMath::Abs(particle->Pt()-track->Pt())>1e-6) {
      printf("ERROR: particles from stack and mcevent are not the same!!!");
      exit(0);
    }
    
    if (track->Charge() == 0) continue;

    if (particle->GetMother(0) == iJpsi) continue; // remove j/psi daugthers
    //    if (particle->Pt()>4.) continue;
    Double_t deltay = (y>0.?y - particle->Y():particle->Y()-y);
    if (deltay<1.5 || deltay>5.0) continue;
    Double_t dphi = phi - particle->Phi();
    if (dphi >  1.5*TMath::Pi()) dphi -= TMath::TwoPi();
    if (dphi < -0.5*TMath::Pi()) dphi += TMath::TwoPi();
    fHistPhiVsPt->Fill(pT,dphi);
    if (particle->Pt()>=0.5) fHistPhiVsPt05->Fill(pT,dphi);
  }
}

//________________________________________________________________________
void AliAnalysisTaskPureMC::Terminate(const Option_t *)
{
  
  //fOutputContainer->Print(); // Will crash on GRID
}

