//
//  AliAnalysisTaskMyMuonTree_AOD.cxx
//
//
//  Created by Marc Arène on 7/31/13.
//
//

#include "AliAnalysisTaskMyMuonTree_AOD.h"
#include "TChain.h"
#include "TTree.h"

#include "TLorentzVector.h"
#include "TMath.h"

#include "TH1I.h"
#include "TList.h"

#include <TClonesArray.h>

#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"

#include "AliAnalysisDataContainer.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisTaskSE.h"

#include "AliAODCluster.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"

#include "AliAODVZERO.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"
#include "AliCounterCollection.h"
#include "AliEventplane.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMuonTrackCuts.h"

ClassImp(MyEventLight);
ClassImp(TrackletLight);
ClassImp(CorrelationLight);
ClassImp(DimuonLight);
ClassImp(MCDimuonLight);

int buildno = 0;
int newlyremoved = 0;

// TClonesArray *MyEventLight::fgTracks = 0;
// TClonesArray *MyEventLight::fgDimuonLights = 0;
// TClonesArray *MyEventLight::fgMCDimuonLights = 0;

// simple task to fill a tree for dimuon analysis
//
// author: Hongyan Yang @ CEA, Saclay
// email: Hongyan.Yang@cea.fr
// 17.01.2011

//________________________________________________________________________
AliAnalysisTaskMyMuonTree_AOD::AliAnalysisTaskMyMuonTree_AOD()
    : AliAnalysisTaskSE("tree Analysis Task"), fMC(kFALSE), fHistEvents(0x0),
      fMyMuonTree(0x0), fOutputList(0x0), fEventCounters(0x0), fEvent(0x0),
      fMuonTrackCuts(0x0) {
  //
  // constructor
  //
  Int_t islot = 1;
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  DefineOutput(islot++, TList::Class());
  DefineOutput(islot++, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskMyMuonTree_AOD::AliAnalysisTaskMyMuonTree_AOD(const char *name)
    : AliAnalysisTaskSE(name), fMC(kFALSE), fHistEvents(0x0), fMyMuonTree(0x0),
      fOutputList(0x0), fEventCounters(0x0), fEvent(0x0), fMuonTrackCuts(0x0) {
  //
  // non-default constructor
  // Define input and output slots here
  //
  Int_t islot = 1;
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  DefineOutput(islot++, TList::Class());
  DefineOutput(islot++, TTree::Class());

  fEvent = new MyEventLight();
}

//________________________________________________________________________
AliAnalysisTaskMyMuonTree_AOD::AliAnalysisTaskMyMuonTree_AOD(
    const char *name, AliMuonTrackCuts *cuts, Bool_t bMC)
    : AliAnalysisTaskSE(name), fMC(bMC), fHistEvents(0x0), fMyMuonTree(0x0),
      fOutputList(0x0), fEventCounters(0x0), fEvent(0x0), fMuonTrackCuts(cuts) {
  //
  // non-default constructor
  // Define input and output slots here
  //
  Int_t islot = 1;
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  DefineOutput(islot++, TList::Class());
  DefineOutput(islot++, TTree::Class());

  fEvent = new MyEventLight();
}

//________________________________________________________________________
AliAnalysisTaskMyMuonTree_AOD::AliAnalysisTaskMyMuonTree_AOD(
    const AliAnalysisTaskMyMuonTree_AOD &ref)
    : AliAnalysisTaskSE(ref), fMC(ref.fMC), fHistEvents(ref.fHistEvents),
      fMyMuonTree(ref.fMyMuonTree), fOutputList(ref.fOutputList),
      fEventCounters(ref.fEventCounters), fEvent(ref.fEvent),
      fMuonTrackCuts(ref.fMuonTrackCuts) {
  //
  // Copy Constructor
  //

  ref.Copy(*this);
}

//________________________________________________________________________
AliAnalysisTaskMyMuonTree_AOD &AliAnalysisTaskMyMuonTree_AOD::operator=(
    const AliAnalysisTaskMyMuonTree_AOD &ref) {
  //
  // assignment operator
  //

  if (this == &ref)
    return *this;
  AliAnalysisTaskSE::operator=(ref);
  fMC = ref.fMC;
  fHistEvents = ref.fHistEvents;
  fMyMuonTree = ref.fMyMuonTree;
  fOutputList = ref.fOutputList;
  fEventCounters = ref.fEventCounters;
  fEvent = ref.fEvent;
  fMuonTrackCuts = ref.fMuonTrackCuts;

  return *this;
}
//________________________________________________________________________
AliAnalysisTaskMyMuonTree_AOD::~AliAnalysisTaskMyMuonTree_AOD() {
  //
  //  default destructor
  //

  if (fHistEvents)
    delete fHistEvents;
  if (fMyMuonTree)
    delete fMyMuonTree;
  if (fOutputList)
    delete fOutputList;
  if (fEventCounters)
    delete fEventCounters;
  if (fEvent)
    delete fEvent;

  Printf("analysis done\n");
}

//________________________________________________________________________
void AliAnalysisTaskMyMuonTree_AOD::UserCreateOutputObjects() {
  //
  // Create output tree
  // Called once
  //

  TString listName = GetOutputSlot(1)->GetContainer()->GetName();
  TString listTitle = "jpsi in Pb+Pb collisions";

  Int_t islot = 1;

  OpenFile(islot++);
  if (!fOutputList)
    fOutputList = new TList();
  fOutputList->SetName(listName);
  fHistEvents =
      new TH1I("histEvents", "Number of Events in the Analysis", 2, -0.5, 1.5);
  fOutputList->Add(fHistEvents);

  fEventCounters = new AliCounterCollection("eventCounters");
  fEventCounters->AddRubric(
      "trgtype", "C0V0M/C0VSC/C0V0H/CV0H7/CMID7/CINT7/CMUL7/CMLL7/CMSL7/CMSH7");
  fEventCounters->AddRubric("trgcluster", "MUFAST/CENT/CENTNOPMD");
  fEventCounters->AddRubric("evt", "ps/nps");
  fEventCounters->AddRubric("v0cent", 101);
  fEventCounters->AddRubric("run", 137);
  fEventCounters->Init();
  fOutputList->Add(fEventCounters);

  OpenFile(islot++);
  if (!fMyMuonTree)
    fMyMuonTree = new TTree("MyMuonTree", "muon track tree");
  fMyMuonTree->Branch("event", &fEvent, 32000, 2);
  //  fMyMuonTree->Print();

  islot = 1;
  PostData(islot++, fOutputList);
  PostData(islot++, fMyMuonTree);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskMyMuonTree_AOD::UserExec(Option_t *) {
  // Main loop
  // Called for each event

  if (!fInputEvent) {
    AliError("Reconstructed Event not available");
    return;
  }

  AliAODEvent *fAOD = dynamic_cast<AliAODEvent *>(InputEvent());

  if (!fAOD) {
    AliError("Error: AOD Event required for AOD Analysis!!!");
    return;
  }

  // ---- CHOICE OF ANALYSIS ----
  if (fMC) {

    if (!(fAOD->FindListObject(AliAODMCParticle::StdBranchName()))) {
      AliError("No MC Array, but MC Data required");
      return;
    }
      Printf("LOLILOL1 - MC Analysis\n");
  }
  Int_t runNumber = fAOD->GetRunNumber();
  Bool_t isPhysSelected =
      (((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()
                                     ->GetInputEventHandler()))
           ->IsEventSelected() &
       (AliVEvent::kINT7 | AliVEvent::kINT7inMUON | AliVEvent::kMUSH7 |
        AliVEvent::kMUS7 | AliVEvent::kMUU7 | AliVEvent::kMUL7));
  Bool_t isPhysSelectedMUFAST =
      (((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()
                                     ->GetInputEventHandler()))
           ->IsEventSelected() &
       (AliVEvent::kINT7inMUON | AliVEvent::kMUSH7 | AliVEvent::kMUS7 |
        AliVEvent::kMUU7 | AliVEvent::kMUL7));
  Bool_t isPhysSelectedMUL =
      (((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()
                                     ->GetInputEventHandler()))
           ->IsEventSelected() &
       (AliVEvent::kMUU7));
  Bool_t isPhysSelectedINT7CENT =
    (((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()
                                   ->GetInputEventHandler()))
         ->IsEventSelected() &
     (AliVEvent::kINT7));
    
    Bool_t isPhysSelectedINT7MUFAST =
    (((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()
                                   ->GetInputEventHandler()))
         ->IsEventSelected() &
     (AliVEvent::kINT7inMUON));
  // Useful trigger class names composed of
  // Class

  const Char_t *trigNames[6] = {
      "CINT7-B-NOPF-MUFAST", "CMSL7-B-NOPF-MUFAST", "CMSH7-B-NOPF-MUFAST",
      "CMUL7-B-NOPF-MUFAST", "CMLL7-B-NOPF-MUFAST", "CINT7-B-NOPF-CENT"}; // >>> CENTNOTRD

  TString firedTrigClasses = fAOD->GetFiredTriggerClasses();

  Bool_t triggerCINT7 = (firedTrigClasses.Contains(trigNames[0])) ? 1 : 0;
  Bool_t triggerCMSL7 = (firedTrigClasses.Contains(trigNames[1])) ? 1 : 0;
  Bool_t triggerCMSH7 = (firedTrigClasses.Contains(trigNames[2])) ? 1 : 0;
  Bool_t triggerCMUL7 = (firedTrigClasses.Contains(trigNames[3])) ? 1 : 0;
  Bool_t triggerCMLL7 = (firedTrigClasses.Contains(trigNames[4])) ? 1 : 0;
  Bool_t triggerCINT7CENT = (firedTrigClasses.Contains(trigNames[5])) ? 1 : 0;

  Bool_t clusterMUFAST = 0;
  Bool_t clusterCENTNOTRD = 0;

  AliMultSelection *MultSelection =
      (AliMultSelection *)fAOD->FindListObject("MultSelection");
  //	Printf("MultSelection %p \n",MultSelection);
  Int_t v0centp =
      TMath::Min(Int_t(MultSelection->GetMultiplicityPercentile("V0M")), 100);

  Bool_t isTriggerSelected =
      (triggerCINT7 || triggerCMSL7 || triggerCMUL7 || triggerCMLL7) ? 1 : 0;
  //	Bool_t isTriggerSelected = (triggerCMULB || triggerCMLLB) ? 1 : 0;

  if (triggerCINT7CENT) {
      Printf("LOLILOL-triggerCINT7CENT is true\n");
    fEventCounters->Count(
        Form("evt:%s/trgtype:%s/trgcluster:%s/run:%d/v0cent:%d",
             (((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()
                                            ->GetInputEventHandler()))
                  ->IsEventSelected() &
              (AliVEvent::kINT7))
                 ? "ps"
                 : "nps",
             "CINT7", "CENT", runNumber, v0centp));
    clusterCENTNOTRD = 1;
  }
  if (triggerCINT7) {
      Printf("LOLILOL-triggerCINT7 is true\n");
    fEventCounters->Count(
        Form("evt:%s/trgtype:%s/trgcluster:%s/run:%d/v0cent:%d",
             (((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()
                                            ->GetInputEventHandler()))
                  ->IsEventSelected() &
              (AliVEvent::kINT7inMUON))
                 ? "ps"
                 : "nps",
             "CINT7", "MUFAST", runNumber, v0centp));
    clusterMUFAST = 1;
  }
  if (triggerCMUL7) {
      Printf("LOLILOL-triggerCMUL7 is true\n");
    fEventCounters->Count(
        Form("evt:%s/trgtype:%s/trgcluster:%s/run:%d/v0cent:%d",
             (((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()
                                            ->GetInputEventHandler()))
                  ->IsEventSelected() &
              (AliVEvent::kMUU7))
                 ? "ps"
                 : "nps",
             "CMUL7", "MUFAST", runNumber, v0centp));
    clusterMUFAST = 1;
  }
  if (triggerCMLL7) {
      Printf("LOLILOL-triggerCMLL7 is true\n");
    fEventCounters->Count(
        Form("evt:%s/trgtype:%s/trgcluster:%s/run:%d/v0cent:%d",
             (((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()
                                            ->GetInputEventHandler()))
                  ->IsEventSelected() &
              (AliVEvent::kMUL7))
                 ? "ps"
                 : "nps",
             "CMLL7", "MUFAST", runNumber, v0centp));
    clusterMUFAST = 1;
  }
  if (triggerCMSL7) {
      Printf("LOLILOL-triggerCMSL7 is true\n");
    fEventCounters->Count(
        Form("evt:%s/trgtype:%s/trgcluster:%s/run:%d/v0cent:%d",
             (((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()
                                            ->GetInputEventHandler()))
                  ->IsEventSelected() &
              (AliVEvent::kMUS7))
                 ? "ps"
                 : "nps",
             "CMSL7", "MUFAST", runNumber, v0centp));
    clusterMUFAST = 1;
  }
  if (triggerCMSH7) {
      Printf("LOLILOL-triggerCMSH7 is true\n");
    fEventCounters->Count(
        Form("evt:%s/trgtype:%s/trgcluster:%s/run:%d/v0cent:%d",
             (((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()
                                            ->GetInputEventHandler()))
                  ->IsEventSelected() &
              (AliVEvent::kMUSH7))
                 ? "ps"
                 : "nps",
             "CMSH7", "MUFAST", runNumber, v0centp));
    clusterMUFAST = 1;
  }

  //	UInt_t isEventSelected =
  //((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  //	//	Bool_t isPhysSelected =
  //(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()
  //& ( AliVEvent::kINT7 | AliVEvent::kMUSH7 | AliVEvent::kMUS7 |
  // AliVEvent::kMUU7 | AliVEvent::kMUL7 )); 	Bool_t isPhysSelected =
  //(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()
  //& AliVEvent::kAny);
  //  if(!isPhysSelected){
  //    AliInfo("AliInfo: not passing physics selection");
  //    //    return;
  //  }
  //  if(!isPhysSelectedMUFAST){
  //    AliInfo("AliInfo: not passing physics selection for MUFAST triggers");
  //    //    return;
  //  }
  //	TString sisPhysSelected = isPhysSelected ? "yes" : "no";
  //	cout << "isEventSelected " << isEventSelected << " isPhysSelected " <<
  // sisPhysSelected << endl;

//    Printf("isPhysSelectedMUFAST: %d",isPhysSelectedMUFAST );
//    Printf("isPhysSelected: %d",isPhysSelected );
//    Printf("triggerCMUL7: %d",triggerCMUL7 );
//    Printf("clusterMUFAST: %d",clusterMUFAST );
//    Printf("Number of events built: %d",buildno );
//    Printf("newlyremoved: %d",newlyremoved );
  if (!isPhysSelectedINT7CENT) //isPhysSelectedINT7CENT isPhysSelectedMUL <<<<
  {
    //  Printf("Event !isPhysSelectedINT7CENT\n");
    fHistEvents->Fill(0);
  }
  else
  {
    fHistEvents->Fill(1);
  //  Printf("Event isPhysSelectedINT7CENT\n");
  }

 //   if (isPhysSelectedINT7CENT) {
// if (isPhysSelectedMUL && triggerCMUL7 && clusterMUFAST) {//isPhysSelectedMUFAST && triggerCMUL7 && clusterMUFAST) {
//if (isPhysSelectedINT7CENT && triggerCINT7CENT) { // >>> CINT7 ou CINT7CENT ?
    // Printf("LOLILOL will build fEvent and fill muon tree\n");
     buildno++;
    fEvent->Build(fAOD, isPhysSelectedINT7CENT, fMuonTrackCuts, fMC); //>>>>> isPhysSelectedMUL isPhysSelectedMUFAST or isPhysSelectedINT7CENT ou INT7MUFAST
    fMyMuonTree->Fill();
 // }

  Int_t islot = 1;
  PostData(islot++, fOutputList);
  PostData(islot++, fMyMuonTree);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskMyMuonTree_AOD::Terminate(Option_t *) {
  //
  // Called once at the end of the query
  //
  //	fEventCounters->Print();
  //	Printf("CMUL7/MUFAST");
  //	fEventCounters->Print("run","evt:PS/trgtype:CMUL7/trgcluster:MUFAST");
  //	fEventCounters->PrintSum("evt:PS/trgtype:CMUL7/trgcluster:MUFAST");
  //	Printf("CINT7/MUFAST");
  //	fEventCounters->Print("run","evt:PS/trgtype:CINT7/trgcluster:MUFAST");
  //	fEventCounters->PrintSum("evt:PS/trgtype:CINT7/trgcluster:MUFAST");
  //	Printf("CINT7/CENTNOPMD");
  //	fEventCounters->Print("run","evt:PS/trgtype:CINT7/trgcluster:CENTNOPMD");
  //	fEventCounters->PrintSum("evt:PS/trgtype:CINT7/trgcluster:CENTNOPMD");
  PostProcess();
  return;
}

//________________________________________________________________________
void AliAnalysisTaskMyMuonTree_AOD::PostProcess() {
  //
  //  do post analysis
  //

  fOutputList = dynamic_cast<TList *>(GetOutputData(1));

  if (!fOutputList) {
    printf("ERROR: fOutputList not available\n");
    return;
  }

  Printf("analysis finished!\n");
}

void AliAnalysisTaskMyMuonTree_AOD::NotifyRun() {
  /// Set run number for cuts
  //	fMuonTrackCuts->SetRun(fCurrentRunNumber); // Run number needed to get
  // parameters from OADB
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
  //	fMuonTrackCuts->SetPassNumber(1);
  fMuonTrackCuts->SetRun(
      (AliInputEventHandler
           *)(AliAnalysisManager::GetAnalysisManager()
                  ->GetInputEventHandler())); // Run number needed to get
                                              // parameters from OADB
}

//___
TrackletLight::TrackletLight()
    : TObject(), fEta(0), fPhi(0), fDPhi(0), fLabel(0) {}

//___
MCDimuonLight::MCDimuonLight()
    : TObject(), fPDG(0), fInvMass(0), fPt(0), fEta(0), fY(0), fPhi(0),
      fLabel(0) {}

//___
DimuonLight::DimuonLight()
    : TObject(), fCharge(0), fInvMass(0), fPt(0), fEta(0), fY(0), fPhi(0),
      fLabel(0) {}

//___
CorrelationLight::CorrelationLight()
    : TObject(), fDeltaPhi(0), fDeltaEta(0),
      fLabel(0) {}

//_______________________________________________________________________
MyEventLight::MyEventLight()
    : TObject(), fNTracklets(0), fNDimuons(0), fNMCDimuons(0), fNCorrelations(0), fRunNumber(0),
    fNPileupVtx(0), fIsPileupFromSPDMultBins(0),fIsPileupClustVsTkl(0) ,fIsPileupReimplemented(0), fVertexZ(0), fVertexXSPD(0), fVertexYSPD(0), fVertexZSPD(0), fVertexSigmaZ(0), fSPDVertexSigmaZ(0), fVertexNC(0), fSPDVertexNC(0), fIsFromVertexerZ(0), fCentralityV0M(0), fCentralityTKL(0), fCentralityZNAC(0), fCentralitySPDTracklets(0), fCentralitySPDClusters(0), fV0MValue(0), fADMValue(0), fADAValue(0), fADCValue(0), fTKLValue(0), fZNACValue(0), fSPDTrackletsValue(0), fSPDClustersValue(0), fPassPhysicsSelection(0), fV0AValue(0), fV0CValue(0), fZNAValue(0), fZNCValue(0), fCentralityV0A(0), fCentralityV0C(0), fCentralityADA(0), fCentralityADC(0), fCentralityADM(0), fCentralityZNA(0), fCentralityZNC(0), fTriggerAnalysis(new AliTriggerAnalysis()),
      fTracklets(new TClonesArray("TrackletLight", 0)),
      fDimuons(new TClonesArray("DimuonLight", 0)),
      fCorrelations(new TClonesArray("CorrelationLight", 0)),
      fMCDimuons(new TClonesArray("MCDimuonLight", 0)) {

  //
  // Create a event object
  //

  //  if(!fgTracks) fgTracks = new TClonesArray("Track", 20);
  //  fTracks = fgTracks;
  //
  //  if(!fgDimuonLights) fgDimuonLights = new TClonesArray("DimuonLight", 400);
  //  fDimuonLights = fgDimuonLights;
  //
  //	if(!fgMCDimuonLights) fgMCDimuonLights = new
  // TClonesArray("MCDimuonLight", 2);
  //  fMCDimuonLights = fgMCDimuonLights;
}

//_______________________________________________________________________
MyEventLight::~MyEventLight() {
  //
  //  destructor
  //

  Reset();
}

//_______________________________________________________________________
void MyEventLight::Build(AliAODEvent *aod, Bool_t isPhysSelected,
                         AliMuonTrackCuts *trackCuts, Bool_t bMC) {
  //
  // build one event
  //
  // save current object count
  Int_t ObjectNumber = TProcessID::GetObjectCount();
  Clear();
    
   // printf("BUILD CALLED: Event starts building\n");

  fNTracklets = fNDimuons = fNMCDimuons = fNCorrelations = 0;
  //  Double_t sigma[3];
  fRunNumber = aod->GetRunNumber();

  fNPileupVtx = aod->GetNumberOfPileupVerticesSPD(); //FIXME ok


  const AliAODVertex *primVtx = aod->GetPrimaryVertex();
  const AliAODVertex *SPDVtx = aod->GetPrimaryVertexSPD();
    Float_t vtxErr[3];
    Double_t covMatrix[6];
    Double_t covMatrixSPD[6];
  fVertexX = primVtx->GetX();
  fVertexY = primVtx->GetY();
  fVertexZ = primVtx->GetZ();
    fVertexXSPD = SPDVtx->GetX();
    fVertexYSPD = SPDVtx->GetY();
    fVertexZSPD = SPDVtx->GetZ();
    primVtx->GetCovarianceMatrix(covMatrix);
    SPDVtx->GetCovarianceMatrix(covMatrixSPD);
 // primVtx->GetSigmaXYZ(vtxErr);
    fVertexSigmaZ = sqrt(covMatrix[5]);
    fSPDVertexSigmaZ = sqrt(covMatrixSPD[5]);
  fVertexNC = primVtx->GetNContributors();
    fSPDVertexNC = SPDVtx->GetNContributors();
    fIsFromVertexerZ = SPDVtx->IsFromVertexerZ();
//    primVtx->GetSigmaXYZ(sigma);
//    fVertexResZ = sigma[2];
//  fVertexNC = -1;

   // printf("PASS PHYSICS SELECTION:\n");
    if (!isPhysSelected){
    fPassPhysicsSelection = 0; // not passing physics selection
   // printf("NO\n");
    }
    else{
    fPassPhysicsSelection = 1; // pass physics selection
   // printf("YES\n");
    }

  const Char_t *trigNames[6] = {"CINT7-B-NOPF-MUFAST", "CMSL7-B-NOPF-MUFAST",
                                "CMSH7-B-NOPF-MUFAST", "CMUL7-B-NOPF-MUFAST",
                                "CMLL7-B-NOPF-MUFAST", "CINT7-B-NOPF-CENT"};

  // fFiredTriggerClass = aod->GetFiredTriggerClasses();

  fTriggerCINT7 = (fFiredTriggerClass.Contains(trigNames[0])) ? 1 : 0;
  fTriggerCMSL7 = (fFiredTriggerClass.Contains(trigNames[1])) ? 1 : 0;
  fTriggerCMSH7 = (fFiredTriggerClass.Contains(trigNames[2])) ? 1 : 0;
  fTriggerCMUL7 = (fFiredTriggerClass.Contains(trigNames[3])) ? 1 : 0;
  fTriggerCMLL7 = (fFiredTriggerClass.Contains(trigNames[4])) ? 1 : 0;
  fTriggerCINT7CENT = (fFiredTriggerClass.Contains(trigNames[5])) ? 1 : 0;

  fClusterMUFAST = 1;
  fClusterCENT = 0;
  fClusterCENTNOPMD = 0;

  //	if(fTriggerCINTB || fTriggerCMSLB || fTriggerCMSHB || fTriggerCMULB ||
  // fTriggerCMLLB)
  //  if(fTriggerCMUL7 || fTriggerCMLL7)
//<<<<< if (fTriggerCMUL7)
   // printf("PASS TRIGGER CINT7CENT:\n");
    if (fTriggerCINT7CENT){
    fPassTriggerSelection = 1;
   // printf("NO\n");
    }
    else{
    fPassTriggerSelection = 0;
  //  printf("YES\n");
    }

  AliMultSelection *MultSelection =
      (AliMultSelection *)aod->FindListObject("MultSelection");
  //	Printf("MultSelection %p \n",MultSelection);
  fCentralityV0M = MultSelection->GetMultiplicityPercentile("V0M");
    fCentralityV0A = MultSelection->GetMultiplicityPercentile("V0A");
    fCentralityV0C = MultSelection->GetMultiplicityPercentile("V0C");
    fCentralityADM = 0; //MultSelection->GetMultiplicityPercentile("ADM");
    fCentralityADA = 0; //MultSelection->GetMultiplicityPercentile("ADA");
    fCentralityADC = 0; //MultSelection->GetMultiplicityPercentile("ADC");
    fCentralityZNA = 0; //MultSelection->GetMultiplicityPercentile("ZNApp");
    fCentralityZNC = 0; //MultSelection->GetMultiplicityPercentile("ZNCpp");
    fCentralityZNAC = 0; //MultSelection->GetMultiplicityPercentile("ZNACpp");
    fCentralitySPDTracklets = MultSelection->GetMultiplicityPercentile("SPDTracklets"); // >>>>>>>> ,true);
    fCentralitySPDClusters = MultSelection->GetMultiplicityPercentile("SPDClusters");
    AliMultEstimator *AliMuEst = MultSelection->GetEstimator("SPDTracklets");
    fSPDTrackletsValue = AliMuEst->GetValue();
    AliMultEstimator *AliMuEst2 = MultSelection->GetEstimator("SPDClusters");
    fSPDClustersValue = aod->GetNumberOfITSClusters(0)+aod->GetNumberOfITSClusters(1); //aod->GetNumberOfITSClusters(0)+aod->GetNumberOfITSClusters(1)
   // AliMultEstimator *AliMuEst5 = MultSelection->GetEstimator("ZNCpp");
    fZNCValue = 0;//AliMuEst5->GetValue();
   // AliMultEstimator *AliMuEst6 = MultSelection->GetEstimator("ZNApp");
    fZNAValue = 0;//AliMuEst6->GetValue();
   // AliMultEstimator *AliMuEst7 = MultSelection->GetEstimator("ZNACpp");
    fZNACValue = 0;//AliMuEst7->GetValue();
    AliMultEstimator *AliMuEst8 = MultSelection->GetEstimator("V0C");
    fV0CValue = AliMuEst8->GetValue();
    AliMultEstimator *AliMuEst9 = MultSelection->GetEstimator("V0A");
    fV0AValue = AliMuEst9->GetValue();
    AliMultEstimator *AliMuEst10 = MultSelection->GetEstimator("V0M");
    fV0MValue = AliMuEst10->GetValue();
  //  AliMultEstimator *AliMuEst11 = MultSelection->GetEstimator("ADC");
    fADCValue = 0;//AliMuEst11->GetValue();
   // AliMultEstimator *AliMuEst12 = MultSelection->GetEstimator("ADA");
    fADAValue = 0;//AliMuEst12->GetValue();
 //   AliMultEstimator *AliMuEst13 = MultSelection->GetEstimator("ADM");
    fADMValue = 0;//AliMuEst13->GetValue();
//    fCentralitySPDClusters = MultSelection->GetMultiplicityPercentile("SPDClusters");
//    fCentralityTKLvsV0M = MultSelection->GetMultiplicityPercentile("TKLvsV0M");
    
    

  AliAODTrack *track;
  TrackletLight *tracklet;

  Int_t nMuonPos = 0;
  Int_t nMuonNeg = 0;

  AliAODTracklets *trklets = (AliAODTracklets *)aod->GetTracklets();
  printf("trklets %p has %d \n", trklets, trklets->GetNumberOfTracklets());
  if (!trklets)
    AliFatal("AliAODTracklets not found");
  trklets->Print("t");
 // fNTracklets = trklets->GetNumberOfTracklets();
  for (Int_t itrklt = 0; itrklt < trklets->GetNumberOfTracklets(); itrklt++) {
//    if (!IsSelectedTracklet(trklets, itrklt))
//      continue;
      if(kFALSE){ //<<
     tracklet = AddTracklet();
     tracklet->fEta = (Float_t)-TMath::Log(TMath::Tan(trklets->GetTheta(itrklt) / 2));
     tracklet->fPhi = (Float_t)trklets->GetPhi(itrklt);
     tracklet->fDPhi = (Float_t)trklets->GetDeltaPhi(itrklt);
     tracklet->fLabel = (Float_t)trklets->GetLabel(itrklt, 0);
      }
  }
    Int_t nTracklets = trklets->GetNumberOfTracklets();
    
    fIsPileupFromSPDMultBins= aod->IsPileupFromSPDInMultBins();
    
    if(nTracklets<40){
        fIsPileupFromSPDMultBins = aod->IsPileupFromSPD(3,0.8);
    }
    else{
        fIsPileupFromSPDMultBins = aod->IsPileupFromSPD(5,0.8);
    } //<<<
    
   // printf("SPDClusters=%i nTracklets=%i et %i ",int(fSPDClustersValue),int(nTracklets),int(100+(4*nTracklets)));
//    if(fSPDClustersValue>(100+(4*fNTracklets))){
//        fIsPileupClustVsTkl = 1;
//       // printf("This is PU mi boi");
//    }
//    else{
//        fIsPileupClustVsTkl = 0;
//      //  printf("Not PU");
//    }
    
//    AliOADBTriggerAnalysis *triggerAnalysisOADB = new AliOADBTriggerAnalysis();
//    fTriggerAnalysis->SetParameters(triggerAnalysisOADB);
//    
//    double fV0MOnVsOfA = fTriggerAnalysis->GetV0MOnVsOfA();
//    double fV0MOnVsOfB = fTriggerAnalysis->GetV0MOnVsOfB();
//    
//    int isV0MOnVsOfPileup = 0;
//    {
//      AliVVZERO* vzero = aod->GetVZEROData();
//      if (!vzero) {
//          isV0MOnVsOfPileup = 3;
//      }
//      if (!vzero->TestBit(AliVVZERO::kTriggerChargeBitsFilled)){
//          isV0MOnVsOfPileup = 3;
//      }
//      // V0A0 excluded from online V0A charge sum => excluding also from offline sum for consistency
//      Float_t on = vzero->GetTriggerChargeA()+vzero->GetTriggerChargeC();
//      Float_t of = vzero->GetMTotV0A()-vzero->GetMRingV0A(0)+vzero->GetMTotV0C();
//    if(isV0MOnVsOfPileup == 0){
//        if(on < fV0MOnVsOfA + fV0MOnVsOfB*of){
//            isV0MOnVsOfPileup = 1;
//        }
//        else{
//            isV0MOnVsOfPileup = 0;
//        }
//    }
//    else if(isV0MOnVsOfPileup == 3){
//        isV0MOnVsOfPileup = 0;
//    }
//    }
//    
//    double fSPDOnVsOfA = fTriggerAnalysis->GetSPDOnVsOfA();
//    double fSPDOnVsOfB = fTriggerAnalysis->GetSPDOnVsOfB();
//    
//    int isSPDOnVsOfPileup = 0;
//    {
//      AliVMultiplicity* mult = aod->GetMultiplicity();
//      if (!mult) {
//        isSPDOnVsOfPileup = 3;
//      }
//      TBits onMap = mult->GetFastOrFiredChips();
//      TBits ofMap = mult->GetFiredChipMap();
//      Int_t on = onMap.CountBits(0);
//      Int_t of = ofMap.CountBits(0);
//    if(isSPDOnVsOfPileup == 0){
//        if(on < fSPDOnVsOfA + fSPDOnVsOfB*of){
//            isSPDOnVsOfPileup = 1;
//        }
//        else{
//            isSPDOnVsOfPileup = 0;
//        }
//    }
//    else if(isSPDOnVsOfPileup == 3){
//        isSPDOnVsOfPileup = 0;
//    }
//    }
//    
//    double fVIRBBAflags = fTriggerAnalysis->GetVIRBBAflags();
//    double fVIRBBCflags = fTriggerAnalysis->GetVIRBBCflags();
//    double fVIRBGAflags = fTriggerAnalysis->GetVIRBGAflags();
//    double fVIRBGCflags = fTriggerAnalysis->GetVIRBGCflags();
//    double fNBCsPast = fTriggerAnalysis->GetNBCsPast();
//    double fNBCsFuture = fTriggerAnalysis->GetNBCsFuture();
//    int isV0PFPileup = 0;
//    {
//      AliVVZERO* vzero = aod->GetVZEROData();
//      if (!vzero) {
//        isV0PFPileup = 3;
//      }
//      if (!vzero->TestBit(AliVVZERO::kPastFutureFlagsFilled)){
//        isV0PFPileup = 3;
//      }
//    if(isV0PFPileup == 0){
//      Bool_t vir[21] = {0};
//      UChar_t bcMod4 = aod->GetBunchCrossNumber()%4;
//
//      for (Int_t bc=0;bc<=20;bc++) {
//        UChar_t nBBA=0;
//        UChar_t nBBC=0;
//        UChar_t nBGA=0;
//        UChar_t nBGC=0;
//        if (fVIRBBAflags<33) for (Int_t i=0;i<32;i++) nBBA+=vzero->GetPFBBFlag(i+32,bc);
//        if (fVIRBBCflags<33) for (Int_t i=0;i<32;i++) nBBC+=vzero->GetPFBBFlag(i   ,bc);
//        if (fVIRBGAflags<33) for (Int_t i=0;i<32;i++) nBGA+=vzero->GetPFBGFlag(i+32,bc);
//        if (fVIRBGCflags<33) for (Int_t i=0;i<32;i++) nBGC+=vzero->GetPFBGFlag(i   ,bc);
//        vir[bc] |= nBBA>=fVIRBBAflags;
//        vir[bc] |= nBBC>=fVIRBBCflags;
//        vir[bc] |= nBGA>=fVIRBGAflags;
//        vir[bc] |= nBGC>=fVIRBGCflags;
//      }
//      
//      // clock index is counting from future to past
//      Int_t bcMin = 10 - fNBCsFuture + bcMod4;
//      Int_t bcMax = 10 + fNBCsPast   + bcMod4;
//      for (Int_t bc=bcMin;bc<=bcMax;bc++) {
//        if (bc==10) continue; // skip current bc
//        if (bc < 0) continue;
//        if (bc >20) continue;
//        if (vir[bc]){
//            isV0PFPileup = 1;
//        }
//      }
//    }
//    else if(isV0PFPileup == 3){
//        isV0PFPileup = 0;
//    }
//    }
//    
//    if(isV0MOnVsOfPileup == 1 || isSPDOnVsOfPileup == 1 || isV0PFPileup == 1){
//        fIsPileupReimplemented = 1;
//    }
//    else{
//        fIsPileupReimplemented = 0;
//    }
//    
//    
//    delete triggerAnalysisOADB;
    
  //  for(Int_t itrack = 0; itrack<(aod->GetNumberOfTracks()) ; itrack ++){
  //    track = (AliAODTrack*)aod->GetTrack(itrack);
  //    if(!track) {
  //      Printf("track not available\n");
  //      continue;
  //    }
  //
  //    if(!IsSelected(track) || !track->IsMuonTrack()) continue;
  //
  //    if (track->Charge()<0) {
  //      nMuonNeg++;
  //    } else if (track->Charge()>0) {
  //      nMuonPos++;
  //    } else {
  //      Printf("What the heck?!\n");
  //    }
  //
  //    real_track = AddTrack();
  //
  //    //    Printf("adding the next muon track: %d -- to the tree\n",
  //    itrack);
  //
  //    real_track->fTriggerData = track->GetMatchTrigger();
  //    //        real_track->fTrackerData = track->ContainTrackerData() ? 1 :
  //    0; real_track->fCharge = track->Charge(); real_track->fPx =
  //    track->Px(); real_track->fPy = track->Py(); real_track->fPz =
  //    track->Pz(); real_track->fPt = track->Pt(); real_track->fP =
  //    track->P(); real_track->fEta = track->Eta(); real_track->fPhi =
  //    track->Phi(); real_track->fRabs = track->GetRAtAbsorberEnd();
  //    real_track->fChi2 = track->GetChi2MatchTrigger();
  //    //        real_track->fNclusters = track->GetNClusters();
  //    real_track->fSelMask = trackCuts->GetSelectionMask(track);
  //    real_track->fLabel = track->GetLabel();
  //  } // track loop

  //	if (nMuonNeg>1 && nMuonPos >1) {

  AliAODTrack *track1;
  AliAODTrack *track2;
  DimuonLight *dimuon;
  CorrelationLight *correlation;
  TLorentzVector muon1L, muon2L, dimuonL;

  TLorentzVector muon1LRotNeg, muon2LRotNeg, dimuonLRotNeg;
  TLorentzVector muon1LRotPos, muon2LRotPos, dimuonLRotPos;

  // dimuon case
  for (Int_t itrack1 = 0; itrack1 < (aod->GetNumberOfTracks()); itrack1++) {

    track1 = (AliAODTrack *)aod->GetTrack(itrack1);

    if (!track1) {
      Printf("track 1 not available\n");
      continue;
    }

    if (!track1->IsMuonTrack())
      continue;
    // first muon lorentz vector
//      if ( ! trackCuts->IsSelected(track1) )
//      {
//          continue;}
    muon1L.SetPxPyPzE(track1->Px(), track1->Py(), track1->Pz(), track1->E());

    // loop over second muon, and construct jpsi candidates
    for (Int_t itrack2 = itrack1 + 1; itrack2 < aod->GetNumberOfTracks();
         itrack2++) {

          // second muon lorentz vector
          track2 = (AliAODTrack *)aod->GetTrack(itrack2);
          if (!track2) {
            Printf("track 2 not available\n");
            continue;
          }

        if (!track2->IsMuonTrack()){
            Printf("track 2 not muon\n");
            continue;
        }
        
//        if ( ! trackCuts->IsSelected(track2) ){
//        continue;}

          muon2L.SetPxPyPzE(track2->Px(), track2->Py(), track2->Pz(), track2->E());
        if (!trackCuts->IsSelected(track1) || !trackCuts->IsSelected(track2)){
            Printf("One of 2 tracks is not selected\n");
            newlyremoved++;
            continue;
        }
        // find the J/psi candidate's lorentz vector
          dimuonL = muon1L + muon2L;
          if (dimuonL.M() < 1.0)
            continue;
        
        if(kFALSE){//<<
          dimuon = AddDimuon();
          dimuon->fCharge = track1->Charge() + track2->Charge();
          dimuon->fInvMass = dimuonL.M();
          dimuon->fPt = dimuonL.Pt();
          dimuon->fEta = dimuonL.Eta();
          dimuon->fY = dimuonL.Rapidity();
          dimuon->fPhi = TMath::Pi() + TMath::ATan2(-dimuonL.Py(), -dimuonL.Px()); // dimuonL.Phi();
          dimuon->fLabel = track1->GetLabel() * track2->GetLabel();
        }
        
        printf("trklets %p has %d \n", trklets, trklets->GetNumberOfTracklets());
        if (!trklets)
          AliFatal("AliAODTracklets not found");
        trklets->Print("t");
        for (Int_t itrklt = 0; itrklt < trklets->GetNumberOfTracklets(); itrklt++) {
//          if (!IsSelectedTracklet(trklets, itrklt))
//            continue;
            if(kFALSE){//<<
           correlation = AddCorrelation();
           correlation->fDeltaPhi = (Float_t)trklets->GetPhi(itrklt) - dimuon->fPhi;
           correlation->fDeltaEta = (Float_t)-TMath::Log(TMath::Tan(trklets->GetTheta(itrklt) / 2)) - dimuonL.Eta();
           correlation->fLabel = track1->GetLabel() * track2->GetLabel(); //same label as dimuon
            
            //Save info on the Dimuon (only one branch could suffice for analysis)
            correlation->fDimuonInvMass = dimuonL.M();
            correlation->fDimuonCharge = track1->Charge() + track2->Charge();
            correlation->fDimuonPt = dimuonL.Pt();
            correlation->fDimuonEta = dimuonL.Eta();
            correlation->fDimuonY = dimuonL.Rapidity();
            correlation->fDimuonPhi = TMath::Pi() + TMath::ATan2(-dimuonL.Py(), -dimuonL.Px()); // dimuonL.Phi();
            correlation->fTrackletEta = (Float_t)-TMath::Log(TMath::Tan(trklets->GetTheta(itrklt) / 2));
            correlation->fTrackletDPhi = (Float_t)trklets->GetDeltaPhi(itrklt);
            }
        }

    } // second muon
  }   // first muon

  MCDimuonLight *mcdimuon;

  //  printf("bMC %s \n",bMC ? "kTRUE" : "kFALSE");
  if (bMC) {
    // Get MC particles
    TClonesArray *mcParticles = static_cast<TClonesArray *>(
        aod->FindListObject(AliAODMCParticle::StdBranchName()));

    // ##########################################
    // loop on MC particles
    // ##########################################
    Int_t nTracks = mcParticles->GetEntries();
    AliAODMCParticle *mcTrack = 0x0;

    for (Int_t itrack = 0; itrack < nTracks; itrack++) {
      //      printf("looping %i \n",itrack);

      mcTrack = dynamic_cast<AliAODMCParticle *>(mcParticles->At(itrack));

      //      printf("mcTrack %p \n",mcTrack);
      if (!mcTrack)
        continue;
      // TParticle *mcPart = stack->Particle(itrack);
      //      printf("mcTrack->GetMother() %i \n",mcTrack->GetMother());
      // only primaries
        if ((mcTrack->IsPhysicalPrimary()) && (mcTrack->Charge()!=0)) { //(TMath::Abs(mcTrack->PdgCode())==211 || TMath::Abs(mcTrack->PdgCode())==321 || TMath::Abs(mcTrack->PdgCode())==2212)
     // if ((mcTrack->GetMother() == -1)) {
        mcdimuon = AddMCDimuon();
        mcdimuon->fPDG = mcTrack->PdgCode();
        mcdimuon->fInvMass = mcTrack->M();
        mcdimuon->fPt = mcTrack->Pt();
        mcdimuon->fEta = mcTrack->Eta();
        mcdimuon->fY = mcTrack->Y();
        mcdimuon->fPhi = mcTrack->Phi();
        mcdimuon->fCharge = mcTrack->Charge();
        mcdimuon->fLabel = mcTrack->GetLabel();
        //        printf("fNMCDimuonLights %i out of %i
        //        mcParticles\n",fNMCDimuonLights,nTracks);
      }
    }
  }

  //	}
  TProcessID::SetObjectCount(ObjectNumber);
}

//______________________________________________________________________________
Bool_t MyEventLight::IsSelected(AliAODTrack *track) {

  /*
   Bool_t isFake = ! (((AliAODTrack*)track)->ContainTriggerData() &&
   ((AliAODTrack*)track)->ContainTrackerData()) ;
   // protection: remove fake track
   if(isFake)  return kFALSE;
   */

  // cut on eta
  if (TMath::Abs(track->Eta()) > 4.0 || TMath::Abs(track->Eta()) < 2.5)
    return kFALSE;
//  // cut on muon Pt
// if (track->Pt() < 0.5)   // >>>>>>
//      return kFALSE;
  // cut on R at the end plane of the absorber: in cm
  if (TMath::Abs(track->GetRAtAbsorberEnd()) > 89.5 ||
      TMath::Abs(track->GetRAtAbsorberEnd()) < 17.6)
    return kFALSE;
  // cut on match trigger
  if (track->GetMatchTrigger() < 2)
    return kFALSE;

  return kTRUE;
}

//______________________________________________________________________________
Bool_t MyEventLight::IsSelectedTracklet(AliAODTracklets *trklets,
                                        Int_t itrklt) {
  printf("tracklet %d, theta = %f Tan: %f Log: %f Eta: %f Phi: %f\n",
         trklets->GetTheta(itrklt), TMath::Tan(trklets->GetTheta(itrklt)),
         TMath::Log(TMath::Tan(trklets->GetTheta(itrklt) / 2)),
         -TMath::Log(TMath::Tan(trklets->GetTheta(itrklt) / 2)),
         trklets->GetPhi(itrklt));
  Float_t eta = -TMath::Log(TMath::Tan(trklets->GetTheta(itrklt) / 2));
  printf("traclet %d, eta = %f\n");
  // cut on eta
  if (eta > 1.0 || eta < -1.0)
    return kFALSE;

  return kTRUE;
}

//______________________________________________________________________________
TrackletLight *MyEventLight::AddTracklet() {
  // Add a new track to the list of tracks for this event.

  TClonesArray &tracklets = *fTracklets;
  TrackletLight *tracklet = new (tracklets[fNTracklets++]) TrackletLight();
  return tracklet;
}

//______________________________________________________________________________
DimuonLight *MyEventLight::AddDimuon() {
  // Add a new dimuon to the list of dimuons for this event.

  TClonesArray &dimuons = *fDimuons;
  DimuonLight *dimuon = new (dimuons[fNDimuons++]) DimuonLight();
  return dimuon;
}

//______________________________________________________________________________
MCDimuonLight *MyEventLight::AddMCDimuon() {
  // Add a new dimuon to the list of dimuons for this event.

  TClonesArray &mcdimuons = *fMCDimuons;
  MCDimuonLight *mcdimuon = new (mcdimuons[fNMCDimuons++]) MCDimuonLight();
  return mcdimuon;
}

//______________________________________________________________________________
CorrelationLight *MyEventLight::AddCorrelation() {
  // Add a new correlation to the list of correlations for this event.

  TClonesArray &correlations = *fCorrelations;
  CorrelationLight *correlation = new (correlations[fNCorrelations++]) CorrelationLight();
  return correlation;
}

//_______________________________________________________________________
void MyEventLight::Clear(Option_t *option) {
  //
  // clear object
  //
  fTracklets->Clear(option);
  fNTracklets = 0;
  fDimuons->Clear(option);
  fNDimuons = 0;
  fMCDimuons->Clear(option);
  fNMCDimuons = 0;
  fCorrelations->Clear(option);
  fNCorrelations = 0;
  if (fTracklets)
    fTracklets->Delete();
  if (fDimuons)
    fDimuons->Delete();
  if (fMCDimuons)
    fMCDimuons->Delete();
  if (fCorrelations)
  fCorrelations->Delete();
}

//_______________________________________________________________________
void MyEventLight::Reset(Option_t *) {

  if (fTracklets)
    fTracklets->Delete();
  if (fDimuons)
    fDimuons->Delete();
  if (fMCDimuons)
    fMCDimuons->Delete();
  if (fCorrelations)
    fCorrelations->Delete();
  //  delete fgTracks; fgTracks = 0;
  //  delete fgDimuonLights; fgDimuonLights = 0;
  //	delete fgMCDimuonLights; fgMCDimuonLights = 0;
}
