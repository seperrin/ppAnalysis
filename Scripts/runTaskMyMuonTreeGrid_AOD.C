#ifdef __INTELLISENSE__
#include "./CreateAlienHandler_AOD.C"
#endif

#if defined(__CLING__)
#include "./CreateAlienHandler_AOD.C"
// Tell  ROOT where to find AliRoot headers
R__ADD_INCLUDE_PATH($ALICE_ROOT)
//#include <ANALYSIS/macros/train/AddESDHandler.C>

// Tell ROOT where to find AliPhysics headers
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
//#include <OADB/macros/AddTaskPhysicsSelection.C>
#elif !(defined(__CLING__) || defined(__CINT__)) || defined(__ROOTCLING__) ||  \
    defined(__ROOTCINT__)
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliMultSelectionTask.h"
#include "AliMuonTrackCuts.h"
#include "AliPhysicsSelectionTask.h"
#include <TInterpreter.h>
#include <TMacro.h>
#include <TROOT.h>
#endif

#include "./AliAnalysisTaskMyMuonTree_AOD.h"

void runTaskMyMuonTreeGrid_AOD(Int_t runno = 0, Bool_t bSig = kFALSE) {
  // Use AliRoot includes to compile our task
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

  Bool_t isMC = kFALSE;
  Bool_t isMMC = kFALSE;

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("muonAnalysis");
// Create and configure the alien handler plugin
#if !(defined(__CLING__)) || defined(__CINT__)
  gROOT->LoadMacro("CreateAlienHandler_AOD.C"); //>>>
#endif
  AliAnalysisGrid *alienHandler = CreateAlienHandler(runno, bSig);
  if (!alienHandler)
    return;
  // Connect plug-in to the analysis manager
  mgr->SetGridHandler(alienHandler);

  AliAODInputHandler *aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);

  // physics selection task
//  TMacro physseladd(gSystem->ExpandPathName(
//      "$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"));
//  AliPhysicsSelectionTask *physseltask =
//      reinterpret_cast<AliPhysicsSelectionTask *>(physseladd.Exec());
   AliPhysicsSelectionTask *physseltask = reinterpret_cast<AliPhysicsSelectionTask *>(gInterpreter->ProcessLine(Form(".x %s(%d,%d)", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"),kFALSE,kTRUE))); //kFALSE, kTRUE (not MC, remove PU)

  // centrality task
  TMacro multseladd(gSystem->ExpandPathName(
      "$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));
  AliMultSelectionTask *taskCentrality =
      reinterpret_cast<AliMultSelectionTask *>(multseladd.Exec()); // user mode:
  taskCentrality->SetSelectedTriggerClass(AliVEvent::kINT7); // default but //kAny pour CMUL //kINT7 pour CINT7 >>>
  //  configurable taskCentrality->SetUseDefaultCalib(kTRUE);
//taskCentrality->SetAlternateOADBFullManualBypass("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/data/OADB-LHC18q.root");

AliMuonTrackCuts *MuonTrackCuts =
      new AliMuonTrackCuts("StandardMuonTrackCuts", "StandardMuonTrackCuts"); //
  MuonTrackCuts->SetAllowDefaultParams(kTRUE);
  MuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuPdca);
    MuonTrackCuts->Print();
    AliOADBMuonTrackCutsParam MuonTrackCutsParam = NULL;
    MuonTrackCutsParam = MuonTrackCuts->GetMuonTrackCutsParam();
    MuonTrackCutsParam.Print();
        
        
  // MuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kEta |
  // AliMuonTrackCuts::kThetaAbs | AliMuonTrackCuts::kPdca |
  // AliMuonTrackCuts::kMatchApt); // Optionally set the cuts you want to apply
  // (if different from the default ones shown here)

  AliAnalysisTaskMyMuonTree_AOD *muonTask = 0x0;
  // muon task
#if !defined(__CINT__) || defined(__CLING__)
  gInterpreter->LoadMacro("AliAnalysisTaskMyMuonTree_AOD.cxx++g"); //>>>
  muonTask = reinterpret_cast<AliAnalysisTaskMyMuonTree_AOD *>(
      gInterpreter->ExecuteMacro("AddTaskMyMuonTreeGrid_AOD.C(\"muonTask\")")); //>>>
  muonTask->SetMC(isMMC);
  muonTask->SetMuonTrackCuts(MuonTrackCuts);
#else
  gROOT->LoadMacro("AliAnalysisTaskMyMuonTree_AOD.cxx++g"); //>>>
  muonTask =
      new AliAnalysisTaskMyMuonTree_AOD("muonTask", MuonTrackCuts, isMMC);
//  gSystem->Load("AliAnalysisTaskMyMuonTree_AOD_cxx.so");
#endif

  mgr->SetDebugLevel(AliLog::kInfo);
  // mgr->SetDebugLevel(AliLog::kFatal);
  // mgr->SetDebugLevel(10);

  if (!mgr->InitAnalysis())
    return;
  // mgr->PrintStatus();
  mgr->StartAnalysis("grid");
};
