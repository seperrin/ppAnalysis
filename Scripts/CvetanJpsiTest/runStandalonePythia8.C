// runmcgenThermalModel.C
// =====================
// This macro runs jet analysis on pythia events with AliMCGen framework
//
// Author: M. Verweij
 
const Bool_t   saveManager         = kFALSE;//kFALSE;//

void runStandalonePythia8(Bool_t prompt = kTRUE, Long64_t nEvents = 10000)
{
//  gSystem->AddIncludePath("-I. -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
    
    // Use AliRoot includes to compile our task
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

	AliAnalysisManager *mgr = new AliAnalysisManager("MCGenPythia8");

	// Handlers
	AliDummyHandler *dumH = new AliDummyHandler();
	//mgr->SetInputEventHandler(dumH);

	// AliDummyHandler *dumHM = static_cast<AliDummyHandler*>mgr->GetInputEventHandler();
	// if(dumHM) Printf("Found dummy handler");
	// if(!dumHM->GetEvent()) Printf("handler has no event");
	AliESDEvent *esdE = new AliESDEvent();
	esdE->CreateStdContent();
	AliESDVertex *vtx = new AliESDVertex(0.,0.,100);
	vtx->SetName("VertexTracks");
	vtx->SetTitle("VertexTracks");
	esdE->SetPrimaryVertexTracks(vtx);
	if(esdE->GetPrimaryVertex()) Printf("vtx set");
	dumH->SetEvent(esdE);
	// if(dumHM->GetEvent()) Printf("handler has event");
	mgr->SetInputEventHandler(dumH);

	AliMCGenHandler* mcInputHandler = new AliMCGenHandler();  
	mgr->SetMCtruthEventHandler(mcInputHandler);

	// Generator
	// gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddMCGenHijing.C");
	// AliGenerator* gener = AddMCGenHijing();
	// mcInputHandler->SetGenerator(gener);
	// mcInputHandler->SetSeedMode(2);
	gROOT->LoadMacro("Add_MCGenPythia8_TuneX.C");
	AliGenerator* gener = Add_MCGenPythia8_TuneX(8000., 5, kTRUE, 1, prompt?2:3, 0., 0.);

	mcInputHandler->SetGenerator(gener);
	mcInputHandler->SetSeedMode(2);

  // Analysis Tasks
	gROOT->LoadMacro("AliAnalysisTaskPureMC.cxx+g");
	gROOT->LoadMacro("AddTask_PureMC.C");
	AliAnalysisTask *taskA = AddTask_PureMC(prompt);		

	// Run analysis
	if(saveManager) {
		TFile *fileMgr = new TFile("AnalysisManager.root","RECREATE");
		mgr->Write();

		fileMgr->Write();
		fileMgr->Close();

		runUseSavedManager("AnalysisManager.root",nEvents);
	}
	else {
		mgr->InitAnalysis();
		mgr->PrintStatus();
		mgr->EventLoop(nEvents);
	}
	(AliPythia8::Instance())->PrintStatistics();

}

void runUseSavedManager(TString fileMgr = "AnalysisManager.root", Int_t nEvents = 1000, Bool_t useGrid = kFALSE) {

  Printf("run saved");

  LoadLibs();

  TFile *f = new TFile(fileMgr.Data());
  AliAnalysisManager *mgr = dynamic_cast<AliAnalysisManager*>f->Get("MCGenPythia8");
  if(!mgr) {
    Printf("Did not find manager");
    return;
  }

  // AliDummyHandler *dumH = static_cast<AliDummyHandler*>mgr->GetInputEventHandler();
  // AliESDEvent *esdE = new AliESDEvent();
  // esdE->CreateStdContent();
  // AliESDVertex *vtx = new AliESDVertex(0.,0.,100);
  // vtx->SetName("VertexTracks");
  // vtx->SetTitle("VertexTracks");
  // esdE->SetPrimaryVertexTracks(vtx);
  // if(esdE->GetPrimaryVertex()) Printf("vtx set");
  // dumH->SetEvent(esdE);

  mgr->SetDebugLevel(11);
  mgr->InitAnalysis();
  mgr->PrintStatus();
  mgr->EventLoop(nEvents);
}


