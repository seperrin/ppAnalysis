#ifdef __INTELLISENSE__
#include "../ReadIntegers.C"
#endif

#if defined(__CLING__)
#include "../ReadIntegers.C"
// Tell  ROOT where to find AliRoot headers
R__ADD_INCLUDE_PATH($ALICE_ROOT)
//#include <ANALYSIS/macros/train/AddESDHandler.C>

// Tell ROOT where to find AliPhysics headers
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
//#include <OADB/macros/AddTaskPhysicsSelection.C>
#elif !(defined(__CLING__) || defined(__CINT__)) || defined(__ROOTCLING__) ||  \
    defined(__ROOTCINT__)
#include <iostream>

#include <TInterpreter.h>
#include <TROOT.h>

#include "AliAODInputHandler.h"
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliMuonTrackCuts.h"
#endif

AliAnalysisGrid *CreateAlienHandler(Int_t runno = 0, Bool_t bSig = kFALSE,
                                    Bool_t iMUONReFit = kTRUE) {

  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init
  // then
  //   source //gclient_env_$UID in the current shell.
  //   if (!AliAnalysisGrid::CreateToken()) return NULL;
  AliAnalysisAlien *plugin = new AliAnalysisAlien();

  // Overwrite all generated files, datasets and output results from a previous
  // session
  plugin->SetOverwriteMode();
  // Set the run mode (can be "full", "test", "offline", "submit" or
  // "terminate") plugin->SetRunMode("submit");
  // plugin->SetRunMode("offline");
// plugin->SetRunMode("test");
  plugin->SetRunMode("full");
 // plugin->SetRunMode("terminate");
  plugin->SetMergeViaJDL(kFALSE); // send jobs to merge the output
                                 // Set versions of used packages

  plugin->SetAPIVersion("V1.1x");
  plugin->SetAliPhysicsVersion("vAN-20210607-1"); // vAN-20210417-1 >>>

#if !(defined(__CLING__)) || defined(__CINT__)
  gROOT->LoadMacro("../ReadIntegers.C");
  //  gSystem->Load("../ReadIntegers_C.so");
#endif
  std::vector<int> fRunList; // input run list

  ReadIntegers("../../RunLists/RunList_Group1_LHC16hMCB.txt", fRunList, kTRUE);
  // Declare input data to be processed.

  // Method 1: Create automatically XML collections using alien 'find' command.
  // Define production directory LFN
  //plugin->SetGridDataDir("/alice/data/2016/LHC16h");
    plugin->SetGridDataDir("/alice/sim/2017/LHC17f5");
  // Set data search pattern
  //  plugin->SetDataPattern("muon_calo_pass3/AOD225/*AliAOD.root");
  // plugin->SetDataPattern("muon_calo_pass1/*AliAOD.Muons.root");
  plugin->SetDataPattern("AOD235/*AliAOD.root");
 //   plugin->SetDataPattern("muon_calo_pass2/AOD/*AliAOD.Muons.root");
  //	plugin->SetDataPattern("AliMUONESDs.root");  // simulated, tags
  // not used 	if (!bSig) plugin->SetDataPattern(Form("p40/*ESDs.root"));
  // // real data check reco pass and data base directory 	else
  //	  plugin->SetDataPattern(Form("p40/*ESDsSignal.root")); // real
  // data check reco pass and data base directory

  if (runno) {
    plugin->SetRunPrefix(""); // real data 000
    //   plugin->SetDataPattern("*tag.root");  // Use ESD tags (same
    //   applies for AOD's)
    // ...then add run numbers to be considered
    //   plugin->AddRunNumber(125020);    // simulated
    plugin->AddRunNumber(runno); // real data
  } else {
    plugin->SetRunPrefix("");//000

    for (std::vector<int>::const_iterator it = fRunList.begin();
         it != fRunList.end(); ++it) {
      int runnumber = *it;
      std::cout << "Run " << runnumber << std::endl;

      plugin->AddRunNumber(Form("%d", runnumber));
    }
  }
  plugin->SetNrunsPerMaster(1);
  plugin->SetOutputToRunNo();
  // Method 2: Declare existing data files (raw collections, xml collections,
  // root file) If no path mentioned data is supposed to be in the work
  // directory (see SetGridWorkingDir()) XML collections added via this method
  // can be combined with the first method if the content is compatible (using
  // or not tags)
  //   plugin->AddDataFile("tag.xml");
  //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");

  // Define alien work directory where all files will be copied. Relative to
  // alien $HOME.
  plugin->SetGridWorkingDir("ppwrk/Manu16hMCINTnospe_PhySelTRUE_SelNtkl1_Vertexer_CorrelVeryLoose");
  //	if (!bSig)
  //		plugin->SetGridWorkingDir("data/2011/LHC11h/pass2embupsi/emb");
  //	else
  //		plugin->SetGridWorkingDir("data/2011/LHC11h/pass2embupsi/sig");
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir(
      "output"); // In this case will be $HOME/work/.../output

  plugin->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

  // Declare the analysis source files names separated by blancs. To be compiled
  // runtime using ACLiC on the worker nodes.
  plugin->SetAnalysisSource("AliAnalysisTaskMyMuonTree_AOD.cxx");

  // Declare all libraries (other than the default ones for the framework. These
  // will be loaded by the generated analysis macro. Add all extra files (task
  // .cxx/.h) here.
  plugin->SetAdditionalLibs(
      "AliAnalysisTaskMyMuonTree_AOD.h AliAnalysisTaskMyMuonTree_AOD.cxx");

  // No need for output file names. Procedure is automatic.
  //   plugin->SetOutputFiles("MuonTask.137161.root");
  //   plugin->SetDefaultOutputs();
  // No need define the files to be archived. Note that this is handled
  // automatically by the plugin.
  //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
  // Set a name for the generated analysis macro (default MyAnalysis.C) Make
  // this unique !
  plugin->SetAnalysisMacro("AMMT.C");
  // Optionally set maximum number of input files/subjob (default 100, put 0 to
  // ignore). The optimum for an analysis is correlated with the run time -
  // count few hours TTL per job, not minutes !
  //  plugin->SetSplitMaxInputFileNumber(300);
  plugin->SetSplitMaxInputFileNumber(30);
  // Optionally modify the executable name (default analysis.sh)
  plugin->SetExecutable("AMMT.sh");
  // Optionally set number of failed jobs that will trigger killing waiting
  // sub-jobs.
  //   plugin->SetMaxInitFailed(5);
  // Set number of files to be read in test mode
  plugin->SetNtestFiles(1);
  // Optionally resubmit threshold.
  plugin->SetMasterResubmitThreshold(99);
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(10800);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName("AMMT.jdl");
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);
  plugin->SetNrunsPerMaster(200);
  // Optionally modify split mode (default 'se')
  plugin->SetSplitMode("se");
  // Optionally Checking possibility to copy files to your AliEn home
  // directory...
  //   plugin->SetCheckCopy(kFALSE);
  plugin->SetKeepLogs(kFALSE);
  plugin->SetMaxMergeStages(1);
  return plugin;
}
