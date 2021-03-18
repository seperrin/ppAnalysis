//
//  AliAnalysisTaskMyMuonTree_AOD.h
//
//
//  Created by Marc Arène on 7/31/13.
//
//

#ifndef _AliAnalysisTaskMyMuonTree_AOD_h
#define _AliAnalysisTaskMyMuonTree_AOD_h

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class TH1I;
class TTree;
class TList;
class TObjArray;
class TClonesArray;
class AliAnalysisManager;
class AliAnalysisTaskSE;
class AliAODEvent;
class AliCounterCollection;
class AliMuonTrackCuts;

#ifndef MYEVENT_H
#define MYEVENT_H

////////////////////////////////
// build an event   //
///////////////////////////////

class AliAODEvent;
class AliAODTrack;
class AliAODTracklets;
class TClonesArray;
class TLorentzVector;
class TrackletLight;
class DimuonLight;
class MCDimuonLight;
class CorrelationLight;

//___________________________________________________________
class MyEventLight : public TObject {
public:
  Int_t fNTracklets;
  Int_t fNDimuons;
  Int_t fNMCDimuons;
  Int_t fNCorrelations;
  Int_t fRunNumber;

  Int_t fNPileupVtx;
    Bool_t fIsPileupFromSPDMultBins;

  Float_t fVertexX; //! vertex X
  Float_t fVertexY; //! vertex Y
  Float_t fVertexZ; // vertex Z
    Float_t fVertexXSPD; // vertex X
    Float_t fVertexYSPD; // vertex Y
    Float_t fVertexZSPD; // vertex Z
  Float_t fVertexSigmaZ; // vertex Z
    Float_t fSPDVertexSigmaZ;
  Int_t fVertexNC;
    Int_t fSPDVertexNC;
//  Double_t fVertexResZ; //Resolution on Z

  Int_t fTriggerCINT7; //! trigger V0C and V0A
  Int_t fTriggerCMSL7; //! trigger V0C and V0A
  Int_t fTriggerCMSH7; //! trigger V0C and V0A
  Int_t fTriggerCMUL7; //! trigger diMuon Unlike sign Low Pt
  Int_t fTriggerCMLL7; //! trigger diMuon Like sign Low Pt
  Int_t fTriggerCINT7CENT; //! CINT7CENT

  Int_t fClusterMUFAST;    //!
  Int_t fClusterCENT;      //!
  Int_t fClusterCENTNOPMD; //!

  TString fFiredTriggerClass; //! fired trigger classes

  Int_t fPassPhysicsSelection; //! event selection criteria
  // 0: events not selected by physics selection;
  // 1: events selected by physics selection;
  Int_t fPassTriggerSelection; //! event selection criteria
  // 0: events not passing trigger selection
  // 1: events passing trigger selection

  Float_t fCentralityV0M;
//  Float_t fCentralityTKL;
//    Float_t fCentralityCL0;
//    Float_t fCentralityCL1;
    Float_t fCentralitySPDTracklets;
    Float_t fSPDTrackletsValue;
    Float_t fSPDClustersValue;
//    Float_t fCentralitySPDClusters;
//  Float_t fCentralityTKLvsV0M;

  TClonesArray *fTracklets; // array of all tracks
  TClonesArray *fDimuons;   // array of all dimuons
  TClonesArray *fMCDimuons; //-> array of all dimuons mc
                            //
                            //  static TClonesArray *fgTracks;
                            //  static TClonesArray *fgDimuons;
                            //	static TClonesArray *fgMCDimuons;
  TClonesArray *fCorrelations; // array of all correlations

public:
  MyEventLight();
  virtual ~MyEventLight();

  void Build(AliAODEvent *aod, Bool_t isSelected, AliMuonTrackCuts *cuts,
             Bool_t bMC = kFALSE);
  void Clear(Option_t *option = "");
  void Reset(Option_t *option = "");

  TrackletLight *AddTracklet();
  DimuonLight *AddDimuon();
  MCDimuonLight *AddMCDimuon();
  CorrelationLight *AddCorrelation();

  Bool_t IsSelected(AliAODTrack *track);
  Bool_t IsSelectedTracklet(AliAODTracklets *tracklet, Int_t itrklt);

  ClassDef(MyEventLight, 2)
};

#endif

//___________________________________________________________
class TrackletLight : public TObject {
public:
  Float_t fEta;  // eta
  Float_t fPhi;  // phi
  Float_t fDPhi; // phi
  Int_t fLabel;  // mc label

public:
  TrackletLight();
  virtual ~TrackletLight() {}

  ClassDef(TrackletLight, 1)
};

//___________________________________________________________
class DimuonLight : public TObject {
public:
  Short_t fCharge;  // charge of the muon
  Float_t fInvMass; // invariant mass of the muon pair
  Float_t fPt;      // pt component of the momentum
  Float_t fEta;     // eta
  Float_t fY;       // y
  Float_t fPhi;     // phi
  Int_t fLabel;     // mc label

public:
  DimuonLight();
  virtual ~DimuonLight() {}

  ClassDef(DimuonLight, 1)
};

//___________________________________________________________
class CorrelationLight : public TObject {
public:
  Float_t fDeltaPhi;  // DeltaPhi of the objects
  Float_t fDeltaEta; // DeltaEta of the objects
  Int_t fLabel;     // mc label
    Short_t fDimuonCharge;  // charge of the dimuon
    Float_t fDimuonInvMass; // invariant mass of the muon pair
    Float_t fDimuonPt;      // pt component of the momentum
    Float_t fDimuonEta;     // eta
    Float_t fDimuonY;       // y
    Float_t fDimuonPhi;     // phi
    Float_t fTrackletEta;  // eta
    Float_t fTrackletDPhi;  // DeltaPhi

public:
  CorrelationLight();
  virtual ~CorrelationLight() {}

  ClassDef(CorrelationLight, 1)
};

//___________________________________________________________
class MCDimuonLight : public TObject {
public:
  Int_t fPDG;       // PDG code
  Float_t fInvMass; // invariant mass of the muon pair
  Float_t fPt;      // pt component of the momentum
  Float_t fEta;     // eta
  Float_t fY;       // y
  Float_t fPhi;     // phi
  Int_t fLabel;     // mc label

public:
  MCDimuonLight();
  virtual ~MCDimuonLight() {}

  ClassDef(MCDimuonLight, 1)
};

//___________________________________________________________
class AliAnalysisTaskMyMuonTree_AOD : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskMyMuonTree_AOD();                 // default constructor
  AliAnalysisTaskMyMuonTree_AOD(const char *name); // non-default constructor
  AliAnalysisTaskMyMuonTree_AOD(const char *name, AliMuonTrackCuts *cuts,
                                Bool_t bMC = kFALSE); // non-default constructor
  AliAnalysisTaskMyMuonTree_AOD(
      const AliAnalysisTaskMyMuonTree_AOD &ref); // copy constructor
  AliAnalysisTaskMyMuonTree_AOD &
  operator=(const AliAnalysisTaskMyMuonTree_AOD &ref); // assignment constructor
  virtual ~AliAnalysisTaskMyMuonTree_AOD();            // destructor

  virtual void NotifyRun();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  void PostProcess(); // post analysis
  void SetMC(Bool_t bMC = kFALSE) { fMC = bMC; }
  void SetMuonTrackCuts(AliMuonTrackCuts *cuts) { fMuonTrackCuts = cuts; }

private:
  Bool_t fMC;
  TH1I *fHistEvents;                    //! event counter
  TTree *fMyMuonTree;                   //! muon tree
  TList *fOutputList;                   //! output list
  Bool_t fIsOutputTree;                 // switch to have tree output or not
  AliCounterCollection *fEventCounters; //! event counters
  MyEventLight *fEvent;
  AliMuonTrackCuts *fMuonTrackCuts;

  ClassDef(AliAnalysisTaskMyMuonTree_AOD, 4);
};
#endif
