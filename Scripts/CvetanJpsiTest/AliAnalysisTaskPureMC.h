#ifndef ALIANLYSISTASKPUREMC_cxx
#define ALIANLYSISTASKPUREMC_cxx

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

class AliAnalysisTaskPureMC : public AliAnalysisTaskSE {
  public:

    AliAnalysisTaskPureMC();
    AliAnalysisTaskPureMC(const char *name);
    virtual ~AliAnalysisTaskPureMC();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *);
    virtual void   Terminate(const Option_t*);

    // MC functions
    void SetIsMC(Int_t isMC){fIsMC=isMC;}
    void SetPrompt(Bool_t prompt){fPrompt = prompt;}
    void ProcessMCParticles();
    void CorrelateWithHadrons(Int_t iJpsi, Double_t pT, Double_t y, Double_t phi);
    
  protected:
    AliVEvent*            fInputEvent;                // current event
    AliMCEvent*           fMCEvent;                   // corresponding MC event
    AliStack*             fMCStack;                   // stack belonging to MC event
    TList*                fOutputContainer;           // Output container
    
    TH1F*                 fHistNEvents;               //! number of events histo
    TH1D*                 fHistXSection;              //! xSection
    TH1F*                 fHistPtJpsi;                //! 
    TH1F*                 fHistPtJpsiSel;             //! 
    TH1F*                 fHistYJpsi;                 //!
    TH2F*                 fHistPhiVsPt;               //!
    TH2F*                 fHistPhiVsPt05;             //!

    Int_t                 fIsMC;                      // MC flag
    Bool_t                fPrompt;                    // prompt J/psi or B-> J/psi
    
  private:
    AliAnalysisTaskPureMC(const AliAnalysisTaskPureMC&); // Prevent copy-construction
    AliAnalysisTaskPureMC &operator=(const AliAnalysisTaskPureMC&); // Prevent assignment

    ClassDef(AliAnalysisTaskPureMC, 1);
};

#endif
