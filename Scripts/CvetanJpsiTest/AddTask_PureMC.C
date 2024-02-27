void AddTask_PureMC(Bool_t prompt) {

	// ================== GetAnalysisManager ===============================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTask_GammaPureMC", "No analysis manager found.");
		return ;
	}

	// ================== GetInputEventHandler =============================
	AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

	
	//================================================
	//========= Add task to the ANALYSIS manager =====
	//================================================
	//            find input container

	AliAnalysisTaskPureMC *task = new AliAnalysisTaskPureMC("PureMC");
	task->SetPrompt(prompt);
	
	//connect containers
	AliAnalysisDataContainer *coutput =
	  mgr->CreateContainer(Form("%s",prompt?"Prompt":"FromB"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:PureMC",AliAnalysisManager::GetCommonFileName()));
		
	mgr->AddTask(task);
	mgr->ConnectInput(task,0,cinput);
	mgr->ConnectOutput(task,1,coutput);
	
	return;
	
}
