AliAnalysisTaskMLPIDResponse *AddTaskMLPIDResponse(TChain* inputChain, UInt_t collisionCandidates=AliVEvent::kINT7, Bool_t isMC=kFALSE)
{
// Macro to connect a centrality selection task to an existing analysis manager.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMLPIDResponse", "No analysis manager to connect to.");
    return 0x0;
  }

  // standard with task
  printf("========================================================================================\n");
  printf("MLPIDResponse: Initialising AliAnalysisMLTaskPIDResponse\n");

  AliAnalysisTaskMLPIDResponse *pidTask = new AliAnalysisTaskMLPIDResponse("MLPIDResponseTask");

  pidTask->setInputTChain(inputChain);
  pidTask->saveCollisionCandidates(collisionCandidates);
  pidTask->SelectCollisionCandidates(collisionCandidates);
  pidTask->setIsMC(isMC);
  
  mgr->AddTask(pidTask);
  mgr->ConnectInput(pidTask, 0, mgr->GetCommonInputContainer());

  return pidTask;
}
