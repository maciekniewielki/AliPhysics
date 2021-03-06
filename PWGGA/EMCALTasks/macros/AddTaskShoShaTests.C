// $Id$

AliAnalysisTaskShoShaTests *AddTaskShoShaTests(
  Double_t trgThresh=4.8, 
  Double_t exoticCut=0.97
)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskShoShaTests", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Create the task and configure it.
  //===========================================================================
  AliAnalysisTaskShoShaTests* ana = new  AliAnalysisTaskShoShaTests("");

  Bool_t isMC = (mgr->GetMCtruthEventHandler() != NULL);

  if(!isMC)
    ana->SelectCollisionCandidates( AliVEvent::kEMC1 | AliVEvent::kMB | AliVEvent::kEMC7 | AliVEvent::kINT7);
  

  ana->SetTrigThresh(trgThresh);
  ana->SetExotCut(exoticCut);
  
  mgr->AddTask(ana);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histosShoShaTests", TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s", AliAnalysisManager::GetCommonFileName()));
  
  mgr->ConnectInput  (ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (ana, 1, coutput1 );
   
  return ana;
}
