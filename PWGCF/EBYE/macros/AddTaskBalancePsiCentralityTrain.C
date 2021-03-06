// now in options
//=============================================//
//const char* centralityEstimator = "V0M";
//const char* centralityEstimator = "CL1";
//const char* centralityEstimator = "TRK";
//=============================================//
//Bool_t gRunShuffling = kFALSE;
//Bool_t gRunShuffling = kTRUE;
//=============================================//

//PID config
Bool_t kUseNSigmaPID = kTRUE;
Double_t nSigmaMax = 3.0;
Bool_t kUseBayesianPID = kFALSE;
Double_t gMinAcceptedProbability = 0.7;

//_________________________________________________________//
AliAnalysisTaskBFPsi *AddTaskBalancePsiCentralityTrain(Double_t centrMin=0.,
						       Double_t centrMax=100.,
						       Bool_t gRunShuffling=kFALSE,
						       Bool_t gRunMixing=kTRUE,
						       Bool_t gRunMixingWithEventPlane=kFALSE,
						       TString centralityEstimator="V0M",
						       Double_t vertexZ=10.,
						       Double_t DCAxy=-1,
						       Double_t DCAz=-1,
						       Double_t ptMin=0.3,
						       Double_t ptMax=1.5,
						       Double_t etaMin=-0.8,
						       Double_t etaMax=0.8,
						       Double_t maxTPCchi2 = -1, 
						       Int_t minNClustersTPC = -1,
						       Bool_t kUsePID = kTRUE,
						       Bool_t bResonancesCut = kTRUE,
						       Bool_t bHBTcut = kTRUE,
						       Double_t HBTCutValue = 0.02,
						       Bool_t bConversionCut = kTRUE,
						       Double_t invMassForConversionCut = 0.04,
						       Bool_t bMomentumDifferenceCut = kTRUE,
						       Double_t fQCutMin = 0.0,
						       Int_t AODfilterBit = 128,
						       AliAnalysisTaskBFPsi::etriggerSel triggerSel = AliAnalysisTaskBFPsi::kINT7,
						       TString fileNameBase="AnalysisResults",
						       TString dirNameExtra="",
						       TString fArgEventClass="Centrality",
						       TString analysisTypeUser="AOD",
						       Bool_t bVertexBinning=kTRUE,
						       Double_t sigmaElectronRejection=3,
						       Bool_t electronExclusiveRejection=kFALSE,
						       TString correctionFileName = "",
						       Int_t nCentralityArrayBinsForCorrection = -1,
						       Double_t *gCentralityArrayForCorrections = 0x0,
						       Bool_t gRunEbyE = kFALSE,
						       Bool_t bMomentumOrdering = kTRUE,
						       AliAnalysisTaskBFPsi::eCorrProcedure corrProc = AliAnalysisTaskBFPsi::kNoCorr,
						       Bool_t bUseSameLabelMCCut = kFALSE) {
  // Creates a balance function analysis task and adds it to the analysis manager.
  // Get the pointer to the existing analysis manager via the static access method.
  TString outputFileName(fileNameBase);
  outputFileName.Append(".root");

  TGrid::Connect("alien:",0,0,"t");
  
  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskBF", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //===========================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskBF", "This task requires an input event handler");
    return NULL;
  }
  TString analysisType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())) analysisType = "MC";
  
  // to set the analysis type manually
  if(analysisTypeUser != ""){
    analysisType = analysisTypeUser;
    ::Info("AddTaskBF",Form("Analysis Type manually set to %s",analysisType.Data()));
  }

  // for local changed BF configuration
  //gROOT->LoadMacro("./configBalanceFunctionPsiAnalysis.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGCF/EBYE/macros/configBalanceFunctionPsiAnalysis.C");
  AliBalancePsi *bf  = 0;  // Balance Function object
  AliBalancePsi *bfs = 0;  // shuffled Balance function object
  AliBalancePsi *bfm = 0;  // mixing Balance function object
  AliBalanceEbyE *bfebye = 0;  // EbyE Balance function object

  //maximum Delta eta range
  Double_t deltaEtaMax=TMath::Abs(etaMax-etaMin);

  if (analysisType=="ESD"){
    bf  = GetBalanceFunctionObject("ESD",centralityEstimator,centrMin,centrMax,kFALSE,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,fArgEventClass,deltaEtaMax,bVertexBinning,bMomentumOrdering);
    if(gRunShuffling) bfs = GetBalanceFunctionObject("ESD",centralityEstimator,centrMin,centrMax,kTRUE,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,fArgEventClass,deltaEtaMax,bVertexBinning,bMomentumOrdering);
    if(gRunMixing)    bfm = GetBalanceFunctionObject("ESD",centralityEstimator,centrMin,centrMax,kFALSE,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,fArgEventClass,deltaEtaMax,bVertexBinning,bMomentumOrdering);
  }
  else if (analysisType=="AOD"){
    bf  = GetBalanceFunctionObject("AOD",centralityEstimator,centrMin,centrMax,kFALSE,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,fArgEventClass,deltaEtaMax,bVertexBinning,bMomentumOrdering);
    if(gRunShuffling) bfs = GetBalanceFunctionObject("AOD",centralityEstimator,centrMin,centrMax,kTRUE,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,fArgEventClass,deltaEtaMax,bVertexBinning,bMomentumOrdering);
    if(gRunMixing)    bfm = GetBalanceFunctionObject("AOD",centralityEstimator,centrMin,centrMax,kFALSE,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,fArgEventClass,deltaEtaMax,bVertexBinning,bMomentumOrdering);
    if(gRunEbyE)      bfebye = GetBalanceFunctionEbyEObject(bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,deltaEtaMax,bMomentumOrdering);
  }
  else if (analysisType=="MC"){
    bf  = GetBalanceFunctionObject("MC",centralityEstimator,centrMin,centrMax,kFALSE,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,fArgEventClass,deltaEtaMax,bVertexBinning,bMomentumOrdering);
    if(gRunShuffling) bfs = GetBalanceFunctionObject("MC",centralityEstimator,centrMin,centrMax,kTRUE,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,fArgEventClass,deltaEtaMax,bVertexBinning,bMomentumOrdering);
    if(gRunMixing)    bfm = GetBalanceFunctionObject("MC",centralityEstimator,centrMin,centrMax,kFALSE,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,fArgEventClass,deltaEtaMax,bVertexBinning,bMomentumOrdering);
  }
  else if (analysisType=="MCAOD"){
    bf  = GetBalanceFunctionObject("MCAOD",centralityEstimator,centrMin,centrMax,kFALSE,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,fArgEventClass,deltaEtaMax,bVertexBinning,bMomentumOrdering);
    if(gRunShuffling) bfs = GetBalanceFunctionObject("MCAOD",centralityEstimator,centrMin,centrMax,kTRUE,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,fArgEventClass,deltaEtaMax,bVertexBinning,bMomentumOrdering);
    if(gRunMixing)    bfm = GetBalanceFunctionObject("MCAOD",centralityEstimator,centrMin,centrMax,kFALSE,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,fArgEventClass,deltaEtaMax,bVertexBinning,bMomentumOrdering);
  }
  else if (analysisType=="MCAODrec"){
    bf  = GetBalanceFunctionObject("MCAODrec",centralityEstimator,centrMin,centrMax,kFALSE,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,fArgEventClass,deltaEtaMax,bVertexBinning,bMomentumOrdering, bUseSameLabelMCCut);
    if(gRunShuffling) bfs = GetBalanceFunctionObject("MCAODrec",centralityEstimator,centrMin,centrMax,kTRUE,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,fArgEventClass,deltaEtaMax,bVertexBinning,bMomentumOrdering, bUseSameLabelMCCut);
    if(gRunMixing)    bfm = GetBalanceFunctionObject("MCAODrec",centralityEstimator,centrMin,centrMax,kFALSE,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,fArgEventClass,deltaEtaMax,bVertexBinning,bMomentumOrdering, bUseSameLabelMCCut);
  }
  else if (analysisType=="AODnano"){
    bf  = GetBalanceFunctionObject("AODnano",centralityEstimator,centrMin,centrMax,kFALSE,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,fArgEventClass,deltaEtaMax,bVertexBinning,bMomentumOrdering);
    if(gRunShuffling) bfs = GetBalanceFunctionObject("AODnano",centralityEstimator,centrMin,centrMax,kTRUE,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,fArgEventClass,deltaEtaMax,bVertexBinning,bMomentumOrdering);
    if(gRunMixing)    bfm = GetBalanceFunctionObject("AODnano",centralityEstimator,centrMin,centrMax,kFALSE,bResonancesCut,bHBTcut,HBTCutValue,bConversionCut,invMassForConversionCut,bMomentumDifferenceCut,fQCutMin,fArgEventClass,deltaEtaMax,bVertexBinning,bMomentumOrdering);
  }
  else{
    ::Error("AddTaskBF", "analysis type NOT known.");
    return NULL;
  }

  // Create the task, add it to manager and configure it.
  //===========================================================================
  AliAnalysisTaskBFPsi *taskBF = new AliAnalysisTaskBFPsi(Form("TaskBFPsi_%.0f-%.0f_Bit%d_%s%s",centrMin,centrMax,AODfilterBit,centralityEstimator.Data(),dirNameExtra.Data()));
  
  //Event characteristics scheme
  taskBF->SetEventClass(fArgEventClass);
  //taskBF->SetCustomBinning("centralityVertex:0,80");
  //taskBF->SetCustomBinning("multiplicity:0,260");
  
  if(fArgEventClass == "Multiplicity") {
    if(analysisType == "MC")
      taskBF->SetMultiplicityRange(centrMin,centrMax);
    else {
    taskBF->SetPercentileRange(centrMin,centrMax);
    taskBF->SetMultiplicityEstimator(centralityEstimator);
    cout<<"Multiplicity estimator "<<centralityEstimator.Data()<<endl;
    }
  }
  else if(fArgEventClass == "Centrality") {
    if(analysisType == "MC")
      taskBF->SetImpactParameterRange(centrMin,centrMax);
    else {
      taskBF->SetPercentileRange(centrMin,centrMax);
      //taskBF->SetCentralityPercentileRange(centrMin,centrMax);
      // centrality estimator (default = V0M)
      taskBF->SetCentralityEstimator(centralityEstimator);
      cout<<"Centrality estimator "<<centralityEstimator.Data()<<endl;
    }
  }

  
  //+++++++++++++++++++++

  taskBF->SetAnalysisObject(bf);
  if(gRunShuffling) taskBF->SetShufflingObject(bfs);
  if(gRunMixing){
    taskBF->SetMixingObject(bfm);
    taskBF->SetMixingTracks(50000);
    if(gRunMixingWithEventPlane){
      taskBF->SetMixingWithEventPlane(gRunMixingWithEventPlane);
    }
  }
  if(gRunEbyE) taskBF->SetEbyEObject(bfebye);

  if(analysisType == "ESD") {
    AliESDtrackCuts *trackCuts = GetTrackCutsObject(ptMin,ptMax,etaMin,etaMax,maxTPCchi2,DCAxy,DCAz,minNClustersTPC);
    taskBF->SetAnalysisCutObject(trackCuts);
    if(kUsePID) {
      if(kUseBayesianPID)
	taskBF->SetUseBayesianPID(gMinAcceptedProbability);
      else if(kUseNSigmaPID)
	taskBF->SetUseNSigmaPID(nSigmaMax);
      taskBF->SetParticleOfInterest(AliAnalysisTaskBFPsi::kPion);
      taskBF->SetDetectorUsedForPID(AliAnalysisTaskBFPsi::kTOFpid);
    }
  }
  else if(analysisType == "AOD" || analysisType == "AODnano") {
    // pt and eta cut (pt_min, pt_max, eta_min, eta_max)
    taskBF->SetAODtrackCutBit(AODfilterBit);
    taskBF->SetKinematicsCutsAOD(ptMin,ptMax,etaMin,etaMax);

    // set extra DCA cuts (-1 no extra cut)
    taskBF->SetExtraDCACutsAOD(DCAxy,DCAz);

    // set extra TPC chi2 / nr of clusters cut
    taskBF->SetExtraTPCCutsAOD(maxTPCchi2, minNClustersTPC);

    // electron rejection (so far only for AOD), <0 --> no rejection
    if(sigmaElectronRejection > 0){
      if(electronExclusiveRejection) taskBF->SetElectronOnlyRejection(sigmaElectronRejection); // no other particle in nsigma 
      else                           taskBF->SetElectronRejection(sigmaElectronRejection); // check only if electrons in nsigma
    }

    //++++++++++++++++//
    /* if(kUsePID) {
      if(kUseBayesianPID)
	taskBF->SetUseBayesianPID(gMinAcceptedProbability);
      else if(kUseNSigmaPID)
	taskBF->SetUseNSigmaPID(nSigmaMax);
      taskBF->SetParticleOfInterest(AliAnalysisTaskBFPsi::kKaon);
      taskBF->SetDetectorUsedForPID(AliAnalysisTaskBFPsi::kTPCTOF); //TOFpid,TPCpid
      }*/
    //++++++++++++++++//

  }
  else if(analysisType == "MC") {
    taskBF->SetKinematicsCutsAOD(ptMin,ptMax,etaMin,etaMax);
  }
  else if(analysisType == "MCAOD") {
    // pt and eta cut (pt_min, pt_max, eta_min, eta_max)
    taskBF->SetAODtrackCutBit(AODfilterBit);
    taskBF->SetKinematicsCutsAOD(ptMin,ptMax,etaMin,etaMax);    
  }
  else if(analysisType == "MCAODrec") {     //++++++++++++++++
    // pt and eta cut (pt_min, pt_max, eta_min, eta_max)
    taskBF->SetAODtrackCutBit(AODfilterBit);
    taskBF->SetKinematicsCutsAOD(ptMin,ptMax,etaMin,etaMax); 

    // set extra DCA cuts (-1 no extra cut)
    taskBF->SetExtraDCACutsAOD(DCAxy,DCAz);

    // set extra TPC chi2 / nr of clusters cut
    taskBF->SetExtraTPCCutsAOD(maxTPCchi2, minNClustersTPC);

    // electron rejection (so far only for AOD), <0 --> no rejection
    if(sigmaElectronRejection > 0){
      if(electronExclusiveRejection) taskBF->SetElectronOnlyRejection(sigmaElectronRejection); // no other particle in nsigma 
      else                           taskBF->SetElectronRejection(sigmaElectronRejection); // check only if electrons in nsigma
    }
  }//++++++++++++++++

  // offline trigger selection (AliVEvent.h)
  // taskBF->UseOfflineTrigger(); // NOT used (selection is done with the AliAnalysisTaskSE::SelectCollisionCandidates()) 
  // with this only selected events are analyzed (first 2 bins in event QA histogram are the same))
  // documentation in https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PWG1EvSelDocumentation

  
  if(triggerSel == AliAnalysisTaskBFPsi::kCentral) taskBF->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
  else if(triggerSel == AliAnalysisTaskBFPsi::kMB) taskBF->SelectCollisionCandidates(AliVEvent::kMB);
  else if(triggerSel == AliAnalysisTaskBFPsi::kINT7) taskBF->SelectCollisionCandidates(AliVEvent::kINT7);

  // centrality estimator (default = V0M)
  taskBF->SetCentralityEstimator(centralityEstimator);
  
  // vertex cut (x,y,z)
  //taskBF->SetVertexDiamond(3.,3.,vertexZ);

  taskBF->SetCorrectionProcedure(corrProc);


  //++++++++++++++++++++++
  // Efficiency + Contamination corrections
  // If correctionFileName = "", do not use corrections
  // if(corrProc == AliAnalysisTaskBFPsi::kMCCorr)
  // taskBF->SetInputCorrection(Form("$ALICE_PHYSICS/PWGCF/EBYE/BalanceFunctions/Corrections/%s",correctionFileName.Data()),nCentralityArrayBinsForCorrection,gCentralityArrayForCorrections);
  
  /*else if (corrProc == AliAnalysisTaskBFPsi::kDataDrivCorr){

    TFile* fNUAFile = TFile::Open(nuaCorrFileName.Data(),"READ");
    TFile* fNUEFile = TFile::Open(nueCorrFileName.Data(),"READ");
    
    if(!fNUAFile) {
      printf(" *** ERROR: NUA file not found! **EXIT** ");
    } 
    TList* fListNUA = dynamic_cast<TList*>(fNUAFile->Get("fListNUA"));
    if(fListNUA)
      taskBF->SetInputListForNUACorr(fListNUA);
    else
      printf(" *** ERROR: NUA List not found! **EXIT**");
    
    
    if(!fNUEFile) {
      printf(" *** ERROR: NUE file not found! **EXIT** ");
    } 
    TList* fListNUE = dynamic_cast<TList*>(fNUEFile->Get("fListNUE"));
    if(fListNUE)
      taskBF->SetInputListForNUECorr(fListNUE);
    else
      printf(" *** ERROR: NUE List not found! **EXIT**");    
  }
  */
  
  //bf->PrintAnalysisSettings();
  mgr->AddTask(taskBF);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGCFEbyE.outputBalanceFunctionPsiAnalysis";
  AliAnalysisDataContainer *coutQA = mgr->CreateContainer(Form("listQAPsi_%.0f-%.0f_Bit%d_%s%s",centrMin,centrMax,AODfilterBit,centralityEstimator.Data(),dirNameExtra.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  AliAnalysisDataContainer *coutBF = mgr->CreateContainer(Form("listBFPsi_%.0f-%.0f_Bit%d_%s%s",centrMin,centrMax,AODfilterBit,centralityEstimator.Data(),dirNameExtra.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  if(gRunShuffling) AliAnalysisDataContainer *coutBFS = mgr->CreateContainer(Form("listBFPsiShuffled_%.0f-%.0f_Bit%d_%s%s",centrMin,centrMax,AODfilterBit,centralityEstimator.Data(),dirNameExtra.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  if(gRunMixing) AliAnalysisDataContainer *coutBFM = mgr->CreateContainer(Form("listBFPsiMixed_%.0f-%.0f_Bit%d_%s%s",centrMin,centrMax,AODfilterBit,centralityEstimator.Data(),dirNameExtra.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  if(kUsePID || sigmaElectronRejection > 0) AliAnalysisDataContainer *coutQAPID = mgr->CreateContainer(Form("listQAPIDPsi_%.0f-%.0f_Bit%d_%s%s",centrMin,centrMax,AODfilterBit,centralityEstimator.Data(),dirNameExtra.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());

  mgr->ConnectInput(taskBF, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskBF, 1, coutQA);
  mgr->ConnectOutput(taskBF, 2, coutBF);
  if(gRunShuffling) mgr->ConnectOutput(taskBF, 3, coutBFS);
  if(gRunMixing) mgr->ConnectOutput(taskBF, 4, coutBFM);
  if(kUsePID||sigmaElectronRejection > 0) mgr->ConnectOutput(taskBF, 5, coutQAPID);
  //if((kUsePID && analysisType == "AOD")||sigmaElectronRejection > 0) mgr->ConnectOutput(taskBF, 5, coutQAPID);
  //if((kUsePID && analysisType == "ESD")||sigmaElectronRejection > 0) mgr->ConnectOutput(taskBF, 5, coutQAPID);

  return taskBF;
}
