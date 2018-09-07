//batch program extracting data from AOD files from range and saving them as a TTree
//in file specified by the user

#ifdef __CLING__
R__ADD_INCLUDE_PATH($ALICE_ROOT)
#include <ANALYSIS/macros/AddTaskPIDResponse.C>

// Tell ROOT where to find AliPhysics headers
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <OADB/macros/AddTaskPhysicsSelection.C>

#include "AliAnalysisTaskPreprocessML.h"
#include "AliAnalysisTaskPreprocessML.cxx"
#endif

void launchMLPIDPreprocessing(const char* filenames, const char* treeName, bool isMC, UInt_t collisionCandidates, const char* outfilename="PreprocessedMLAOD.root", int Nevents=20000000) {
	TStopwatch timer;
	timer.Start();

    TGrid::Connect("alien://");
	//Setting up required packages
	gSystem->Load("libTree.so");
	gSystem->Load("libGeom.so");
	gSystem->Load("libVMC.so");
	gSystem->Load("libPhysics.so");
	gSystem->Load("libSTEERBase.so");
	gSystem->Load("libESD.so");
	gSystem->Load("libAOD.so");
	gSystem->Load("libANALYSIS.so");
	gSystem->Load("libANALYSISalice.so");
	gSystem->Load("libEVGEN.so");
	gSystem->Load("libpythia6_4_25.so");
	gROOT->ProcessLine(".include $ALICE_ROOT/include");

	gSystem->AddIncludePath("-I$ALICE_ROOT/include");
	gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");

	//add files to analysis
	TChain *chain = new TChain(treeName);
    TObjArray *ty = TString(filenames).Tokenize(":");
    for (Int_t i = 0; i < ty->GetEntries(); i++)
        chain->Add(((TObjString *)(ty->At(i)))->String().Data());


    //char* path = isMC ? "276437" : "000276437/pass1";
	//char pname[2000];
	//for( int i = start_range; i < stop_range; ++i) {
		//sprintf(pname, "%s/Ali%d.root", path, i);
		//chain->Add(pname);
	//}

	// Make the analysis manager
	AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
	AliAODInputHandler* aodH = new AliAODInputHandler;
	mgr->SetInputEventHandler(aodH);

	//AddTaskPIDResponse - get labels from Monte Carlo data
    //ROOT5
	//gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");

    if (!isMC)
    {
        //Physics Selection
         AliPhysicsSelectionTask *test = AddTaskPhysicsSelection(isMC, kTRUE);
        //ROOT5
        //gROOT->ProcessLine(".L $ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
        //AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(isMC,kTRUE);
    }
	AddTaskPIDResponse(isMC);

	//Efficiency task - extracting data from tracks
	//gROOT->LoadMacro("AliAnalysisTaskPreprocessML.cxx+g");
	AliAnalysisTaskPreprocessML*myTask = new AliAnalysisTaskPreprocessML("MyTask");
	myTask->SelectCollisionCandidates(collisionCandidates);
	if( !myTask )
		exit(-1);
	mgr->AddTask(myTask);

	// Create containers for input/output
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("MyTree",
			TList::Class(),AliAnalysisManager::kOutputContainer,outfilename);

	//connect them to future analysis
	mgr->ConnectInput(myTask,0,cinput);
	mgr->ConnectOutput(myTask,1,coutput2);

	if ( !mgr->InitAnalysis() )
		return;
	mgr->PrintStatus();

	//start analysis
	mgr->StartAnalysis("local",chain,Nevents);

	timer.Stop();
	timer.Print();

	//save result in a file specified by user
	//extracting TTree from TList first
	TFile *ifile = new TFile(outfilename);
	TList *lists;
	ifile->GetObject("MyTree", lists);
	TFile *ofile = new TFile(outfilename, "RECREATE");

	//iterate over every object in TList and save it
	TIter next(lists);
	TObject *obj;
	while (obj = next())
	{
		ofile->cd();
		obj->Write();
	}
	delete ofile;
    
    const char* argv[] = {"-i PreprocessedMLAOD.root", "-g" };
    TPython::ExecScript("classifier.py", 2, argv);

}
