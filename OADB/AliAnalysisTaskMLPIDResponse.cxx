#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <cmath>
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TObjArray.h"
#include "TString.h"
#include "TParticle.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"

#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"

#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliInputEventHandler.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliPIDResponse.h"

#include "AliAODMLpidUtil.h"
#include "AliAODHeader.h"

#include "AliAnalysisTaskMLPIDResponse.h"


class AliAnalysisUtils;

ClassImp(AliAnalysisTaskMLPIDResponse)



AliAnalysisTaskMLPIDResponse::AliAnalysisTaskMLPIDResponse() :
	AliAnalysisTaskSE(), pidUtil(0), trackPipe(0), pidPipe(0)
{
}

AliAnalysisTaskMLPIDResponse::AliAnalysisTaskMLPIDResponse(const Char_t *partName) :
	AliAnalysisTaskSE(partName), pidUtil(0), trackPipe(0), pidPipe(0)
{
    //std::cout<<"Constructing MLPIDResponse"<<std::endl;
}


AliAnalysisTaskMLPIDResponse::~AliAnalysisTaskMLPIDResponse()
{
}

void AliAnalysisTaskMLPIDResponse::UserCreateOutputObjects()
{
	//********** PID ****************
    
    //std::cout<<"AliAnalysisTaskMLPIDResponse::UserCreateOutputObjects"<<std::endl;

    binaryCommunication = true;

    if (binaryCommunication)
    {
        std::cout<<"Starting the MLPIDResponse task with binary communication"<<std::endl;
        trackPipe = new std::ofstream("MLPIDTrackPipe", std::ios::binary);
        pidPipe = new std::ifstream("MLPIDProbabilityPipe", std::ios::binary);
    }
    else
    {
        std::cout<<"Starting the MLPIDResponse task with txt/ascii communication"<<std::endl;
        trackPipe = new std::ofstream("MLPIDTrackPipe");
        pidPipe = new std::ifstream("MLPIDProbabilityPipe");
    }

    // Ask the classifier for the number of classes
    // And for the classes themselves
    std::cout<<"Asking for the classes"<<std::endl;
    if (binaryCommunication)
    {
        Int_t data = 0;
        trackPipe->write(reinterpret_cast<char*>(&data), 4);
        trackPipe->flush();
        std::cout<<"Asked for the classes"<<std::endl;
        
        Int_t nClasses = 0;
        pidPipe->read(reinterpret_cast<char*>(&nClasses), 4);
        std::cout<<"Just read the number of classes: "<<nClasses<<std::endl;
        
        possiblePdgCodes.resize(nClasses);
        pidPipe->read(reinterpret_cast<char*>(&possiblePdgCodes[0]), nClasses*4);
        std::cout<<"Just read all the classes"<<std::endl;
    }
    else
    {
        (*trackPipe)<<0<<std::endl;
        trackPipe->flush();

        Int_t nClasses;
        Int_t currentClass;
        (*pidPipe)>>nClasses;
        possiblePdgCodes.resize(nClasses);

        for(Int_t i=0; i<nClasses; i++)
        {
            (*pidPipe)>>currentClass;
            possiblePdgCodes[i] = currentClass;
        }
    }
    std::cout<<"Got "<<possiblePdgCodes.size()<<" classes:"<<std::endl;
    for (Int_t i=0; i<possiblePdgCodes.size(); i++)
        std::cout<<possiblePdgCodes[i]<<" ";
    std::cout<<std::endl;


    //Translating TChain to list of files


    const Int_t numberOfTypes = 5;
    Int_t particleTypeToRead[] = {11, 13, 211, 321, 2212};
    Int_t eventID = 0;
    Int_t trackID = 0;

	pidUtil = new AliAODMLpidUtil();
    //const char* argv[] = {"-i PreprocessedMLAOD.root", "-g" };
    //TPython::ExecScript("classifier.py", 2, argv);


	AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
	AliAODInputHandler* inputHandler = (AliAODInputHandler*) (man->GetInputEventHandler());
    inputHandler->SetMLpidUtil(pidUtil);
	fpidResponse = inputHandler->GetPIDResponse();
}

bool AliAnalysisTaskMLPIDResponse::isTrackValid(double eta, double pt, bool covxyz)
{
	//typical cuts (selection criteria)
	if(eta < -0.8 || eta > 0.8) //pseudorapidity range of TPC
		return false;

	if (pt < 0.2 || pt > 20) //transverse momentum
		return false;

	if ( !covxyz )
		return false;

	return true;
}

void AliAnalysisTaskMLPIDResponse::predictTracksPID(AliAODEvent* aodEvent)
{
    // Clear the responses from previous event
    pidUtil->clearResponses();

    // Used in binary
    std::vector<float> binTracksData;
    std::vector<float> trackProbabilities;
    // Used in txt
    std::stringstream txtTracksData;

    std::vector<Int_t> validTrackIDs;

    const Int_t n_features = 25;
    float trackData[n_features];
    const Double_t speedLight = 2.99792458e-2;
    
	for (Int_t iTracks = 0; iTracks < aodEvent->GetNumberOfTracks(); iTracks++) {
		//get track
		AliAODTrack *track = (AliAODTrack*)aodEvent->GetTrack(iTracks);
		if (!track)
			continue;

		UInt_t filterBit = 96;
		if(!track->TestFilterBit(filterBit))
			continue;

		//check if particle was correctly measured
		double cov[21];
		if( !isTrackValid( track->Eta(), track->Pt(), track->GetCovarianceXYZPxPyPz(cov) ) )
			continue;
        
		//check particle nsigma status
		nSigmaTOFPi = -9999;
		nSigmaTOFK = -9999;
		nSigmaTOFP = -9999;
		nSigmaTOFe = -9999;
		Int_t statusTOFPi =
			fpidResponse->NumberOfSigmas(AliPIDResponse::kTOF,track,AliPID::kPion, nSigmaTOFPi);
		Int_t statusTOFK =
			fpidResponse->NumberOfSigmas(AliPIDResponse::kTOF,track,AliPID::kKaon, nSigmaTOFK);
		Int_t statusTOFP =
			fpidResponse->NumberOfSigmas(AliPIDResponse::kTOF,track,AliPID::kProton, nSigmaTOFP);
		Int_t statusTOFe =
			fpidResponse->NumberOfSigmas(AliPIDResponse::kTOF,track,AliPID::kElectron, nSigmaTOFe);
		if((statusTOFPi != AliPIDResponse::kDetPidOk)
				|| (statusTOFK != AliPIDResponse::kDetPidOk)
				|| (statusTOFP != AliPIDResponse::kDetPidOk)
				|| (statusTOFe != AliPIDResponse::kDetPidOk))
			continue;

		nSigmaTPCPi = -9999;
		nSigmaTPCK = -9999;
		nSigmaTPCP = -9999;
		nSigmaTPCe = -9999;
		Int_t statusTPCPi =
			fpidResponse->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kPion, nSigmaTPCPi);
		Int_t statusTPCK =
			fpidResponse->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kKaon, nSigmaTPCK);
		Int_t statusTPCP =
			fpidResponse->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kProton, nSigmaTPCP);
		Int_t statusTPCe =
			fpidResponse->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kElectron, nSigmaTPCe);
		if((statusTPCPi != AliPIDResponse::kDetPidOk)
				|| (statusTPCK != AliPIDResponse::kDetPidOk)
				|| (statusTPCP != AliPIDResponse::kDetPidOk)
				|| (statusTPCe != AliPIDResponse::kDetPidOk))
		    continue;

        float tTofSig = track->GetTOFsignal(); //TOF signal
        double TOFTime = tTofSig - fAODpidUtil->GetTOFResponse().GetStartTime(track->P());
        if (TOFTime == 0.0)
            continue;

        
        double velocity = 0.0;
        velocity = track->GetIntegratedLength() / TOFTime;
        double tbeta = velocity / speedLight;

        validTrackIDs.push_back(track->GetID());
	    //fill track array with values
        trackData[0] = track->GetTPCNcls();
        trackData[1] = track->GetTPCsignal();
        trackData[2] = track->P();
        trackData[3] = track->Pt();
        trackData[4] = track->Px();
        trackData[5] = track->Py();
        trackData[6] = track->Pz();
        trackData[7] = nSigmaTOFPi;
        trackData[8] = nSigmaTOFK;
        trackData[9] = nSigmaTOFP;
        trackData[10] = nSigmaTOFe;
        trackData[11] = nSigmaTPCPi;
        trackData[12] = nSigmaTPCK;
        trackData[13] = nSigmaTPCP;
        trackData[14] = nSigmaTPCe;
        trackData[15] = tbeta;
        trackData[16] = cov[2];
        trackData[17] = cov[5];
        trackData[18] = cov[9];
        trackData[19] = cov[13];
        trackData[20] = cov[14];
        trackData[21] = cov[17];
        trackData[22] = cov[18];
        trackData[23] = cov[19];
        trackData[24] = cov[20];
        
        if (binaryCommunication)
            for (Int_t j=0; j<n_features; j++)
                binTracksData.push_back(trackData[j]);
        else
        {
            for (Int_t j=0; j<n_features-1;j++)
                txtTracksData<<trackData[j]<<" ";
            txtTracksData<<trackData[n_features-1]<<std::endl;
        }

	  }
    Int_t n_tracks = validTrackIDs.size();

    if (n_tracks<=0)
        return;

    //std::cout<<"Sending "<<n_tracks<<" to the classifier"<<std::endl;
    if (binaryCommunication)
    {
        // send
        trackPipe->write(reinterpret_cast<char*>(&n_tracks), 4);
        trackPipe->write(reinterpret_cast<char*>(&binTracksData[0]), binTracksData.size()*4); 
        trackPipe->flush();

        // receive
        Int_t floatsToReceive = n_tracks * possiblePdgCodes.size();
        trackProbabilities.resize(floatsToReceive);
        pidPipe->read(reinterpret_cast<char*>(&trackProbabilities[0]), floatsToReceive * 4);
    }
    else
    {
        (*trackPipe)<<n_tracks<<"\n";
        (*trackPipe)<<n_features<<"\n";
        (*trackPipe)<<txtTracksData.rdbuf();
        trackPipe->flush();
    }
    


    for(Int_t i=0; i<n_tracks; i++)
    {
        float maxProb = 0;
        float currentProb = 0;
        Int_t maxProbPdg = 0;
        Int_t currentPdg = 0;

        AliMLPIDResponse* resp = new AliMLPIDResponse;
        for(Int_t j=0; j<possiblePdgCodes.size(); j++)
        {
            currentPdg = possiblePdgCodes[j];

            if (binaryCommunication)
                currentProb = trackProbabilities[i*possiblePdgCodes.size() + j];
            else
                (*pidPipe)>>currentProb;

            resp->probabilities[currentPdg] = currentProb;
            if (currentProb > maxProb)
            {
                maxProb = currentProb;
                maxProbPdg = currentPdg;
            }
        }
        resp->predictedPDG = maxProbPdg;
        pidUtil->addPIDResponse(validTrackIDs[i], resp);
    }
}

//function iterating over every available track to get its info
//and save it into TTree
void AliAnalysisTaskMLPIDResponse::UserExec(Option_t *)
{
	AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler *>(
			AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
	AliAODEvent *fAOD = aodH->GetEvent();
	fAODpidUtil = aodH->GetAODpidUtil();

	//get event
	AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(InputEvent());
	if (!aodEvent) return;
	AliAODHeader *fAODheader = (AliAODHeader*)aodEvent->GetHeader();
	Double_t mult = fAODheader->GetRefMultiplicity();
	if(mult<0) return;

	// EVENT SELECTION ********************
	const AliAODVertex* vertex =(AliAODVertex*) aodEvent->GetPrimaryVertex();
    float fV1[3];
	vertex->GetPosition(fV1); //position of primary vertex
	if (!vertex || vertex->GetNContributors()<=0) return;

	AliAnalysisUtils *anaUtil=new AliAnalysisUtils();

	Bool_t fMVPlp = kFALSE;
	Bool_t fisPileUp = kTRUE;
	Int_t fMinPlpContribMV = 0;
	Int_t fMinPlpContribSPD = 3;

	if(fMVPlp)
		anaUtil->SetUseMVPlpSelection(kTRUE);
	else
		anaUtil->SetUseMVPlpSelection(kFALSE);

	if(fMinPlpContribMV)
		anaUtil->SetMinPlpContribMV(fMinPlpContribMV);

	if(fMinPlpContribSPD)
		anaUtil->SetMinPlpContribSPD(fMinPlpContribSPD);

	if(fisPileUp)
		if(anaUtil->IsPileUpEvent(aodEvent))
			return;

	delete anaUtil;

	Float_t zvtx = vertex->GetZ();
	if (TMath::Abs(zvtx) > 10)
		return;

    predictTracksPID(aodEvent);
}

void AliAnalysisTaskMLPIDResponse::Terminate(Option_t *)
{
    if (binaryCommunication)
    {
        Int_t endFlag = -1;
        trackPipe->write(reinterpret_cast<char*>(&endFlag), 4);
    }
    else
        (*trackPipe)<<-1<<std::endl;

    trackPipe->close();
    pidPipe->close();

}

