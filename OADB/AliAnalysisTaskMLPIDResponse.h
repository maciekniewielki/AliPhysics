#ifndef __ALIANALYSISTASKMLPIDRESPONSE__
#define __ALIANALYSISTASKMLPIDRESPONSE__

#define MULTBINS 1
#define PARTTYPES 5

#include "TChain.h"
#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"
#include <fstream>
#include <sstream>
#include "AliAnalysisUtils.h"
#include "AliAODMLpidUtil.h"
#include "AliAODpidUtil.h"
#include "AliPIDResponse.h"

//main class used for iterating over particle collisions ( events )
//and saving all data into designated TTree container
class AliAnalysisTaskMLPIDResponse : public AliAnalysisTaskSE{
public:
	AliAnalysisTaskMLPIDResponse();
	AliAnalysisTaskMLPIDResponse(const Char_t *partName);
	virtual ~AliAnalysisTaskMLPIDResponse();

	//prepare output objects
	virtual void UserCreateOutputObjects();

	//main function used for iterating over collisions
	virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);

    void setInputTChain(TChain* chain){inputChain=chain;}
    void setIsMC(bool isMC){this->isMC=isMC;}
    void saveCollisionCandidates(UInt_t collCand){collisionCandidates=collCand;}

    bool isTrackValid(double eta, double pt, bool covxyz);
    void predictTracksPID(AliAODEvent* aodEvent);



private:
	AliAnalysisTaskMLPIDResponse(const AliAnalysisTaskMLPIDResponse &); // copy constructor
	AliAnalysisTaskMLPIDResponse &operator=(const AliAnalysisTaskMLPIDResponse &); // operator=
    
    TChain* inputChain;
    bool isMC;
    UInt_t collisionCandidates;
    bool binaryCommunication;

	AliPIDResponse *fpidResponse;
    AliAODMLpidUtil* pidUtil;
    AliAODpidUtil *fAODpidUtil;

    std::ofstream* trackPipe;
    std::ifstream* pidPipe;
    std::vector<int> possiblePdgCodes;
    double nSigmaTOFPi;
    double nSigmaTOFK;
    double nSigmaTOFP;
    double nSigmaTOFe;
    double nSigmaTPCPi;
    double nSigmaTPCK;
    double nSigmaTPCP;
    double nSigmaTPCe;

	ClassDef(AliAnalysisTaskMLPIDResponse, 0);
};

#endif // __ALIANALYSISTASKMLPIDRESPONSE__
