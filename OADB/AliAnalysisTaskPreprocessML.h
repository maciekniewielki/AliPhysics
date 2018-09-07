#ifndef __ALIANALYSISTASKPREPROCESSML__
#define __ALIANALYSISTASKPREPROCESSML__

#define MULTBINS 1
#define PARTTYPES 5

#include "TFile.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliAODTrack.h"
#include "AliPIDResponse.h"
#include "AliAODpidUtil.h"
#include "THnSparse.h"
#include <fstream>

class AliAnalysisUtils;

//main class used for iterating over particle collisions ( events )
//and saving all data into designated TTree container
class AliAnalysisTaskPreprocessML : public AliAnalysisTaskSE{
public:
	AliAnalysisTaskPreprocessML();
	AliAnalysisTaskPreprocessML(const Char_t *partName);
	virtual ~AliAnalysisTaskPreprocessML();

	//prepare output objects
	virtual void UserCreateOutputObjects();

	//main function used for iterating over collisions
	virtual void UserExec(Option_t *option);

	//tracks criteria
	bool isTrackValid(double eta, double pt, bool covxyz);
	bool IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP, float TOFtime);
	bool IsPionNSigma(float mom, float nsigmaTPCPi, float nSigmaTOFPi, float TOFtime);
	bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK, float TOFtime);

	//load attributes of track to designated variables
	void loadTrackInfo(AliAODTrack *track);

private:
	AliAnalysisTaskPreprocessML(const AliAnalysisTaskPreprocessML &); // copy constructor
	AliAnalysisTaskPreprocessML &operator=(const AliAnalysisTaskPreprocessML &); // operator=

	//output list holding resulting TTree
	TList *fHistoList;

	//PID objects
	AliPIDResponse *fpidResponse;
	AliAODpidUtil	*fAODpidUtil;

	//container holding all attributes of observations
	TTree *treeOutput;

	//TTree branches containing tracks attributes
	//and traditional classification results
	int PDGCode;
	double TPCNcls;
	double TPCsignal;
	double ped;
	double Pt;
	double Px;
	double Py;
	double Pz;
	double nSigmaTOFPi;
	double nSigmaTOFK;
	double nSigmaTOFP;
	double nSigmaTOFe;
	double nSigmaTPCPi;
	double nSigmaTPCK;
	double nSigmaTPCP;
	double nSigmaTPCe;
	double tbeta;
	double cov0;
	double cov1;
	double cov2;
	double cov3;
	double cov4;
	double cov5;
	double cov6;
	double cov7;
	double cov8;
	double cov9;
	double cov10;
	double cov11;
	double cov12;
	double cov13;
	double cov14;
	double cov15;
	double cov16;
	double cov17;
	double cov18;
	double cov19;
	double cov20;
	int trackID;
	long eventID;
	//int isPion;
	//int isKaon;
	//int isProton;
    double fV1[3];
	ClassDef(AliAnalysisTaskPreprocessML, 0);
};

#endif // __ALIANALYSISTASKPREPROCESSML__
