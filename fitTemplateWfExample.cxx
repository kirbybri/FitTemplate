#include "FitTemplateWf.hxx"
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

void fitTemplateWfExample(){

        TCanvas* c0 = new TCanvas("c0", "c0",1400,800);

	//vectors to hold samples, sample valid flag
	std::vector<double> tempVal;
	std::vector<double> dataVal;
	std::vector<bool> dataQuality;

        ifstream inFile;
        inFile.open("atb_tempWf.txt");
        if (!inFile) {
        	std::cout << "Unable to open file" << std::endl;
        	exit(1); // terminate with error
        }
        double samp;
    	while (inFile >> samp) {
		tempVal.push_back( samp );
        }

        ifstream inFile_roi;
        inFile_roi.open("atb_roi.txt");
        if (!inFile_roi) {
        	std::cout << "Unable to open file" << std::endl;
        	exit(1); // terminate with error
        }
    	while (inFile_roi >> samp) {
		if( dataVal.size() > 9 ) continue;
		dataVal.push_back(samp);
		dataQuality.push_back(1);
        }

	FitTemplateWf *tempFit = new FitTemplateWf();
	tempFit->setSampleError(0.7);
	std::cout << "SAMPLE ERROR " << tempFit->sampleErr << std::endl;
        tempFit->addTemplate(tempVal,1./30.);
	tempFit->addData(dataVal, dataQuality);

        TGraph* gFit_temp = new TGraph();
	TGraph* gFit_data = new TGraph();

	for(unsigned int s = 0 ; s < tempFit->tempData->samp.size() ; s++ ){
		samp = tempFit->tempData->samp[s] + 1046.82;
		gFit_temp->SetPoint( s, s*tempFit->tempData->sampPeriod , samp);
	}

	for(unsigned int s = 0 ; s < tempFit->fitData_vals.size() ; s++ ){
		samp = tempFit->fitData_vals[s];
		gFit_data->SetPoint( s, s , samp);
	}

	tempFit->initVals.push_back(1046.82);
	tempFit->initVals.push_back(0.5);
	tempFit->initVals.push_back(0.9);
	tempFit->initErrs.push_back(0.5);
	tempFit->initErrs.push_back(0.1);
	tempFit->initErrs.push_back(0.1);

	tempFit->fixFitVars.push_back(0);

	tempFit->showOutput = 1;
	tempFit->doFit();

	double fit_offset = tempFit->fitVals[0];
	double fit_pulseStartTime = tempFit->fitVals[1];
	double fit_amp = tempFit->fitVals[2];

        TGraph* gFit_result = new TGraph();
	for(unsigned int s = 0 ; s < tempFit->fitData_vals.size() ; s++ ){
		for( unsigned int frac = 0 ; frac < 100 ; frac++ ){
			double time = s + frac/100.;
			double fitVal = 0;
			tempFit->tempData->getSignalValue(time, fit_offset, fit_pulseStartTime, fit_amp,fitVal);
			gFit_result->SetPoint( gFit_result->GetN(), time , fitVal);
		}
	}

	gFit_temp->SetLineColor(kBlue);
	gFit_data->SetMarkerStyle(21);
        gFit_result->SetLineColor(kRed);

	c0->Clear();
	//gFit_temp->Draw("AL");
	gFit_data->Draw("AP");
        gFit_result->Draw("L");
	c0->Update();
}
