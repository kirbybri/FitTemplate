//compile independently with: g++ -std=c++11 -o doTemplateWfFit doTemplateWfFit.cxx `root-config --cflags --glibs` -lMinuit
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
using namespace std;

#include "TROOT.h"
#include "TMath.h"
#include "TApplication.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TImage.h"

#include "FitTemplateWf.hxx"

//global TApplication object declared here for simplicity
TApplication *theApp;

class Analyze {
  public:
  Analyze();
  ~Analyze();
  void doAnalysis();
  void printResults();

  //Files

  //data objects
  TCanvas* c0;
  TGraph* gFit_result;
  TGraph* gFit_data;
  TImage *img;

  //Fit objects
  FitTemplateWf *tempFit;
	
  //constants

  //variables
  int status;
  double fit_offset;
  double fit_pulseStartTime;
  double fit_amp;
  std::vector<double> tempVal;
  std::vector<double> dataVal;
  std::vector<bool> dataQuality;
  std::vector<double> fitConfig;

  //histograms
};

Analyze::Analyze(){
  status = 1;
  fit_offset = 0;
  fit_pulseStartTime = 0;
  fit_amp = 0;
  //data objects
  c0 = new TCanvas("c0", "c0",1200,600);
  gFit_result = new TGraph();
  gFit_data = new TGraph();
  img = TImage::Create();
  //fit objects	
  tempFit = new FitTemplateWf();
}

Analyze::~Analyze(){
  delete c0;
  delete tempFit;
  delete gFit_result;
  delete gFit_data;
  delete img;
}

void Analyze::doAnalysis(){
  //load template file
  ifstream inFile;
  inFile.open("atb_tempWf.txt");
  if (!inFile) {
    std::cout << "Unable to open file" << std::endl;
    return; // terminate with error
  }
  double samp;
  while (inFile >> samp)
    tempVal.push_back( samp );

  //load ROI file
  ifstream inFile_roi;
  inFile_roi.open("atb_roi.txt");
  if (!inFile_roi) {
    std::cout << "Unable to open file" << std::endl;
    return; // terminate with error
  }
  while (inFile_roi >> samp) {
    dataVal.push_back(samp);
    dataQuality.push_back(1);
  }

  //load ROI file
  ifstream inFile_config;
  inFile_config.open("atb_fitConfig.txt");
  if (!inFile_config) {
    std::cout << "Unable to open file" << std::endl;
    std::cout << status << "\t" << fit_offset << "\t" << fit_pulseStartTime << "\t" << fit_amp << std::endl;
    exit(1); // terminate with error
  }
  while (inFile_config >> samp)
    fitConfig.push_back(samp);

  //load data into fit object
  tempFit->addTemplate(tempVal,fitConfig[0]); //template vector, sample period
  tempFit->addData(dataVal, dataQuality);

  //set sample uncertainty
  tempFit->setSampleError(fitConfig[1]);

  //set initial conditions
  tempFit->initVals.push_back(fitConfig[2]); //init offset
  tempFit->initVals.push_back(fitConfig[3]); //init start time
  tempFit->initVals.push_back(fitConfig[4]); //init amp
  tempFit->initErrs.push_back(fitConfig[5]); //init offset err
  tempFit->initErrs.push_back(fitConfig[6]); //init start time err
  tempFit->initErrs.push_back(fitConfig[7]); //init amp err

  //fix fit parameter
  tempFit->fixFitVars.push_back(0);

  //do the fit
  tempFit->showOutput = 0;
  tempFit->doFit();

  //get fit results
  status = tempFit->status;
  fit_offset = tempFit->fitVals[0];
  fit_pulseStartTime = tempFit->fitVals[1];
  fit_amp = tempFit->fitVals[2];

  return;

  //draw result
  for(unsigned int s = 0 ; s < tempFit->fitData_vals.size() ; s++ ){
    samp = tempFit->fitData_vals[s];
    gFit_data->SetPoint( s, s , samp);
  }

  for(unsigned int s = 0 ; s < tempFit->fitData_vals.size() ; s++ ){
    for( unsigned int frac = 0 ; frac < 100 ; frac++ ){
      double time = s + frac/100.;
      double fitVal = 0;
      tempFit->tempData->getSignalValue(time, fit_offset, fit_pulseStartTime, fit_amp,fitVal);
      gFit_result->SetPoint( gFit_result->GetN(), time , fitVal);
    }
  }


  gFit_data->SetMarkerStyle(21);
  gFit_result->SetLineColor(kRed);

  c0->Clear();
  gFit_data->Draw("AP");
  gFit_result->Draw("L");
  c0->Modified();
  c0->Update();

  img->FromPad(c0);
  img->WriteImage("atb_fitImage.png");

  usleep(1000*1000);
  return;
}

void Analyze::printResults(){
  std::cout << status << "\t" << fit_offset << "\t" << fit_pulseStartTime << "\t" << fit_amp << std::endl;
}

void process() {
  Analyze ana;
  ana.doAnalysis();
  ana.printResults();
  return;
}

int main(int argc, char *argv[]){
  if(argc!=1){
    cout<<"Usage: doTemplateWfFit"<<endl;
    return 0;
  }

  //define ROOT application object
  //theApp = new TApplication("App", &argc, argv);
  process(); 

  //return 1;
  gSystem->Exit(0);
}
