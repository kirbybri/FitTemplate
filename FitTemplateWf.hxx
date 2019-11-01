#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "TMinuit.h"

///////////TEMPLATE DATA CLASS/////////////////
class TemplateData {
private:

public:
  TemplateData();
  ~TemplateData();
  void addTemplate(const std::vector<double>& val, double period);
  void getSignalValue(double time, double offset, double pulseStartTime, double amp,double& simVal);
  std::vector<double> samp;
  double sampPeriod;
};

TemplateData::TemplateData(){
  sampPeriod = 0.01; //fractions of sample
  samp.clear();
}

TemplateData::~TemplateData(){
}

void TemplateData::addTemplate(const std::vector<double>& val, double period){
  samp.clear();
  samp = val;
  sampPeriod = period;
  return;
}

//get response function value quickly from vector interpolation
void TemplateData::getSignalValue(double time, double offset, double templateStartTime, double amp,double& simVal){
  simVal = offset; //simVal defaults to baseline if interpolation fails

  double sampTime = 0;
  if( time > templateStartTime )
   sampTime = time - templateStartTime;

  //determine position of time value in signal vector
  double templateSampTime = sampTime/sampPeriod; //convert req time into sample # within template array
  unsigned int templateSampNum = floor(templateSampTime); //get actual template array element #
  if( templateSampNum >= samp.size() - 1 ) //if req time exceeds template array, simVal is just the baseline
    return;
  
  //do linear interpolation
  double sigVal = samp[templateSampNum] + ( samp[templateSampNum+1] - samp[templateSampNum] )*(templateSampTime - templateSampNum);
  //scale by amplitude factor, this assumes template is pedestal subtracted
  sigVal = amp*sigVal;
  simVal += sigVal;
  return;
}

//Stupid TMinuit global variables/functions
static void ffer_fitFuncML(int& npar, double* gout, double& result, double par[], int flg);
void ffer_calcLnL(double par[], double& result);
double ffer_sampleErr;
std::vector<double> *ffer_fitData_vals;
std::vector<bool> *ffer_fitData_quality;
TemplateData *ffer_tempData;

class FitTemplateWf {
private:

public:

  FitTemplateWf();
  ~FitTemplateWf();

  void clearData();
  void addTemplate(const std::vector<double>& val, double period);
  void addData(const std::vector<double>& val, const std::vector<bool>& valQuality);
  void doFit();
  void setSampleError(double err);

  bool showOutput;
  double sampleErr;
  int status;
  TemplateData *tempData;
  std::vector<double> fitData_vals;
  std::vector<bool> fitData_quality;
  std::vector<double> initVals;
  std::vector<double> initErrs;
  std::vector<double> fitVals;
  std::vector<double> fitValErrs;
  std::vector<unsigned int> fixFitVars;
};

using namespace std;

FitTemplateWf::FitTemplateWf(){
  //initial values
  showOutput = 0;
  sampleErr = 1.;
  status = -1;
  tempData = new TemplateData();
}

FitTemplateWf::~FitTemplateWf(){
  delete tempData;
}

void FitTemplateWf::clearData(){
  fitData_vals.clear();
  fitData_quality.clear();
  tempData->samp.clear();
  initVals.clear();
  initErrs.clear();
  fitVals.clear();
  fitValErrs.clear();
  fixFitVars.clear();
}

void FitTemplateWf::setSampleError(double err){
  sampleErr = 1.;
  if(err <= 0. )
    return;
  sampleErr = err;
  return;
}

void FitTemplateWf::addTemplate(const std::vector<double>& val, double period){
  tempData->addTemplate(val,period);
  return;
}

void FitTemplateWf::addData(const std::vector<double>& val, const std::vector<bool>& valQuality){
  fitData_vals.clear();
  fitData_quality.clear();
  fitData_vals = val;
  fitData_quality = valQuality;
  return;
}

//wrapper function for TMinuit
void FitTemplateWf::doFit(){
  //sanity checks
  if( fitData_vals.size() == 0 || tempData->samp.size() == 0 || initVals.size() == 0 || initErrs.size() == 0 ){
    std::cout << "Invalid # of data or template vectors" << std::endl;
    return;
  }

  //give the global variables to the class objects (very lame)
  ffer_sampleErr = sampleErr;
  ffer_fitData_vals = &fitData_vals;
  ffer_fitData_quality = &fitData_quality;
  ffer_tempData = tempData;

  //initialize variables
  status = -1;
  unsigned int numParameters = 3;
  TMinuit *minimizer = new TMinuit(numParameters);

  //Set print level , -1 = suppress, 0 = info
  minimizer->SetPrintLevel(-1);
  if( showOutput == 1 )
    minimizer->SetPrintLevel(3);

  //define fit parameters
  minimizer->SetFCN(ffer_fitFuncML);
  minimizer->DefineParameter(0, "Offset", initVals[0], initErrs[0],0,0);
  minimizer->DefineParameter(1, "Time", initVals[1], initErrs[1],0,0);
  minimizer->DefineParameter(2, "Amp", initVals[2], initErrs[2],0,0);

  //optionally fix parameters
  for( unsigned int i = 0 ; i <  fixFitVars.size() ; i++ ){
    if( fixFitVars.at(i) < numParameters )
      minimizer->FixParameter( fixFitVars.at(i) );
  }

  //Set Minuit flags
  Double_t arglist[10];
  arglist[0] = 0.5; //what is this
  Int_t ierflg = 0;
  minimizer->mnexcm("SET ERR", arglist ,1,ierflg);  //command, arguments, # arguments, error flag

  //Do MIGRAD minimization
  Double_t tmp[1];
  tmp[0] = 100000;
  Int_t err;
  minimizer->mnexcm("MIG", tmp ,1,err);
  //minimizer->mnexcm("HES", tmp ,1,err);
  status = err;

  fitVals.clear();
  fitValErrs.clear();
  for(unsigned int i = 0 ; i < numParameters ; i++ ){
    double fitVal, fitValErr;
    minimizer->GetParameter(i, fitVal, fitValErr);
    fitVals.push_back( fitVal );
    fitValErrs.push_back( fitValErr );
  }

  delete minimizer; //memory leak?
  return;
}

//likelihood calc - note not included in class
void calcLnL(double par[], double& result){
  double diffSq = 0;
  double norm = 1./(ffer_sampleErr)/(ffer_sampleErr);

  //loop over data vector elements
  for(unsigned int num = 0 ; num < ffer_fitData_vals->size() ; num++ ){
    //skip over invalid data elements
    if( ffer_fitData_quality->operator[](num) == 0 )
      continue;

    //create fit hypothesis
    double simVal = 0;
    ffer_tempData->getSignalValue(num, par[0], par[1], par[2],simVal);
           
    //do lnL calc
    double dataVal = ffer_fitData_vals->operator[](num);
    diffSq = diffSq + (dataVal - simVal)*(dataVal - simVal)*norm; //gauss err assumed, noise is correlated but assume small
  }//end of element loop
  
  //calculate value to minimize
  result = -0.5*diffSq;
  return;
}

//Fit wrapper function - used by Minuit - has to be static void, annoying
static void ffer_fitFuncML(int& npar, double* gout, double& result, double par[], int flg){
  calcLnL(par, result);
  result = -1.*result;//Minuit is minimizing result ie maximizing LnL
  return;
}
