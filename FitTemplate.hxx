#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "TMinuit.h"

///////////STRUCTS FOR PULSE AND FIT CLASS/////////////////
struct FitTemplate_data{
  std::vector<double> val;
  std::vector<bool> valQuality;
};

struct FitTemplate_template{
  std::string name;
  std::vector<double> val;
};

///////////FIT CLASS/////////////////

//Stupid TMinuit global variables/functions
static void ffer_fitFuncML(int& npar, double* gout, double& result, double par[], int flg);
void ffer_calcLnL(double par[], double& result);
double ffer_sampleErr;
std::vector<FitTemplate_template> *ffer_templates;
std::vector<FitTemplate_data> *ffer_data;

class FitTemplate {
private:

public:

  FitTemplate();
  ~FitTemplate();

  void clearData();
  void addTemplate(const std::vector<double>& val, std::string name = "default");
  void addData(const std::vector<double>& val, const std::vector<bool>& valQuality);
  void doFit();
  void setSampleError(double err);

  bool showOutput;
  double sampleErr;
  double minLnL;
  int status;
  std::vector<FitTemplate_template> fitTemplates;
  std::vector<FitTemplate_data> fitData;
  std::vector<double> initVals;
  std::vector<double> initErrs;
  std::vector<double> fitVals;
  std::vector<double> fitValErrs;
  std::vector<unsigned int> fixFitVars;
};

using namespace std;

FitTemplate::FitTemplate(){
  //initial values
  showOutput = 0;
  sampleErr = 1.;
  minLnL = 0;
  status = -1;
}

FitTemplate::~FitTemplate(){
}

void FitTemplate::clearData(){
  fitTemplates.clear();
  fitData.clear();
  initVals.clear();
  initErrs.clear();
  fitVals.clear();
  fitValErrs.clear();
  fixFitVars.clear();
}

void FitTemplate::setSampleError(double err){
  sampleErr = 1.;
  if(err <= 0. )
    return;
  sampleErr = err;
  return;
}

void FitTemplate::addTemplate(const std::vector<double>& val, std::string name = "default"){
  FitTemplate_template fTemp;
  fTemp.val = val;
  fTemp.name = name;
  fitTemplates.push_back(fTemp);
  return;
}

void FitTemplate::addData(const std::vector<double>& val, const std::vector<bool>& valQuality){
  FitTemplate_data fTemp;
  fTemp.val = val;
  fTemp.valQuality = valQuality;
  fitData.push_back(fTemp);
  return;
}

//wrapper function for TMinuit
void FitTemplate::doFit(){
  //sanity checks
  if( fitData.size() == 0 || fitTemplates.size() == 0 || initVals.size() == 0 || initErrs.size() == 0 ){
    std::cout << "Invalid # of data or template vectors" << std::endl;
    return;
  }
  if( fitTemplates.size() != initVals.size() ){
    std::cout << "Invalid # of data or template vectors" << std::endl;
    return;
  }
  if( fitTemplates.size() != initErrs.size() ){
    std::cout << "Invalid # of data or template vectors" << std::endl;
    return;
  }

  //give the global variables to the class objects (very lame)
  ffer_sampleErr = sampleErr;
  ffer_templates = &fitTemplates;
  ffer_data = &fitData;

  //initialize variables
  status = -1;
  unsigned int numParameters = fitTemplates.size();
  //std::unique_ptr<TMinuit> minimizer (new TMinuit(numParameters) );
  TMinuit *minimizer = new TMinuit(numParameters);

  //Set print level , -1 = suppress, 0 = info
  minimizer->SetPrintLevel(-1);
  if( showOutput == 1 )
    minimizer->SetPrintLevel(3);

  //define fit parameters
  minimizer->SetFCN(ffer_fitFuncML);
  for(unsigned int parNum = 0 ; parNum < fitTemplates.size() ; parNum++ ){
    char name[200];
    memset(name,0,sizeof(char)*100 );
    sprintf(name,"ParNum_%i",parNum);
    minimizer->DefineParameter(parNum, fitTemplates[parNum].name.c_str(), initVals[parNum], initErrs[parNum],0,0);
  }

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

  //loop over data vectors
  for(unsigned int dataNum = 0 ; dataNum < ffer_data->size() ; dataNum++ ){
    //loop over data vector elements
    for(unsigned int num = 0 ; num < ffer_data->operator[](dataNum).val.size() ; num++ ){
      //skip over data elements
      if( ffer_data->operator[](dataNum).valQuality[num] == 0 )
        continue;

      //create fit hypothesis
      double fitVal = 0;
      //std::cout << "Element # " << num << std::endl;
      for(unsigned int tempNum = 0 ; tempNum < ffer_templates->size() ; tempNum++ ){
 	//std::cout << "\ttemp # " << tempNum << "\tval " << ffer_templates->operator[](tempNum).val[num] << "\tpar " << par[tempNum] << std::endl;
        fitVal += ffer_templates->operator[](tempNum).val[num]*par[tempNum];
      }
           
      //do lnL calc
      double dataVal = ffer_data->operator[](dataNum).val[num];
      diffSq = diffSq + (dataVal - fitVal)*(dataVal - fitVal)*norm; //gauss err assumed, noise is correlated but assume small
      //std::cout << "num " << num << "\tdata " << dataVal << "\tfit " << fitVal << "\tdiff " << diffSq<< std::endl;
    }//end of element loop
  }//end of data loop

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
