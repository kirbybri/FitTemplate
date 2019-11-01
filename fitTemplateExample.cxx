#include "FitTemplate.hxx"

void fitTemplateExample(){
	//declare fit object, ROOT objects
	FitTemplate *tempFit = new FitTemplate();
	TGraph* gData = new TGraph();
	TGraph* gFit = new TGraph();
	TCanvas* c0 = new TCanvas("c0", "c0",1400,800);
	TRandom3* rand = new TRandom3(0);
	tempFit->clearData();

	//define simulated data
	//vectors to hold samples, sample valid flag
	std::vector<double> dataVal;
	std::vector<bool> dataQuality;
        dataVal.push_back(130);
        dataQuality.push_back(1);
        dataVal.push_back(6);
        dataQuality.push_back(1);
        dataVal.push_back(1);
        dataQuality.push_back(1);
        dataVal.push_back(15);
        dataQuality.push_back(1);
        dataVal.push_back(4);
        dataQuality.push_back(1);
        dataVal.push_back(9);
        dataQuality.push_back(1);
        dataVal.push_back(7);
        dataQuality.push_back(1);
        tempFit->addData(dataVal, dataQuality);

	//define templates
	std::vector<double> tempVal;
        tempVal.push_back(6);
        tempVal.push_back(0.5);
        tempVal.push_back(0.34);
        tempVal.push_back(0.18);
        tempVal.push_back(0.1);
        tempVal.push_back(0.034);
        tempVal.push_back(0.23);
        tempFit->addTemplate(tempVal, "Almond");

	tempVal.clear();
        tempVal.push_back(4);
        tempVal.push_back(0);
        tempVal.push_back(0);
        tempVal.push_back(0);
        tempVal.push_back(0);
        tempVal.push_back(0);
        tempVal.push_back(0.8);
        tempFit->addTemplate(tempVal, "EggWhite");

	tempVal.clear();
        tempVal.push_back(2.75);
        tempVal.push_back(0);
        tempVal.push_back(0);
        tempVal.push_back(0.75);
        tempVal.push_back(0.075);
        tempVal.push_back(0.625);
        tempVal.push_back(0.025);
        tempFit->addTemplate(tempVal, "Dates");

	tempVal.clear();
        tempVal.push_back(3);
        tempVal.push_back(0.2);
        tempVal.push_back(0);
        tempVal.push_back(0.6);
        tempVal.push_back(0.4);
        tempVal.push_back(0.0);
        tempVal.push_back(0.2);
        tempFit->addTemplate(tempVal, "Cocoa");

	//get initial values for fitter
	tempFit->initVals.push_back(10);
	tempFit->initErrs.push_back(0.1);

	tempFit->initVals.push_back(10);
	tempFit->initErrs.push_back(0.1);

	tempFit->initVals.push_back(10);
	tempFit->initErrs.push_back(0.1);

	tempFit->initVals.push_back(10);
	tempFit->initErrs.push_back(0.1);

	//do fit
	tempFit->setSampleError(0.5);
	tempFit->showOutput = 1;
	tempFit->doFit();

	//load fit data into 
	TH1F *hFitData = new TH1F("hFitData","",7,-0.5,7-0.5);
	for(unsigned int s = 0 ; s < tempFit->fitTemplates.size() ; s++ ){
		for( unsigned int t = 0 ; t < tempFit->fitTemplates[s].val.size() ; t++ ){
			hFitData->Fill(t, tempFit->fitTemplates[s].val[t]*tempFit->fitVals[s] );
		}
	}

	//load data into tgraph
	gData->Set(0);
	for(unsigned int s = 0 ; s < dataVal.size() ; s++ ){
		gData->SetPoint( gData->GetN(), s , dataVal[s]);
	}

	//draw graphs
	gData->SetMarkerStyle(21);
	gData->SetMarkerColor(kBlack);
	gFit->SetLineColor(kBlue);
	gData->GetXaxis()->SetTitle("Sample Time");
	gData->GetYaxis()->SetTitle("Sample Value");

        auto legend = new TLegend(0.6,0.7,0.9,0.9);
        legend->AddEntry(gData,"Data","p");
        legend->AddEntry(gFit,"Fit","l");

	c0->Clear();
	gData->Draw("AP");
	//gFit->Draw("LP");
	//gSim->Draw("LP");
	hFitData->Draw("samehist");
        legend->Draw("same");
	c0->Update();

}
