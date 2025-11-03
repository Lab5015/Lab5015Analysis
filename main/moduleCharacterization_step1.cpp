#include "interface/TOFHIRThresholdZero.h"
#include "interface/AnalysisUtils.h"
#include "interface/FitUtils.h"
#include "interface/SetTDRStyle.h"
#include "interface/Na22SpectrumAnalyzer.h"
#include "interface/Na22SpectrumAnalyzerSingleBar.h"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <time.h>
#include <stdio.h>
#include <sys/stat.h>
#include <algorithm>
#include <iterator>
#include <dirent.h>

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRandom3.h"

#include <cmath>


int main(int argc, char** argv) {
    setTDRStyle();
    float cpu[2]{0}, mem[2]={0}, vsz[2]={0}, rss[2]={0};  
  

    if( argc < 2 ) {
       std::cout << ">>> moduleCharacterization_step1::usage:   " << argv[0] << " configFile.cfg" << std::endl;
       return -1;
    }
  
  
    //--- parse the config file
    CfgManager opts;
    opts.ParseConfigFile(argv[1]);
    int debugMode = 0;
    if( argc > 2 ) debugMode = atoi(argv[2]);
  

    //--- open files and make the tree chain
    std::string inputDir = opts.GetOpt<std::string>("Input.inputDir");
    std::string fileBaseName = opts.GetOpt<std::string>("Input.fileBaseName");
    std::string runs = opts.GetOpt<std::string>("Input.runs");
    int maxEntries = opts.GetOpt<int>("Input.maxEntries");
    int usePedestals = opts.GetOpt<int>("Input.usePedestals");
    std::string source = opts.GetOpt<std::string>("Input.sourceName");
    int useTrackInfo = opts.GetOpt<int>("Input.useTrackInfo");
	int useMCP = opts.GetOpt<int>("Input.useMCP");
    float my_step1 = opts.GetOpt<float>("Input.vov") ;

    int DUTasic = opts.GetOpt<int>("Channels.DUTasic");
    int REFasic = opts.GetOpt<int>("Channels.REFasic");
    std::string inputDirH4 = opts.GetOpt<std::string>("Input.inputDirH4");
    std::string discCalibrationFile = opts.GetOpt<std::string>("Input.discCalibration");
    TOFHIRThresholdZero thrZero(discCalibrationFile,1);

    TChain* tree = new TChain("data","data");
	TChain* treeH4 = new TChain("h4","h4");
  
  
    std::stringstream ss(runs); 
    std::string token;
    while( std::getline(ss,token,',') ) {
        std::stringstream ss2(token);
        std::string token2;
        int runMin = -1;
        int runMax = -1;
        while( std::getline(ss2,token2,'-') ) {
           if( runMin != -1 && runMax == -1 ) runMax = atoi(token2.c_str());
           if( runMin == -1 ) runMin = atoi(token2.c_str());
        }
        if( runMax == -1 ) runMax = runMin;
    
        for(int run = runMin; run <= runMax; ++run) {
            // -- analyze only spills at a chosen OV to speed up analysis
            // - list of files in run folder
            DIR *dir_ptr;
            struct dirent *diread;
            std::vector<std::string> filenames;
            std::string directory_path = Form("%s/%04d/",inputDir.c_str(),run);
            std::cout << directory_path.c_str()<<std::endl;

            if ((dir_ptr = opendir(directory_path.c_str())) != nullptr) {
                while ((diread = readdir(dir_ptr)) != nullptr) {
	                //std::cout << diread->d_name << std::endl;
	                std::string fname(diread->d_name);
	                filenames.push_back(fname);
	            }
	            closedir(dir_ptr);    
            }	

            for (auto fname: filenames) {
	            //std::cout << fname.c_str() << std::endl;;

	            if (fname == ".") continue;
	            if (fname == "..") continue;

	            // -- check if Vov selected
	            bool addFile = true;
	            TFile *f = TFile::Open((directory_path+fname).c_str());
	            TTree *tmpTree = f->Get<TTree>("data");
	            float step1;
	            tmpTree->SetBranchAddress("step1",&step1);
	            tmpTree->GetEntry(0);
	            if (my_step1 > 0 && step1!=my_step1) addFile = false;
	            delete tmpTree;
	            f->Close();
	
	            if (addFile){
	                std::cout << ">>> step1 = " << step1 << " --> Adding file: " << fname.c_str()<< std::endl;
	                tree->Add((directory_path+fname).c_str());
	            }
            }
			if(useMCP){
				std::string fileNameH4 = Form("%s/%04d.root",inputDirH4.c_str(),run);
			    std::cout <<">>> from MCP directory  --> Adding file: " << fileNameH4.c_str() << std::endl;
			    treeH4->Add(fileNameH4.c_str());
			}
            

      
            //
            //std::string fileName;
            ////if( !usePedestals ) fileName = Form("%s/%s%04d_*e.root",inputDir.c_str(),fileBaseName.c_str(),run); // pc-mtd-mib01
            //if( !usePedestals ) fileName = Form("%s/%04d/*_e.root",inputDir.c_str(),run); // pc-mtd-tb01 
            ////if( !usePedestals ) fileName = Form("%s/%s%05d_*e.root",inputDir.c_str(),fileBaseName.c_str(),run); // cmslpc
            //else                fileName = Form("%s/%04d/*ped_e.root",inputDir.c_str(),run);
            //std::cout << ">>> Adding file " << fileName << std::endl;
            //tree -> Add(fileName.c_str());
            //


           /* 
			std::cout<<Form("/data1/cmsdaq/tofhir2/h8/raw/%04d/",run)<<std::endl;
            struct stat t_stat;
            //stat(Form("/data/TOFHIR2/raw/run%04d.rawf",run), &t_stat);
            stat(Form("/data1/cmsdaq/tofhir2/h8/raw/%04d/",run), &t_stat);
            //stat(Form("/eos/uscms/store/group/cmstestbeam/2023_03_cmstiming_BTL/TOFHIR/RawData/%s%05d.rawf/",fileBaseName.c_str(),run), &t_stat); // cmslpc
            struct tm * timeinfo = localtime(&t_stat.st_mtime);
            std::cout << "Time and date of raw file of run" << run << ": " << asctime(timeinfo); */
        }
    }
    //--- define channels (read mapping from the configuration file)
    std::vector<unsigned int> channelMapping = opts.GetOpt<std::vector<unsigned int> >("Channels.channelMapping");
  
    int chL[16];
    int chR[16];
  
    for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar){

       chL[iBar] = channelMapping[iBar*2+0]+32*DUTasic;
       chR[iBar] = channelMapping[iBar*2+1]+32*DUTasic;


        //if(opts.GetOpt<int>("Channels.array")==0){
        //    chL[iBar] = channelMapping[iBar*2+0];
        //    chR[iBar] = channelMapping[iBar*2+1];
        //}
        //if(opts.GetOpt<int>("Channels.array")==1){
        //    chL[iBar] = channelMapping[iBar*2+0]+64;
        //    chR[iBar] = channelMapping[iBar*2+1]+64;
        //}
        std::cout << "Bar: " << iBar << "   chL: "<< chL[iBar] << "    chR: " <<chR[iBar] <<std::endl;
    }
  
    //--- define branches
    float step1, step2;
    int channelIdx[256];
    std::vector<unsigned short> *qfine = 0;
    std::vector<float> *tot = 0;
    std::vector<float> *energy = 0;
    std::vector<long long> *time = 0;
    std::vector<float>* qT1 = 0;
    std::vector<unsigned short>* t1fine = 0;
  
    int nhits;
	int iev_H4DAQ;
    int iev_TOFHIR;
    float x, y;
	unsigned long long t_TOFHIR,t_H4DAQ;

	int MCP_H4;
    int CFD_H4;
    int CLK_P_H4;
    int CLK_M_H4;
    int CLK_C_H4;
    float amp_max_H4[100];
    float time_H4[100];
    float fit_time_H4[100];
  
    tree -> SetBranchStatus("*",0);
    tree -> SetBranchStatus("step1",  1); tree -> SetBranchAddress("step1",  &step1);
    tree -> SetBranchStatus("step2",  1); tree -> SetBranchAddress("step2",  &step2);
  
    tree -> SetBranchStatus("channelIdx",  1); tree -> SetBranchAddress("channelIdx",  channelIdx);
  
    tree -> SetBranchStatus("qfine",  1); tree -> SetBranchAddress("qfine",   &qfine);  
    tree -> SetBranchStatus("tot",    1); tree -> SetBranchAddress("tot",       &tot);
    tree -> SetBranchStatus("energy", 1); tree -> SetBranchAddress("energy", &energy);
    tree -> SetBranchStatus("time",   1); tree -> SetBranchAddress("time",     &time);
  
    tree -> SetBranchStatus("qT1",       1); tree -> SetBranchAddress("qT1",      &qT1);
    tree -> SetBranchStatus("t1fine",    1); tree -> SetBranchAddress("t1fine",   &t1fine);
  
    if ( !opts.GetOpt<std::string>("Input.sourceName").compare("TB") &&  useTrackInfo ){
        tree -> SetBranchStatus("nhits_WC", 1); tree -> SetBranchAddress("nhits_WC",  &nhits);
        tree -> SetBranchStatus("x_WC", 1);     tree -> SetBranchAddress("x_WC",          &x);
        tree -> SetBranchStatus("y_WC", 1);     tree -> SetBranchAddress("y_WC",          &y);
    }

	if ( !opts.GetOpt<std::string>("Input.sourceName").compare("TB") && useMCP  ){
        tree -> SetBranchStatus("iev_H4DAQ", 1);  tree -> SetBranchAddress("iev_H4DAQ",  &iev_H4DAQ);
        tree -> SetBranchStatus("iev_TOFHIR", 1); tree -> SetBranchAddress("iev_TOFHIR",  &iev_TOFHIR);
        tree -> SetBranchStatus("t_H4DAQ", 1);    tree -> SetBranchAddress("t_H4DAQ",  &t_H4DAQ);
        tree -> SetBranchStatus("t_TOFHIR", 1);   tree -> SetBranchAddress("t_TOFHIR",  &t_TOFHIR);
      }
  
    if ( !opts.GetOpt<std::string>("Input.sourceName").compare("TB") && useMCP){ 
      //from H4DAQ ntuple
        treeH4 -> SetBranchStatus("MCP",1);      treeH4 -> SetBranchAddress("MCP",&MCP_H4);
        treeH4 -> SetBranchStatus("CFD",1);      treeH4 -> SetBranchAddress("CFD",&CFD_H4);
        treeH4 -> SetBranchStatus("CLK_P",1);    treeH4 -> SetBranchAddress("CLK_P",&CLK_P_H4);
        treeH4 -> SetBranchStatus("CLK_M",1);    treeH4 -> SetBranchAddress("CLK_M",&CLK_M_H4);
        treeH4 -> SetBranchStatus("CLK_C",1);    treeH4 -> SetBranchAddress("CLK_C",&CLK_C_H4);
        treeH4 -> SetBranchStatus("amp_max",1);  treeH4 -> SetBranchAddress("amp_max",amp_max_H4);
        treeH4 -> SetBranchStatus("time",1);     treeH4 -> SetBranchAddress("time",time_H4);
        treeH4 -> SetBranchStatus("fit_time",1); treeH4 -> SetBranchAddress("fit_time",fit_time_H4);
    }

    //--- get plot settings
    std::vector<float> Vov = opts.GetOpt<std::vector<float> >("Plots.Vov");
    std::vector<int> energyBins = opts.GetOpt<std::vector<int> >("Plots.energyBins");
    std::vector<int> energyMins = opts.GetOpt<std::vector<int> >("Plots.energyMins");
    std::vector<int> energyMaxs = opts.GetOpt<std::vector<int> >("Plots.energyMaxs");
  
    std::map<float,int> map_energyBins;
    std::map<float,int> map_energyMins;
    std::map<float,int> map_energyMaxs;
    for(unsigned int ii = 0; ii < Vov.size(); ++ii){
        map_energyBins[Vov[ii]] = energyBins[ii];
        map_energyMins[Vov[ii]] = energyMins[ii];
        map_energyMaxs[Vov[ii]] = energyMaxs[ii];
    }
  
    //float vetoEnergyThreshold = opts.GetOpt<float>("Plots.vetoEnergyThreshold");   


    // -- read minimum energy for each bar from file
    int vetoOtherBars = opts.GetOpt<int>("Cuts.vetoOtherBars");
    std::string minEnergiesFileName = opts.GetOpt<std::string>("Cuts.minEnergiesFileName");
    std::map < std::pair<int, float>, float> minE; 
    if (minEnergiesFileName != "") {
        std::ifstream minEnergiesFile;
        minEnergiesFile.open(minEnergiesFileName);
        std::string line;
        int bar;
        float ov;
        float value;
        while ( minEnergiesFile.good() ){
            getline(minEnergiesFile, line);
            std::istringstream ss(line);
            ss >> bar >> ov >> value; 
            minE[std::make_pair(bar,ov)] = value; 
            std::cout<< bar <<  "   " << ov << "  " << minE[std::make_pair(bar,ov)] <<std::endl;
        }
    }
    else{
        for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar){
            for(unsigned int ii = 0; ii < Vov.size(); ++ii){
	            minE[std::make_pair(iBar, Vov[ii])] = map_energyMins[Vov[ii]];
            }
        }
    }
  

    //--- define histograms
    std::string outFileName = opts.GetOpt<std::string>("Output.outFileNameStep1");
    TFile* outFile = TFile::Open(Form("%s",outFileName.c_str()),"RECREATE");
    outFile -> cd();
  
    std::map<int,TTree*> outTrees;

    std::map<int,TH1F*> h1_qfineL;
    std::map<int,TH1F*> h1_qfineR;  
    std::map<int,TH1F*> h1_totL;
    std::map<int,TH1F*> h1_totR;
    std::map<int,TH1F*> h1_energyL;
    std::map<int,TH1F*> h1_energyR;
    std::map<int,TH1F*> h1_energyLR;
    std::map<int,TH1F*> h1_energyLR_ext;
	std::map<int,TH1F*> h1_energyLR_PRE;
	std::map<int,TH1F*> h1_energyLR_POST;
	std::map<int,TH1F*> h1_energyLR_PREPOST;
	std::map<int,TH1F*> h1_energyLR_PRE_sum;
	std::map<int,TH1F*> h1_energyLR_POST_sum;
	std::map<int,TH1F*> h1_energyLR_PREPOST_sum;
	std::map<int,TH1F*> h1_energyL_POST;
	std::map<int,TH1F*> h1_energyR_POST;
	std::map<int,TH1F*> h1_energyL_PRE;
	std::map<int,TH1F*> h1_energyR_PRE;
	std::map<int,TH1F*> h1_energyL_POST_sum;
	std::map<int,TH1F*> h1_energyR_POST_sum;
	std::map<int,TH1F*> h1_energyL_PRE_sum;
	std::map<int,TH1F*> h1_energyR_PRE_sum;
	std::map<int,TH1F*> h1_counter;
    std::map<int,TH2F*> h2_xy_double_noCoin;
	std::map<int,TH2F*> h2_xy_double;
    std::map<int,TH2F*> h2_xy_REF;
	std::map<int,TH2F*> h2_xy_REF_cut;

	std::map<int,TH2F*> h2_tREF_vs_eREF;
	std::map<int, TProfile*> p1_tREF_vs_eREF;
	std::map<int,TH2F*> h2_tREF_vs_eREF2;
	std::map<int, TProfile*> p1_tREF_vs_eREF2;

    std::map<int,TH1F*> h1_energy_LR_PRE_triple;
	std::map<int,TH1F*> h1_energy_LR_POST_triple;
    std::map<int,TCanvas*> c;
    std::map<int,std::vector<float>*> rangesLR;
    std::map<int,bool> acceptEvent;


    std::map<int,TH1F*> h1_events_type;
    std::map<int,TH1F*> h1_events_type_REF;
    std::map<int,TH1F*> h1_events_type_3;	

	int nActiveBarsArray_ext;

    // -- Coincidence pre loop
    if( !opts.GetOpt<std::string>("Coincidence.status").compare("yes") &&
        ( !opts.GetOpt<std::string>("Input.sourceName").compare("Na22SingleBar") ||
	     !opts.GetOpt<std::string>("Input.sourceName").compare("Na22") ||
	     !opts.GetOpt<std::string>("Input.sourceName").compare("TB") ||
	     !opts.GetOpt<std::string>("Input.sourceName").compare("keepAll") ) ) 
	  {
        float energyL_ext;
        float energyR_ext;
        int chL_ext = opts.GetOpt<float>("Coincidence.chL");//NB: gli passo direttamente da cfg il ch barra 8 +64
        int chR_ext = opts.GetOpt<float>("Coincidence.chR");//NB: gli passo direttamente da cfg il ch barra 8 +64
      
        int nEntries = tree->GetEntries();
        if( maxEntries > 0 ) nEntries = maxEntries;
        for(int entry = 0; entry < nEntries; ++entry){
	        tree -> GetEntry(entry);
	        acceptEvent[entry] = false;
	        if( entry%200000 == 0 ){
	            std::cout << "\n>>> external bar loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << std::endl;
	            TrackProcess(cpu, mem, vsz, rss);
	        }
      
	        float Vov = step1;
            float vth1 = float(int(step2/10000)-1);
            float vth2 = float(int((step2-10000*(vth1+1))/100.)-1);
            float vth = 0.;
	        if(!opts.GetOpt<std::string>("Input.vth").compare("vth1"))  { vth = vth1;}
	        if(!opts.GetOpt<std::string>("Input.vth").compare("vth2"))  { vth = vth2;}

	        // select only one OV
	        if (my_step1 > 0  && my_step1 != step1) continue;

	        if (channelIdx[chL_ext] <0 || channelIdx[chR_ext] <0) continue;
	

	        // --- calculate energy sum/Nbars for module - useful to remove showering events/cross talk
	        if (!opts.GetOpt<std::string>("Input.sourceName").compare("TB")){

	            auto it = std::find( channelMapping.begin(), channelMapping.end(), chL_ext);       	
	            int thisBar = 0;
	            if (it != channelMapping.end()) thisBar = int((it - channelMapping.begin()))/2;

	            float energySumArray = 0.; 
	            int   nActiveBarsArray = 0;
	  
	            // check adjacent bars
	            /*for( int iBar = thisBar - 2; iBar < thisBar+3; ++iBar) {                                                                                                                            
	                if (iBar<0 || iBar>15 ) continue;                                                                                                                                                   
	                if (iBar == thisBar) continue;
	                //int chL_iext = channelMapping[iBar*2+0];// array0 for coincidence is hard coded... - to be fixed
	                //int chR_iext = channelMapping[iBar*2+1];// array0 for coincidence is hard coded... - to be fixed
	                int chL_iext = channelMapping[iBar*2+0]+64;// array1 for coincidence is hard coded... - to be fixed
	                int chR_iext = channelMapping[iBar*2+1]+64;// array1 for coincidence is hard coded... - to be fixed
	                float energyL_iext = (*energy)[channelIdx[chL_iext]];              
	                float energyR_iext = (*energy)[channelIdx[chR_iext]]; 
	                float totL_iext    = 0.001*(*tot)[channelIdx[chL_iext]];              
	                float totR_iext    = 0.001*(*tot)[channelIdx[chR_iext]]; 
	                if ( totL_iext > 0 && totL_iext < 100 && totR_iext > 0 && totR_iext < 100   ){
	                    float energyMean=(energyL_iext+energyR_iext)/2;
	                    if (energyMean>200 && energyMean < 1024){
		                    energySumArray+=energyMean;
		                    nActiveBarsArray+=1;
	                    }
	                }
	            }
	            if (nActiveBarsArray > 0 ) continue; // veto signals in the adjacent bars
	            */


	            for(int iBar = 0; iBar < int(channelMapping.size())/2; ++iBar) {
		        int chL_iext = channelMapping[iBar*2+0] + REFasic*32;
                        int chR_iext = channelMapping[iBar*2+1] + REFasic*32;
    
	               // int chL_iext = channelMapping[iBar*2+0];// module under test is array1, coincidence channel in array0 
	               // int chR_iext = channelMapping[iBar*2+1];// module under test is array1, coincidence channel in array0 
	               // if (opts.GetOpt<int>("Channels.array")==0) {
	               //     int chL_iext = channelMapping[iBar*2+0]+64;// module under test is array0, coincidence channel in array1
	               //     int chR_iext = channelMapping[iBar*2+1]+64;// module under test is array0, coincidence channel in array1
	               // }
	                float energyL_iext = (*energy)[channelIdx[chL_iext]];              
	                float energyR_iext = (*energy)[channelIdx[chR_iext]]; 
	                float totL_iext    = 0.001*(*tot)[channelIdx[chL_iext]];              
	                float totR_iext    = 0.001*(*tot)[channelIdx[chR_iext]]; 
	                if ( totL_iext > -10 && totL_iext < 100 && totR_iext > -10 && totR_iext < 100   ){
	                    float energyMean=(energyL_iext+energyR_iext)/2;
	                    if (energyMean>0){
		                    energySumArray+=energyMean;
		                    nActiveBarsArray+=1;
	                    }
	                }
	            }
	            //if (nActiveBarsArray > 5 ) continue;
	            if (nActiveBarsArray > 3 ) continue;
				nActiveBarsArray_ext = nActiveBarsArray;
	        }
	
	        energyL_ext = (*energy)[channelIdx[chL_ext]];
	        energyR_ext = (*energy)[channelIdx[chR_ext]];
	
	        int index( (10000*int(Vov*100.)) + (100*vth) + 99 );
	
	        //--- create histograms, if needed
	        if( h1_energyLR_ext[index] == NULL ){
	            c[index] = new TCanvas(Form("c1_Vov%.2f_th%02.0f",Vov,vth), Form("c1_Vov%.2f_th%02.0f",Vov,vth));
	            c[index] -> cd();
	            c[index] ->SetLogy();
	            h1_energyLR_ext[index] = new TH1F(Form("h1_energy_external_barL-R_Vov%.2f_th%02.0f",Vov,vth),"",map_energyBins[Vov],map_energyMins[Vov],map_energyMaxs[Vov]);
			    h1_events_type_REF[index] = new TH1F(Form("h1_events_type_REF_Vov%.2f_th%02.0f",Vov,vth),"",16,0.5,16.5);
				if(useTrackInfo){
					h2_xy_REF[index] = new TH2F(Form("h2_trackInfo_REF_Vov%.2f_th%02.0f",Vov,vth),Form("h2_trackInfo_REF_Vov%.2f_th%02.0f",Vov,vth),100,-100,100,100,-100,100);
				}
	        }
	
	        acceptEvent[entry] = true;
	        h1_events_type_REF[index] ->Fill(nActiveBarsArray_ext);
	        h1_energyLR_ext[index] -> Fill(0.5*(energyL_ext + energyR_ext));
            
			float totL_ext    = 0.001*(*tot)[channelIdx[chL_ext]];              
	        float totR_ext    = 0.001*(*tot)[channelIdx[chR_ext]];
			if(useTrackInfo && nhits>0 && x>-100 && y>-100){
				if(totL_ext >-10 && totL_ext < 100 && totR_ext >-10 && totR_ext < 100){
				    h2_xy_REF[index] -> Fill(x,y);
			    } 
		    }  
	
        }
    
        std::cout << std::endl;
      
        if (!opts.GetOpt<std::string>("Input.sourceName").compare("Na22SingleBar") || !opts.GetOpt<std::string>("Input.sourceName").compare("Na22")){		
	        for( auto index : h1_energyLR_ext){
	            rangesLR[index.first] = new std::vector<float>;
	            rangesLR[index.first]->push_back(30);
	            rangesLR[index.first]->push_back(950);
	  
	            if(!opts.GetOpt<std::string>("Input.sourceName").compare("Na22SingleBar")) Na22SpectrumAnalyzerSingleBar(index.second,rangesLR[index.first]);
	            if(!opts.GetOpt<std::string>("Input.sourceName").compare("Na22")) Na22SpectrumAnalyzer(index.second,rangesLR[index.first]);
	        }
        }
      
        if (!opts.GetOpt<std::string>("Input.sourceName").compare("TB")){
	        for( auto index : h1_energyLR_ext){
	            rangesLR[index.first] = new std::vector<float>;


	            float Vov = float ((int(index.first /10000))/100.);
	            float vth1 = float(int((index.first-Vov*10000*100)/100.));
	            float vth2 = float(int((step2-10000*(vth1+1))/100.)-1);
                float vth = 0;
                if(!opts.GetOpt<std::string>("Input.vth").compare("vth1"))  { vth = vth1;}
                if(!opts.GetOpt<std::string>("Input.vth").compare("vth2"))  { vth = vth2;}
	  
	            index.second->GetXaxis()->SetRangeUser(200,900);

	            float max = index.second->GetBinCenter(index.second->GetMaximumBin());
	            index.second->GetXaxis()->SetRangeUser(0,1024);
	  
	            TF1* f_pre = new TF1(Form("fit_energy_coincBar_Vov%.2f_vth1_%02.0f",Vov,vth), "[0]*TMath::Landau(x,[1],[2])", 0, 1000.); 
	            f_pre -> SetRange(max*0.85, max*1.4);
	            f_pre -> SetLineColor(kBlack);
                    f_pre -> SetLineWidth(2);
                    f_pre -> SetParameters(index.second->Integral(index.second->GetMaximumBin(), index.second->GetNbinsX())/10, max, 0.1*max);
	            f_pre -> SetParLimits(1, 0, 9999);
	            f_pre -> SetParLimits(2, 0, 9999);
	            index.second->Fit(f_pre, "QRS+");      
	  
	            if (f_pre->GetParameter(1)>20) rangesLR[index.first] -> push_back( 0.7*f_pre->GetParameter(1));
	            else   rangesLR[index.first] -> push_back( 20 );
	  			//rangesLR[index.first] -> push_back( 350 );
	  			//rangesLR[index.first] -> push_back( 600 );
	            rangesLR[index.first] -> push_back( 600 );
	  
	            std::cout << "Vov = " << Vov << "  vth1 = " << vth1 << "   vth2 = " << vth2 
		                  << "    Coincidence bar - energy range:  " << rangesLR[index.first]->at(0) << " - " << rangesLR[index.first]->at(1)<< std::endl;
	        }
        }
    }
  
  
  
    //------------------------
    //--- 1st loop over events
  
    ModuleEventClass anEvent;
    // DUT 
    unsigned short qfineL[16];
    unsigned short qfineR[16];    
    float totL[16];
    float totR[16];
    long long timeL[16];
    long long timeR[16];
    unsigned short t1fineL[16]; 
    unsigned short t1fineR[16]; 
    float qT1L[16]; 
    float qT1R[16]; 
    float energyL[16];
    float energyR[16];
    // REF
	long long timeL_ext;
	long long timeR_ext;
	float energyL_ext;
	float energyR_ext;
	float totL_ext;              
	float totR_ext;
	unsigned short t1fineL_ext; 
    unsigned short t1fineR_ext; 
  
	//--- my counters for debugging
	int mio_counter=0;
	int pre_counter=0;
	int post_counter=0;
	int prepost_counter=0;
	int single_counter=0;
	int bar0_counter=0;
	int bar15_counter=0;
	int multi_counter=0;

	int mio_counter7=0;
	int pre_counter7=0;
	int post_counter7=0;
	int prepost_counter7=0;
	int single_counter7=0;

	int mio_counter11=0;
	int pre_counter11=0;
	int post_counter11=0;
	int prepost_counter11=0;
	int single_counter11=0;
    
	int mio_counter15=0;
	int pre_counter15=0;
	int post_counter15=0;
	int prepost_counter15=0;
	int single_counter15=0;
    
	int pre_post_diff=0;
	int pre_post_diff_old=0;
	bool print= false;
	int xy_counter=0;
	int xy_counter0=0;

	int c_entry7=0;
	int c_entry11=0;
	int c_entry15=0;
	int count_ref=0;
	//--- 
	int nEntries = tree->GetEntries();
    if( maxEntries > 0 ) nEntries = maxEntries;
    for(int entry = 0; entry < nEntries; ++entry) {
		
        tree -> GetEntry(entry);
        if( entry%200000 == 0 ) {
	        std::cout << "\n>>> 1st loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << std::endl;
	        TrackProcess(cpu, mem, vsz, rss);
        }

        if (useTrackInfo && nhits > 0 &&  (x < -100 || y < -100 ) ) continue;
		if( useTrackInfo && nhits > 0 && x>=-100 && y>=-100) xy_counter0++; 

		//Get corresponding event in H4DAQ ntuple
	    if (useMCP && iev_H4DAQ>-1){
	        std::cout << "Getting H4 " << iev_H4DAQ << "," << iev_TOFHIR << "," << t_H4DAQ << "," << t_TOFHIR << std::endl;
	        treeH4->GetEntry(iev_H4DAQ);
	        std::cout << "time MCP " << time_H4[1] << std::endl;
	    }
        
        float Vov = step1;
        float vth1 = float(int(step2/10000)-1);
        float vth2 = int((step2-10000*(vth1+1))/100.)-1;
        float vth = 0;
        std::string vthMode = opts.GetOpt<std::string>("Input.vth");
        if(!opts.GetOpt<std::string>("Input.vth").compare("vth1"))  { vth = vth1;}
        if(!opts.GetOpt<std::string>("Input.vth").compare("vth2"))  { vth = vth2;}
        // float vthe = float(int((step2-10000*vth1-step2-100*vth2)/1)-1);
        if(vth==7){c_entry7++;}
		else if (vth==11){c_entry11++;}
		else {c_entry15++;}

        // select only one OV
        if (my_step1 > 0  && my_step1 != step1) continue;
        // --- check coincidence with another channel 
        if(!opts.GetOpt<std::string>("Coincidence.status").compare("yes")) {
	        if(!acceptEvent[entry] ) continue;

	        int chL_ext = opts.GetOpt<int>("Coincidence.chL");
	        int chR_ext = opts.GetOpt<int>("Coincidence.chR");
	        /*float*/ energyL_ext = (*energy)[channelIdx[chL_ext]];
	        /*float*/ energyR_ext = (*energy)[channelIdx[chR_ext]];
	
	        int label = (10000*int(Vov*100.)) + (100*vth) + 99;
	        int eBin = opts.GetOpt<int>("Coincidence.peak511eBin");
	        float avEn = 0.5 * ( energyL_ext + energyR_ext);
	        if ( (!opts.GetOpt<std::string>("Input.sourceName").compare("Na22SingleBar") || !opts.GetOpt<std::string>("Input.sourceName").compare("Na22")) &&  avEn > rangesLR[label]-> at(eBin)) {
	            continue;
	        }
	
	        if ( (!opts.GetOpt<std::string>("Input.sourceName").compare("TB")) && ( avEn < rangesLR[label]-> at(0) || avEn > rangesLR[label]-> at(1) ) ) {
	            continue;
	        } 

			timeL_ext = (*time)[channelIdx[chL_ext]];
			timeR_ext = (*time)[channelIdx[chR_ext]];
			t1fineL_ext = (*t1fine)[channelIdx[chL_ext]];
			t1fineR_ext = (*t1fine)[channelIdx[chR_ext]];
			int indexx( (10000*int(Vov*100.)) + (100*vth) + 99 );
			if(useTrackInfo && h2_xy_REF_cut[indexx] == NULL) h2_xy_REF_cut[indexx] = new TH2F(Form("h2_trackInfo_REF_cuts_Vov%.2f_th%02.0f",Vov,vth),Form("h2_trackInfo_REF_cuts_Vov%.2f_th%02.0f",Vov,vth),100,-100,100,100,-100,100);
			totL_ext    = 0.001*(*tot)[channelIdx[chL_ext]];              
	        totR_ext    = 0.001*(*tot)[channelIdx[chR_ext]];
			if(useTrackInfo && nhits>0 && x>-100 && y>-100){
				if(totL_ext >-10 && totL_ext < 100 && totR_ext >-10 && totR_ext < 100){
				    h2_xy_REF_cut[indexx] -> Fill(x,y);
			    } 
		    }
			
        }
		if(print){std::cout<<"\n----------Entry: "<<entry<<"-----------"<<std::endl;}
        for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar){
			if (channelIdx[chL[iBar]] >=0 && channelIdx[chR[iBar]] >=0){
				qfineL[iBar]=(*qfine)[channelIdx[chL[iBar]]];
				qfineR[iBar]=(*qfine)[channelIdx[chR[iBar]]];
				totL[iBar]=0.001*(*tot)[channelIdx[chL[iBar]]];
				totR[iBar]=0.001*(*tot)[channelIdx[chR[iBar]]];
				energyL[iBar]=(*energy)[channelIdx[chL[iBar]]];
				energyR[iBar]=(*energy)[channelIdx[chR[iBar]]];
				timeL[iBar]=(*time)[channelIdx[chL[iBar]]];
				timeR[iBar]=(*time)[channelIdx[chR[iBar]]];
				t1fineL[iBar]=(*t1fine)[channelIdx[chL[iBar]]];
				t1fineR[iBar]=(*t1fine)[channelIdx[chR[iBar]]];
				if(print){
				  std::cout<<"Entry: "<<entry<<" Barra: "<<iBar<<" chL["<<iBar<<"]= "<<chL[iBar]<<" channelIdx["<<chL[iBar]<<"]= "<<channelIdx[chL[iBar]]<<"\nEnergyL: "<<energyL[iBar]<<std::endl;
				}

			}
			else{
				qfineL[iBar]=-10;
				qfineR[iBar]=-10;
				totL[iBar]=-9999;
				totR[iBar]=-9999;
				energyL[iBar]=-10;
				energyR[iBar]=-10;
				timeL[iBar]=-10;
				timeR[iBar]=-10;
				t1fineL[iBar]=-10;
				t1fineR[iBar]=-10;
			} 
			if(print){
				std::cout<<"\ntotL["<<iBar<<"]: "<<totL[iBar]<<"\t   totR["<<iBar<<"]: "<<totR[iBar]<<std::endl;
			    std::cout<<"energyL["<<iBar<<"]: "<<energyL[iBar]<<"\t   energyR["<<iBar<<"]: "<<energyR[iBar]<<std::endl;
			}
			if(print){
				  std::cout<<"Entry: "<<entry<<" Barra: "<<iBar<<"\nchL["<<iBar<<"]= "<<chL[iBar]<<" channelIdx["<<chL[iBar]<<"]= "<<channelIdx[chL[iBar]]<<"  totL["<<iBar<<"]: "<<totL[iBar]<<"  energyL["<<iBar<<"]: "<<energyL[iBar]<<std::endl;
				}
        }// end loop over bars
        if(print) std::cout<<" "<<std::endl;
	    

        int maxEn=0;
        int maxBar=0;

        float energySumArray = 0;
        int   nActiveBarsArray = 0;
        int nBarsVeto[16];


	    int index2( int(Vov*10000) + vth );
		if(h1_events_type[index2]==NULL){
			h1_events_type[index2] = new TH1F(Form("h1_events_type_Vov%.2f_th%02.0f",Vov,vth),"",16,0.5,16.5);
			h1_events_type_3[index2] = new TH1F(Form("h1_events_type3_Vov%.2f_th%02.0f",Vov,vth),"",3,0.5,3.5);
		}

        for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar) {
	        nBarsVeto[iBar] = 0;
	
	        if (totL[iBar]>-10 && totR[iBar]>-10 && totL[iBar]<100 && totR[iBar]<100) {
	  
	            float energyMean=(energyL[iBar]+energyR[iBar])/2;
	            if (energyL[iBar]>0 && energyR[iBar]>0 && energyMean > 0){
	           	    //if (energyMean > 0){
	                energySumArray+=energyMean;
	                nActiveBarsArray+=1;
	            }
	    

	            // check energy in adjacent bars
	            for (int jBar = int(iBar) - 2; jBar < int(iBar) + 3; ++jBar){
	                if (jBar == int(iBar)) continue;
	                if (jBar < 0 || jBar > 15 ) continue;
	                if (totL[jBar]<-10 || totL[jBar]>100) continue;
	                if (totR[jBar]<-10 || totR[jBar]>100) continue;
	                float en = (energyL[jBar]+energyR[jBar])/2;
	                //if ( en > minE[std::make_pair(jBar, Vov)] && en<1024 ){
	                if ( en > minE[std::make_pair(jBar, Vov)] && minE[std::make_pair(jBar, Vov)]>1 && en<1024 ){
		                nBarsVeto[iBar]+=1;
	                }
	            }
	    
	            // find max bar
	            if(energyMean>maxEn){
	                maxEn = energyMean;
	                maxBar = iBar;
	            }
	        }
        }// end loop over bars
		h1_events_type[index2]->Fill(nActiveBarsArray);
		h1_events_type_3[index2]->Fill(nActiveBarsArray);
        
        double mean_en_sum=0;
		int active_counter=0;
	    int dead_counter=0;
		int loop_counter=0;
		pre_post_diff_old=abs(pre_counter+bar15_counter-post_counter-bar0_counter);
		float temp_energy1=0;
		float temp_energy2=0;
		int temp_entry1=-1;
		int temp_entry2=-1;
        for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar) {
            if (totL[iBar]>-10 && totR[iBar]>-10 && totL[iBar]<100 && totR[iBar]<100){  
                if(print){std::cout<<"entro in "<<iBar<<std::endl;}
	            int index( (10000*int(Vov*100.)) + (100*vth) + iBar );
	            //--- create histograms, if needed
	            if( h1_totL[index] == NULL ) {
	                h1_qfineL[index] = new TH1F(Form("h1_qfine_bar%02dL_Vov%.2f_th%02.0f",iBar,Vov,vth),"",512,-0.5,511.5);
	                h1_qfineR[index] = new TH1F(Form("h1_qfine_bar%02dR_Vov%.2f_th%02.0f",iBar,Vov,vth),"",512,-0.5,511.5);
	  
	                h1_totL[index] = new TH1F(Form("h1_tot_bar%02dL_Vov%.2f_th%02.0f",iBar,Vov,vth),"",400,-5.,35.);
	                h1_totR[index] = new TH1F(Form("h1_tot_bar%02dR_Vov%.2f_th%02.0f",iBar,Vov,vth),"",400,-5.,35.);
	  
	                h1_energyL[index] = new TH1F(Form("h1_energy_bar%02dL_Vov%.2f_th%02.0f",iBar,Vov,vth),"",map_energyBins[Vov],map_energyMins[Vov],map_energyMaxs[Vov]);
	                h1_energyR[index] = new TH1F(Form("h1_energy_bar%02dR_Vov%.2f_th%02.0f",iBar,Vov,vth),"",map_energyBins[Vov],map_energyMins[Vov],map_energyMaxs[Vov]);

					h1_energyL_PRE_sum[index] = new TH1F(Form("h1_energy_PRE_bar%02dL_Vov%.2f_th%02.0f",iBar,Vov,vth),"",1012,0,2024);
					h1_energyR_PRE_sum[index] = new TH1F(Form("h1_energy_PRE_bar%02dR_Vov%.2f_th%02.0f",iBar,Vov,vth),"",1012,0,2024);
					h1_energyL_POST_sum[index] = new TH1F(Form("h1_energy_POST_bar%02dL_Vov%.2f_th%02.0f",iBar,Vov,vth),"",1012,0,2024);
					h1_energyR_POST_sum[index] = new TH1F(Form("h1_energy_POST_bar%02dR_Vov%.2f_th%02.0f",iBar,Vov,vth),"",1012,0,2024);
					h1_energyL_PRE[index] = new TH1F(Form("h1_energy_PRE_ONLY_bar%02dL_Vov%.2f_th%02.0f",iBar,Vov,vth),"",1012,0,2024);
					h1_energyR_PRE[index] = new TH1F(Form("h1_energy_PRE_ONLY_bar%02dR_Vov%.2f_th%02.0f",iBar,Vov,vth),"",1012,0,2024);
					h1_energyL_POST[index] = new TH1F(Form("h1_energy_POST_ONLY_bar%02dL_Vov%.2f_th%02.0f",iBar,Vov,vth),"",1012,0,2024);
					h1_energyR_POST[index] = new TH1F(Form("h1_energy_POST_ONLY_bar%02dR_Vov%.2f_th%02.0f",iBar,Vov,vth),"",1012,0,2024);
	  
	                outTrees[index] = new TTree(Form("data_bar%02dL-R_Vov%.2f_th%02.0f",iBar,Vov,vth),Form("data_bar%02dL-R_Vov%.2f_th%02.0f",iBar,Vov,vth));
	                outTrees[index] -> Branch("event",&anEvent);

	                //cambio range per asse x lasciando la stessa larghezza bin(2)
	                h1_energyLR[index] = new TH1F(Form("h1_energy_bar%02dL-R_Vov%.2f_th%02.0f",iBar,Vov,vth),"",1012,0,2024);
                    // for bar 7 (as example) PRE means 7+6 , POST means 7+8, PREPOST means 6+7+8 
					h1_energyLR_PRE_sum[index]=  new TH1F(Form("h1_energy_PRE_bar%02dL-R_Vov%.2f_th%02.0f",iBar,Vov,vth),"",1012,0,2024);
					h1_energyLR_POST_sum[index]=  new TH1F(Form("h1_energy_POST_bar%02dL-R_Vov%.2f_th%02.0f",iBar,Vov,vth),"",1012,0,2024);
					h1_energyLR_PREPOST_sum[index]= new TH1F(Form("h1_energy_PREPOST_bar%02dL-R_Vov%.2f_th%02.0f",iBar,Vov,vth),"",1012,0,2024);
                    // this one refers to the energy seen by a single bar in case of PRE;POST and PREPOST events
					h1_energyLR_PRE[index]=  new TH1F(Form("h1_energy_PRE_ONLY_bar%02dL-R_Vov%.2f_th%02.0f",iBar,Vov,vth),"",1012,0,2024);
					h1_energyLR_POST[index]=  new TH1F(Form("h1_energy_POST_ONLY_bar%02dL-R_Vov%.2f_th%02.0f",iBar,Vov,vth),"",1012,0,2024);
					h1_energyLR_PREPOST[index]= new TH1F(Form("h1_energy_PREPOST_ONLY_bar%02dL-R_Vov%.2f_th%02.0f",iBar,Vov,vth),"",1012,0,2024);
					// histos for trackInfo
					h2_xy_double_noCoin[index] = new TH2F(Form("h2_xy_double_noCoin_bar%02d_Vov%.2f_th%02.0f",iBar,Vov,vth),"",100,-100,100,100,-100,100);
					h2_xy_double[index] = new TH2F(Form("h2_xy_double_bar%02d_Vov%.2f_th%02.0f",iBar,Vov,vth),"",100,-100,100,100,-100,100);
	            }
				if(h1_counter[index2]==NULL){
					// events counter (single double triple)
					h1_counter[index2] = new TH1F(Form("h1_counter_Vov%.2f_th%02.0f",Vov,vth),"",3,0.5,3.5);
					// time vs energy of external bar
					h2_tREF_vs_eREF[index2] = new TH2F(Form("h2_tREF_L_vs_eREF_Vov%.2f_th%02.0f",Vov,vth),"",100,100,450,300,-500,500);
					p1_tREF_vs_eREF[index2] = new TProfile(Form("p1_tREF_L_vs_eREF_Vov%.2f_th%02.0f",Vov,vth),"",50,100,450);
					h2_tREF_vs_eREF2[index2] = new TH2F(Form("h2_tREF_R_vs_eREF_Vov%.2f_th%02.0f",Vov,vth),"",100,100,450,300,-500,500);
					p1_tREF_vs_eREF2[index2] = new TProfile(Form("p1_tREF_R_vs_eREF_Vov%.2f_th%02.0f",Vov,vth),"",50,100,450);
				}

				if(h1_energy_LR_PRE_triple[index]== NULL){
					// in case of triple for bar 7 (as example) events PRE_triple refers to the energy seen by bar 6, POST_triple bar 8
					h1_energy_LR_PRE_triple[index] = new TH1F(Form("h1_energy_LR_PRE_triple_bar%02d_Vov%.2f_th%02.0f",iBar,Vov,vth),"",1012,0,2024);
					h1_energy_LR_POST_triple[index] = new TH1F(Form("h1_energy_LR_POST_triple_bar%02d_Vov%.2f_th%02.0f",iBar,Vov,vth),"",1012,0,2024);
				}
				
	
	            //--- fill histograms for each bar for Co60 & TB analysis
	            if( !opts.GetOpt<std::string>("Input.sourceName").compare("Co60") ||
	                !opts.GetOpt<std::string>("Input.sourceName").compare("Co60SumPeak") ||
	                !opts.GetOpt<std::string>("Input.sourceName").compare("TB") ) {

	                if( totL[iBar] <= -10. || totR[iBar] <= -10. ) continue;
	                if( totL[iBar] >= 50. ||  totR[iBar] >= 50.) continue;
	                if( ( thrZero.GetThresholdZero(chL[iBar],vthMode) + vth) > 63. ) continue;
                    if( ( thrZero.GetThresholdZero(chR[iBar],vthMode) + vth) > 63. ) continue;
                    if(energyL[iBar]<0 || energyR[iBar]<0 ) continue;
					//if (!opts.GetOpt<std::string>("Input.sourceName").compare("TB") && (energySumArray > 900 || nActiveBarsArray > 5)) continue; // to remove showering events
	                //if (!opts.GetOpt<std::string>("Input.sourceName").compare("TB") && (vetoOtherBars && nBarsVeto[iBar] > 0)) continue; // to remove showering events/ cross talk
	                //if (!opts.GetOpt<std::string>("Input.sourceName").compare("TB") && (vetoOtherBars && nActiveBarsArray > 5)) continue; // to remove showering events
	                //if (!opts.GetOpt<std::string>("Input.sourceName").compare("TB") && (vetoOtherBars && nActiveBarsArray > 3)) continue; // to remove showering events
	                int maxActiveBars = 3;
	                if (!opts.GetOpt<std::string>("Input.sourceName").compare("TB") && (vetoOtherBars && nActiveBarsArray > maxActiveBars)) continue; // to remove showering events
                    active_counter++;
                    int mio_check=0;//controlla se è un evento singolo, pre, post o pre-post e somma energy hist
					mean_en_sum=0;
					double timeDiff;
					temp_entry1=entry;
					anEvent.nClusters=-1;
					if((iBar > 0) && (iBar < 15) ){
						// triple events
						if( nActiveBarsArray==3 && (totL[iBar-1]>-10 && totR[iBar-1]>-10 && totL[iBar-1]<50. && totR[iBar-1]<50.)&&
							(totL[iBar+1]>-10 && totR[iBar+1]>-10 && totL[iBar+1]<50. && totR[iBar+1]<50.) &&
						    (energyL[iBar-1]>0 && energyR[iBar-1]>0 && energyL[iBar+1]>0 && energyR[iBar+1]>0 ) ){
								mio_check=3;//pre-post
								prepost_counter++;
								if(vth==7) {prepost_counter7++;}
								else if(vth==11) {prepost_counter11++;}
								else {prepost_counter15++;}
								if(print) {std::cout<<"\ntriple coincidence (n."<<prepost_counter<<") for entry: "<<entry<<" and bar: "<<iBar<<"\nmean_en_sum:  ";}
								for(int jBar=int(iBar)-1; jBar<int(iBar)+2; jBar++){
									mean_en_sum+=0.5*(energyL[jBar]+energyR[jBar]);
									if(print){std::cout<<" + "<<0.5*(energyL[jBar]+energyR[jBar]);}
								}
								if(print){std::cout<<"\nenergy wrote on histo: "<<mean_en_sum<<std::endl;}
								anEvent.energySum=mean_en_sum;
								anEvent.nClusters=3;
								h1_energyLR_PREPOST_sum[index] -> Fill(mean_en_sum);
								mean_en_sum=0;
								h1_energyLR_PREPOST[index] -> Fill(0.5*(energyL[iBar]+energyR[iBar]));
								h1_energy_LR_PRE_triple[index] -> Fill(0.5*(energyL[iBar-1]+energyR[iBar-1]));
								h1_energy_LR_POST_triple[index] -> Fill(0.5*(energyL[iBar+1]+energyR[iBar+1]));
							}	
						else if((iBar!=14)&&(totL[iBar+1]>-10 && totR[iBar+1]>-10 && totL[iBar+1]<50. && totR[iBar+1]<50.) && (totL[iBar+2]>-10 && totR[iBar+2]>-10 && totL[iBar+2]<50. && totR[iBar+2]<50.) &&
						        (energyL[iBar+1]>0 && energyR[iBar+1]>0 && energyL[iBar+2]>0 && energyR[iBar+2]>0 )) {
							loop_counter++;
							 continue;}
						//double events (post)	 
						else if(nActiveBarsArray==2 && totL[iBar+1]>-10 && totR[iBar+1]>-10 && totL[iBar+1]<50. && totR[iBar+1]<50. && (energyL[iBar+1]>0 && energyR[iBar+1]>0) ) { 
							mio_check=-1;//post	
							post_counter++;
							if(vth==7) {post_counter7++;}
							else if(vth==11) {post_counter11++;}
							else {post_counter15++;}
							if(print){std::cout<<"\ndouble post coincidence (n."<<post_counter<<") for entry: "<<entry<<" and bar: "<<iBar<<"\nmean_en_sum:  ";}
							for(int jBar=int(iBar); jBar<int(iBar)+2; jBar++){
								mean_en_sum+=0.5*(energyL[jBar]+energyR[jBar]);
								if(print){std::cout<<" + "<<0.5*(energyL[jBar]+energyR[jBar]);}
							}
							if(print){std::cout<<"\n energy wrote on histo: "<<mean_en_sum<<std::endl;}
							anEvent.energySum=mean_en_sum;
							anEvent.nClusters=2;
							h1_energyLR_POST_sum[index] -> Fill(mean_en_sum);
		                    h1_energyL_POST[index] -> Fill(energyL[iBar]);
							h1_energyR_POST[index] -> Fill(energyR[iBar]);
							h1_energyL_POST_sum[index] -> Fill(energyL[iBar]+energyL[iBar+1]);
							h1_energyR_POST_sum[index] -> Fill(energyR[iBar]+energyR[iBar+1]);
							mean_en_sum=0;
							h1_energyLR_POST[index] -> Fill(0.5*(energyL[iBar]+energyR[iBar]));
							timeDiff=0.5*(timeL[iBar]+timeR[iBar] -(timeL[iBar+1]+timeR[iBar+1]) );
							temp_energy1=energyL_ext;
							if(temp_energy2!=temp_energy1){
								h2_tREF_vs_eREF[index2] -> Fill(0.5*(energyL_ext+energyR_ext),0.5*(timeL_ext+timeR_ext) - timeL_ext);
							    p1_tREF_vs_eREF[index2] -> Fill(0.5*(energyL_ext+energyR_ext),0.5*(timeL_ext+timeR_ext) - timeL_ext);
								h2_tREF_vs_eREF2[index2] -> Fill(0.5*(energyL_ext+energyR_ext),0.5*(timeL_ext+timeR_ext) - timeR_ext);
							    p1_tREF_vs_eREF2[index2] -> Fill(0.5*(energyL_ext+energyR_ext),0.5*(timeL_ext+timeR_ext) - timeR_ext);
								temp_energy2=temp_energy1;
							}
							if(!opts.GetOpt<std::string>("Coincidence.status").compare("no") && useTrackInfo && nhits>0 && x>-100 && y>-100) {
								if(iBar==6){
									h2_xy_double_noCoin[index]-> Fill(x,y);
								}
							}
							if(!opts.GetOpt<std::string>("Coincidence.status").compare("yes") && useTrackInfo && nhits>0 && x>-100 && y>-100) {
								if(iBar==6){
									h2_xy_double[index]-> Fill(x,y);
								}
							}
						}
						else if((iBar!=1)&&(totL[iBar-2]>-10 && totR[iBar-2]>-10 && totL[iBar-2]<50. && totR[iBar-2]<50.)&&(totL[iBar-1]>-10 && totR[iBar-1]>-10 && totL[iBar-1]<50. && totR[iBar-1]<50.) && 
						        (energyL[iBar-2]>0 && energyR[iBar-2]>0 && energyL[iBar-1]>0 && energyR[iBar-1]>0 )) {
							loop_counter++;
							continue;}
						// double events (pre)
						else if(nActiveBarsArray==2 && totL[iBar-1]>-10 && totR[iBar-1]>-10 && totL[iBar-1]<50 && totR[iBar-1]<50 && energyL[iBar-1]>0 && energyR[iBar-1]>0 ) { 
							mio_check=2;//pre
							pre_counter++;
							if(vth==7) {pre_counter7++;}
							else if(vth==11) {pre_counter11++;}
							else {pre_counter15++;}
							if(print){std::cout<<"\ndouble pre coincidence (n."<<pre_counter<<") for entry: "<<entry<<" and bar: "<<iBar<<"\nmean_en_sum:  ";}
							for(int jBar=int(iBar)-1; jBar<int(iBar)+1; jBar++){
								mean_en_sum+=0.5*(energyL[jBar]+energyR[jBar]);
								if(print){std::cout<<" + "<<0.5*(energyL[jBar]+energyR[jBar]);}
							}
							if(print){std::cout<<"\nenergy wrote on histo: "<<mean_en_sum<<std::endl;}
							anEvent.energySum=mean_en_sum;
							anEvent.nClusters=4; // to distinguish from post events, i should rename the label to anEvent.Type  
							h1_energyLR_PRE_sum[index] -> Fill(mean_en_sum);
							h1_energyL_PRE[index] -> Fill(energyL[iBar]);
							h1_energyR_PRE[index] -> Fill(energyR[iBar]);
							h1_energyL_PRE_sum[index] -> Fill(energyL[iBar]+energyL[iBar-1]);
						    h1_energyR_PRE_sum[index] -> Fill(energyR[iBar]+energyR[iBar-1]);
							mean_en_sum=0;
							h1_energyLR_PRE[index] -> Fill(0.5*(energyL[iBar]+energyR[iBar]));
							timeDiff=0.5*(timeL[iBar-1]+timeR[iBar-1] -(timeL[iBar]+timeR[iBar]) );
						} 	
					}
					// check on the bars 0 and 15 for double events
					if(iBar==0 && (totL[iBar+1]>-10 && totR[iBar+1]>-10 && totL[iBar+1]<50 && totR[iBar+1]<50) && (totL[iBar+2]>-10 && totR[iBar+2]>-10 && totL[iBar+2]<50 && totR[iBar+2]<50) 
				    	&& energyL[iBar+1]>0 && energyR[iBar+1]>0 && energyL[iBar+2]>0 && energyR[iBar+2]>0) {
						loop_counter++;
						continue;}
					else if(nActiveBarsArray==2 && iBar==0 && (totL[iBar+1]>-10 && totR[iBar+1]>-10 && totL[iBar+1]<50 && totR[iBar+1]<50)  && energyL[iBar+1]>0 && energyR[iBar+1]>0 ) {
						mio_check=-1;//post
						post_counter++;
						if(vth==7) {post_counter7++;}
						else if(vth==11) {post_counter11++;}
						else {post_counter15++;}
						if(print){std::cout<<"\ndouble post 0 coincidence (n."<<bar0_counter<<") for entry: "<<entry<<" and bar: "<<iBar<<"\nmean_en_sum:  ";}
						for(int jBar=int(iBar); jBar<int(iBar)+2; jBar++){
							mean_en_sum+=0.5*(energyL[jBar]+energyR[jBar]);
							if(print){std::cout<<" + "<<0.5*(energyL[jBar]+energyR[jBar]);}
						}
						if(print){std::cout<<"\nenergy wrote on histo: "<<mean_en_sum<<std::endl;}
						anEvent.energySum=mean_en_sum;
						anEvent.nClusters=2;
						h1_energyLR_POST_sum[index] -> Fill(mean_en_sum);
						h1_energyL_POST[index] -> Fill(energyL[iBar]);
						h1_energyR_POST[index] -> Fill(energyR[iBar]);
						h1_energyL_POST_sum[index] -> Fill(energyL[iBar]+energyL[iBar+1]);
						h1_energyR_POST_sum[index] -> Fill(energyR[iBar]+energyR[iBar+1]);
						mean_en_sum=0;
						h1_energyLR_POST[index] -> Fill(0.5*(energyL[iBar]+energyR[iBar]));
						timeDiff=0.5*(timeL[iBar]+timeR[iBar] -(timeL[iBar+1]+timeR[iBar+1]) );
					}
					if(iBar==15 && (totL[iBar-1]>-10 && totR[iBar-1]>-10 && totL[iBar-1]<50 && totR[iBar-1]<50) &&  (totL[iBar-2]>-10 && totR[iBar-2]>-10 && totL[iBar-2]<50 && totR[iBar-2]<50) &&
					    (energyL[iBar-1]>0 && energyR[iBar-1]>0 && energyL[iBar-2]>0 && energyR[iBar-2]>0)  ) {
						loop_counter++; 
						continue;}
					else if(nActiveBarsArray==2 && iBar==15 && (totL[iBar-1]>-10 && totR[iBar-1]>-10 && totL[iBar-1]<50 && totR[iBar-1]<50) && energyL[iBar-1]>0 && energyR[iBar-1]>0 ) {
						mio_check=2;//pre
						pre_counter++;
						if(vth==7) {pre_counter7++;}
						else if(vth==11) {pre_counter11++;}
						else {pre_counter15++;}
						if(print){std::cout<<"\ndouble pre 15 coincidence (n."<<bar15_counter<<") for entry: "<<entry<<" and bar: "<<iBar<<"\nmean_en_sum:  ";}
						for(int jBar=int(iBar)-1; jBar<int(iBar)+1; jBar++){
							mean_en_sum+=0.5*(energyL[jBar]+energyR[jBar]);
							if(print){std::cout<<" + "<<0.5*(energyL[jBar]+energyR[jBar]);}
						}
						if(print){std::cout<<"\nenergy wrote on histo: "<<mean_en_sum<<std::endl;}
						anEvent.energySum=mean_en_sum;
						anEvent.nClusters=4; // to distinguish from post events, i should rename the label to anEvent.Type
						h1_energyLR_PRE_sum[index] -> Fill(mean_en_sum);
						h1_energyLR_PRE[index] -> Fill(0.5*(energyL[iBar]+energyR[iBar]));
						h1_energyL_PRE[index] -> Fill(energyL[iBar]);
						h1_energyR_PRE[index] -> Fill(energyR[iBar]);
						h1_energyL_PRE_sum[index] -> Fill(energyL[iBar]+energyL[iBar-1]);
						h1_energyR_PRE_sum[index] -> Fill(energyR[iBar]+energyR[iBar-1]);
						timeDiff=0.5*(timeL[iBar]+timeR[iBar] -(timeL[iBar-1]+timeR[iBar-1]) );
					}
                    //single event
					if(nActiveBarsArray==1 && (energyL[iBar]>0 && energyR[iBar]>0) ) {
							single_counter++;
							mio_check=1;
							if(vth==7) {single_counter7++;}
							else if(vth==11) {single_counter11++;}
							else {single_counter15++;}
							h1_energyLR[index] -> Fill(0.5*(energyL[iBar]+energyR[iBar]));
							anEvent.energySum=0.5*(energyL[iBar]+energyR[iBar]);
							anEvent.nClusters=1;
							if(print){std::cout<<"\nsingle coincidence (n."<<single_counter<<") for entry: "<<entry<<" and bar: "<<iBar<<"\n mean_en_sum: "<<0.5*(energyL[iBar]+energyR[iBar])<<std::endl;}
							h1_qfineL[index] -> Fill( qfineL[iBar] );
	                		h1_totL[index] -> Fill( totL[iBar]  );
	                		h1_energyL[index] -> Fill( energyL[iBar] );
	  
	                		h1_qfineR[index] -> Fill( qfineR[iBar] );
	                		h1_totR[index] -> Fill( totR[iBar] );
	                		h1_energyR[index] -> Fill( energyR[iBar] );
				    }

					h1_counter[index2]-> Fill(mio_check);
                    
	                anEvent.barID = iBar;
	                anEvent.Vov = Vov;
	                anEvent.vth1 = vth;
	                anEvent.energyL = energyL[iBar];
	                anEvent.energyR = energyR[iBar];
	                anEvent.totL = totL[iBar];
	                anEvent.totR = totR[iBar];
	                anEvent.timeL = timeL[iBar];
	                anEvent.timeR = timeR[iBar];
	                anEvent.t1fineL = t1fineL[iBar];
	                anEvent.t1fineR = t1fineR[iBar];
					

					anEvent.timeL_pre = -10;
					anEvent.timeR_pre = -10;
					anEvent.timeL_post = -10;
					anEvent.timeR_post = -10;

					anEvent.energyL_pre = -10;
					anEvent.energyR_pre= -10;
					anEvent.energyL_post = -10;
					anEvent.energyR_post= -10;
					
					anEvent.timeDiff= -99999;

					anEvent.t1fineL_post = -10;
	                anEvent.t1fineR_post = -10;
					anEvent.t1fineL_pre = -10;
	                anEvent.t1fineR_pre = -10;



					if(mio_check==3){
						anEvent.timeL_pre = timeL[iBar-1];
						anEvent.timeR_pre = timeR[iBar-1];
						anEvent.timeL_post = timeL[iBar+1];
						anEvent.timeR_post = timeR[iBar+1];

						anEvent.energyL_pre = energyL[iBar-1];
						anEvent.energyR_pre= energyR[iBar-1];
						anEvent.energyL_post = energyL[iBar+1];
						anEvent.energyR_post= energyR[iBar+1];

						anEvent.t1fineL_pre = t1fineL[iBar-1];
	                    anEvent.t1fineR_pre = t1fineR[iBar-1];
						anEvent.t1fineL_post = t1fineL[iBar+1];
	                    anEvent.t1fineR_post = t1fineR[iBar+1];
					}

					else if(mio_check==2){
						anEvent.timeL_pre = timeL[iBar-1];
						anEvent.timeR_pre = timeR[iBar-1];

						anEvent.energyL_pre = energyL[iBar-1];
						anEvent.energyR_pre = energyR[iBar-1];

						anEvent.t1fineL_pre = t1fineL[iBar-1];
	                    anEvent.t1fineR_pre = t1fineR[iBar-1];

						anEvent.timeDiff = timeDiff;
					}

					else if(mio_check==-1){
						anEvent.timeL_post = timeL[iBar+1];
						anEvent.timeR_post = timeR[iBar+1];

						anEvent.energyL_post = energyL[iBar+1];
						anEvent.energyR_post= energyR[iBar+1];

						anEvent.t1fineL_post = t1fineL[iBar+1];
	                    anEvent.t1fineR_post = t1fineR[iBar+1];

						anEvent.timeDiff = timeDiff ;
					}

	                if(useTrackInfo){
	                    anEvent.nhits = nhits;
	                    anEvent.x = x;
	                    anEvent.y = y;
	                }
	                else{
	                    anEvent.nhits = -1;
	                    anEvent.x = -999.;
	                    anEvent.y = -999.;
	                }

					// variable for REF 
                    if(mio_check!=0 && totL_ext>-10 && totR_ext>-10 && totL_ext<100 && totR_ext<100 && energyL_ext>0 && energyR_ext>0){
						anEvent.timeL_ext = timeL_ext;
					    anEvent.timeR_ext = timeR_ext;
					    anEvent.energyL_ext = energyL_ext;
					    anEvent.energyR_ext = energyR_ext;
					    anEvent.t1fineL_ext = t1fineL_ext;
					    anEvent.t1fineR_ext = t1fineR_ext;
						count_ref++;
					}
					else {
						anEvent.timeL_ext = -10;
					    anEvent.timeR_ext = -10;
					    anEvent.energyL_ext = -10;
					    anEvent.energyR_ext = -10;
					    anEvent.t1fineL_ext = -10;
					    anEvent.t1fineR_ext = -10;
					}

					if(useMCP && iev_H4DAQ>-1){
						std::cout << "Chek pre anEvent "<< std::endl;
		                anEvent.amp_MCP = amp_max_H4[MCP_H4];
		                anEvent.t_MCP = fit_time_H4[MCP_H4];//+CFD
			            anEvent.t_CFD_MCP = time_H4[MCP_H4+CFD_H4];
		                anEvent.t_CLK_P = fit_time_H4[CLK_P_H4];//prima CLK_P_H4
		                anEvent.t_CLK_M = fit_time_H4[CLK_M_H4];// 
			            //anEvent.t_CLK_C = fit_time_H4[CLK_C_H4];
						std::cout << "Chek post anEvent "<< std::endl;
		            }
		            else{
		                anEvent.amp_MCP = -999;
		                anEvent.t_MCP = -999;
		                anEvent.t_CLK_P = -999;
		                anEvent.t_CLK_M = -999;
		            }
    

	                outTrees[index] -> Fill();
					/*if(0.5*(anEvent.timeL+anEvent.timeR) - 0.5*(anEvent.timeL_post+anEvent.timeR_post) < 0.0001){
						std::cout<<"\nEntry: "<<entry<<" counter: "<<mio_counter<<" nHits: "<<anEvent.nClusters
						         <<"\nenergyL["<<iBar-1<<"]: "<<anEvent.energyL_pre<<"\tenergyR["<<iBar-1<<"]: "<<anEvent.energyR_pre
						         <<"\nenergyL["<<iBar<<"]: "<<anEvent.energyL<<"\tenergyR["<<iBar<<"]: "<<anEvent.energyR
								 <<"\nenergyL["<<iBar+1<<"]: "<<anEvent.energyL_post<<"\tenergyR["<<iBar+1<<"]: "<<anEvent.energyR_post
								 <<"\ntimeL["<<iBar-1<<"]: "<<anEvent.timeL_pre<<"\ttimeR["<<iBar-1<<"]: "<<anEvent.timeR_pre
						         <<"\ntimeL["<<iBar<<"]: "<<anEvent.timeL<<"\ttimeR["<<iBar<<"]: "<<anEvent.timeR
								 <<"\ntimeL["<<iBar+1<<"]: "<<anEvent.timeL_post<<"\ttimeR["<<iBar+1<<"]: "<<anEvent.timeR_post
								 <<"\ntime difference: "<< anEvent.timeDiff<<std::endl;
					}*/

                    
	            }	  
	        }
			else {dead_counter++;}
			loop_counter++;
     	}// -- end loop over bars
        // some print used for debugging
		pre_post_diff=abs(pre_counter+bar15_counter-post_counter-bar0_counter);
		//std::cout<<"\n iteration"<<mio_counter<<") new: "<<pre_post_diff<<"\n  post counter: "<<post_counter<<"  post 0 counter: "<<bar0_counter<<"\n  pre counter: "<<pre_counter<<"  pre 15 counter: "<<bar15_counter<<std::endl;
		if(pre_post_diff > pre_post_diff_old){
			std::cout<<"\n"<<mio_counter<<") strange event spotted for entry: "<<entry<<"\n  post counter: "<<post_counter<<"  post 0 counter: "<<bar0_counter<<"\n  pre counter: "<<pre_counter<<"  pre 15 counter: "<<bar15_counter<<"\n"<<pre_post_diff<<" "<<pre_post_diff_old<<std::endl;
		}
        if(loop_counter!=active_counter+dead_counter){std::cout<<" problem in loop bover bar for entry: "<<entry<<"  active bars: "<<active_counter<<"  dead bars: "<<dead_counter<<"  loop over "<<loop_counter<<" bars"<<std::endl;}
		else{
			if(print){std::cout<<"\n\nfor entry: "<<entry<<"  active bars: "<<active_counter<<"  dead bars: "<<dead_counter<<std::endl; } 
		}
		if(nActiveBarsArray>3){multi_counter++;}
        // --- for Na22 or Laser analysis use only the bar with max energy to remove cross-talk between adjacent bars
        if( !opts.GetOpt<std::string>("Input.sourceName").compare("Na22") ||
	        !opts.GetOpt<std::string>("Input.sourceName").compare("Na22SingleBar") ||
	        !opts.GetOpt<std::string>("Input.sourceName").compare("Laser") ||
	        !opts.GetOpt<std::string>("Input.sourceName").compare("keepAll") ) {
	        int index( (10000*int(Vov*100.)) + (100*vth) + maxBar );
	
	        if( totL[maxBar] <= -10. || totR[maxBar] <= -10. ) continue;
	        if( totL[maxBar] >= 50. ||  totR[maxBar] >= 50.) continue;
	        if( ( thrZero.GetThresholdZero(chL[maxBar],vthMode) + vth) > 63. ) continue;
            if( ( thrZero.GetThresholdZero(chR[maxBar],vthMode) + vth) > 63. ) continue;
	
	        //--- fill histograms
	        h1_qfineL[index] -> Fill( qfineL[maxBar] );
	        h1_totL[index] -> Fill( totL[maxBar] );
	        h1_energyL[index] -> Fill( energyL[maxBar] );
	
	        h1_qfineR[index] -> Fill( qfineR[maxBar] );
	        h1_totR[index] -> Fill( totR[maxBar] );
	        h1_energyR[index] -> Fill( energyR[maxBar] );
	  
	        h1_energyLR[index] -> Fill(0.5*(energyL[maxBar]+energyR[maxBar]));
	
	        anEvent.barID = maxBar;
	        anEvent.Vov = Vov;
	        anEvent.vth1 = vth;
	        anEvent.energyL = energyL[maxBar];
	        anEvent.energyR = energyR[maxBar];
	        anEvent.totL = totL[maxBar];
	        anEvent.totR = totR[maxBar];
	        anEvent.timeL = timeL[maxBar];
	        anEvent.timeR = timeR[maxBar];
	        anEvent.t1fineL = t1fineL[maxBar];
	        anEvent.t1fineR = t1fineR[maxBar];
	
	        if(useTrackInfo){
	            anEvent.nhits = nhits;
	            anEvent.x = x;
	            anEvent.y = y;
	        }
	        else{
	            anEvent.nhits = -1;
	            anEvent.x = -999.;
	            anEvent.y = -999.;
	        }
	        outTrees[index] -> Fill();

			
        }
		if(dead_counter==loop_counter) continue;
		mio_counter++;
		if(vth==7) {mio_counter7++;}
		else if(vth==11) {mio_counter11++;}
		else {mio_counter15++;}
		/*if(useTrackInfo && nhits>0 && (x!=-999 && y!=-999)){
			std::cout<<"entry: "<<entry<<" counter: "<<mio_counter<<"\nx: "<<x<<"\ty: "<<y<<" nhits: "<< nhits<<std::endl;
			xy_counter++;
			std::cout<<"xy_counter: "<<xy_counter<<std::endl;
		}*/
    } // --- end loop over events
    
	
	
    int extrass=(mio_counter-(pre_counter+prepost_counter+single_counter));
	int extrass7=(mio_counter7-(pre_counter7+prepost_counter7+single_counter7));
	int extrass11=(mio_counter11-(pre_counter11+prepost_counter11+single_counter11));
	int extrass15=(mio_counter15-(pre_counter15+prepost_counter15+single_counter15));

    std::cout<<"\nRun "<<runs<<" Started with "<<nEntries<<" events.  Accepted "<<mio_counter<<" events (" << 100.*mio_counter/nEntries << "%)."<<std::endl; 
	std::cout<<"Of which \n"<<single_counter<<" (" << 100.*single_counter/mio_counter << "%) single events\n"<< post_counter<<" ("
	         << 100.*post_counter/mio_counter <<"%) [post] \t"<< pre_counter<<" ("<<100.*pre_counter/mio_counter <<"%) [pre]   double events\n"/*<<bar0_counter<<" [post 0] \t"<<bar15_counter<<" [pre 15]   double events\ndifference between all pre and post events: "<<pre_post_diff<<"\n"*/
			 <<prepost_counter<<" (" << 100.*prepost_counter/mio_counter << "%) triple events"<<std::endl;
    
	std::cout<<"\nExtra events = "<<extrass<<" (" << 100.* extrass/mio_counter << "%)\n"<<std::endl;
	std::cout<<"\nEvents in reference module: "<<count_ref<<" ("<<100.* count_ref/mio_counter<<"%)\n"<<std::endl;
    /*
	std::cout<<"events with th 7: "<<mio_counter7<<" (" << 100.*mio_counter7/mio_counter << "%).\n"<<single_counter7<<" (" << 100.*single_counter7/mio_counter7 << "%) single events\n"<< post_counter7<<" ("
	<< 100.*post_counter7/mio_counter7 <<"%) [post] \t"<< pre_counter7<<" ("<<100.*pre_counter7/mio_counter7 <<"%) [pre]   double events\n"
	<<prepost_counter7<<" (" << 100.*prepost_counter7/mio_counter7 << "%) triple events"<<std::endl;
    std::cout<<"\nExtra events = "<<extrass7<<" (" << 100.* extrass7/mio_counter7 << "%)\n"<<std::endl;

	std::cout<<"events with th 11: "<<mio_counter11<<" (" << 100.*mio_counter11/mio_counter << "%).\n"<<single_counter11<<" (" << 100.*single_counter11/mio_counter11 << "%) single events\n"<< post_counter11<<" ("
	<< 100.*post_counter11/mio_counter11 <<"%) [post] \t"<< pre_counter11<<" ("<<100.*pre_counter11/mio_counter11 <<"%) [pre]   double events\n"
	<<prepost_counter11<<" (" << 100.*prepost_counter11/mio_counter11 << "%) triple events"<<std::endl;
    std::cout<<"\nExtra events = "<<extrass11<<" (" << 100.* extrass11/mio_counter11 << "%)\n"<<std::endl;

	std::cout<<"events with th 15: "<<mio_counter15<<" (" << 100.*mio_counter15/mio_counter << "%).\n"<<single_counter15<<" (" << 100.*single_counter15/mio_counter15 << "%) single events\n"<< post_counter15<<" ("
	<< 100.*post_counter15/mio_counter15 <<"%) [post] \t"<< pre_counter15<<" ("<<100.*pre_counter15/mio_counter15 <<"%) [pre]   double events\n"
	<<prepost_counter15<<" (" << 100.*prepost_counter15/mio_counter15 << "%) triple events"<<std::endl;
    std::cout<<"\nExtra events = "<<extrass15<<" (" << 100.* extrass15/mio_counter15 << "%)\n"<<std::endl;

    std::cout<<"Threshold 7 started with "<<c_entry7<<" events.  Accepted "<<mio_counter7<<" events (" << 100.*mio_counter7/c_entry7 << "%)."<<std::endl;
	std::cout<<"Threshold 11 started with "<<c_entry11<<" events.  Accepted "<<mio_counter11<<" events (" << 100.*mio_counter11/c_entry11 << "%)."<<std::endl;
	std::cout<<"Threshold 15 started with "<<c_entry15<<" events.  Accepted "<<mio_counter15<<" events (" << 100.*mio_counter15/c_entry15 << "%).\n"<<std::endl; 
    */
    std::cout<<"Total number of events with more than 3 active bars: "<<multi_counter<<" (" << 100.*multi_counter/mio_counter << "%).\n"<<std::endl;
    if(useTrackInfo){
		std::cout<<"Total number of events with xy info: "<<xy_counter0<<" (" << 100.*xy_counter0/nEntries << "%).\n"<<std::endl;
		std::cout<<"Total number of accepted events with xy info: "<<xy_counter<<" (" << 100.*xy_counter/mio_counter << "%).\n"<<std::endl;
	}
	
    int bytes = outFile -> Write();
    std::cout << "============================================"  << std::endl;
    std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
    std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
    std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
    std::cout << "============================================"  << std::endl;
}
