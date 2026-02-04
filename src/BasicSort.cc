#include <pthread.h>
#include <time.h>

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cctype>

#include <libpixie/pixie.hh>


struct CalibrationParameter {
    int channelID;
    double intercept;
    double slope;
};



  std::vector<std::vector<double>> readCalibrationParameters(const std::string &filename) {
    std::vector<std::vector<double>> parameters;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return parameters;
    }

    std::string line;
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        int channelID;
        double intercept, slope;
        
        // Skip non-numeric characters, extract three numbers per line
        if (iss >> channelID >> intercept >> slope) {
            if (channelID >= parameters.size()) {
                parameters.resize(channelID + 1);
            }
            parameters[channelID] = {intercept, slope};
	    std::cout<<intercept<<'\n';
        } else {
            std::cerr << "Warning: Unable to parse line: " << line << std::endl;
        }
    }

    file.close();
    return parameters;
}
int main(int argc, const char **argv) {
  if (argc < 3) { std::cout << "Useage: ./BasicSort input_files.txt output.root" << std::endl; return 0; }
  //set options for reader
  std::vector<std::string> dataPaths;

  // float a;
  // while (myfile >> a)
  // {
  //     printf("%f ", a);
  // }

   

   std::string filename = "calibration.csv"; // Replace with your file path
    std::vector<std::vector<double>> parameters = readCalibrationParameters(filename);
    std::cout<<"size: "<<parameters.size()<<'\n';
    // Print out the parameters to verify
    for (int i = 0; i < parameters.size(); ++i) {
        if (!parameters[i].empty()) {
            std::cout << "Channel ID: " << i
                      << ", Intercept: " << parameters[i][0]
                      << ", Slope: " << parameters[i][1] << std::endl;
        }
    }
  
    // for (const auto &param : parameters) {
    //     std::cout << "Channel ID: " << param.channelID
    //               << ", Intercept: " << param.intercept
    //               << ", Slope: " << param.slope << std::endl;
    // }
  //read files to sort
  std::ifstream files(argv[1]);
  std::string path;
  while (files >> path) {
    std::cout << path << std::endl; 
    dataPaths.push_back(path);
  }    
  files.close();

  for (int i = 0; i < dataPaths.size(); i++) {
    std::cout << dataPaths[i] << std::endl;
  }
  
  float coincWindow = 1500; //in ns
  double dither = 0.0;
  
  std::string defPath = "experimentdefinition.def";
  
  std::cout << argv[2] << std::endl;
  TFile file(argv[2], "recreate");
  //ROOT histogram definitions
  // "raw" histograms - on sub events before GaGG logic
  TH2F *raw_e = new TH2F("raw_e", "Raw Energy; Energy; Channel ID", 8192, 0, 8192, 16*13, 0, 16*13);  
  TH1I *raw_pu = new TH1I("raw_pu", "Raw Pileup; Channel ID; Counts", 16*13, 0, 16*13);
  TH1I *raw_dpu = new TH1I("raw_dpu", "Raw Pileup (doubles); Channel ID; Counts", 16*13, 0, 16*13);
  TH2I *raw_pu2 = new TH2I("raw_pu2", "Raw Pileup; Channel ID; Channel ID", 16*13, 0, 16*13, 16*13, 0, 16*13);
  TH1I *raw_or = new TH1I("raw_or", "Raw Out of Range; Channel ID; Counts", 16*13, 0, 16*13);
  TH1I *raw_tot = new TH1I("raw_tot", "Raw Total Counts; Channel ID; Counts", 16*13, 0, 16*13);
  TH2I *raw_mult = new TH2I("raw_mult", "Raw Subevent Multiplicity; Wall Time (s); Multiplicity", 4096, 0, 4096, 50, 0, 50);
  TH2F *raw_seqtdiff = new TH2F("raw_seqtdiff", "Raw sequential time difference; Time (us); Channel ID", 4000, 0, 400, 13*16, 0, 13*16);
  TH1F *raw_tot_seqtdiff = new TH1F("raw_tot_seqtdiff", "Raw sequential time difference; Time (us); Counts", 4000, 0, 400);
  TH2F *rates = new TH2F("rates", "Rates; Wall Time (s); Channel ID", 4096, 0, 4096, 16*13, 0, 16*13);  
  TH1F *rates_ep = new TH1F("rates_ep", "Rates; Wall Time (s); Counts/s", 4096, 0, 4096);
  TH1F *rates_sp = new TH1F("rates_sp", "Rates; Wall Time (s); Counts/s", 4096, 0, 4096);
  TH2F *hitpat = new TH2F("hitpat", "Hit pattern; ID1; ID2", 13*16, 0, 13*16, 13*16, 0, 13*16);
  TH2F *hitpat_db = new TH2F("hitpat_db", "Hit pattern; ID1; ID2", 13*16, 0, 13*16, 13*16, 0, 13*16);  
  TH2F *traces = new TH2F("traces", "Traces; Tick; Trace No.", 500, 0, 500, 1000, 0, 1000);  
  TH2F *trap_e = new TH2F("trap_e", "Trace Energy; Energy; Channel ID", 8192, 0, 16*4096, 16*13, 0, 16*13);
  TH2F *qdcE = new TH2F("qdcE", "QDC Energy; Energy; Channel ID", 8192, 0, 8192, 16*13, 0, 16*13);
  TH2F *qdcE_cfd = new TH2F("qdcE_cfd_trigger", "QDC Energy With CFD ; Energy; Channel ID", 8192, 0, 8192, 16*13, 0, 16*13);
  TH2F *qdcE_not_cfd = new TH2F("qdcE_not_cfd", "QDC Energy Without CFD; Energy; Channel ID", 8192, 0, 8192, 16*13, 0, 16*13);
  TH2D *qdcE_Vs_tdiff[16][16];
  TH2D *qdcE_Vs_tdiff_cfd[16][16];
  TH2D *qdcE_Vs_tdiff_not_cfd[16][16];
  TH1D* labr_tdiff[16];
  double centroids[16]={
    0,
-5.03292,
1.64954,
5.80224,
-8.54909,
0.131311,
0.491325,
1.06116,
-4.26671,
-2.72089,
-8.711,
-3.11793,
-39.216,
-17.9703,
-2.22478
  };

    for(int i=0; i<16; i++){
      labr_tdiff[i]=new TH1D(Form("tdiff_0_%d",i), Form("tdiff_0_%d",i), 3000, -150, 150);
      for(int j=i+1; j<16; j++){
	qdcE_Vs_tdiff[i][j]= new TH2D(Form("qdcE_Vs_tdiff_%d_%d", i, j), Form("QDC Energy Vs TDiff_LaBr3_%d_%d; TLaBr3_%d-TLaBr3_%d [ns]; LaBr3_%d QDC Energy ", i, j, i, j), 300, -150, 150, 8000, 0, 80000);
	qdcE_Vs_tdiff_cfd[i][j]	= new TH2D(Form("qdcE_Vs_tdiff_with_cfd_%d_%d", i, j), Form("QDC Energy Vs TDiff_LaBr3 With cfd_%d_%d; TLaBr3_%d-TLaBr3_%d [ns]; LaBr3_%d QDC Energy ", i, j, i, j), 300, -150, 150, 8000, 0, 80000);
	qdcE_Vs_tdiff_not_cfd[i][j]= new TH2D(Form("qdcE_Vs_tdiff_without_cfd_%d_%d", i, j), Form("QDC Energy Vs TDiff_LaBr3 Without cfd_%d_%d; TLaBr3_%d-TLaBr3_%d [ns]; LaBr3_%d QDC Energy ", i, j, i, j), 300, -150, 150, 8000, 0, 80000);
      }
    }
  int trace_chan = 1;
  int trace_mod = 12;
  int tracect = 0;
  
  long long unsigned int last_time[13*16] = {0};
  long long unsigned int last_time_tot = 0;
  
  int ds = 0;
  for (int file_num = 0; file_num < dataPaths.size(); ++file_num) {
    std::string listPath = dataPaths.at(file_num);
    PIXIE::Reader reader;
    //set up reader for file(s)
    std::cout << std::endl;
    std::cout << "Opening " << listPath << " for sorting" << std::endl;

    reader.openfile(listPath); //potential fix of changing open to openfile
    reader.definition.open(defPath);
    reader.definition.read();
    
    std::cout<< "coincidence window is " << reader.set_coinc(coincWindow) << std::endl;
    
    
    time_t now;
    time(&now);
    std::cout << now << std::endl;

    //THIS IS THE CODE TIME SENT TO TEST!!!!!
    //check if someone else is using disk
    while (true) {
      struct stat buffer;
      if (stat(".ioLock", &buffer)==0) {
	usleep(300000);
	continue;
      }
      else {
	break;
      }
    }
    std::ofstream lockFile(".ioLock");
    //lockFile.close();
    std::cout << "Populating buffer with entire file..." << std::endl;
    time_t pre_buffer;
    time(&pre_buffer);
    reader.loadbuffer();
    time_t post_buffer;
    time(&post_buffer);
    std::cout << "\nBuffer populated in " << post_buffer-pre_buffer << " s" << std::endl;
    std::remove(".ioLock");

    //TIM CODE BLOCK ENDS HERE!!!!!!!
    
    reader.start(); //this is for timing
    while (true) {
      ++ds;
      //std::cout << "Are we getting to this point " << std::endl;
      int retval = reader.read();
      //std::cout << retval << std::endl;
      if (retval) { 
        std::cout << "Exiting with retval " << retval << std::endl;
        break;
      }
      //std::cout << "Are we getting to this point " << std::endl;
      int eventIndx = reader.eventCtr-1;
      PIXIE::Event *e = &(reader.events[eventIndx]);
      //fill histograms
      //std::cout << "Are we getting to this point " << std::endl;
      int nMeas = e->nMeas;
      raw_mult->Fill(reader.measurements[e->fMeasurements[0]].eventTime/(3276.8*1e9), nMeas);
      long double TimeLaBr3[16]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}; 
      long double QDCELaBr3[16]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
      int CFDflag[16]={-20,20,-20,-20,-20,-20,-20,-20,-20,-20,-20,-20,-20,-20,-20,-20};
      bool pileup[2]={false, false};
      for (int i=0; i<nMeas; ++i) {
        auto &meas = reader.measurements[e->fMeasurements[i]];
        int chan = meas.channelNumber;
        int mod = meas.slotID-2;
        int tracelength = meas.traceLength;
        double energy = meas.eventEnergy;
        raw_tot->Fill(mod*16+chan);
	// std::cout<<mod*16+chan<<'\n';
	// if( parameters.at(mod*16+chan).channelID==mod*16+chan){
	int index = mod * 16 + chan;
	//std::cout << "index = " << index << std::endl;
	//std::cout << "parameters[" << index << "][0] = " << parameters[index][0] << ", parameters[" <<index << "][1] = " << parameters[index][1] << std::endl;
        raw_e->Fill( parameters[index][0]+energy*parameters[index][1], index);
 	// }

// // Ensure the index is within bounds
// if (index >= 0 && index < parameters.size()) {
//     // Ensure the channel ID matches
//     if (parameters.at(index).channelID == index) {
//         // Check if energy is properly initialized (assuming energy is a double)
//         if (raw_e) {  // Ensure raw_e is not nullptr
//             raw_e->Fill(parameters.at(index).intercept + energy * parameters.at(index).slope, index);
//         } else {
//             std::cerr << "Error: raw_e is not initialized." << std::endl;
//         }
//     } else {
//         std::cerr << "Channel ID mismatch for index " << index << std::endl;
//     }
// } else {
//     std::cerr << "Index out of bounds: " << index << std::endl;
// }
	trap_e->Fill(meas.trace_meas[3].datum, mod*16+chan);
        if (meas.finishCode) {
          raw_pu->Fill(mod*16+chan);
        }
        if (meas.outOfRange) {
          raw_or->Fill(mod*16+chan);
        }
        raw_seqtdiff->Fill((double)(meas.eventTime - last_time[mod*16+chan])/3276.8/1000.0, mod*16 + chan);
        raw_tot_seqtdiff->Fill((double)(meas.eventTime - last_time_tot)/3276.8/1000.);
        last_time[mod*16+chan] = meas.eventTime;
        last_time_tot = meas.eventTime;
        rates->Fill(meas.eventTime/(3276.8*1e9), mod*16+chan);

        if (tracect < 1000) {
        if (chan == trace_chan && mod == trace_mod) {
          tracect += 1;
          for (int k=0; k<tracelength; ++k) {
            traces->SetBinContent(k+1,tracect, meas.trace[k]);
          }
        }

        }

        if (mod == 12) {
          //std::cout << chan << "   " << meas.trace_meas[0].datum<< "   " << meas.trace_meas[3].datum << std::endl;
	  // trap_e->Fill(meas.trace_meas[3].datum, mod*16+chan);
	  qdcE->Fill((meas.QDCSums[1]-meas.QDCSums[0])/50, mod*16+chan);
	  if(meas.CFDForce==0)
	    qdcE_cfd->Fill((meas.QDCSums[1]-meas.QDCSums[0])/50, mod*16+chan);
	  if(meas.CFDForce!=0)
	    qdcE_not_cfd->Fill((meas.QDCSums[1]-meas.QDCSums[0])/50, mod*16+chan);

	  TimeLaBr3[chan]=meas.eventTime/3276.8;
	  QDCELaBr3[chan]=meas.QDCSums[1]-meas.QDCSums[0];
	  CFDflag[chan]=meas.CFDForce;
	  // if(chan==0){
	  //      TimeLaBr3[0]=meas.eventTime/3276.8;
	  //      QDCELaBr3[0]=meas.QDCSums[1]-meas.QDCSums[0];
	  //      CFDflag[0]=meas.CFDForce;
	  //      if(meas.finishCode)
	  // 	 pileup[0]=true;
	  //    }
	  //      if(chan==1){
	  //      TimeLaBr3[1]=meas.eventTime/3276.8;
	  //      QDCELaBr3[1]=meas.QDCSums[1]-meas.QDCSums[0];
	  //      CFDflag[1]=meas.CFDForce;
	  //      if(e->pileups)
	  // 	 pileup[1]=true;
	  //    }
	     
	  // qdc_slow->Fill(meas.qdcSlow.datum, mod*16+chan);
	  // qdc_tot->Fill(meas.qdcTotal.datum, mod*16+chan);
	  // qdc_pid->Fill(meas.qdcPID.datum, mod*16+chan);
        }
        for (int j=i+1; j<reader.events[eventIndx].nMeas; ++j) {
          auto &meas2 = reader.measurements[e->fMeasurements[j]];
          int chan2 = meas2.channelNumber;
          int mod2 = meas2.slotID-2;
          hitpat->Fill(mod*16+chan, mod2*16+chan2);

          if (meas.finishCode || meas2.finishCode) {
            raw_pu2->Fill(mod*16+chan, mod2*16+chan2);
            raw_pu2->Fill(mod2*16+chan2, mod*16+chan);
          }
        }
        if (reader.events[eventIndx].nMeas > 1) {
          if (meas.finishCode) {
            raw_dpu->Fill(mod*16+chan);
          }
        }
      }
      for(int l=0; l<16; l++){
	 if(TimeLaBr3[l]!=-1 && TimeLaBr3[0]!=-1 && TimeLaBr3[l]>10 &&TimeLaBr3[0]>10 && CFDflag[l]==0 &&CFDflag[0]==0)
	labr_tdiff[l]->Fill(TimeLaBr3[0]-TimeLaBr3[l]-centroids[l]);
      	for(int m=l+1; m<16; m++){
      	  if(TimeLaBr3[l]!=-1 && TimeLaBr3[m]!=-1 && TimeLaBr3[m]>10 &&TimeLaBr3[l]>10){
      	    qdcE_Vs_tdiff[l][m]->Fill(TimeLaBr3[l]+centroids[l]-TimeLaBr3[m]-centroids[m], QDCELaBr3[m]);
      	    if(CFDflag[l]==0 &&CFDflag[m]==0)
      	      qdcE_Vs_tdiff_cfd[l][m]->Fill(TimeLaBr3[l]+centroids[l]-TimeLaBr3[m]-centroids[m], QDCELaBr3[m]);
      	    if(CFDflag[l]!=0 && CFDflag[m]!=0)
      	      qdcE_Vs_tdiff_not_cfd[l][m]->Fill(TimeLaBr3[l]+centroids[l]-TimeLaBr3[m]-centroids[m], QDCELaBr3[m]);
      	  }
      	}
      }

      for(int l=0; l<16; l++){
        TimeLaBr3[l]=-1;
	QDCELaBr3[l]=-1;
	CFDflag[l]=-20;
      }
      if (ds%10000 == 0) {
        reader.printUpdate();
      }
    }
    time_t finaltime;
    time(&finaltime);
    std::cout << "Total fill time = " << finaltime-reader.starttime << " s " << std::endl;
    reader.printSummary();\
    //moved code block
    std::cout << "closing file" << std::endl;
    reader.closefile();
    std::cout << "clearing buffer" << std::endl;
    reader.clearbuffer();
    //moved code block
  }
  

  //SECOND TEST CODE BLOCK!!!!!!!!
  //std::cout << "closing file" << std::endl;
  //reader.closefile();
  //std::cout << "clearing buffer" << std::endl;
  //reader.clearbuffer();
  file.Write();
  file.Close();

  //END OF SECOND CODE BLOCK!!!!!!!!!!
  
  return 0;
}
