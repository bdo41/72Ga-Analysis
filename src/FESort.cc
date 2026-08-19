#include<time.h>
#include<string>
#include<cmath>

//Adding these to help with testing calibrations
#include <map>
#include <sstream>
#include <vector>

#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TF1.h"
//Testing Calibration block

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCutG.h"

#include <libpixie/pixie.hh>

#include <libClarionTrinity/Event.hh>
#include <libClarionTrinity/GAGG.hh>
#include <libClarionTrinity/Trinity.hh>
#include <libClarionTrinity/Clover.hh>
#include <libClarionTrinity/ClarionHit.hh>
#include <libClarionTrinity/Clarion.hh>

#include "SortFunctions.hh"

//global clean condition
bool Clean(const ClarionTrinity::Clover *clover) {
  return !(clover->nNonPrompt > 0 || clover->nSuppress > 0);
}

////Another block to test calibrations 
// -------------------- Calibration-check helpers (linear only) --------------------

struct CalRow {
  int ybin;
  double ycenter;
  double ch;
  double cherr;
  double Eexp;
};

struct CalCoef {
  double p0 = 0.0; // intercept
  double p1 = 1.0; // slope
  bool ok = false;
};

struct CloverMapResult {
  bool ok = false;
  char clover = '?';   // 'A','B',...
  int xtal = -1;       // 0..3 within clover
  int crystal_id = -1; // ybin-1
  int ybin = -1;
};

// Your ybin->clover ranges
static CloverMapResult MapYbinToClover(int ybin) {
  CloverMapResult r;
  r.ybin = ybin;
  r.crystal_id = ybin - 1;

  struct Range { char c; int lo; int hi; };
  static const Range ranges[] = {
    {'A',  1,  4},
    {'B',  6,  9},
    {'C', 11, 14},
    {'D', 17, 20},
    {'E', 22, 25},
    {'G', 33, 36},
    {'H', 38, 41},
    {'I', 43, 46},
    {'J', 49, 52},
  };

  for (const auto& rr : ranges) {
    if (ybin >= rr.lo && ybin <= rr.hi) {
      r.ok = true;
      r.clover = rr.c;
      r.xtal = ybin - rr.lo; // 0..3
      return r;
    }
  }
  return r; // ok=false if ybin not in any clover range
}

//This is here to read a small BEGe calibration file so that it can calibrate the BEGe detectors. I'm hoping that I can add this into a similar slot at clovers later
struct BEGeCalibration {
    double intercept = 0.0;
    double slope = 1.0;
    bool valid = false;
};

static bool ReadBEGeCalibration(
    const std::string& filename,
    BEGeCalibration cal[3])
{
    std::ifstream input(filename);

    if (!input.is_open()) {
        std::cerr << "ERROR: could not open BEGe calibration file: "
                  << filename << "\n";
        return false;
    }

    std::string line;
    int rowsLoaded = 0;

    while (std::getline(input, line)) {

        // Ignore blank lines
        if (line.empty()) {
            continue;
        }

        // Ignore comments, including comments with leading spaces
        const std::size_t firstCharacter =
            line.find_first_not_of(" \t\r\n");

        if (firstCharacter == std::string::npos) {
            continue;
        }

        if (line[firstCharacter] == '#') {
            continue;
        }

        std::istringstream stream(line);

        int detector = 0;
        int crystal = 0;
        double intercept = 0.0;
        double slope = 1.0;
        double shift = 0.0;

        // Expected format:
        // detector crystal intercept slope shift
        if (!(stream >> detector
                     >> crystal
                     >> intercept
                     >> slope
                     >> shift)) {

            std::cerr << "WARNING: could not parse BEGe calibration row:\n"
                      << line << "\n";
            continue;
        }

        if (detector < 1 || detector > 2) {
            std::cerr << "WARNING: unsupported BEGe detector number "
                      << detector << "\n";
            continue;
        }

        // Each BEGe is treated as one detector, crystal 0
        if (crystal != 0) {
            std::cerr << "WARNING: BEGe detector " << detector
                      << " has unexpected crystal number "
                      << crystal << "\n";
        }

        cal[detector].intercept = intercept;
        cal[detector].slope = slope;
        cal[detector].valid = true;

        ++rowsLoaded;

        std::cout << "Loaded BEGe " << detector
                  << " calibration: intercept=" << intercept
                  << ", slope=" << slope << "\n";
    }

    std::cout << "Loaded " << rowsLoaded
              << " BEGe calibration row(s) from "
              << filename << "\n";

    if (!cal[1].valid) {
        std::cerr << "ERROR: no valid calibration found for BEGe E1\n";
    }
    return cal[1].valid;
}

// Clover letter -> numeric index used in your cal-coefficient file
// Based on your file having clovers 1..5,7..10 and placeholders for 6,11,12.
static int CloverLetterToNumber(char c) {
  switch (c) {
    case 'A': return 1;
    case 'B': return 2;
    case 'C': return 3;
    case 'D': return 4;
    case 'E': return 5;
    case 'G': return 7;
    case 'H': return 8;
    case 'I': return 9;
    case 'J': return 10;
    default:  return -1;
  }
}

// Read centroid TSV: ybin ycenter channel chan_err energy_keV
static bool ReadEu152CentroidsTSV(const std::string& path,
                                 std::map<int, std::vector<CalRow>>& byYbin)
{
  byYbin.clear();
  std::ifstream fin(path);
  if (!fin) {
    std::cerr << "ERROR: cannot open centroid TSV: " << path << "\n";
    return false;
  }

  std::string line;
  bool first = true;
  while (std::getline(fin, line)) {
    if (line.empty()) continue;
    if (first) { first = false; continue; } // skip header row

    std::istringstream iss(line);
    CalRow r{};
    if (!(iss >> r.ybin >> r.ycenter >> r.ch >> r.cherr >> r.Eexp)) continue;
    byYbin[r.ybin].push_back(r);
  }
  return !byYbin.empty();
}

// Read cal coefficients: clover crystal intercept slope [maybe extra column(s)]
// We'll read first 4 tokens and ignore anything else on the line.
static bool ReadLinearCalCoefFile(const std::string& path,
                                 std::map<std::pair<int,int>, CalCoef>& coef)
{
  coef.clear();
  std::ifstream fin(path);
  if (!fin) {
    std::cerr << "ERROR: cannot open cal coef file: " << path << "\n";
    return false;
  }

  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    std::istringstream iss(line);
    int cl = -1, xt = -1;
    double p0 = 0.0, p1 = 1.0;

    if (!(iss >> cl >> xt >> p0 >> p1)) continue;

    CalCoef c;
    c.p0 = p0;
    c.p1 = p1;
    c.ok = true;
    coef[{cl, xt}] = c;
  }
  return !coef.empty();
}

// Make and write one residual canvas per crystal into out->/cal_check/
// ?E = (p0+p1*ch) - E_expected
static void MakeCalResidualPlots_Linear(TFile* out,
                                       const std::string& centroidTSV,
                                       const std::string& calfile,
                                       bool onlyMappedYbins = true)
{
  std::map<int, std::vector<CalRow>> byYbin;
  if (!ReadEu152CentroidsTSV(centroidTSV, byYbin)) return;

  std::map<std::pair<int,int>, CalCoef> coef;
  if (!ReadLinearCalCoefFile(calfile, coef)) return;

  TDirectory* saveDir = gDirectory;
  TDirectory* d = out->GetDirectory("cal_check");
  if (!d) d = out->mkdir("cal_check");
  d->cd();

  int made = 0;

  for (auto& kv : byYbin) {
    const int ybin = kv.first;
    auto& rows = kv.second;
    if ((int)rows.size() < 2) continue;

    auto info = MapYbinToClover(ybin);
    if (onlyMappedYbins && !info.ok) continue;

    int cloverNum = info.ok ? CloverLetterToNumber(info.clover) : -1;
    int xtalNum   = info.ok ? info.xtal : -1;
    if (onlyMappedYbins && (cloverNum < 0 || xtalNum < 0)) continue;

    auto it = coef.find({cloverNum, xtalNum});
    if (it == coef.end() || !it->second.ok) {
      std::cerr << "WARN: no coef for ybin " << ybin
                << " (clover " << info.clover << "=" << cloverNum
                << ", xtal " << xtalNum << ")\n";
      continue;
    }
    const auto& cc = it->second;

    // Build graph ?E(E)
    auto* g = new TGraphErrors((int)rows.size());
    g->SetName(Form("g_cal_res_ybin%d", ybin));
    g->SetMarkerStyle(20);

    // --- NEW: graph for "calibration line + peaks" plot
    // Points are (channel centroid, expected energy)
    auto* gCal = new TGraphErrors((int)rows.size());
    gCal->SetName(Form("g_cal_linepoints_ybin%d", ybin));
    gCal->SetMarkerStyle(20);

    auto* gEE = new TGraphErrors((int)rows.size());
    gEE->SetName(Form("g_Ecal_vs_Eexp_ybin%d", ybin));
    gEE->SetMarkerStyle(20);



    double Emin = 1e9, Emax = -1e9;
    double dEmin = 1e9, dEmax = -1e9;

    // --- NEW: residual histogram for this crystal (ybin)
    const double hLo = -2.0;   // keV, change if you want wider/narrower
    const double hHi =  2.0;   // keV
    const int    hBins = 80;

    TH1D* hRes = new TH1D(Form("h_residual_ybin%d", ybin),
			  Form("Residuals ybin %d;#DeltaE = E_{cal} - E_{exp} (keV);Counts", ybin),
			  hBins, hLo, hHi);


    for (int i = 0; i < (int)rows.size(); i++) {
      const double Eexp = rows[i].Eexp;
      const double ch   = rows[i].ch;
      const double dch  = rows[i].cherr;

      const double Ecal = cc.p0 + cc.p1 * ch;
      const double dE   = Ecal - Eexp;
      hRes->Fill(dE);//NEW

      const double dEerr = std::abs(cc.p1) * dch; // linear propagation

      g->SetPoint(i, Eexp, dE);
      g->SetPointError(i, 0.0, dEerr);

      gCal->SetPoint(i, ch, Eexp);
      gCal->SetPointError(i, dch, 0.0);

      gEE->SetPoint(i, Eexp, Ecal);
      gEE->SetPointError(i, 0.0, dEerr);  // x error 0, y error from centroid uncertainty



      Emin = std::min(Emin, Eexp);
      Emax = std::max(Emax, Eexp);
      dEmin = std::min(dEmin, dE - dEerr);
      dEmax = std::max(dEmax, dE + dEerr);
    }

    // Nice label
    std::string tag;
    if (info.ok) {
      tag = Form("Clover %c crystal %d (ybin %d, crystalID %d)  p0=%.5f p1=%.5f",
                 info.clover, info.xtal, info.ybin, info.crystal_id, cc.p0, cc.p1);
    } else {
      tag = Form("ybin %d (crystalID %d)  p0=%.5f p1=%.5f", ybin, ybin-1, cc.p0, cc.p1);
    }

    TCanvas* c = new TCanvas(Form("cal_residual_ybin%d", ybin),
                             tag.c_str(), 900, 450);

    g->SetTitle(Form("%s;Energy (keV);#DeltaE = E_{cal} - E_{exp} (keV)", tag.c_str()));
    g->Draw("AP");

    // force a sensible y-range
    double pad = 0.15 * std::max(1e-6, (dEmax - dEmin));
    g->GetYaxis()->SetRangeUser(dEmin - pad, dEmax + pad);

    // 0-line
    TLine* l0 = new TLine(Emin, 0.0, Emax, 0.0);
    l0->SetLineColor(kRed);
    l0->SetLineStyle(2);
    l0->SetLineWidth(2);
    l0->Draw("same");

    c->Write();
    g->Write();
    hRes->Write();   // NEW

    //ANOTHER TESTING CAL PLOT
    // --- NEW: "calibration line + peaks" canvas
    TCanvas* c2 = new TCanvas(Form("cal_linepoints_ybin%d", ybin),
			      Form("%s (line+points)", tag.c_str()), 900, 450);

    gCal->SetTitle(Form("%s;Centroid channel;Expected energy (keV)", tag.c_str()));
    gCal->Draw("AP");

    // Determine channel range for the red calibration line
    double chMin = 1e99, chMax = -1e99;
    for (int i = 0; i < gCal->GetN(); i++) {
      double x, y;
      gCal->GetPoint(i, x, y);
      chMin = std::min(chMin, x);
      chMax = std::max(chMax, x);
    }

    // Red calibration line: E = p0 + p1*ch
    TLine* calLine = new TLine(chMin, cc.p0 + cc.p1*chMin,
			       chMax, cc.p0 + cc.p1*chMax);
    calLine->SetLineColor(kRed);
    calLine->SetLineWidth(2);
    calLine->Draw("same");

    // Optional legend
    TLegend* leg = new TLegend(0.12, 0.75, 0.35, 0.88);
    leg->AddEntry(gCal, "Eu-152 peaks (expected E)", "p");
    leg->AddEntry(calLine, "Calibration: E = p0 + p1*ch", "l");
    leg->Draw();

    c2->Write();
    gCal->Write();

    //YET ANOTHER TESTING PLOT
    TCanvas* cEE = new TCanvas(Form("cal_Ecal_vs_Eexp_ybin%d", ybin),
                           Form("%s (Ecal vs Eexp)", tag.c_str()), 900, 450);

    gEE->SetTitle(Form("%s;Expected energy E_{exp} (keV);Calibrated energy E_{cal} (keV)", tag.c_str()));
    gEE->Draw("AP");

    // Determine plot range
    double Emin2 = 1e99, Emax2 = -1e99;
    for (int i = 0; i < gEE->GetN(); i++) {
      double x, y;
      gEE->GetPoint(i, x, y);
      Emin2 = std::min(Emin2, x);
      Emax2 = std::max(Emax2, x);
    }

    // Ideal line y = x
    TLine* yEqx = new TLine(Emin2, Emin2, Emax2, Emax2);
    yEqx->SetLineColor(kRed);
    yEqx->SetLineWidth(2);
    yEqx->Draw("same");

    cEE->Write();
    gEE->Write();

    delete yEqx;
    delete cEE;
    delete gEE;


    delete leg;
    delete calLine;
    delete c2;


    delete hRes;     // NEW
    delete l0;
    delete c;
    delete g;
    delete gCal; //NEW
    made++;
  }

  std::cout << "CalCheck: made " << made << " residual canvases in /cal_check/\n";
  saveDir->cd();
}
///END OF Calibration Block


int main(int argc, char **argv) {
  double pi = 3.1415926535;
  //process options and command line arguments
  int opt;
  int verbose = 0;
  int projectileDetected = 0;
  int makeCalPlots = 0;
  while ((opt = getopt(argc, argv, "vpC")) != -1) {
  switch (opt) {
  case 'v':
    verbose = 1;
    break;
  case 'p':
    projectileDetected = 1;
    break;
  case 'C':
    makeCalPlots = 1;
    break;
    default:
      abort();
  }
  }

  if (argc - optind < 4) {
    std::cout << "usage ./FESort [file list] [output file] [write option] [clarion_cal] "  << std::endl;
    return 0;
  }
  
  std::string files(argv[optind]); 
  std::cout << "                         Files from : " << ANSI_COLOR_GREEN << files << ANSI_COLOR_RESET << std::endl;
  std::string outfilename(argv[optind+1]);
  std::cout << "                      Histograms to : " << ANSI_COLOR_GREEN << outfilename << ANSI_COLOR_RESET << std::endl;
  std::string fileopt(argv[optind+2]);
  std::string calfile(argv[optind+3]); 


  ////THIS IS FOR THE CALIBRATION TESTING FEB 4
  std::string centroidTSV = "";
  std::string coefFile    = "";

  if (makeCalPlots) {
    if (argc - optind < 5) {
      std::cout << "usage: ./FESort -C [file list] [output file] [write option] "
	"[clarion_cal] [eu152_peaks.tsv]\n";
      return 0;
    }
    centroidTSV = argv[optind+4];
  }
  ///CALIBRATION TESTING 


  //data files
  std::vector<std::string> dataPaths;
  //establish if the file passed is a text file containing paths to evt.to files, or an evt.to file itself
  std::string suffix = files.substr(files.find_last_of(".")+1, files.size() - files.find_last_of(".")-1);
  if (suffix == "to" || suffix == "evt") {
    dataPaths.push_back(files); 
  }
  else if ( suffix=="txt" || suffix=="dat") {
  std::ifstream inFiles(files);
  std::string path;
  while (inFiles >> path) {
    dataPaths.push_back(path);
  }
  }
  else {
    std::cout << "Suffix for file (list) unclear: please use *.evt, *.to for listmode data file, or *.txt, *.dat for a list of paths to data files" << std::endl;
    return 0;
  }

  //experimental definition, event-building parameters
  std::string defPath = "experimentdefinition.def";
  PIXIE::Experiment_Definition definition;
  float coincWindow = 1500; //in ns
  if (definition.open(defPath) < 0) { return -1; }
  if (definition.read() < 0) { return -1; }
  std::cout << "Experimental definition loaded from : " << ANSI_COLOR_GREEN << defPath << ANSI_COLOR_RESET << std::endl;
  if (verbose) { definition.print(); }

  //trinity definition, angles, configuration
  ClarionTrinity::Trinity trinity("trinity.conf");
  if (trinity.ReadAngleMap("trinity.map") < 0) { return -1; }      
  std::cout << "              Trinity configuration : " << ANSI_COLOR_GREEN << "trinity.conf" << ANSI_COLOR_RESET << std::endl;
  std::cout << "                                    : " << ANSI_COLOR_GREEN << "trinity.map"  << ANSI_COLOR_RESET << std::endl;
  if (verbose) { trinity.PrintConf(); }
  trinity.conf.rejectPU = true;
  trinity.conf.rejectOOR = true;
  
  //for comissioning: all QDCs
  trinity.conf.SetPIDType(ClarionTrinity::PIDType::kQDC);

  if (trinity.conf.ReadQDCParams("QDCparams.txt") < 0) { return -1; }
  //clean particle threshold
  double tdiff_thresh = 100.0; //ns

  // Prompt gamma-gamma coincidence gate.
  // Tune this value later using a gamma-gamma time-difference spectrum if needed.
  const double gg_tdiff_thresh = 200.0; // ns

  auto gg_dt_ns = [](const auto* first, const auto* second) {
    return ((long long int)first->timestamp -
	    (long long int)second->timestamp) / 3276.8;
  };

  auto same_clover_pair_index = [](int xtal1, int xtal2) {
    if (xtal1 > xtal2) {
      int tmp = xtal1;
      xtal1 = xtal2;
      xtal2 = tmp;
    }

    if (xtal1 == 0 && xtal2 == 1) return 0;
    if (xtal1 == 0 && xtal2 == 2) return 1;
    if (xtal1 == 0 && xtal2 == 3) return 2;
    if (xtal1 == 1 && xtal2 == 2) return 3;
    if (xtal1 == 1 && xtal2 == 3) return 4;
    if (xtal1 == 2 && xtal2 == 3) return 5;

    return -1;
  };

  //auto gg_prompt = [gg_tdiff_thresh, gg_dt_ns](const auto* first, const auto* second) {
  // const double dt =
  // ((long long int)first->timestamp -
  //  (long long int)second->timestamp) / 3276.8;
  //Changing gg_prompt to match with gg_dt_ns for comparison
  auto gg_prompt = [gg_tdiff_thresh, gg_dt_ns](const auto* first, const auto* second) {
    const double dt = gg_dt_ns(first, second);
    return std::abs(dt) <= gg_tdiff_thresh;
  };

  int peak[2] = {75, 94};
  int tail[2] = {105,159};

  //clarion definition, angles, calibration, configuration, efficiency
  ClarionTrinity::Clarion clarion("clarion.conf");
  if (clarion.ReadAngleMap("clarion.map") < 0) { return -1; }
  if (clarion.conf.ReadCal(calfile) < 0) { return -1; }
  clarion.conf.ReadCTCal("clar_ct.cal");
  clarion.conf.BGOVetoTime = 500.0; //ns
  clarion.conf.EnThresh = 50.0; //keV
  clarion.conf.AddBackTDiff = 200.0; //ns
  clarion.conf.RawUpThresh = 60000; //ADC units
  clarion.conf.RawBGOThresh = 20; //ADC units

  if (clarion.conf.ReadEffIDCal("cal/clar_effi_ab.cal") < 0 ) { return -1; } //calibration for each detector individually
  if (clarion.conf.ReadEfficiencyCal("cal/clar_eff_ab.cal") < 0 ) { return -1; } //calibration for whole array  

  std::cout << "              Clarion configuration : " << ANSI_COLOR_GREEN << "clarion.conf"        << ANSI_COLOR_RESET << std::endl;
  std::cout << "                                    : " << ANSI_COLOR_GREEN << "clarion.map"         << ANSI_COLOR_RESET << std::endl;
  std::cout << "                                    : " << ANSI_COLOR_GREEN << "clarion.cal"         << ANSI_COLOR_RESET << std::endl;
  std::cout << "                                    : " << ANSI_COLOR_GREEN << "cal/clar_eff_ab.cal" << ANSI_COLOR_RESET << std::endl;

  if (verbose) {
  std::cout << "Efficiency at: " << std::endl;
  std::cout << "121: " << clarion.conf.Efficiency(121.0) << std::endl;
  std::cout << "159: " << clarion.conf.Efficiency(159.0) << std::endl;
  std::cout << "500: " << clarion.conf.Efficiency(500.0) << std::endl;
  std::cout << "1000: " << clarion.conf.Efficiency(1000.0) << std::endl;
  std::cout << "1092: " << clarion.conf.Efficiency(1092.0) << std::endl;
  std::cout << "1332: " << clarion.conf.Efficiency(1332.0) << std::endl;
  std::cout << "1500: " << clarion.conf.Efficiency(1500.0) << std::endl;
  std::cout << "1550: " << clarion.conf.Efficiency(1550.0) << std::endl;
  std::cout << "1650: " << clarion.conf.Efficiency(1650.0) << std::endl;
  std::cout << "1750: " << clarion.conf.Efficiency(1750.0) << std::endl;
  }

  if (verbose) { clarion.PrintConf(); }

  // BEGe setup This is to find the BEGe calibration file in FE_sorted
  std::string begeCalFile =
    "/data/FSU/rawdata/FSUexp/MSU_Ga72/FE_Sorted/133BaCal.cal";

  BEGeCalibration begeCal[3];

  if (!ReadBEGeCalibration(begeCalFile, begeCal)) {
    std::cerr << "ERROR: failed to read BEGe calibration file: "
              << begeCalFile << "\n";
    return -1;
  }

  std::cout << "                 BEGe calibration : "
	    << begeCalFile << "\n";

  ClarionTrinity::Event event(clarion, trinity);

  //ROOT histogram definitions
  TFile file(outfilename.c_str(), fileopt.c_str());
  
  //check option -> we should never want a read-only option
  if (file.GetOption() == "READ") { std::cout << "Output histogram file opened in read-only mode" << std::endl; return -1; }

 // clarion histograms
  //changing the binning of this first histo to test cal.
  TH2I *clar_hits = get_TH2I(&file, "clar_hits", "Energy; Energy; Clover Channel ID", 8192, 0, 4096, 16*4, 0, 16*4);
  TH2I *clar_hitsclean = get_TH2I(&file, "clar_hitsclean", "Clean Energy; Energy; Clover Channel ID", 8192, 0, 4096, 16*4, 0, 16*4);
  TH2I *clar_bgotdiff = get_TH2I(&file, "clar_bgotdiff", "BGO Time Difference; Time (ns); Clover Channel ID", 400, -2000, 2000, 16*4, 0, 16*4);
  TH2I *clar_e = get_TH2I(&file, "clar_e", "Energy; Energy; Clover Channel ID", 8192, 0, 4096, 16*4, 0, 16*4);  
  TH2I *clar_edirty = get_TH2I(&file, "clar_edirty", "Energy; Energy; Clover Channel ID", 8192, 0, 4096, 16*4, 0, 16*4);
  TH2D *clar_e_eff = get_TH2D(&file, "clar_e_eff", "Energy; Energy; Clover Channel ID", 8192, 0, 4096, 16*4, 0, 16*4);
  //Chat Testing Gamma Spectra
  TH1I *pg_e_prompt_cut = get_TH1I(gDirectory,
				   "pg_e_prompt_cut", "Gamma energy (PID+prompt);Energy (keV);Counts",
				   10000, 0, 10000);
  // NEW: clean addback gamma singles with no PID or timing gate
  TH1I *g_e_ungated = get_TH1I(
			       gDirectory,
			       "g_e_ungated",
			       "Ungated clean gamma singles;Energy (keV);Counts",
			       10000, 0, 10000
			       );
  //Calibration Hist
  // Global/summed spectrum (singles)
  TH1F *h_calib_en_all   = new TH1F("calib_en_all",   "Clarion singles (all crystals);Energy [keV];Counts", 8192, 0, 4096);

  // "raw" histograms - on sub events before GaGG logic
  TH2I *raw_e = get_TH2I(&file, "raw_e", "Raw Energy; Energy; Channel ID", 4096, 0, 8*4096, 16*13, 0, 16*13); 
  // adding this hist for making calibrations using fit click tawfik
  //TH2D *raw_e_cal = get_TH2D(&file, "raw_e_cal", "Raw Energy For Cal; Energy; Channel ID", 4096, 0, 8*4096, 16*4, 0, 16*4);

  /// Standalone BEGe detector spectra
  TH1I *bege_e1_raw = get_TH1I(
			       &file,
			       "bege_e1_raw",
			       "BEGe E1 Raw Spectrum;ADC Channel;Counts",
			       32768, 0, 32768
			       );

  TH1I *bege_e2_raw = get_TH1I(
			       &file,
			       "bege_e2_raw",
			       "BEGe E2 Raw Spectrum;ADC Channel;Counts",
			       32768, 0, 32768
			       );

  TH1I *bege_e1_raw_clean = get_TH1I(
				     &file,
				     "bege_e1_raw_clean",
				     "BEGe E1 Raw Spectrum, no pileup or out-of-range;ADC Channel;Counts",
				     32768, 0, 32768
				     );

  TH1I *bege_e2_raw_clean = get_TH1I(
				     &file,
				     "bege_e2_raw_clean",
				     "BEGe E2 Raw Spectrum, no pileup or out-of-range;ADC Channel;Counts",
				     32768, 0, 32768
				     );

  TH1I *bege_e1_cal = get_TH1I(
			       &file,
			       "bege_e1_cal",
			       "BEGe E1 Calibrated Spectrum;Energy (keV);Counts",
			       8192, 0, 4096
			       );

  TH1I *bege_e1_cal_pid = get_TH1I(
				   &file,
				   "bege_e1_cal_pid",
				   "BEGe E1 Calibrated PID-gated Spectrum;Energy (keV);Counts",
				   8192, 0, 4096
				   );

  // Raw spectra for every channel in the LaBr module, slot S14
  TH2I *labr_raw = get_TH2I(
			    &file,
			    "labr_raw",
			    "LaBr raw spectra, slot S14;Raw ADC Energy;S14 Channel",
			    8192, 0, 32768,
			    16, 0, 16
			    );

  TH2I *labr_raw_clean = get_TH2I(
				  &file,
				  "labr_raw_clean",
				  "LaBr raw spectra, slot S14, clean;Raw ADC Energy;S14 Channel",
				  8192, 0, 32768,
				  16, 0, 16
				  );
  //Testing for 1d for the labr detectors (NEED TO FILL THESE LATER ON)
  TH1I *labr_single_raw = get_TH1I(
			       &file,
			       "labr_single_raw",
			       "LaBr Raw Spectrum;ADC Channel;Counts",
			       32768, 0, 32768
			       );

  TH1I *labr_single_raw_clean = get_TH1I(
				     &file,
				     "labr_single_raw_clean",
				     "LaBr Raw Spectrum, no pileup or out-of-range;ADC Channel;Counts",
				     32768, 0, 32768
				     );


  TH1I *raw_pu = get_TH1I(&file, "raw_pu", "Raw Pileup; Channel ID; Counts", 16*13, 0, 16*13);
  TH1I *raw_or = get_TH1I(&file, "raw_or", "Raw Out of Range; Channel ID; Counts", 16*13, 0, 16*13);
  TH2I *raw_mult = get_TH2I(&file, "raw_mult", "Raw Subevent Multiplicity; Wall Time (s); Multiplicity", 4096, 0, 57600, 50, 0, 50);
  TH2I *raw_seqtdiff = get_TH2I(&file, "raw_seqtdiff", "Raw sequential time difference; Time (us); Channel ID", 4000, 0, 400, 13*16, 0, 13*16);
  TH2I *raw_seqtdiff_clean = get_TH2I(&file, "raw_seqtdiff_clean", "Raw sequential time difference; Time (us); Channel ID", 4000, 0, 400, 13*16, 0, 13*16);
  TH2I *raw_rates = get_TH2I(&file, "raw_rates", "Rates; Wall Time (s); Channel ID", 4096, 0, 57600, 16*13, 0, 16*13);  
  TH2I *raw_hitpat = get_TH2I(&file, "raw_hitpat", "Hit pattern; ID1; ID2", 13*16, 0, 13*16, 13*16, 0, 13*16);
  TH2I *raw_hitpat_db = get_TH2I(&file, "raw_hitpat_db", "Hit pattern; ID1; ID2", 13*16, 0, 13*16, 13*16, 0, 13*16);
  TH2I *raw_e_wt = get_TH2I(&file, "raw_e_wt", "Raw Energy vs Wall Time; Wall Time (s); Raw Energy", 4096, 0, 57600, 4096, 0, 8*4096);
  
  //Hist to look at count rate
  TH1D *rate_all = new TH1D("rate_all",
			    "All events vs elapsed time;Elapsed Time (s);Counts / s",
			    57600, 0, 57600);

  TH1I* h_clover_xtal_sum[17] = {nullptr};   // 1..16
  TH1I* h_clover_xtal_sum_clean[17] = {nullptr}; // optional: suppressed-clean version
  
  //This is for the anti-gated run by run test
  TH1I* h_bg_sum = new TH1I("h_bg_sum",
  "Background-like summed gamma singles;Energy (keV);Counts",
  4096, 0, 4096);

  //Histograms to look at clover sums and basic gamma gamma
  for (int cid = 1; cid <= 16; ++cid) {
    h_clover_xtal_sum[cid] =
      new TH1I(Form("clover%02d_xtal_sum", cid),
	       Form("Clover %d crystal-hit singles (no addback);Energy [keV];Counts", cid),
	       8192, 0, 4096);

    h_clover_xtal_sum_clean[cid] =
      new TH1I(Form("clover%02d_xtal_sum_clean", cid),
	       Form("Clover %d crystal-hit singles (Suppress==false);Energy [keV];Counts", cid),
	       8192, 0, 4096);
  }

  TH2I* gg_sameclover[17] = {nullptr};       // 1..16
  TH2I* gg_sameclover_clean[17] = {nullptr}; // optional

  for (int cid = 1; cid <= 16; ++cid) {
    gg_sameclover[cid] =
      new TH2I(Form("gg_sameclover%02d", cid),
	       Form("#gamma-#gamma same clover %d (crystal hits);E1 [keV];E2 [keV]", cid),
	       2048, 0, 4096, 2048, 0, 4096);

    gg_sameclover_clean[cid] =
      new TH2I(Form("gg_sameclover%02d_clean", cid),
	       Form("#gamma-#gamma same clover %d (Suppress==false);E1 [keV];E2 [keV]", cid),
	       2048, 0, 4096, 2048, 0, 4096);
  }

  // PID-gated versions THESE ARE FOR USING THE CUTS ON THE GG AND CLOVER HISTS
  TH1I* h_clover_xtal_sum_pid[17]        = {nullptr};  // 1..16
  TH1I* h_clover_xtal_sum_clean_pid[17]  = {nullptr};  // 1..16

  TH2I* gg_sameclover_pid[17]            = {nullptr};  // 1..16
  TH2I* gg_sameclover_pid_clean[17]      = {nullptr};  // 1..16

  for (int cid = 1; cid <= 16; ++cid) {
    h_clover_xtal_sum_pid[cid] =
      new TH1I(Form("clover%02d_xtal_sum_pid", cid),
	       Form("Clover %d crystal-hit singles (PID gated);Energy [keV];Counts", cid),
	       8192, 0, 4096);

    h_clover_xtal_sum_clean_pid[cid] =
      new TH1I(Form("clover%02d_xtal_sum_clean_pid", cid),
	       Form("Clover %d crystal-hit singles (Suppress==false, PID gated);Energy [keV];Counts", cid),
	       8192, 0, 4096);

    gg_sameclover_pid[cid] =
      new TH2I(Form("gg_sameclover%02d_pid", cid),
	       Form("#gamma-#gamma same clover %d (PID gated);E1 [keV];E2 [keV]", cid),
	       2048, 0, 4096, 2048, 0, 4096);

    gg_sameclover_pid_clean[cid] =
      new TH2I(Form("gg_sameclover%02d_clean_pid", cid),
	       Form("#gamma-#gamma same clover %d (Suppress==false, PID gated);E1 [keV];E2 [keV]", cid),
	       2048, 0, 4096, 2048, 0, 4096);
  }
  
  //THESE ARE THE TEST HISTS FOR A CLOVER-CLOVER GG
  TH2I* gg_cloverclover = new TH2I(
				   "gg_cloverclover",
				   "#gamma-#gamma different clovers (all clover pairs);E1 [keV];E2 [keV]",
				   2048, 0, 4096,
				   2048, 0, 4096
				   );

  TH2I* gg_cloverclover_clean = new TH2I(
					 "gg_cloverclover_clean",
					 "#gamma-#gamma different clovers (Suppress==false);E1 [keV];E2 [keV]",
					 2048, 0, 4096,
					 2048, 0, 4096
					 );

  TH2I* gg_cloverclover_pid = new TH2I(
				       "gg_cloverclover_pid",
				       "#gamma-#gamma different clovers (PID gated);E1 [keV];E2 [keV]",
				       2048, 0, 4096,
				       2048, 0, 4096
				       );

  TH2I* gg_cloverclover_clean_pid = new TH2I(
					     "gg_cloverclover_clean_pid",
					     "#gamma-#gamma different clovers (Suppress==false, PID gated);E1 [keV];E2 [keV]",
					     2048, 0, 4096,
					     2048, 0, 4096
					     );



  // "trin" histograms - after forming GaGG detectors 
  TH2I *trin_gainmatch = get_TH2I(&file, "trin_gainmatch", "GaGG gain match; Ratio; Trinity ID", 512, 0.5, 2.5, 600, 1, 601);
  TH2I *trin_e = get_TH2I(&file, "trin_e", "GaGG energy; Energy; Trinity ID", 4096, 0, 8*8192, 600, 1, 601);
  TH2I *trin_rates = get_TH2I(&file, "trin_rates", "GaGG rates; Wall Time (s); Trinity ID", 4096, 0, 32768, 600, 1, 601);
  TH2I *trin_hitpat = get_TH2I(&file, "trin_hitpat", "Gagg ID; GaggID; Counts", 600, 1, 601, 600, 1, 601);
  TH2I *trin_pu = get_TH2I(&file, "trin_pu", "GaGG pileups; Energy; Trinity ID", 4096, 0, 8*8192, 600, 1, 601);
  TH2I *trin_tdiff = get_TH2I(&file, "trin_tdiff", "GaGG self time difference; Time (ns); Trinity ID", 200, -1000, 1000, 600, 1, 601);
  TH2I *trin_seqtdiff = get_TH2I(&file, "trin_seqtdiff", "GaGG sequential time difference; Time (us); Trinity ID", 4000, 0, 400, 600, 1, 601);
  TH2I *trin_backdiff1 = get_TH2I(&file, "trin_backdiff1", "Background Difference 1; Difference; Trinity ID", 1024, -512, 512, 600, 1, 601);  
  TH2I *trin_backdiff2 = get_TH2I(&file, "trin_backdiff2", "Background Difference 2; Difference; Trinity ID", 1024, -512, 512, 600, 1, 601);  
  TH2I *trin_z_wt = get_TH2I(&file, "trin_z_wt", "Trinity zero degree vs Wall time; Wal Time (s); Energy", 4096, 0, 36000, 4096, 0, 8*8192);
  TH2I *trin_feat = get_TH2I(&file, "trin_feat", "Feature distinguishing; Ratio; Counts", 1024, 0, 1, 600, 1, 601);


  TH1I *part_dmult = get_TH1I(&file, "part_dmult", "Dirty particle multiplicity; Multipliticy;", 100, 0, 100);
  TH1I *part_mult = get_TH1I(&file, "part_mult", "Clean particle multiplicity; Multipliticy;", 100, 0, 100);
  TH1I *part_himult = get_TH1I(&file, "part_himult", "Valid particle multiplicity; Multipliticy;", 100, 0, 100);
  TH1I *part_combmult = get_TH1I(&file, "part_combmult", "Combination Particle Multiplicity; Multipliticy;", 1000, 0, 1000);

  TH2I *part_dmult_id = get_TH2I(&file, "part_dmult_id", "Dirty particle multiplicity; Multipliticy; Trinity ID", 10, 0, 10, 600, 0, 601);
  TH2I *part_mult_id = get_TH2I(&file, "part_mult_id", "Clean particle multiplicity; Multipliticy; Trinity ID", 10, 0, 10, 600, 0, 601);
  TH2I *part_himult_id = get_TH2I(&file, "part_himult_id", "Valid particle multiplicity; Multipliticy; Trinity ID", 10, 0, 10, 600, 0, 601);

  TH1I *part_singles= get_TH1I(&file, "part_singles", "Valid particle singles;", 600, 0, 600);
  TH2I *part_singles_wt = get_TH2I(&file, "part_singles_wt", "Valid particle singles vs Wall Time; Wall Time (s); GaGG ID;", 4096, 0, 57600, 600, 0, 600);
  TH2I *part_en_wt = get_TH2I(&file, "part_en_wt", "Particle energy (201) vs wall time; Wall Time (s); Particle energy", 4096, 0, 36000, 4096, 0, 8*8192); 
  TH1I *part_en[6];

  // particle-gamma
  TH2I *pg_enr_p = get_TH2I(&file, "pg_enr_p", "Particle-gamma prompt Energy vs Theta; Energy; Ring", 8192, 0, 4096, 6, 0, 6);
  TH2I *pg_enr_np = get_TH2I(&file, "pg_enr_np", "Particle-gamma non-prompt Energy vs Theta; Energy; Ring", 8192, 0, 4096, 6, 0, 6);
  
  TH2I *pg_enwt = get_TH2I(&file, "pg_enwt", "Particle-gamma Energy vs Wall Time; Wall Time (s); Energy (keV)", 8192, 0, 21600, 4096, 0, 2048);

  TH2I *pg_nden_trin_p = get_TH2I(&file, "pg_nden_trin_p", "Particle-gamma Energy vs GAGG ID (prompt); Energy (keV); GAGG ID", 4096, 0, 4096, 600, 1, 601);
  TH2I *pg_nden_trin_lnp = get_TH2I(&file, "pg_nden_trin_lnp", "Particle-gamma Energy vs GAGG ID (low non prompt); Energy (keV); GAGG ID", 4096, 0, 4096, 600, 1, 601);
  TH2I *pg_nden_trin_hnp = get_TH2I(&file, "pg_nden_trin_hnp", "Particle-gamma Energy vs GAGG ID (high non prompt); Energy (keV); GAGG ID", 4096, 0, 4096, 600, 1, 601);

  TH2I *pg_clar_trin_p = get_TH2I(&file, "pg_clar_trin_p", "Crystal ID vs GAGG ID (gated, prompt); Energy (keV); GAGG ID", 16*4, 0, 16*4, 600, 1, 601);
  TH2I *pg_clar_trin_np = get_TH2I(&file, "pg_clar_trin_np", "Crystal ID vs GAGG ID (gated, non-prompt); Energy (keV); GAGG ID", 16*4, 0, 16*4, 600, 1, 601);

  TH2I *pgg = get_TH2I(&file, "pgg", "Particle-Gamma-Gamma; Energy (keV); Energy (keV)", 4096, 0, 4096, 4096, 0, 4096);

  TH2I *pg_tdiff[6];
  TH2I *pg_ent[6];
  TH2I *pg_tent[6];
  TH2I *pg_ndent[6];
  TH2D *pg_ent_eff[6];
  TH2D *pg_tent_eff[6];
  TH2D *pg_mb_ent_eff[6][10];

  float hpge_angs[4] = {48.5, 90.0, 131.4, 149.4};
  int rings[6] = {0, 1, 1, 1, 1, 1};
  int trin_xtls[6] = {1, 8, 10, 14, 16, 16};

  //THIS IS TESTING ADDBACK 
  TH1I* h_addback_all = new TH1I(
				 "h_addback_all",
				 "All clover addback gammas;Energy [keV];Counts",
				 8192, 0, 4096
				 );

  TH1I* h_noaddback_all = new TH1I(
				   "h_noaddback_all",
				   "All clover crystal-hit singles, no addback;Energy [keV];Counts",
				   8192, 0, 4096
				   );

  // Addback gamma spectrum when TRINITY does not fire
  TH1I* h_addback_no_trinity = new TH1I(
					"h_addback_no_trinity",
					"Clover addback gammas with no TRINITY hit;Energy [keV];Counts",
					8192, 0, 4096
					);

  // Optional sanity-check spectrum: addback gammas when TRINITY does fire
  TH1I* h_addback_with_trinity = new TH1I(
					  "h_addback_with_trinity",
					  "Clover addback gammas with at least one TRINITY hit;Energy [keV];Counts",
					  8192, 0, 4096
					  );

  TH1I* h_addback_no_trinity_clean = new TH1I(
					      "h_addback_no_trinity_clean",
					      "Clean clover addback gammas with no TRINITY hit;Energy [keV];Counts",
					      8192, 0, 4096
					      );

  TH2I* gg_cloverclover_ab_no_trinity = new TH2I(
						 "gg_cloverclover_ab_no_trinity",
						 "#gamma-#gamma clover-clover addback (no TRINITY);E1 [keV];E2 [keV]",
						 2048, 0, 4096,
						 2048, 0, 4096
						 );
  //doubled binning for Dr. Haring-Kaye GNUscope
  TH2I* gg_cloverclover_ab_clean_no_trinity = new TH2I(
						       "gg_cloverclover_ab_clean_no_trinity",
						       "#gamma-#gamma clover-clover addback clean (no TRINITY);E1 [keV];E2 [keV]",
						       2048, 0, 4096,
						       2048, 0, 4096
						       );

  //More ADDBACK Histograms
  TH2I* gg_cloverclover_ab = new TH2I(
				      "gg_cloverclover_ab",
				      "#gamma-#gamma clover-clover addback;E1 [keV];E2 [keV]",
				      2048, 0, 4096,
				      2048, 0, 4096
				      );

  TH2I* gg_cloverclover_ab_pid = new TH2I(
					  "gg_cloverclover_ab_pid",
					  "#gamma-#gamma clover-clover addback (PID gated);E1 [keV];E2 [keV]",
					  2048, 0, 4096,
					  2048, 0, 4096
					  );
  //doubled binning for GNUscope
  TH2I* gg_cloverclover_ab_clean = new TH2I(
					    "gg_cloverclover_ab_clean",
					    "#gamma-#gamma clover-clover addback clean;E1 [keV];E2 [keV]",
					    2048, 0, 4096,
					    2048, 0, 4096
					    );
  //moving the binning around here to see if speed increases when trying to run files
  TH2I* gg_cloverclover_ab_clean_pid = new TH2I(
						"gg_cloverclover_ab_clean_pid",
						"#gamma-#gamma clover-clover addback clean (PID gated);E1 [keV];E2 [keV]",
						8192, 0, 4096,
						8192, 0, 4096
						);

  TH1D* gg_sameclover_tdiff_all = new TH1D(
					   "gg_sameclover_tdiff_all",
					   "#gamma-#gamma time difference within same clover;#Delta t_{#gamma#gamma} [ns];Counts",
					   2000, -1000, 1000
					   );

  TH1D* gg_sameclover_tdiff_pid = new TH1D(
					   "gg_sameclover_tdiff_pid",
					   "#gamma-#gamma time difference within same clover, PID gated;#Delta t_{#gamma#gamma} [ns];Counts",
					   2000, -1000, 1000
					   );

  TH2D* gg_sameclover_e_vs_tdiff = new TH2D(
					    "gg_sameclover_e_vs_tdiff",
					    "#gamma energy vs #gamma-#gamma time difference within same clover;#Delta t_{#gamma#gamma} [ns];E_{#gamma} [keV]",
					    1000, -1000, 1000,
					    1024, 0, 4096
					    );

  TH2D* gg_sameclover_pair_vs_tdiff = new TH2D(
					       "gg_sameclover_pair_vs_tdiff",
					       "Same-clover crystal-pair timing;#Delta t_{#gamma#gamma} [ns];Crystal pair",
					       1000, -1000, 1000,
					       6, 0, 6
					       );

  //file.mkdir("ac_phi");
  for (int r = 1; r<6; ++r) {
    if (rings[r] == 0) { continue; }
    //Changed the binning and upper limit of this histogram as it is overflowing. Is this energy the same as peak->sum
    part_en[r] = get_TH1I(&file, ("part_en_r"+std::to_string(r)), "Valid particle singles energy;", 8192, 0, 65536);

    pg_tdiff[r] = get_TH2I(&file, ("pg_tdiff_r"+std::to_string(r)), 
        ("Particle-gamma time difference (ring "+std::to_string(r)+"); Time (ns); Crystal ID; Counts"), 400, -2000, 2000, 16*4, 0, 16*4);

    pg_ent[r] = get_TH2I(&file, ("pg_ent_r"+std::to_string(r)),
        ("Particle-gamma energy vs time (ring "+std::to_string(r)+"); Time Difference (ns); Energy (keV)"),
        400, -2000, 2000, 8192, 0, 4096) ;

    pg_tent[r] = get_TH2I(&file, ("pg_tent_r"+std::to_string(r)),
        ("Particle-gamma energy vs time (ring "+std::to_string(r)+"); Time Difference (ns); Energy (keV)"),
        400, -2000, 2000, 8192, 0, 4096) ;

    pg_ndent[r] = get_TH2I(&file, ("pg_ndent_r"+std::to_string(r)),
        ("Particle-gamma energy vs time (ring "+std::to_string(r)+", no doppler correction); Time Difference (ns); Energy (keV)"),
        400, -2000, 2000, 4096, 0, 2048) ;

    pg_ent_eff[r] = get_TH2D(&file, ("pg_ent_r"+std::to_string(r)+"_eff"),
        ("Particle-gamma energy vs time (ring "+std::to_string(r)+"); Time Difference (ns); Energy (keV)"),
        400, -2000, 2000, 8192, 0, 4096);

    pg_tent_eff[r] = get_TH2D(&file, ("pg_tent_r"+std::to_string(r)+"_eff"),
			      ("Particle-gamma energy vs time (ring "+std::to_string(r)+"); Time Difference (ns); Energy (keV)"),
			      400, -2000, 2000, 8192, 0, 4096);
  }

  //if ( !file.cd("trin_pid") ) {
    //std::cout << "Couldn't get into trin_pid, making it now" << std::endl;
  //file.mkdir("trin_pid");
  //file.cd("trin_pid");
  //}
  
  if(file.GetDirectory("trin_pid") == nullptr){
    file.mkdir("trin_pid");
    }
  file.cd("trin_pid");

  //TRIN histograms
  TH2I *trin_pid[600]; //most of these will be null
  //TH2I *trin_pid_cut[600]; //most of these will be null
  TH2I *trin_ggpid[600]; //most of these will be null
  TH2I *trin_evspc[600]; //most of these will be null
  TH2I *trin_evsr[600]; //most of these will be null
  TH2I *trin_gevsr[600]; //most of these will be null
  TH2I *trin_enen[600]; //most of these will be null
  TH2I *trin_pepe[600]; //most of these will be null
  TH2I *trin_tata[600]; //most of these will be null
  for (int i=0; i<600; ++i) {
    trin_pid[i] =  NULL;
    //trin_pid_cut[i] = NULL;
    trin_ggpid[i] =  NULL;
    trin_evsr[i] = NULL;
    trin_enen[i] = NULL;
    trin_pepe[i] = NULL;
    trin_tata[i] = NULL;
    trin_gevsr[i] = NULL;
  }
  // zero-degree
  int upper_energy = 10*4096;
  int nbins = 4096;
  for (int r = 0; r < 6; ++r) {
    if (r > 0 && rings[r] == 0) { continue; }
    for (int i = r*100 + 1; i <= r*100+trin_xtls[r]; ++i) {
      trin_pid[i] = get_TH2I(gDirectory, ("trin_pid_"+std::to_string(i)),
			     ("Particle ID "+std::to_string(i)+";Tail Sum; Peak Sum"),
			     4096, 0, upper_energy, 4096, 0, upper_energy);
      //trin_pid_cut[i] = get_TH2I(gDirectory, ("trin_pid_cut_"+std::to_string(i)),
      // ("PID (gated) "+std::to_string(i)+";Tail Sum; Peak Sum"),
      // 4096, 0, upper_energy, 4096, 0, upper_energy);
    //trin_ggpid[i] = get_TH2I(gDirectory, ("trin_ggpid_"+std::to_string(i)),
    //    ("Gamma Gated Particle ID "+std::to_string(i)+";Tail sum; Peak sum"),
    //    4096, 0, upper_energy, 4096, 0, upper_energy);
    trin_evsr[i] = get_TH2I(gDirectory, ("trin_evsr_"+std::to_string(i)),
        ("Energy vs Ratio "+std::to_string(i)+"; Ratio; Energy"),
        4096, 0, 2, 4096, 0, 2*upper_energy);
    trin_gevsr[i] = get_TH2I(gDirectory, ("trin_gevsr_"+std::to_string(i)), ("Energy vs Ratio "+std::to_string(i)+"; Ratio; Energy"),
        4096, 0, 2, 4096, 0, 2*upper_energy);
    //trin_evspc[i] = get_TH2I(gDirectory, ("trin_evspc_"+std::to_string(i)),
    //    ("Energy vs PC0 "+std::to_string(i)+"; Ratio; Energy"),
    //    4096, -0.5,0.5, 4096, 0, 2*upper_energy);
    //trin_enen[i] = get_TH2I(gDirectory, ("trin_enen_"+std::to_string(i)),
    //    ("Energy vs Energy "+std::to_string(i)+"; Energy 1; Energy 2"),
    //    512, 0, 2*upper_energy, 512, 0, 2*upper_energy);
    //trin_tata[i] = get_TH2I(gDirectory, ("trin_tata_"+std::to_string(i)),
    //    ("Tail vs Tail "+std::to_string(i)+"; Tail 1; Tail 2"),
    //    512, 0, upper_energy, 512, 0, upper_energy);
    //trin_pepe[i] = get_TH2I(gDirectory, ("trin_pepe_"+std::to_string(i)),
    //    ("Peak vs Peak "+std::to_string(i)+"; Peak 1; Peak 2"),
    //    512, 0, upper_energy, 512, 0, upper_energy);
  }

  }

  file.cd("../");

  //trace writeout
  //for debug purposes
  int nTraces = 0;
  int nSuperPulse = 0;
  int nTotSuper = 100000;
  TFile *tracefile = new TFile("tracefile.root", "recreate"); 
  TH1D* super_pulse = new TH1D("super_pulse", "Super Pulse", 200, 0, 200);
  TH2D* traces = new TH2D("traces", "traces", 200, 0, 200, 4096, 0, 4096);

  float pg_p[2] = {-86.0, 254.0};  //this is a narrow gate for better SNR with high energy gammas
  float pg_hnp[2] = {300.0, 1400.0};
  float pg_lnp[2] = {-1400.0, -100.0};

  long long unsigned int last_trin_time[600] = {0}; 
  long long unsigned int last_time[13*16] = {0};
  long long unsigned int last_time_clean[13*16] = {0};

  int last_pRing = -1;
  float last_pPhi = -1;

  /*
  TFile *cutfile = NULL;
  cutfile = new TFile("CoulexCuts.root");
  if (cutfile->IsZombie()) {    
    std::cout << "CoulexCuts.root does not exist" << std::endl;
    return -1;
  }
  std::cout << "loaded CoulexCuts.root" << std::endl;
*/

  TCutG* cuts[600];
  int nCuts = 0;
  int cut_ids[600];
  for (int i=0; i<600; ++i) {
    cuts[i] = NULL;
    cut_ids[i] = -1;
  }
  read_cuts("Cuts.cut", &cuts[0], &cut_ids[0], nCuts); 
  for (int i = 0; i < 600; ++i) {
    if (cuts[i]) {
      printf("Loaded cut at ID %d with %d points\n", i, cuts[i]->GetN());
    }
  }

  std::cout << nCuts << " cuts loaded from file" << std::endl;

  std::cout << "Sorting " << dataPaths.size() << " files now" << std::endl;

  //Also for the counter hist
  bool firstTimeSeen = false;
  double t0_seconds = 0.0;
  /////////

  for (int fn = 0; fn<dataPaths.size(); ++fn) {
    std::string listPath = dataPaths.at(fn); 
    int ds = 0;
    PIXIE::Reader reader;
    //set up reader for file(s)
    std::cout << std::endl;
    std::cout << "Opening " << listPath << " for sorting (" << fn+1 << "/" << dataPaths.size() << ")" << std::endl;
    if (reader.openfile(listPath) < 0) { std::cout << "listmode file not found" << std::endl; return -1; } //once again changing open to openfile

    reader.definition.open(defPath);
    reader.definition.read();

    reader.set_coinc(coincWindow);
    std::cout << "coincidence window is " << coincWindow << std::endl;

    reader.RejectSCPU = true;

    time_t now;
    time(&now);
    std::cout << now << std::endl;
    std::cout << "Starting the code block now" << std::endl;

    //THIS IS THE CODE TIM SENT TO TEST!!!!!
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

    //End of Tim's Code Block

    std::cout << "Stopping the code block now" << std::endl;

    reader.start(); //this is for timing
    int nValid[6] = {0}; 
    //int counter = 0;

    while (true) {
      ++ds;
      
      //std::cout << "counter = " << counter <<std::endl;
      //++counter;
      int retval = reader.read();
      if (retval) { 
        std::cout << "Exiting with retval " << retval << std::endl;
        break;
      }

      int eventIndx = reader.eventCtr-1;
      PIXIE::Event *e = &(reader.events[eventIndx]);
      event.Set(reader, eventIndx);  //this sets up Clarion and Trinity objects
      
      //Once again more counts testing/////
      double eventTimeSeconds =
	reader.measurements[e->fMeasurements[0]].eventTime / (3276.8 * 1e9);

      if (!firstTimeSeen) {
	t0_seconds = eventTimeSeconds;
	firstTimeSeen = true;
      }

      double elapsedTime = eventTimeSeconds - t0_seconds;
      if (event.clarion.nClovers > 0) {
	rate_all->Fill(elapsedTime);
      }
      ///////////
      
      //fill histograms
      int nMeas = reader.events[eventIndx].nMeas;
      raw_mult->Fill(reader.measurements[e->fMeasurements[0]].eventTime/(3276.8*1e9), nMeas);
      //This should allow for events in the BEGe beyond just raw
      struct BEGeEventHit {
	int detector = 0;
	double rawEnergy = 0.0;
	unsigned long long timestamp = 0;
	bool clean = false;
      };

      std::vector<BEGeEventHit> begeHits;


      for (int i=0; i<reader.events[eventIndx].nMeas; ++i) {
        auto &meas = reader.measurements[e->fMeasurements[i]];
        int chan = meas.channelNumber;
        int mod = meas.slotID-2;
	int ID   = mod * 16 + chan;
        int tracelength = meas.traceLength;
        double energy = meas.eventEnergy;
        raw_e->Fill(energy, mod*16+chan);

	// Detector map:
	// S2 Ch.15 = global ID 15 = BEGe E1
	// S3 Ch.15 = global ID 31 = BEGe E2
        if (ID == 15 || ID == 31) {
	  const int begeDetector = (ID == 15) ? 1 : 2;
	  const bool begeClean = !meas.finishCode && !meas.outOfRange;

	  if (begeDetector == 1) {
	    bege_e1_raw->Fill(energy);

	    if (begeClean) {
	      bege_e1_raw_clean->Fill(energy);
	    }
	  }
	  else {
	    bege_e2_raw->Fill(energy);

	    if (begeClean) {
	      bege_e2_raw_clean->Fill(energy);
	    }
	  }

	  BEGeEventHit hit;
	  hit.detector = begeDetector;
	  hit.rawEnergy = energy;
	  hit.timestamp = meas.eventTime;
	  hit.clean = begeClean;

	  begeHits.push_back(hit);
	}

	// All channels belonging to slot S14
	if (meas.slotID == 14) {
	  labr_raw->Fill(energy, chan);

	  if (!meas.finishCode && !meas.outOfRange) {
            labr_raw_clean->Fill(energy, chan);
	  }
	}

	// filling calibration histogram
	//raw_e_cal->Fill(energy, mod*16+chan);
        if (mod*16 + chan == 8) {
          raw_e_wt->Fill(meas.eventTime/(3276.8 * 1e9), energy);
        }

        if (meas.finishCode) {
          raw_pu->Fill(mod*16+chan);
        }
        if (meas.outOfRange) {
          raw_or->Fill(mod*16+chan);
        }

        raw_seqtdiff->Fill((double)(meas.eventTime - last_time[mod*16+chan])/3276.8/1000.0, mod*16 + chan);
        last_time[mod*16+chan] = meas.eventTime;
        if (energy > 30.0 && energy < 60000) {
          if (!meas.finishCode && !meas.outOfRange) {
            raw_seqtdiff_clean->Fill((double)(meas.eventTime - last_time_clean[mod*16+chan])/3276.8/1000.0, mod*16 + chan);
            last_time_clean[mod*16+chan] = meas.eventTime;
          } 
        }
        raw_rates->Fill(meas.eventTime/(3276.8*1e9), mod*16+chan);

        for (int j=i+1; j<reader.events[eventIndx].nMeas; ++j) {
          auto &meas2 = reader.measurements[e->fMeasurements[j]];
          int chan2 = meas2.channelNumber;
          int mod2 = meas2.slotID-2;
          raw_hitpat->Fill(mod*16+chan, mod2*16+chan2);          
        }
      }

      //THIS IS TO MAKE SURE THE PASS_PID IS KNOWN BEFORE IT'S USED. THIS IS FOR THE GG GATED HIST
      // ----------------------------
      // NEW: event-level PID gate (compute BEFORE gamma fills)
      // ----------------------------
      int cutMult_pid = 0;

      for (int i = 0; i < event.trinity.nParts; ++i) {
	auto part = &(event.trinity.parts[i]);

	// keep these consistent with how you define "good" particles elsewhere
	if (!part->clean) continue;

	if (cuts[part->GaggID] == NULL) continue;
	if (!cuts[part->GaggID]->IsInside((float)part->tail, (float)part->peak)) continue;

	++cutMult_pid;
      }

      // matches your "if (cutMult != 1) continue;" style strictness
      bool pass_pid = (cutMult_pid == 1);
      bool anti_pid = (cutMult_pid == 0);

      // True event-level no-particle condition.
      // This means the event builder found zero TRINITY particle objects in this event.
      bool no_trinity = (event.trinity.nParts == 0);

      // Optional opposite condition, useful for testing
      bool has_trinity = (event.trinity.nParts > 0);

      //This block should calibrate the BEGe
      for (const auto& hit : begeHits) {
	if (!hit.clean) {
	  continue;
	}

	const auto& cal = begeCal[hit.detector];

	if (!cal.valid) {
	  continue;
	}

	const double calibratedEnergy =
	  cal.intercept + cal.slope * hit.rawEnergy;

	if (calibratedEnergy < 0.0 ||
	    calibratedEnergy >= 4096.0) {
	  continue;
	}

	if (hit.detector == 1) {
	  bege_e1_cal->Fill(calibratedEnergy);

	  if (pass_pid) {
            bege_e1_cal_pid->Fill(calibratedEnergy);
	  }
	}
      }

      //hits, no addback
      /*for (int ih=0; ih<event.clarion.nHits; ++ih) {
        ClarionTrinity::ClarionHit *hit = event.clarion.hits[ih];
        clar_hits->Fill(hit->Energy, (hit->CloverID - 1)*4 + hit->CrystalID);

        if (hit->Suppress == false) {
          clar_hitsclean->Fill(hit->Energy, (hit->CloverID - 1)*4 + hit->CrystalID);
        }
	}*/

      // ----------------------------
      // hits, no addback
      // ----------------------------
      std::vector<const ClarionTrinity::ClarionHit*> hits_by_clover[17];
      std::vector<const ClarionTrinity::ClarionHit*> hits_by_clover_clean[17]; // optional

      for (int ih = 0; ih < event.clarion.nHits; ++ih) {
	ClarionTrinity::ClarionHit* hit = event.clarion.hits[ih];

	int cid = hit->CloverID;
	int xid = (hit->CloverID - 1)*4 + hit->CrystalID;

	// existing
	clar_hits->Fill(hit->Energy, xid);

	// NEW: global no-addback singles spectrum
	if (hit->Energy >= 30 && hit->Energy <= 4096) {
	  h_noaddback_all->Fill(hit->Energy);
	}

	// NEW: per-clover singles (no addback)
	if (cid >= 1 && cid <= 16) {
	  h_clover_xtal_sum[cid]->Fill(hit->Energy);
	  hits_by_clover[cid].push_back(hit);
	  // PID-gated versions (NEW)
	  if (pass_pid) {
	    h_clover_xtal_sum_pid[cid]->Fill(hit->Energy);
	  }
	}

	// existing + NEW clean variants
	if (hit->Suppress == false) {
	  clar_hitsclean->Fill(hit->Energy, xid);

	  if (cid >= 1 && cid <= 16) {
	    h_clover_xtal_sum_clean[cid]->Fill(hit->Energy);
	    hits_by_clover_clean[cid].push_back(hit);
	    if (pass_pid) {
	      h_clover_xtal_sum_clean_pid[cid]->Fill(hit->Energy);
	    }
	  }
	}
      }

      // ----------------------------
      // gamma-gamma: same clover, crystal-crystal (no addback)
      // ----------------------------
      for (int cid = 1; cid <= 16; ++cid) {
	auto& v = hits_by_clover[cid];
	int n = (int)v.size();
	if (n < 2) continue;

	for (int a = 0; a < n; ++a) {
	  
	  for (int b = a + 1; b < n; ++b) {
	    //Can prob remove the line under this comment
	    if (v[a]->CrystalID == v[b]->CrystalID) continue; // require different crystals

	    //THIS IS TO TEST TIMING WINDOW FOR ADDBACK PURPOSES
	    // Calculate gamma-gamma time difference in ns.
	    double dt = gg_dt_ns(v[a], v[b]);

	    double e1_pre = v[a]->Energy;
	    double e2_pre = v[b]->Energy;

	    int xtal1 = v[a]->CrystalID;
	    int xtal2 = v[b]->CrystalID;
	    int pair_index = same_clover_pair_index(xtal1, xtal2);

	    // ----------------------------
	    // Fill timing diagnostics BEFORE prompt cut
	    // ----------------------------

	    gg_sameclover_tdiff_all->Fill(dt);

	    if (pass_pid) {
	      gg_sameclover_tdiff_pid->Fill(dt);
	    }

	    gg_sameclover_e_vs_tdiff->Fill(dt, e1_pre);
	    gg_sameclover_e_vs_tdiff->Fill(dt, e2_pre);

	    if (pair_index >= 0) {
	      gg_sameclover_pair_vs_tdiff->Fill(dt, pair_index);
	    }

	    // ----------------------------
	    // Existing prompt cut
	    // ----------------------------
	    if (!gg_prompt(v[a], v[b])) continue;             // require prompt timing
	     

	    double e1 = v[a]->Energy;
	    double e2 = v[b]->Energy;

	    gg_sameclover[cid]->Fill(e1, e2);
	    gg_sameclover[cid]->Fill(e2, e1); // symmetric (optional)

	    if (pass_pid) {
	      gg_sameclover_pid[cid]->Fill(e1, e2);
	      gg_sameclover_pid[cid]->Fill(e2, e1);
	    }
	  }
	}
      }

      // Optional: clean gg
      for (int cid = 1; cid <= 16; ++cid) {
	auto& v = hits_by_clover_clean[cid];
	int n = (int)v.size();
	if (n < 2) continue;

	for (int a = 0; a < n; ++a) {
	  for (int b = a + 1; b < n; ++b) {
	    if (v[a]->CrystalID == v[b]->CrystalID) continue;
	    if (!gg_prompt(v[a], v[b])) continue;

	    gg_sameclover_clean[cid]->Fill(v[a]->Energy, v[b]->Energy);
	    gg_sameclover_clean[cid]->Fill(v[b]->Energy, v[a]->Energy);

	    if (pass_pid) {
	      gg_sameclover_pid_clean[cid]->Fill(v[a]->Energy, v[b]->Energy);
	      gg_sameclover_pid_clean[cid]->Fill(v[b]->Energy, v[a]->Energy);
	    }
	  }
	}
      }
      
      // ----------------------------
      // gamma-gamma: different clovers (crystal hits, no addback)
      // ----------------------------
      for (int cid1 = 1; cid1 <= 16; ++cid1) {
	auto& v1 = hits_by_clover[cid1];
	if (v1.empty()) continue;

	for (int cid2 = cid1 + 1; cid2 <= 16; ++cid2) {
	  auto& v2 = hits_by_clover[cid2];
	  if (v2.empty()) continue;

	  for (size_t a = 0; a < v1.size(); ++a) {
	    for (size_t b = 0; b < v2.size(); ++b) {
	      if (!gg_prompt(v1[a], v2[b])) continue;


	      double e1 = v1[a]->Energy;
	      double e2 = v2[b]->Energy;

	      gg_cloverclover->Fill(e1, e2);
	      gg_cloverclover->Fill(e2, e1); // symmetric optional

	      if (pass_pid) {
		gg_cloverclover_pid->Fill(e1, e2);
		gg_cloverclover_pid->Fill(e2, e1);
	      }
	    }
	  }
	}
      }

      // Clean versions (Suppress==false hits only)
      for (int cid1 = 1; cid1 <= 16; ++cid1) {
	auto& v1 = hits_by_clover_clean[cid1];
	if (v1.empty()) continue;

	for (int cid2 = cid1 + 1; cid2 <= 16; ++cid2) {
	  auto& v2 = hits_by_clover_clean[cid2];
	  if (v2.empty()) continue;

	  for (size_t a = 0; a < v1.size(); ++a) {
	    for (size_t b = 0; b < v2.size(); ++b) {

	      if (!gg_prompt(v1[a], v2[b])) continue;

	      double e1 = v1[a]->Energy;
	      double e2 = v2[b]->Energy;

	      gg_cloverclover_clean->Fill(e1, e2);
	      gg_cloverclover_clean->Fill(e2, e1);

	      if (pass_pid) {
		gg_cloverclover_clean_pid->Fill(e1, e2);
		gg_cloverclover_clean_pid->Fill(e2, e1);
	      }
	    }
	  }
	}
      }

      //clovers, addback
      //So addback testing works
      std::vector<const ClarionTrinity::Gamma*> gammas_by_clover[17];
      std::vector<const ClarionTrinity::Gamma*> gammas_by_clover_clean[17];


      for (int iCl=0; iCl<event.clarion.nClovers; ++iCl) {
        auto clover_i = &(event.clarion.clovers[iCl]);
        int ID_i = clover_i->CloverID;

        for (int hi=0; hi<clover_i->nHits; ++hi) {
          for (int bj=0; bj<clover_i->nBGOs; ++bj) {
            double tdiff = (double) ((long long int)clover_i->hits[hi].timestamp - (long long int) clover_i->bgos[bj].timestamp)/3276.8;

            if (clover_i->bgos[bj].RawEnergy < 30 || clover_i->bgos[bj].RawEnergy > 60000) { continue; }
            if (clover_i->hits[hi].Energy < 30 || clover_i->hits[hi].RawEnergy > 60000) { continue; }

            clar_bgotdiff->Fill(tdiff, (ID_i-1)*4 + clover_i->hits[hi].CrystalID);
          }
        }

        for (int iGam=0; iGam<clover_i->nGammas; ++iGam) {
          //not sure if this should be here or not
          //will keep for consistency with Mitch for now MIGHT NOT NEED THIS 
	  // if (iGam > 0) { continue; }

          ClarionTrinity::Gamma *gam_i = &(clover_i->gammas[iGam]);

	  // NEW: addback singles test + save gamma by clover
	  if (ID_i >= 1 && ID_i <= 16) {
	    if (gam_i->Energy >= 30 && gam_i->Energy <= 4096) {
	      h_addback_all->Fill(gam_i->Energy);

	      // Addback spectrum for events where TRINITY did not fire
	      if (no_trinity) {
		h_addback_no_trinity->Fill(gam_i->Energy);
	      }

	      // Optional sanity-check spectrum for events where TRINITY did fire
	      if (has_trinity) {
		h_addback_with_trinity->Fill(gam_i->Energy);
	      }
	      gammas_by_clover[ID_i].push_back(gam_i);
	    }
	  }

          int iMaxEnInd = gam_i->MaxEnInd;
          int idi = (ID_i-1)*4 + clover_i->hits[iMaxEnInd].CrystalID;

          clar_edirty->Fill(gam_i->Energy, idi);

	  // Ungated singles for calibration (works on pure 60Co runs)
	  h_calib_en_all->Fill(gam_i->Energy);

	  // Anti-PID singles for background-line checking
	  if (anti_pid) {
	    h_bg_sum->Fill(gam_i->Energy);
	  }

          if (!Clean(clover_i)) { continue; }

	  // NEW: no PID cut and no particle-gamma timing cut
	  g_e_ungated->Fill(gam_i->Energy);

	  // Clean no-TRINITY addback spectrum
	  if (no_trinity && gam_i->Energy >= 30 && gam_i->Energy <= 4096) {
	    h_addback_no_trinity_clean->Fill(gam_i->Energy);
	  }

	  if (ID_i >= 1 && ID_i <= 16) {
	    if (gam_i->Energy >= 30 && gam_i->Energy <= 4096) {
	      gammas_by_clover_clean[ID_i].push_back(gam_i);
	    }
	  }

          clar_e->Fill(gam_i->Energy, idi);

          double efficiencyID = 0.0;
          if (gam_i->Energy > 50.0) {
            efficiencyID = clarion.conf.Efficiency(gam_i->Energy, ID_i);
          }

          clar_e_eff->Fill(gam_i->Energy, idi, 1.0/efficiencyID);
        }
      }

      // ----------------------------
      // gamma-gamma: different clovers using addback gammas
      // ----------------------------
      for (int cid1 = 1; cid1 <= 16; ++cid1) {
	auto &v1 = gammas_by_clover[cid1];
	if (v1.empty()) continue;

	for (int cid2 = cid1 + 1; cid2 <= 16; ++cid2) {
	  auto &v2 = gammas_by_clover[cid2];
	  if (v2.empty()) continue;

	  for (size_t a = 0; a < v1.size(); ++a) {
	    for (size_t b = 0; b < v2.size(); ++b) {

	      if (!gg_prompt(v1[a], v2[b])) continue;

	      double e1 = v1[a]->Energy;
	      double e2 = v2[b]->Energy;

	      gg_cloverclover_ab->Fill(e1, e2);
	      gg_cloverclover_ab->Fill(e2, e1);  // symmetric fill

	      // NEW: fill only when TRINITY did not fire
	      if (no_trinity) {
		gg_cloverclover_ab_no_trinity->Fill(e1, e2);
		gg_cloverclover_ab_no_trinity->Fill(e2, e1);
	      }

	      if (pass_pid) {
		gg_cloverclover_ab_pid->Fill(e1, e2);
		gg_cloverclover_ab_pid->Fill(e2, e1);
	      }
	    }
	  }
	}
      }
      
      // ----------------------------
      // gamma-gamma: different clovers using CLEAN addback gammas
      // ----------------------------
      for (int cid1 = 1; cid1 <= 16; ++cid1) {
	auto &v1 = gammas_by_clover_clean[cid1];
	if (v1.empty()) continue;

	for (int cid2 = cid1 + 1; cid2 <= 16; ++cid2) {
	  auto &v2 = gammas_by_clover_clean[cid2];
	  if (v2.empty()) continue;

	  for (size_t a = 0; a < v1.size(); ++a) {
	    for (size_t b = 0; b < v2.size(); ++b) {

	      if (!gg_prompt(v1[a], v2[b])) continue;

	      double e1 = v1[a]->Energy;
	      double e2 = v2[b]->Energy;

	      gg_cloverclover_ab_clean->Fill(e1, e2);
	      gg_cloverclover_ab_clean->Fill(e2, e1);

	      // NEW: clean gamma-gamma with no TRINITY
	      if (no_trinity) {
		gg_cloverclover_ab_clean_no_trinity->Fill(e1, e2);
		gg_cloverclover_ab_clean_no_trinity->Fill(e2, e1);
	      }

	      // optional PID-gated clean version
	       if (pass_pid) {
	         gg_cloverclover_ab_clean_pid->Fill(e1, e2);
	         gg_cloverclover_ab_clean_pid->Fill(e2, e1);
	       }
	    }
	  }
	}
      }

      int dirtyMult = 0;  //number of particles not passing clean conditions
      int cleanMult = 0;  //number of particles passing clean conditions
      int cutMult = 0;    //number of (clean) particles passing cut
      int cutInds[100] = {0}; //indexes of valid HI particles
      int iMaxEn = -1;
      double maxEn = 0;
      for (int i=0; i<event.trinity.nParts; ++i) { //-> vector of GaGG objects 
        auto part = &(event.trinity.parts[i]);
        int ID = part->GaggID;       
        if (event.trinity.parts[i].fired1 && event.trinity.parts[i].fired2) {
          float diff1 = (float)part->background1 - (float)part->postbackground1;
          float diff2 = (float)part->background2 - (float)part->postbackground2;

          if (part->energy > 400) {
            trin_backdiff1->Fill((float)part->background1 - (float)part->postbackground1, ID); 
            trin_backdiff2->Fill((float)part->background2 - (float)part->postbackground2, ID); 
          }

          if (part->oor1 || part->oor2) {
            std::cout << "Warning! GAGG " << ID << " has gone out of range" << std::endl;
          }

          if (part->peak < 0) {
            part->clean = 0;
          }

          if (part->energy < 300) {
            part->clean = 0;
          }

          double tdiff = (double)((long long int)(part->time1 - part->time2))/3276.8;
          trin_tdiff->Fill(tdiff, ID);
          if (abs(tdiff) > tdiff_thresh) {
            //std::cout << tdiff << std::endl; 
            part->clean = 0;
          }         
          if (ID == 1) { //special rules for zero degree
            double en1 = part->energy1;
            double en2 = part->energy2;
            if ((abs(en1-en2)/part->energy) > 0.04) {
              part->clean = 0;
            }
          }
        } //double fired end 

        if (event.trinity.parts[i].clean) {
          ++cleanMult;
          //now that we know it's clean, do time alignment and re-calculate peak and tail

          auto &meas1 = reader.measurements[e->fMeasurements[part->iMeas1]];
          auto &meas2 = reader.measurements[e->fMeasurements[part->iMeas2]];

          if (cuts[part->GaggID] == NULL) { continue; }
          if (!cuts[part->GaggID]->IsInside((float)part->tail, (float)part->peak)) { continue; }

          cutInds[cutMult] = i;
          ++cutMult;

          //if there's more than one particle inside the cut, record the max energy one          
          if (part->energy > maxEn) {
            maxEn = part->energy;
            iMaxEn = i;
          }
        }
        else {
          ++dirtyMult;
        }
      } //particle loop end

      // --- NEW: event-level PID gate for gamma histograms ---
      // THIS IS HELPING APPLY GATES TO THE GG HISTOGRAMS
      //bool pass_pid = (cutMult == 1);   // matches your later "if (cutMult != 1) continue;" This is a redefinition that shouldn't be here but keeping it for the time being

      part_dmult->Fill(dirtyMult);
      part_mult->Fill(cleanMult);
      part_himult->Fill(cutMult);
      if ((dirtyMult < 10) && (cleanMult-cutMult < 10 && cutMult < 10)) {
        part_combmult->Fill(dirtyMult + (cleanMult-cutMult)*10 + cutMult*100);
      }

      for (int i=0; i<event.trinity.nParts; ++i) {
        auto part = &(event.trinity.parts[i]);
        int ID = part->GaggID;
        part_dmult_id->Fill(dirtyMult, ID);
        part_mult_id->Fill(cleanMult, ID);
        part_himult_id->Fill(cutMult, ID);

        if (!part->clean) { continue; }

        for (int j=i+1; j<event.trinity.nParts; ++j) {
          auto part_j = event.trinity.parts[j];
          if (!part_j.clean) { continue; }
          trin_hitpat->Fill(ID, part_j.GaggID);
          trin_hitpat->Fill(part_j.GaggID, ID);
        }

        trin_rates->Fill(part->time/(3276.8*1e9), ID);        
        trin_seqtdiff->Fill((double)(part->time - last_trin_time[ID])/(3276800.0), ID);


        //fill Trinity histograms
        last_trin_time[ID] = part->time;
        trin_e->Fill(part->energy, ID);
        if (ID == 1) { //zero-degree detector
          trin_z_wt->Fill(part->time/(3276.8*1e9), part->energy);
        }

        if (part->energy > 400) {
          trin_gainmatch->Fill((float)part->peak1/(float)part->peak2, ID);
        }

        double tdiff = (double)((long long int)(part->time1 - part->time2))/3276.8;
	//std::cout << ID << std::endl;

        if (trin_pid[ID] != NULL) {
          trin_pid[ID]->Fill(part->tail, part->peak);
	  //std::cout << ID << " part->tail " << (double)part->tail << " part->peak " << (double)part->peak << std::endl;
	  
          //std::cout << tdiff << std::endl;
        }
        else {
          std::cout << "Warning: Trinity ID " << ID << " in data stream but PID histogram not defined";
        }

	//if 
        
        if (trin_evsr[ID] != NULL) {
          if (part->peak <= 0) { continue; }
          trin_evsr[ID]->Fill((double)part->tail/(double)part->peak, part->energy);
        }

        //TRACE WRITEOUT//            
        //TRACE WRITEOUT END//
	//ChatGPT Test script
	// Debug: print the first few entries for GaggID 401
	//if (part->GaggID == 401) {
	 //static int printed = 0;  // limit how many we print
	//if (printed < 20) {
	//  std::cout << "[DEBUG] GaggID=" << part->GaggID
	//	      << "  tail=" << part->tail
	//	      << "  peak=" << part->peak << std::endl;
	//  printed++;
	//}
	//}


	

        if (cuts[part->GaggID] == NULL) { continue; }
        if (!cuts[part->GaggID]->IsInside((float)part->tail, (float)part->peak)) { continue; }

	//Filling Cut PID
	//if (trin_pid_cut[ID]) trin_pid_cut[ID]->Fill(part->tail, part->peak);
      }

      for (int i=0; i<cutMult; ++i) {
        auto part = &(event.trinity.parts[cutInds[i]]);

        //if more than one particle inside cut (???)
        if (cutMult != 1) { continue; }
        //even more stringent: only one valid GAGG
        //if (cleanMult > 1) { continue; }

        int ring = (int)(((float)part->GaggID)/100.0);
        ++nValid[ring];
        part_singles->Fill(part->GaggID);
        part_singles_wt->Fill(part->time/(3276.8*1e9), part->GaggID);
        part_en[ring]->Fill(part->energy);
        if (part->GaggID == 416) {
          part_en_wt->Fill(part->time/(3276.8*1e9), part->energy);
        }

        float pTheta = event.trinity.conf.theta[part->GaggID];
        float pPhi = event.trinity.conf.phi[part->GaggID];
	//TEMP DEFINITIONS
	//////
	/////
	/////
	/////
	//float tTheta = 0.0;
	//float tPhi = 0.0;
	//////
	//////
	/////
	/////
	/////

        //PARTICLE-GAMMA
        for (int iCl=0; iCl<event.clarion.nClovers; ++iCl) {
          auto clover_i = &(event.clarion.clovers[iCl]);

          //gamma clean condition
          if (!Clean(clover_i)) { continue; }

          for (int iGam=0; iGam<event.clarion.clovers[iCl].nGammas; ++iGam) {
            ClarionTrinity::Gamma *gam_i = &(clover_i->gammas[iGam]);

            int iMaxEnInd = gam_i->MaxEnInd;
            int CloverID = clover_i->CloverID;
            int CrystalID = clover_i->hits[iMaxEnInd].CrystalID;
            int idi = 4*(CloverID - 1) + CrystalID;

            double tdiff = (double)((long long int)(gam_i->timestamp - part->time))/3276.8;

            pg_ndent[ring]->Fill(tdiff, gam_i->Energy);
	    //testing this line right here
	    pg_tdiff[ring]->Fill(tdiff, idi);

	    

            //get corresponding Ti theta:
            if (pPhi > 360.0) { pPhi = pPhi - 360.0; }
            float gTheta = event.clarion.conf.theta[CloverID][CrystalID];
            float gPhi = event.clarion.conf.phi[CloverID][CrystalID];

            int gInd = -1;
            float gThetaC = event.clarion.conf.theta[CloverID][4];
            for (int gi = 0; gi < 4; ++gi) {
              if (hpge_angs[gi] - 5.0 <=  gThetaC && gThetaC <= hpge_angs[gi] + 5.0) {
                gThetaC = hpge_angs[gi];
                gInd = gi;
              }
            }

            //DOPPLER CORRECTION GOES HERE Jan 28 2026 Uncommenting from here to double e_dop_alt = gam_i->Energy/beta_corr_alt; for addback testing (JAN 28 Changed tTheta to pTheta to test it, will get the real value soon)
            
	    /* double pgcos = std::sin(pTheta*pi/180.0)*std::sin(gTheta*pi/180.0)*std::cos(pPhi*pi/180. - gPhi*pi/180.)
              + std::cos(pTheta*pi/180.0)*std::cos(gTheta*pi/180.0);
            double targ_pgcos = std::sin(tTheta*pi/180.0)*std::sin(gTheta*pi/180.0)*std::cos(tPhi*pi/180. - gPhi*pi/180.)
              + std::cos(tTheta*pi/180.0)*std::cos(gTheta*pi/180.0);

            if (pgcos > 1 || pgcos < -1) { 
              std::cout << "Warning! pgcos(theta) unphysical: " << pgcos  << std::endl;
              std::cout << "    pTheta: " << pTheta << "  pPhi: " << pPhi << std::endl;
              std::cout << "    gTheta: " << gTheta << "  gPhi: " << gPhi << std::endl;
               }
            if (targ_pgcos > 1 || targ_pgcos < -1) { 
              std::cout << "Warning! target pgcos(theta) unphysical: " << targ_pgcos  << std::endl; 
              std::cout << "    tTheta: " << tTheta << "  tPhi: " << tPhi << std::endl;
              std::cout << "    gTheta: " << gTheta << "  gPhi: " << gPhi << std::endl;
            }
            double dPhi = pPhi - gPhi;
            if (dPhi < 0) { dPhi += 360; }
            if (dPhi > 360.0) { dPhi -= 360; }

            double beta = betas[ring];
            double beta_corr = std::sqrt(1.0 - beta * beta)/(1.0 - beta * pgcos);
            double targ_beta = targ_betas[ring];
            double targ_beta_corr = std::sqrt(1.0 - targ_beta * targ_beta)/(1.0 - targ_beta * targ_pgcos);

            double e_dop = gam_i->Energy/beta_corr;
            double targ_e_dop = gam_i->Energy/targ_beta_corr;

            double beta_alt = betas_alt[ring];
            double beta_corr_alt = std::sqrt(1.0 - beta_alt * beta_alt)/(1.0 - beta_alt * pgcos);

            double e_dop_alt = gam_i->Energy/beta_corr_alt; */


            
            if (pg_p[0] < tdiff && tdiff < pg_p[1]) {
	      
	      //Filling cut and tdiff gated spectra
	      pg_e_prompt_cut->Fill(gam_i->Energy);

              for (int jCl = iCl+1; jCl < event.clarion.nClovers; ++jCl) {
                auto clover_j = &(event.clarion.clovers[jCl]);
                if (!Clean(clover_j)) { continue; }
                for (int jGam = 0; jGam <event.clarion.clovers[jCl].nGammas; ++jGam) {
                  ClarionTrinity::Gamma *gam_j = &(clover_j->gammas[jGam]);

                  int jMaxEnInd = gam_j->MaxEnInd;
                  int jCloverID = clover_j->CloverID;
                  int jCrystalID = clover_j->hits[jMaxEnInd].CrystalID;
                  int idj = 4*(jCloverID - 1) + jCrystalID;

                  double tdiff_j = (double)((long long int)(gam_j->timestamp - part->time))/3276.8;


                  if (pg_p[0] < tdiff_j && tdiff_j < pg_p[1]) {
                    float gTheta_j = event.clarion.conf.theta[jCloverID][jCrystalID];
                    float gPhi_j = event.clarion.conf.phi[jCloverID][jCrystalID];
                    // DO DOPPLER CORRECTION HERE uncommented up tp pgg->Fill(e_dop_j, e_dop)
		    /* double pgcos_j = std::sin(pTheta*pi/180.0)*std::sin(gTheta_j*pi/180.0)*std::cos(pPhi*pi/180. - gPhi_j*pi/180.)
                      + std::cos(pTheta*pi/180.0)*std::cos(gTheta_j*pi/180.0);

                    double beta_corr_j = std::sqrt(1.0 - beta * beta)/(1.0 - beta * pgcos_j);
                    double e_dop_j = gam_j->Energy/beta_corr_j;
                    pgg->Fill(e_dop, e_dop_j);
                    pgg->Fill(e_dop_j, e_dop); */
		    
                  }
                }
              }
            }


            if ((pg_hnp[0] < tdiff && tdiff < pg_hnp[1]) || (pg_lnp[0] < tdiff && tdiff < pg_lnp[1])) {
              //non-prompt histograms here
            }

            //USE DOPPLER-CORRECTED VALUES reintroduced pg_ent and pg_ent (jan 28 uncommented pg_ent and pg_tent for both lines below and in if
	    // pg_ent[ring]->Fill(tdiff, e_dop);
	    //pg_tent[ring]->Fill(tdiff, targ_e_dop);
            double efficiency = 0.0;
            double efficiencyID = 0.0;
            if (gam_i->Energy > 50.0) {
              efficiency = clarion.conf.Efficiency(gam_i->Energy);
              efficiencyID = clarion.conf.Efficiency(gam_i->Energy, CloverID);
            }
            if (efficiency > 0.0) {
              //USE DOPPLER-CORRECTED VALUES
              //pg_ent_eff[ring]->Fill(tdiff, e_dop, 1.0/efficiency);
	      // pg_tent_eff[ring]->Fill(tdiff, targ_e_dop, 1.0/efficiency);
            }
          } //gamma loop          
        } //clover loop
      } //cutMult loop

      if (ds%10000 == 0) {
        reader.printUpdate();
      }

    } //reader loop end
    time_t finaltime;
    time(&finaltime);
    std::cout << "Total fill time = " << finaltime-reader.starttime << " s " << std::endl;
    reader.printSummary();    

    std::cout << "Valid GAGG: " << std::endl;
    std::cout << "  R1 : " << nValid[1] << std::endl;
    std::cout << "  R2 : " << nValid[2] << std::endl;
    std::cout << "  R3 : " << nValid[3] << std::endl;
    std::cout << "  R4 : " << nValid[4] << std::endl;
    std::cout << "  R5 : " << nValid[5] << std::endl;
  }  //file loop end

  file.cd();

  tracefile->Purge();  
  tracefile->Write();
  tracefile->Close(); 
  h_calib_en_all->Write();

  //AGAIN CALIBRATION TESTING
  if (makeCalPlots) {
    std::cout << "CalCheck: using calfile=" << calfile
	      << " centroids=" << centroidTSV << "\n";
    MakeCalResidualPlots_Linear(&file, centroidTSV, calfile, true);
}
  ///CALIBRATION TESTING

  file.Purge();
  file.GetDirectory("trin_pid")->Purge();
  file.Write(NULL, TObject::kOverwrite);
  file.Purge();
  file.GetDirectory("trin_pid")->Purge();
  file.ls();
  file.Close();
}

