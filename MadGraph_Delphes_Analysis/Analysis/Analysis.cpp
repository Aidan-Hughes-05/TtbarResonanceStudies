#include "Analysis.h"
#include "Logger.h"
#include "Kinematics.h"
#include "Constants.h"

#include <ROOT/RDataFrame.hxx>
#include <TLorentzVector.h>
#include <TFile.h>
#include <vector>
#include <string>
#include <map>

void PartonAnalysisHistograms(const std::string& fileName, const std::string& outputFileName, float lumiWeight, bool doWeights){
    ROOT::RDataFrame df("mytree", fileName);
    ROOT::RDF::RNode node = df; 
    
    // 1) Build the top container
    // Status = 1 means final state particles
    auto myFilter = [](std::vector<int> pdgid, std::vector<int> status) { 
        std::vector<int> indices{};
        for (int i{} ; i<pdgid.size() ; i++){ 
            if (abs(pdgid.at(i)) == 6 && status.at(i) == 1) indices.push_back(i); 
        } 
        return indices; 
    };

    // Define the dummy column for the weight or apply lumi scaling
    std::string weightColumn = "";
    if (!doWeights){
        node = node.Define("dummyWeight",[](const std::vector<int>& pdgid) { return 1.0; },{"pdgID"});
        weightColumn = "dummyWeight";
    } else {
        node = node.Define("lumiWeight",[&lumiWeight](float weightColumn) { return lumiWeight*weightColumn; },{"eventWeight"});
        weightColumn = "lumiWeight";
    }

    node = node.Define("top_indices",myFilter,{"pdgID","pdgStatus"});

    auto getSize = [](std::vector<int> vec) { return vec.size(); };
    node = node.Define("top_indices_size",getSize,{"top_indices"});

    auto buildParticlesFromIndices = [](std::vector<int>& indices, std::vector<float>& px, std::vector<float>& py, std::vector<float>& pz, std::vector<float>& E){
        std::vector<TLorentzVector> particles;
        for (auto i : indices){
            TLorentzVector p;
            p.SetPxPyPzE(px.at(i), py.at(i), pz.at(i), E.at(i));
            particles.push_back(p);
        }

        // Order by pt
        std::sort(particles.begin(), particles.end(), [](const TLorentzVector& a, const TLorentzVector& b) { return a.Pt() > b.Pt(); });

        return particles;
    };
    node = node.Define("Tops",buildParticlesFromIndices,{"top_indices","px","py","pz","E"});

    // 3) Define jets and leptons pt
    node = node.Define("Top_1_pt","Tops.at(0).Pt()");
    node = node.Define("Top_2_pt","Tops.at(1).Pt()");

    // Some more kinematics
    // 4)
    auto deltaRapidity = [](std::vector<TLorentzVector> particles) { return abs(particles.at(0).Rapidity() - particles.at(1).Rapidity()); };
    // 5)
    auto invariantMassTwoLeading = [](std::vector<TLorentzVector> particles) { return (sqrt( 2.0*particles.at(0).Dot(particles.at(1))) ); };
    // 6)
    auto dilepCentrality = [](std::vector<TLorentzVector> leptons, std::vector<TLorentzVector> jets) { 
        float lepton_xi = (leptons.at(0) + leptons.at(1)).Rapidity();
        float dijet_xi = jets.at(0).Rapidity() + jets.at(1).Rapidity();
        float delta_y = abs(jets.at(0).Rapidity() - jets.at(1).Rapidity());
        float centrality = abs((lepton_xi - 0.5*dijet_xi)/delta_y);
        return centrality;
    };
    // 7)
    auto dileptonPt = [](std::vector<TLorentzVector> leptons) { 
        if (leptons.size() < 2) return 0.0;
        return (leptons.at(0) + leptons.at(1)).Pt(); 
    };
    // 8)
    auto deltaPhi = [](std::vector<TLorentzVector> particles) { return Kinematics::del_phi(particles.at(0).Phi(),particles.at(1).Phi()); };

    // Create histogram container
    std::vector histograms = {node.Histo1D("top_indices", weightColumn)};
    histograms.push_back(node.Histo1D("top_indices_size", weightColumn));
    histograms.push_back(node.Histo1D({"Top_1_pt_basic", "Top_1_pt_basic", 500, 0, 500},"Top_1_pt", weightColumn));
    histograms.push_back(node.Histo1D({"Top_2_pt_basic", "Top_2_pt_basic", 500, 0, 500},"Top_2_pt", weightColumn));

    node = node.Define("deltaPhi_Tops",deltaPhi,{"Tops"});
    node = node.Define("deltaRapidity",deltaRapidity,{"Tops"});
    node = node.Define ("mttbar",invariantMassTwoLeading,{"Tops"});
    node = node.Define("diTopPt",dileptonPt,{"Tops"});

    histograms.push_back(node.Histo1D({"deltaPhi_Tops_basic", "deltaPhi_Tops_basic", 64, 0, 3.2},"deltaPhi_Tops", weightColumn));
    histograms.push_back(node.Histo1D({"deltaRapidity_basic", "deltaRapidity_basic", 100, 0, 10},"deltaRapidity", weightColumn));
    histograms.push_back(node.Histo1D({"mttbar_basic", "mttbar_basic", 1000, 0, 1000},"mttbar", weightColumn));
    histograms.push_back(node.Histo1D({"diTopPt_basic", "diTopPt_basic", 1000,0,1000},"diTopPt", weightColumn));


    // Save histograms
    auto output = TFile::Open(outputFileName.c_str(), "RECREATE");
    output->cd();
    output->mkdir("histograms");
    output->cd("histograms");
    for (auto h : histograms){
        h->Write();
    }
    output->Close();
}

void DelphesAnalysis(const std::string& fileName, const std::string& outputFileName, bool doWeights){
    ROOT::RDataFrame df("Delphes", fileName);
    ROOT::RDF::RNode node = df; 

    // Print column names
    Constants::g_LOG(LogLevel::INFO, "Columns available:");
    std::vector<std::string> columnNames = node.GetColumnNames();
    for (auto name : columnNames){
        std::cout << name << std::endl;
    }

    // Define the dummy column for the weight
    std::string weightColumn = "EvtW";
    float nEvents = 50000;
    float luminosity = 139.0;
    if (!doWeights){
        node = node.Define(weightColumn,[](const std::vector<float>& dummyV) { return 1.0; },{"Event.Weight"});
    } else {
        node = node.Define("SampleWeight",[&nEvents, &luminosity](const ROOT::VecOps::RVec<float>& eventWeight) { return eventWeight[10] * 1000 * luminosity / (nEvents * eventWeight[10]); },{"Event.Weight"});
        node = node.Define(weightColumn,[&nEvents, &luminosity](const ROOT::VecOps::RVec<float>& eventWeight, float sampleWeight) { return eventWeight[0] * sampleWeight; },{"Event.Weight","SampleWeight"});
    }


    // Define functions
    // 1)
    auto getSize = [](ROOT::VecOps::RVec<Float_t> vec) { return vec.size(); };
    // 2)
    auto buildParticlesZeroMass = [](ROOT::VecOps::RVec<Float_t>& pt){
        ROOT::VecOps::RVec<Float_t> m;
        for (std::size_t i{}; i < pt.size(); i++){
            m.push_back(0.0f);
        }
        return m;
    };
    // 3)
    auto buildParticlesFromContainer = [](ROOT::VecOps::RVec<Float_t>& pt, ROOT::VecOps::RVec<Float_t>& eta, ROOT::VecOps::RVec<Float_t>& phi, ROOT::VecOps::RVec<Float_t>& m){
        ROOT::VecOps::RVec<TLorentzVector> particles;
        for (std::size_t i{}; i < pt.size(); i++){
            TLorentzVector p;
            p.SetPtEtaPhiM(pt.at(i), eta.at(i), phi.at(i), m.at(i));
            particles.push_back(p);
        }
        return particles;
    };
    // 4)
    auto deltaPhiMuons = [](ROOT::VecOps::RVec<TLorentzVector> particles) { return Kinematics::del_phi(particles.at(0).Phi(),particles.at(1).Phi()); };
    // 5) 
    auto isPtOrdered = [](ROOT::VecOps::RVec<Float_t> pt) { 
        if (pt.size() < 2) Constants::g_LOG(LogLevel::ERROR, "Function called with less than two particles to compare!");
        return pt.at(0) > pt.at(1); 
    };
    // 6)
    auto numberOfHeavyJets = [](ROOT::VecOps::RVec<unsigned int> btagFlagVec) { 
        int numberBJets = 0;
        for (auto btag : btagFlagVec){
            if (btag == 1) numberBJets++;
        }
        return numberBJets;
    };
    // 7)
    auto deltaRapidity = [](ROOT::VecOps::RVec<TLorentzVector> particles) { return abs(particles.at(0).Rapidity() - particles.at(1).Rapidity()); };
    // 8)
    auto invariantMassTwoLeading = [](ROOT::VecOps::RVec<TLorentzVector> particles) { return (sqrt( 2.0*particles.at(0).Dot(particles.at(1))) ); };
    // 9)
    auto isJetInGap = [](const TLorentzVector& j1,const TLorentzVector& j2, const TLorentzVector& testJet) {
        float delta_y_j1j2 = abs(j1.Rapidity()-j2.Rapidity());
        float delta_y_j1test = abs(j1.Rapidity()-testJet.Rapidity());
        float delta_y_j2test = abs(j2.Rapidity()-testJet.Rapidity());
        return (delta_y_j1test>delta_y_j1j2 || delta_y_j2test>delta_y_j1j2);
    };

    auto nJetsInGap = [&isJetInGap](ROOT::VecOps::RVec<TLorentzVector> jets) { 
        int gapJets = 0;
        TLorentzVector j1 = jets.at(0);
        TLorentzVector j2 = jets.at(1);

        if (jets.size() == 2) return 0;

        for (std::size_t i{2}; i < jets.size(); i++){
            TLorentzVector testJet = jets.at(i);
            if (isJetInGap(j1,j2,testJet)) gapJets++;
        }
        return gapJets;
    };
    // 9) 
    auto ptBalance = [&isJetInGap](ROOT::VecOps::RVec<TLorentzVector> leptons, ROOT::VecOps::RVec<TLorentzVector> jets) { 
        float scalarPtSum{0};
        TLorentzVector vectorPtSum;
        for (std::size_t i{}; i < leptons.size(); i++){
            if ( i > 2) break;
            scalarPtSum += leptons.at(i).Pt();
            vectorPtSum += leptons.at(i);
        }

        // Jets
        for (std::size_t i{}; i < jets.size(); i++){

            if ( i>1 ){
                if ( isJetInGap(jets.at(0),jets.at(1),jets.at(i)) ){
                    scalarPtSum += jets.at(i).Pt();
                    vectorPtSum += jets.at(i);
                }
            } else {
                scalarPtSum += jets.at(i).Pt();
                vectorPtSum += jets.at(i);
            }
        }
        return vectorPtSum.Pt()/scalarPtSum;
    };
    // 10)
    auto dilepCentrality = [](ROOT::VecOps::RVec<TLorentzVector> leptons, ROOT::VecOps::RVec<TLorentzVector> jets) { 
        float lepton_xi = (leptons.at(0) + leptons.at(1)).Rapidity();
        float dijet_xi = jets.at(0).Rapidity() + jets.at(1).Rapidity();
        float delta_y = abs(jets.at(0).Rapidity() - jets.at(1).Rapidity());
        float centrality = abs((lepton_xi - 0.5*dijet_xi)/delta_y);
        return centrality;
    };
    // 11)
    auto dileptonPt = [](ROOT::VecOps::RVec<TLorentzVector> leptons) { 
        if (leptons.size() < 2) return 0.0;
        return (leptons.at(0) + leptons.at(1)).Pt(); 
    };


    // Define new columns
    node = node.Define("Electron_mass",buildParticlesZeroMass,{"Electron.PT"});
    node = node.Define("electrons_p4",buildParticlesFromContainer, {"Electron.PT", "Electron.Eta", "Electron.Phi", "Electron_mass"});
    node = node.Define("jets_p4",buildParticlesFromContainer,{"Jet.PT", "Jet.Eta", "Jet.Phi","Jet.Mass"});    
    node = node.Define("deltaPhi_electrons",deltaPhiMuons,{"electrons_p4"});
    node = node.Define("ElectronIsPtOrdered",isPtOrdered,{"Electron.PT"});
    node = node.Define("JetisPtOrdered",isPtOrdered,{"Jet.PT"});
    node = node.Define("Electron_1_pt","Electron.PT[0]");
    node = node.Define("Electron_2_pt","Electron.PT[1]");
    node = node.Define("Jet_1_pt","Jet.PT[0]");
    node = node.Define("Jet_2_pt","Jet.PT[1]");
    node = node.Define("Bjet_size",numberOfHeavyJets,{"Jet.BTag"});
    node = node.Define("deltaRapidity",deltaRapidity,{"jets_p4"});
    node = node.Define("mjj",invariantMassTwoLeading,{"jets_p4"});
    node = node.Define("mll",invariantMassTwoLeading,{"electrons_p4"});
    node = node.Define("nJetsInGap",nJetsInGap,{"jets_p4"});
    node = node.Define("ptBalance",ptBalance,{"electrons_p4","jets_p4"});
    node = node.Define("dilepCentrality",dilepCentrality,{"electrons_p4","jets_p4"});
    node = node.Define("dileptonPt",dileptonPt,{"electrons_p4"});
    
    
    // Filters
    node = node.Filter("Electron_size == 2");
    node = node.Filter("Muon_size == 0");
    node = node.Filter("Jet_size == 2 || Jet_size == 3");
    node = node.Filter("Electron.Charge[0] != Electron.Charge[1]");

    node = node.Filter("Electron.PT[0] >= 27");
    node = node.Filter("Electron.PT[1] >= 27");
    node = node.Filter("Jet.PT[0] >= 25");
    node = node.Filter("Jet.PT[1] >= 25");
    node = node.Filter("Jet_size == 3 ? Jet.PT[2] >= 25 : 1");
    node = node.Filter("mjj >= 500");
    node = node.Filter("mll > 40");

    std::vector histograms = {node.Histo1D("ElectronIsPtOrdered")};
    histograms.push_back(node.Histo1D({"lep1_pt_basic", "Electron_1_pt_basic", 500, 0, 500},"Electron_1_pt", weightColumn));
    histograms.push_back(node.Histo1D({"lep2_pt_basic", "Electron_2_pt_basic", 500, 0, 500},"Electron_2_pt", weightColumn));
    histograms.push_back(node.Histo1D({"ljet0_pt_basic", "Jet_1_pt_basic", 100, 0, 1000},"Jet_1_pt", weightColumn));
    histograms.push_back(node.Histo1D({"ljet1_pt_basic", "Jet_2_pt_basic", 100, 0, 1000},"Jet_2_pt", weightColumn));   
    histograms.push_back(node.Histo1D({"delta_phi_basic", "deltaPhi_electrons_basic", 64, 0, 3.2},"deltaPhi_electrons", weightColumn));
    histograms.push_back(node.Histo1D({"delta_y_basic", "deltaRapidity_basic", 100, 0, 10},"deltaRapidity", weightColumn));
    histograms.push_back(node.Histo1D({"mass_jj_basic", "mjj_basic", 50, 0, 5000},"mjj", weightColumn));
    histograms.push_back(node.Histo1D({"inv_mass_basic", "mll_basic", 300, 0, 300},"mll", weightColumn));
    histograms.push_back(node.Histo1D({"n_jets_interval_basic", "nJetsInGap_basic", 2, 0, 2},"nJetsInGap", weightColumn));
    histograms.push_back(node.Histo1D({"pt_bal_basic", "ptBalance_basic", 100, 0, 1},"ptBalance", weightColumn));
    histograms.push_back(node.Histo1D({"Z_centrality_basic", "dilepCentrality_basic", 100, 0, 1},"dilepCentrality", weightColumn));
    histograms.push_back(node.Histo1D({"Z_pt_reco_basic", "dileptonPt_basic", 100, 0, 1000},"dileptonPt", weightColumn));

    node = node.Filter("Bjet_size == 0");
    node = node.Filter("Electron.PT[0] >= 50");
    node = node.Filter("Electron.PT[1] >= 40");
    node = node.Filter("Jet.PT[0] >= 75");
    node = node.Filter("Jet.PT[1] >= 70");
    node = node.Filter("deltaRapidity >= 2.0");
    node = node.Filter("mjj >= 1000");
    node = node.Filter("mll > 101 && mll > 81");
    node = node.Filter("nJetsInGap == 0");
    node = node.Filter("ptBalance <= 0.15");
    node = node.Filter("dilepCentrality < 0.5");


    // More Histograms
    histograms.push_back(node.Histo1D("JetisPtOrdered", weightColumn));
    histograms.push_back(node.Histo1D("Electron_size", weightColumn));
    histograms.push_back(node.Histo1D("Jet_size", weightColumn));
    histograms.push_back(node.Histo1D("Bjet_size", weightColumn));
    histograms.push_back(node.Histo1D({"lep1_pt", "lep1_pt", 500, 0, 500},"Electron_1_pt", weightColumn));
    histograms.push_back(node.Histo1D({"lep2_pt", "lep2_pt", 500, 0, 500},"Electron_2_pt", weightColumn));
    histograms.push_back(node.Histo1D({"ljet0_pt", "ljet0_pt", 100, 0, 1000},"Jet_1_pt", weightColumn));
    histograms.push_back(node.Histo1D({"ljet1_pt", "ljet1_pt", 100, 0, 1000},"Jet_2_pt", weightColumn));   
    histograms.push_back(node.Histo1D({"delta_phi", "delta_phi", 64, 0, 3.2},"deltaPhi_electrons", weightColumn));
    histograms.push_back(node.Histo1D({"delta_y", "delta_y", 100, 0, 10},"deltaRapidity", weightColumn));
    histograms.push_back(node.Histo1D({"mass_jj", "mass_jj", 50, 0, 5000},"mjj", weightColumn));
    histograms.push_back(node.Histo1D({"inv_mass", "inv_mass", 300, 0, 300},"mll", weightColumn));
    histograms.push_back(node.Histo1D({"nJetsInGap", "nJetsInGap", 2, 0, 2},"nJetsInGap", weightColumn));
    histograms.push_back(node.Histo1D({"ptBalance", "ptBalance", 100, 0, 1},"ptBalance", weightColumn));
    histograms.push_back(node.Histo1D({"Z_centrality", "Z_centrality", 100, 0, 1},"dilepCentrality", weightColumn));
    histograms.push_back(node.Histo1D({"Z_pt_reco_basic_all", "Z_pt_reco_basic_all", 100, 0, 1000},"dileptonPt", weightColumn));


    // Save histograms
    auto output = TFile::Open(outputFileName.c_str(), "RECREATE");
    output->cd();
    //output->mkdir("histograms");
    //output->cd("histograms");
    for (auto h : histograms){
        h->Write();
    }
    output->Close();
}

float CalculateLumiWeight(const std::string& fileName){
    // Map with the luminosity weights for each process.
    std::map<std::string, float> xSecMap = {
        {"BG.root", 518.5}, // in pb ttbar SM
        {"Interference.root", 0.02899},
        {"Signal.root", 0.07545},
    };
    // Find the correct cross-section based on the file name
    float xSec = 0.0;
    for(const auto& pair : xSecMap) {
        if (fileName.find(pair.first) != std::string::npos) {
            xSec = pair.second;
            break;
        }
    }
    // Check something was found
    if (xSec == 0.0) {
        throw std::runtime_error("Cross-section not found for the given file name: " + fileName);
    }

    // Calculate sum of weights
    ROOT::RDataFrame df("mytree", fileName);
    auto sumWeights = df.Sum("eventWeight").GetValue();
    float lumiWeight = 140.1 * 1000 * xSec / sumWeights;

    std::cout << "Calculated luminosity weight: " << lumiWeight << std::endl;
    return lumiWeight;
}