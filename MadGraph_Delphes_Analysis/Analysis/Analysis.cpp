#include "Analysis.h"
#include "Logger.h"
#include "Kinematics.h"
#include "Constants.h"

#include <ROOT/RDataFrame.hxx>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TChain.h>

#include <vector>
#include <string>
#include <map>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <map>

std::map<std::string, float> getTheoryParametersFromFileName(const std::filesystem::path& filePath){
    std::map<std::string, float> theoryParameters;

    // Get file name:
    std::string fileName = filePath.filename().string();
    if (fileName == "") {
        Constants::g_LOG(LogLevel::ERROR, "Empty file name provided for theory parameter extraction.");
        throw std::runtime_error("");
    }

    theoryParameters["Coupling"] = 0.0f;
    theoryParameters["HMass"] = 0.0f;
    theoryParameters["HWidth"] = -1.0f;
    // 1) Coupling: Find by appearance of "_c" in the file name.
    std::size_t pos = fileName.find("_c");
    if (pos != std::string::npos){
        std::size_t start = pos + 2; // position after "_c"
        std::size_t end = start+3; // Always _cXpY
        std::string couplingString = fileName.substr(start,end);
        // Replace "p" by "."
        std::replace(couplingString.begin(), couplingString.end(), 'p', '.');
        float coupling = std::stof(couplingString);
        theoryParameters["Coupling"] = coupling;
        Constants::g_LOG(LogLevel::INFO, "Extracted coupling parameter: " + std::to_string(coupling) + " from file name: " + fileName);
    } else {
        Constants::g_LOG(LogLevel::INFO, "Coupling parameter not found in file name: " + fileName+ " setting to 0.0");
    }

    // 2) Mass: Find by appearance of "_m" in the file name.
    std::size_t pos_mass = fileName.find("_m");
    if (pos_mass != std::string::npos){
        std::size_t start = pos_mass + 2; // position after "_m"
        std::size_t end = start+3; // Always _mXXX
        std::string massString = fileName.substr(start,end);
        float mass = std::stof(massString);
        theoryParameters["HMass"] = mass;
        Constants::g_LOG(LogLevel::INFO, "Extracted mass parameter: " + std::to_string(mass) + " from file name: " + fileName);
    } else {
        Constants::g_LOG(LogLevel::INFO, "Mass parameter not found in file name: " + fileName+ " setting to 0.0");
    }

    // 3) Width: Find by appearance of "_w" in the file name.
    std::size_t pos_width = fileName.find("_w");
    if (pos_width != std::string::npos){
        std::string after_width_string = fileName.substr(pos_width);
        std::size_t pos_end_width = after_width_string.find("_", 1); // Find the next underscore after "_w"
        std::string widthString = fileName.substr(pos_width+2, pos_width + pos_end_width);
        float width = std::stof(widthString);
        theoryParameters["HWidth"] = width;
        Constants::g_LOG(LogLevel::INFO, "Extracted width parameter: " + std::to_string(width) + " from file name: " + fileName);
    } else {
        if (fileName.find("Signal") != std::string::npos || fileName.find("Interference") != std::string::npos){
            Constants::g_LOG(LogLevel::INFO, "Width parameter not found in file name: " + fileName+ ", but sample seems to be signal or interference. Setting to 2.0!");
            theoryParameters["HWidth"] = 2.0f;
        } else {
            Constants::g_LOG(LogLevel::INFO, "Width parameter not found in file name: " + fileName+ ", background sample setting width to -1.0");
        }
    }

    return theoryParameters;

}


void PartonHistograms(const std::vector<std::filesystem::path>& filesToProcess, const std::string& outputFileName, const MetadataObject& metadataObject, bool doWeights){

    // Build the RDF RDatasetSpec
    ROOT::RDF::Experimental::RDatasetSpec spec;

    // Convert file paths into strings vector for RDF
    std::vector<std::string> fileNames;
    int fileIndex = 1;
    for (const auto& filePath : filesToProcess) {
        // Create the RMetaData object which will be used to create the RSample object.
        ROOT::RDF::Experimental::RMetaData meta;
        auto theoryParameters = getTheoryParametersFromFileName(filePath);
        meta.Add("Coupling", theoryParameters["Coupling"]);
        meta.Add("HMass", theoryParameters["HMass"]);
        meta.Add("HWidth", theoryParameters["HWidth"]);
        meta.Add("lWeight", CalculateLumiWeight(filePath, metadataObject)); 

        // Create the RSample object "iSample" with:
        // sample name "mySampleName",
        // tree name "outputTree1",
        // file name "outputFile.root",
        // and associated metadata information.
        ROOT::RDF::Experimental::RSample iSample("sample_" + std::to_string(fileIndex), "mytree", filePath.string(), meta);
        spec.AddSample(iSample);
        fileIndex++;
    }
        
    
    // Finally, create an RDataFrame from the RDatasetSpec object.
    auto df = ROOT::RDataFrame(spec);
    ROOT::RDF::RNode node = df;
    
    // 1) Build the top container
    // Status = 1 means final state particles
    auto myFilter = [](std::vector<int> pdgid, std::vector<int> status) { 
        std::vector<int> indices{};
        for (std::size_t i{} ; i<pdgid.size() ; i++){ 
            if (abs(pdgid.at(i)) == 6 && status.at(i) == 1) indices.push_back(i); 
        } 
        return indices; 
    };

    // Define the dummy column for the weight or apply lumi scaling
    std::string weightColumn = "";
    if (!doWeights){
        node = node.Define("dummyWeight",[]() { return 1.0; },{"pdgID"});
        weightColumn = "dummyWeight";
    } else { // Return +1 or -1 weight
        node = node.DefinePerSample("lWeight", []([[maybe_unused]] unsigned int slot, const ROOT::RDF::RSampleInfo &id) { return static_cast<float>(id.GetD("lWeight")); });
        node = node.Define("lumiWeight",[](float lWeight, float eventWeight) { return lWeight*eventWeight; },{"lWeight","eventWeight"});
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

    // 3) Define tops pt, eta, phi
    node = node.Define("Top_1_pt","Tops.at(0).Pt()");
    node = node.Define("Top_2_pt","Tops.at(1).Pt()");
    node = node.Define("Top_1_eta","Tops.at(0).Eta()");
    node = node.Define("Top_2_eta","Tops.at(1).Eta()");
    node = node.Define("Top_1_phi","Tops.at(0).Phi()");
    node = node.Define("Top_2_phi","Tops.at(1).Phi()");

    // Some more kinematics
    // 
    auto deltaRapidity = [](std::vector<TLorentzVector> particles) { return abs(particles.at(0).Rapidity() - particles.at(1).Rapidity()); };
    //
    auto invariantMassTwoLeading = [](std::vector<TLorentzVector> particles) { return (sqrt( 2.0*particles.at(0).Dot(particles.at(1))) ); };
    //
    auto dileptonPt = [](std::vector<TLorentzVector> leptons) { 
        if (leptons.size() < 2) return 0.0;
        return (leptons.at(0) + leptons.at(1)).Pt(); 
    };
    //
    auto deltaPhi = [](std::vector<TLorentzVector> particles) { return Kinematics::del_phi(particles.at(0).Phi(),particles.at(1).Phi()); };

    // Create histogram container
    std::vector histograms = {node.Histo1D("top_indices", weightColumn)};
    histograms.push_back(node.Histo1D("top_indices_size", weightColumn));
    histograms.push_back(node.Histo1D({"Top_1_pt_basic", "Top_1_pt_basic", 500, 0, 500},"Top_1_pt", weightColumn));
    histograms.push_back(node.Histo1D({"Top_2_pt_basic", "Top_2_pt_basic", 500, 0, 500},"Top_2_pt", weightColumn));
    histograms.push_back(node.Histo1D({"Top_1_eta_basic", "Top_1_eta_basic", 100, -5, 5},"Top_1_eta", weightColumn));
    histograms.push_back(node.Histo1D({"Top_2_eta_basic", "Top_2_eta_basic", 100, -5, 5},"Top_2_eta", weightColumn));
    histograms.push_back(node.Histo1D({"Top_1_phi_basic", "Top_1_phi_basic", 64, -3.2, 3.2},"Top_1_phi", weightColumn));
    histograms.push_back(node.Histo1D({"Top_2_phi_basic", "Top_2_phi_basic", 64, -3.2, 3.2},"Top_2_phi", weightColumn));

    node = node.Define("deltaPhi_Tops",deltaPhi,{"Tops"});
    node = node.Define("deltaRapidity",deltaRapidity,{"Tops"});
    node = node.Define ("mttbar",invariantMassTwoLeading,{"Tops"});
    node = node.Define("diTopPt",dileptonPt,{"Tops"});

    histograms.push_back(node.Histo1D({"deltaPhi_Tops_basic", "deltaPhi_Tops_basic", 64, 0, 3.2},"deltaPhi_Tops", weightColumn));
    histograms.push_back(node.Histo1D({"deltaRapidity_basic", "deltaRapidity_basic", 100, 0, 10},"deltaRapidity", weightColumn));
    histograms.push_back(node.Histo1D({"mttbar_basic", "mttbar_basic", 1000, 0, 1000},"mttbar", weightColumn));
    histograms.push_back(node.Histo1D({"diTopPt_basic", "diTopPt_basic", 1000,0,1000},"diTopPt", weightColumn));
    histograms.push_back(node.Histo1D({"weight_basic", "weight_basic", 100, -5, 5},weightColumn));

    node = node.Filter("mttbar > 500");
    histograms.push_back(node.Histo1D({"weight_mttbar", "weight_mttbar", 100, -5, 5},weightColumn));
    histograms.push_back(node.Histo1D({"mttbar_final", "mttbar_final", 1200, 0, 1200},"mttbar", weightColumn));
    
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

void PartonNtuples(const std::vector<std::filesystem::path>& filesToProcess, const std::string& outputFileName, const MetadataObject& metadataObject, bool doWeights){
    // Build the RDF RDatasetSpec
    ROOT::RDF::Experimental::RDatasetSpec spec;

    // Convert file paths into strings vector for RDF
    std::vector<std::string> fileNames;
    int fileIndex = 1;
    for (const auto& filePath : filesToProcess) {
        // Create the RMetaData object which will be used to create the RSample object.
        ROOT::RDF::Experimental::RMetaData meta;
        auto theoryParameters = getTheoryParametersFromFileName(filePath);
        meta.Add("Coupling", theoryParameters["Coupling"]);
        meta.Add("HMass", theoryParameters["HMass"]);
        meta.Add("HWidth", theoryParameters["HWidth"]);
        meta.Add("lWeight", CalculateLumiWeight(filePath, metadataObject)); 

        // Create the RSample object "iSample" with:
        // sample name "mySampleName",
        // tree name "outputTree1",
        // file name "outputFile.root",
        // and associated metadata information.
        ROOT::RDF::Experimental::RSample iSample("sample_" + std::to_string(fileIndex), "mytree", filePath.string(), meta);
        spec.AddSample(iSample);
        fileIndex++;
    }
        
    
    // Finally, create an RDataFrame from the RDatasetSpec object.
    auto df = ROOT::RDataFrame(spec);
    ROOT::RDF::RNode node = df;
    // 1) Build the top container
    // Status = 1 means final state particles
    auto myFilter = [](std::vector<int> pdgid, std::vector<int> status) { 
        std::vector<int> indices{};
        for (std::size_t i{} ; i<pdgid.size() ; i++){ 
            if (abs(pdgid.at(i)) == 6 && status.at(i) == 1) indices.push_back(i); 
        } 
        return indices; 
    };

    // Define the dummy column for the weight or apply lumi scaling
    std::string weightColumn = "";
    if (!doWeights){
        node = node.Define("dummyWeight",[]() { return 1.0; },{"pdgID"});
        weightColumn = "dummyWeight";
    } else { // Return +1 or -1 weight
        node = node.DefinePerSample("lWeight", []([[maybe_unused]] unsigned int slot, const ROOT::RDF::RSampleInfo &id) { return static_cast<float>(id.GetD("lWeight")); });
        node = node.Define("lumiWeight",[](float lWeight, float eventWeight) { return lWeight*eventWeight; },{"lWeight","eventWeight"});
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

    // 3) Define tops pt, eta, phi
    node = node.Define("Top_1_pt","Tops.at(0).Pt()");
    node = node.Define("Top_2_pt","Tops.at(1).Pt()");
    node = node.Define("Top_1_eta","Tops.at(0).Eta()");
    node = node.Define("Top_2_eta","Tops.at(1).Eta()");
    node = node.Define("Top_1_phi","Tops.at(0).Phi()");
    node = node.Define("Top_2_phi","Tops.at(1).Phi()");

    // Some more kinematics
    // 
    // auto deltaRapidity = [](std::vector<TLorentzVector> particles) { return abs(particles.at(0).Rapidity() - particles.at(1).Rapidity()); };
    //
    auto invariantMassTwoLeading = [](std::vector<TLorentzVector> particles) { return (sqrt( 2.0*particles.at(0).Dot(particles.at(1))) ); };
    //
    // auto dileptonPt = [](std::vector<TLorentzVector> leptons) { 
    //    if (leptons.size() < 2) return 0.0;
    //    return (leptons.at(0) + leptons.at(1)).Pt(); 
    //};
    //
    // auto deltaPhi = [](std::vector<TLorentzVector> particles) { return Kinematics::del_phi(particles.at(0).Phi(),particles.at(1).Phi()); };

    node = node.Define ("mttbar",invariantMassTwoLeading,{"Tops"});

    node = node.Filter("mttbar > 500");

    std::vector<std::string> columnsToSave = {"Top_1_pt","Top_2_pt","Top_1_eta","Top_2_eta","Top_1_phi","Top_2_phi",weightColumn,"mttbar"};

    // Add Coupling per sample
    node = node.DefinePerSample("Coupling",
                   []([[maybe_unused]] unsigned int slot, const ROOT::RDF::RSampleInfo &id)
                   { return id.GetD("Coupling"); });
    columnsToSave.push_back("Coupling");

    // Add HMass per sample
    node = node.DefinePerSample("HMass",
                   []([[maybe_unused]] unsigned int slot, const ROOT::RDF::RSampleInfo &id)
                   { return id.GetD("HMass"); });
    columnsToSave.push_back("HMass");

    // Add HWidth per sample
    node = node.DefinePerSample("HWidth",
                   []([[maybe_unused]] unsigned int slot, const ROOT::RDF::RSampleInfo &id)
                   { return id.GetD("HWidth"); });
    columnsToSave.push_back("HWidth");

    node.Snapshot("Ttbar",outputFileName,columnsToSave);
}

void DelphesHistograms(const std::string& fileName, const std::string& outputFileName, bool doWeights){
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
        node = node.Define(weightColumn,[]() { return 1.0; },{"Event.Weight"});
    } else {
        node = node.Define("SampleWeight",[&nEvents, &luminosity](const ROOT::VecOps::RVec<float>& eventWeight) { return eventWeight[10] * 1000 * luminosity / (nEvents * eventWeight[10]); },{"Event.Weight"});
        node = node.Define(weightColumn,[&nEvents, &luminosity](const ROOT::VecOps::RVec<float>& eventWeight, float sampleWeight) { return eventWeight[0] * sampleWeight; },{"Event.Weight","SampleWeight"});
    }


    // Define functions
    // 1)
    // auto getSize = [](ROOT::VecOps::RVec<Float_t> vec) { return vec.size(); };
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

std::map<std::string, std::tuple<float, long, float>> FormMetadataObject(std::string pathToMetadataFile){
    // Create the metadata map
    std::map<std::string, std::tuple<float, long, float>> metadataMap;

    // Now open the metadata map file and parse it
    std::filesystem::path metadataFilePath = pathToMetadataFile;
    std::ifstream xSecFile(metadataFilePath);
    if (!xSecFile.is_open()){
        Constants::g_LOG(LogLevel::ERROR, "Could not open cross-section file: " + metadataFilePath.string());
        throw std::runtime_error("");
    }

    // Name, cross-section (pb), number of events, sum of weights
    std::string line;
    while (std::getline(xSecFile, line)){
        // Split by comma
        std::istringstream ss(line);
        std::string token;
        std::vector<std::string> tokens;
        while (std::getline(ss, token, ',')){
            tokens.push_back(token);
        }
        if (tokens.size() != 4){
            Constants::g_LOG(LogLevel::ERROR, "Malformed line in cross-section file: " + line);
            throw std::runtime_error("");
        }
        std::string sampleName = tokens[0];
        float crossSection = std::stof(tokens[1]);
        long nEvents = std::stol(tokens[2]);
        float sumWeights = std::stof(tokens[3]);
        auto metadataTuple = std::make_tuple(crossSection, nEvents, sumWeights);
        metadataMap[sampleName] = metadataTuple;

    }

    xSecFile.close();
    return metadataMap;
}

float CalculateLumiWeight(const std::filesystem::path& fileToProcess, const MetadataObject& metadataObject){

    // Get the file name without the full path
    std::string baseFileName = fileToProcess.filename().string();

    // Remove the last part "_X.root" by first splitting by '_' and then removing the last element
    std::size_t lastUnderscorePos = baseFileName.find_last_of('_');
    if (lastUnderscorePos == std::string::npos){
        Constants::g_LOG(LogLevel::ERROR, "File name does not contain expected '_' character: " + baseFileName);
        throw std::runtime_error("");
    }
    std::string fileNameWithoutLastPart = baseFileName.substr(0, lastUnderscorePos);


    // Find the entry in the metadata object
    const auto& metadata = metadataObject.at(fileNameWithoutLastPart);
    float xSec = std::get<0>(metadata);
    // long nEvents = std::get<1>(metadata);
    float sumWeights = std::get<2>(metadata);

    float lumiWeight = 140.1 * 1000 * xSec / sumWeights;

    return lumiWeight;
}