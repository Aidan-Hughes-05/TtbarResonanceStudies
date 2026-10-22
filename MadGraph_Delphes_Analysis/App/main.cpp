#include "Analysis.h"
#include "Constants.h"

#include <string>
#include <iostream>
#include <filesystem>
#include <vector>

#include <TROOT.h>
#include <TChain.h>


std::vector<std::filesystem::path> getFilteredFiles(const std::string& inputPath, const std::string& filter){
    std::vector<std::filesystem::path> filteredFiles;
    namespace fs = std::filesystem;

    if (fs::is_regular_file(inputPath)){ // One file
        if(filter!="") {
            Constants::g_LOG(LogLevel::ERROR, "Filter cannot be used when a single file is provided as input.");
            throw std::runtime_error("");
        }
        if(fs::path(inputPath).extension() != ".root"){
            Constants::g_LOG(LogLevel::ERROR, "Input file is not a .root file.");
            throw std::runtime_error("");
        }
        filteredFiles.push_back(fs::path(inputPath));
    } else if (fs::is_directory(inputPath)){ // Directory
        Constants::g_LOG(LogLevel::INFO, "Input is a directory.");
        for (const auto& entry : fs::directory_iterator(inputPath)){
            // Only one filter passed.
            if (filter.find(',') == std::string::npos){
                if (entry.path().extension() == ".root" && entry.path().string().find(filter) != std::string::npos){
                    Constants::g_LOG(LogLevel::INFO, "Adding file to processing list: " + entry.path().string());
                    filteredFiles.push_back(entry.path());
                }
            } else { // Multiple filters passed.
                // Split by comma
                std::vector<std::string> filters;
                std::string temp;
                for (char c : filter){
                    if (c == ','){
                        filters.push_back(temp);
                        temp = "";
                    } else {
                        temp += c;
                    }
                }
                // Add last filter
                if (temp != "") filters.push_back(temp);

                // Check each filter
                bool passedAllFilters = true;
                for (const auto& f : filters){
                    if (entry.path().string().find(f) == std::string::npos){
                        passedAllFilters = false;
                        break;
                    }
                }
                if (entry.path().extension() == ".root" && passedAllFilters){
                    Constants::g_LOG(LogLevel::INFO, "Adding file to processing list: " + entry.path().string());
                    filteredFiles.push_back(entry.path());
                }
            }
        }
    } else {
        Constants::g_LOG(LogLevel::ERROR, "Input path is neither a file nor a directory.");
        throw std::runtime_error("");
    }

    if (filteredFiles.empty()){
        Constants::g_LOG(LogLevel::ERROR, "No files found matching the criteria.");
        throw std::runtime_error("");
    }

    return filteredFiles;
}

int main(int argc, char** argv){
    if (argc < 4){
        std::cout << "Usage: " << argv[0] << " <input_file.root OR path to files> <output_file.root> <mode> <comma_separated_file_filter>" << std::endl;
        return 1;
    }

    std::string inputPath = argv[1];
    std::string outputFileName = argv[2];
    std::string mode = argv[3];
    std::string fileFilter = (argc > 4) ? argv[4] : "";

    // Get the list of files to process
    std::vector<std::filesystem::path> filesToProcess = getFilteredFiles(inputPath, fileFilter);

    auto metadataObject = FormMetadataObject("/home/hughesa/dip_hunting/TtbarResonanceStudies/MadGraph_Delphes_Analysis/Metadata/metadata.txt");

    ROOT::EnableImplicitMT(8); // Enable multi-threading with 8 threads

    if (mode=="h") PartonHistograms(filesToProcess, outputFileName, metadataObject, true);
    if (mode=="n") PartonNtuples(filesToProcess, outputFileName, metadataObject, true);

    return 0;
}