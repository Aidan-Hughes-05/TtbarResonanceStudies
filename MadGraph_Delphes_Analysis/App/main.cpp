#include "Analysis.h"
#include <string>
#include <iostream>

int main(int argc, char** argv){
    if (argc != 4){
        std::cout << "Usage: " << argv[0] << " <input_file.root> <output_file.root> <mode>" << std::endl;
        return 1;
    }

    std::string inputFileName = argv[1];
    std::string outputFileName = argv[2];
    std::string mode = argv[3];

    float lumiWeight = CalculateLumiWeight(inputFileName);

    if (mode=="h") PartonAnalysisHistograms(inputFileName, outputFileName, lumiWeight, true);
    // DelphesAnalysis(inputFileName, outputFileName, useLumiWeights);

    return 0;
}