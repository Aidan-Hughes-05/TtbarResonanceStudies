#include "Analysis.h"

#include <string>
#include <iostream>

#include <TROOT.h>

int main(int argc, char** argv){
    if (argc != 4){
        std::cout << "Usage: " << argv[0] << " <input_file.root> <output_file.root> <mode>" << std::endl;
        return 1;
    }

    std::string inputFileName = argv[1];
    std::string outputFileName = argv[2];
    std::string mode = argv[3];

    float lumiWeight = CalculateLumiWeight(inputFileName);


    ROOT::EnableImplicitMT(8); // Enable multi-threading with 8 threads

    if (mode=="h") PartonHistograms(inputFileName, outputFileName, lumiWeight, true);
    if (mode=="n") PartonNtuples(inputFileName, outputFileName, lumiWeight, true);

    return 0;
}