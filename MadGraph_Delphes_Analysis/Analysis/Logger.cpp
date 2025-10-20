#include "Logger.h"
#include "Constants.h"

const std::string Logger::getLogMessage(LogLevel level){
    std::string logMessage;
    if (level == LogLevel::INFO) {
        logMessage = Constants::g_INFO_MESSAGE;
    } else if (level == LogLevel::DEBUG) {
        logMessage = Constants::g_DEBUG_MESSAGE;
    } else if (level == LogLevel::ERROR) {
        logMessage = Constants::g_ERROR_MESSAGE;
    }
    return logMessage;
}