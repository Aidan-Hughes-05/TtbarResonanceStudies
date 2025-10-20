#pragma once
#include "Logger.h"

#include <string>


namespace Constants {
// Error message with colors!
inline const std::string g_ERROR_MESSAGE = "\033[1;31mERROR:\033[0m";
inline const std::string g_DEBUG_MESSAGE = "\033[1;94mDEBUG:\033[0m";
inline const std::string g_INFO_MESSAGE = "\033[1;92mINFO:\033[0m";

// Logger object
extern const Logger g_LOG;
}