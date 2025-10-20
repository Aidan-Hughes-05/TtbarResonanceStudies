#pragma once

#include <iostream>
#include <string>

// Enum for log levels.
enum LogLevel {
    ERROR = 1,
    INFO,
    DEBUG,
};


// Class for a logger object.
class Logger {
    public: 
        // Constructor
        Logger(LogLevel level) {
            m_logLevel = level;
        }

        // Destructor
        ~Logger() {}

        static const std::string getLogMessage(LogLevel level);

        LogLevel getLogLevel() const{
            return m_logLevel;
        }

        // Functions to print the log message.
        template <typename T>
        void operator () (LogLevel level, const T& message) const {
            if (level <= getLogLevel()) std::cout << getLogMessage(level) << " " << message << std::endl;
        }

        template <typename T, typename V>
        void operator () (LogLevel level, const T& message1, const V& message2) const{
            if (level <= getLogLevel()) std::cout << getLogMessage(level) << " " << message1 << message2 << std::endl;
        }

        // Specialization for string literals
        void operator () (LogLevel level, const char* message) const {
            if (level <= getLogLevel()) std::cout << getLogMessage(level) << " " << std::string(message) << std::endl;
        }

    private:
        LogLevel m_logLevel = LogLevel::ERROR;
};