#pragma once
#include "string"
#include "thread"
#include "chrono"

namespace progress{
    typedef std::chrono::time_point<std::chrono::high_resolution_clock> clock;
    std::string getColorCode(int percent);

    class progBar{
        const int barLen = 30;
        std::string prevStr;
        double state;
        double min;
        double max;
        clock startTime;
        bool cont;
        std::thread mainThread;
        void mainLoop();
    public:
        progBar(double min, double max);
        progBar(double num);
        void add(double extra);
        std::thread start();
        void finish();
        void reset(int num);
        void update(double status);
    };
}