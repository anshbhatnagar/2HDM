#include <vector>
#include "progress.hpp"
#include "iostream"
#include "string"
#include "math.h"
#include "thread"
#include "chrono"

progress::progBar::progBar(double min, double max) : min(min), max(max) {
    state = min;
    startTime = std::chrono::high_resolution_clock::now();
    std::string toPrint;
    std::string bar(barLen, '_');
    toPrint = "|" + bar + "| " + std::to_string(100*(state-min)/(max-min)) + "%";
    prevStr = toPrint;
    std::cout << "Progress:" << std::endl;
    std::cout << toPrint << std::flush;
}

void progress::progBar::update(double status) {
    state = status;
    clock now = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime).count();
    double frac = (state-min)/(max-min);
    auto ETA = (double)duration/frac - (double)duration;
    int filled = std::floor(barLen*frac);
    std::string toPrint;
    if (filled > 0 && barLen-filled > 0) {
        std::string bar1(filled, '#');
        std::string bar2(barLen - filled, '_');

        toPrint = "|" + bar1 + bar2 + "| " + std::to_string(100 * frac) + "%; Elapsed Time:" + std::to_string(1e-3*(double)duration/60) + "min; ETA: " + std::to_string(1e-3*ETA/60) + "min";
        std::string clear((int) prevStr.length(), '\b');
        std::cout << clear;
        std::cout << toPrint << std::flush;
        prevStr = toPrint;
    }
}

void progress::progBar::add(double extra) {
    state += extra;
}

progress::progBar::progBar(double num) {
    min = 0;
    state = min;
    max = num;
//    mainThread = std::thread;
}

void progress::progBar::mainLoop() {
    std::string toPrint;
    std::string bar(barLen, '_');
    // toPrint = "|" + bar + "| " + std::to_string(floor(100 * (state - min) / (max - min))) + "%";
    // prevStr = toPrint;
    std::cout << "Progress:" << std::endl;
    // std::cout << toPrint << std::flush;

    cont = true;
    while (cont) {
        try {
            clock now = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime).count();
            double frac = (state-min)/(max-min);
            auto ETA = (double)duration/frac - (double)duration;

            int filled = std::ceil(barLen*frac);
            std::string toPrint;
            std::string bar1(filled, '=');
            std::string bar2(barLen - filled, '_');

            double durationSeconds = 1e-3*(double)duration;
            std::vector<int> hms(3);
            std::vector<std::string> shms(3);
            int durSec = int(durationSeconds);
            hms[0] = durSec / 3600;
            hms[1] = durSec/60 - hms[0]*60;
            hms[2] = durSec - hms[0]*60*60 - hms[1]*60;
            for (int i=0; i<3; i++){
                if (hms[i] < 1) {
                    hms[i] = 0;
                }
                if (hms[i] < 10) {
                    shms[i] = "0" + std::to_string(hms[i]);
                } else {
                    shms[i] = std::to_string(hms[i]);
                }

            }
            std::string timeTaken = shms[0] + ":" + shms[1] + ":" + shms[2];

            durationSeconds = 1e-3*(double)ETA;
            durSec = int(durationSeconds);
            hms[0] = durSec / 3600;
            hms[1] = durSec/60 - hms[0]*60;
            hms[2] = durSec - hms[0]*60*60 - hms[1]*60;
            for (int i=0; i<3; i++){
                if (hms[i] < 1) {
                    hms[i] = 0;
                }
                if (hms[i] < 10) {
                    shms[i] = "0" + std::to_string(hms[i]);
                } else {
                    shms[i] = std::to_string(hms[i]);
                }
            }
            std::string timeETA = shms[0] + ":" + shms[1] + ":" + shms[2];

            int percent = ceil(100 * frac);
            std::string colourCode = getColorCode(percent);

            toPrint = "|\x1B[32m" + bar1 + ">\x1B[0m" + bar2 + "| " + colourCode + std::to_string(percent) + "%\x1B[0m Time: " + timeTaken + " ETA: " + timeETA;
            std::string clear((int) prevStr.length(), '\b');
            std::cout << clear;
            std::cout << toPrint << std::flush;
            prevStr = toPrint;
        } catch (...) {
            printf("\nExceeded length of progress bar! \n");
            cont = false;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds (100));
    }

}

std::thread progress::progBar::start() {
    startTime = std::chrono::high_resolution_clock::now();
    return std::thread(&progBar::mainLoop, this);
}

void progress::progBar::finish() {
    cont = false;
}

void progress::progBar::reset(int num) {
    min = 0;
    state = min;
    max = num;
}

std::string progress::getColorCode(int percent)
{
    std::string colourCode = "\x1B[0m";
    if (percent < 33) {
        colourCode = "\x1B[31m";
    } else if (percent < 66) { 
        colourCode = "\x1B[33m";
    } else {
        colourCode = "\x1B[32m";
    }
    return colourCode;
}
