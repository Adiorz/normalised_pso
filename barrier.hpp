#pragma once

#include <mutex>
#include <condition_variable>

class Barrier {
public:
    explicit Barrier(std::size_t iCount);

    void wait();

private:
    std::mutex mMutex;
    std::condition_variable mCond;
    std::size_t mThreshold;
    std::size_t mCount;
    std::size_t mGeneration;
};
