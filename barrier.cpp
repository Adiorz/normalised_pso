#include "barrier.hpp"

Barrier::Barrier(std::size_t iCount) :
    mThreshold(iCount),
    mCount(iCount),
    mGeneration(0) {
}

void Barrier::wait() {
    std::unique_lock<std::mutex> lLock{mMutex};
    auto lGen = mGeneration;
    if (!--mCount) {
        mGeneration++;
        mCount = mThreshold;
        mCond.notify_all();
    } else {
        mCond.wait(lLock, [this, lGen] { return lGen != mGeneration; });
    }
}
