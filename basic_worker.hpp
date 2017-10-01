#pragma once

#include <stddef.h>

class BasicWorker {
protected:
    size_t id;
public:
    BasicWorker(size_t id) { this->id = id; }
    size_t get_id() const { return id; }
};
