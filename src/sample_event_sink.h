#ifndef SRC_SAMPLE_EVENT_SINK_H_
#define SRC_SAMPLE_EVENT_SINK_H_

#include <functional>
#include <vector>

#include "data_struct.h"

// Called synchronously with one sample, its zero-based index, and the event count.
using SampleEventSink = std::function<void(std::vector<iSS_Hadron> &, int, int)>;

#endif  // SRC_SAMPLE_EVENT_SINK_H_
