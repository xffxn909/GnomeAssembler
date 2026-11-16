#pragma once

#include <cstddef>
#include <string>

namespace gnome {

struct CleaningOptions {
    std::string reads_path;
    std::string output_prefix;
    int minimum_quality = 20;
    std::size_t minimum_read_length = 50;
    std::string adapter_prefix;
    std::string adapter_suffix;
};

} // namespace gnome
