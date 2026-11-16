#pragma once

#include <cstddef>
#include <string>

namespace gnome {

struct FastqRecord {
    std::string name;
    std::string sequence;
    std::string quality;
    std::size_t source_record = 0;
};

} // namespace gnome
