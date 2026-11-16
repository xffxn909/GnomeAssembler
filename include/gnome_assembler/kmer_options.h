#pragma once

#include <cstddef>
#include <cstdint>

namespace gnome {

struct KmerOptions {
    int kmer_length = 21;
    uint32_t min_kmer_abundance = 2;
};

} // namespace gnome
