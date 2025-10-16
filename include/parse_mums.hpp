// parse_mums.hpp
// Common parsers for .mums (text) and .bumbl (binary) files
// Exposes a shared Mum structure usable from multiple compilation units

#pragma once

#include <cstdint>
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <functional>

namespace mumsio {

struct Mum {
    uint32_t length;
    std::vector<int64_t> offsets; // supports partials as -1
    std::vector<bool> strands;
};

// Common binary read helper: read exactly n bytes or throw
inline void readExact(std::istream& in, char* dst, std::streamsize n) {
    in.read(dst, n);
    if (in.gcount() != n) {
        throw std::runtime_error("Unexpected EOF while reading bumbl file");
    }
}

// Parse a .mums text file (length\toffsets\tstrands[\textra...])
// If forbidPartials is true, throws if any offset field is empty.
inline std::vector<Mum> parse_mums(const std::string& path, bool noPartials = true) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Failed to open MUMs file: " + path);
    }

    std::vector<Mum> mums;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string length_str, offsets_str, strands_str;
        iss >> length_str >> offsets_str >> strands_str; // ignore extra fields if present
        if (length_str.empty() || offsets_str.empty() || strands_str.empty()) {
            throw std::runtime_error("Malformed MUMs line: " + line);
        }

        uint32_t len = static_cast<uint32_t>(std::stoul(length_str));
        std::vector<int64_t> row_offsets;
        std::vector<bool> row_strands;

        // parse offsets
        {
            std::stringstream sso(offsets_str);
            std::string tok;
            while (std::getline(sso, tok, ',')) {
                if (tok.empty()) {
                    if (noPartials) {
                        throw std::runtime_error(
                            "Cannot parse partial MUMs: empty offset encountered. Filter to strict MUMs.");
                    }
                    row_offsets.push_back(-1);
                } else {
                    row_offsets.push_back(static_cast<int64_t>(std::stoll(tok)));
                }
            }
        }

        // parse strands
        {
            std::stringstream sss(strands_str);
            std::string tok;
            while (std::getline(sss, tok, ',')) {
                row_strands.push_back(tok == "+");
            }
        }

        if (row_offsets.size() != row_strands.size()) {
            throw std::runtime_error("Offsets and strands column size mismatch in MUMs file");
        }
        Mum m;
        m.length = len;
        m.offsets = std::move(row_offsets);
        m.strands = std::move(row_strands);
        mums.emplace_back(std::move(m));
    }

    return mums;
}

// Parse a .bumbl binary file written by Python utils or mem_finder (harmonized flags)
// Layout (little-endian):
// uint16 flags; uint64 n_seqs; uint64 n_mums; lengths[n_mums] (u32 if bit15 length32==1, else u16);
// starts[n_mums * n_seqs] (int64); strands[ceil(n_mums*n_seqs/8)] (packed bits, little bit order);
// if flags.coll_blocks: uint64 num_blocks; blocks[num_blocks][2] (uint32 pairs)
inline std::vector<Mum> parse_bumbl(const std::string& path, bool noPartials = true) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Failed to open bumbl file: " + path);
    }


    uint16_t flags = 0;
    readExact(in, reinterpret_cast<char*>(&flags), sizeof(flags));
    // bit13=partial, bit14=coll_blocks, bit15=length32
    bool flag_partial    = (flags & static_cast<uint16_t>(1u << 13)) != 0;
    bool flag_collblocks = (flags & static_cast<uint16_t>(1u << 14)) != 0;
    bool flag_length32   = (flags & static_cast<uint16_t>(1u << 15)) != 0;

    uint64_t n_seqs = 0, n_mums = 0;
    readExact(in, reinterpret_cast<char*>(&n_seqs), sizeof(n_seqs));
    readExact(in, reinterpret_cast<char*>(&n_mums), sizeof(n_mums));

    // If noPartials and header indicates partials, we can early reject; otherwise proceed and validate per-start
    if (noPartials && flag_partial) {
        throw std::runtime_error("Cannot parse partial bumbl: header indicates partial MUMs");
    }

    std::vector<uint32_t> lengths(n_mums);

    if (flag_length32) {
        // lengths are uint32
        readExact(in, reinterpret_cast<char*>(lengths.data()), static_cast<std::streamsize>(n_mums * sizeof(uint32_t)));
    } else {
        // lengths are uint16
        std::vector<uint16_t> lens16(n_mums);
        readExact(in, reinterpret_cast<char*>(lens16.data()), static_cast<std::streamsize>(n_mums * sizeof(uint16_t)));
        for (size_t i = 0; i < n_mums; ++i) {
            lengths[i] = static_cast<uint32_t>(lens16[i]);
        }
    }

    // starts: int64, row-major [n_mums][n_seqs]
    const size_t total = static_cast<size_t>(n_mums) * static_cast<size_t>(n_seqs);
    std::vector<int64_t> flat_starts(total);
    if (total > 0) {
        readExact(in, reinterpret_cast<char*>(flat_starts.data()), static_cast<std::streamsize>(total * sizeof(int64_t)));
    }

    // strands: packed bits, MSB-first within each byte (compatible with numpy's default and mem_finder)
    std::vector<uint8_t> packed;
    if (total > 0) {
        const size_t num_bytes = (total + 7) / 8;
        packed.resize(num_bytes);
        readOrThrow(reinterpret_cast<char*>(packed.data()), static_cast<std::streamsize>(num_bytes));
    }

    // optional blocks
    if (flag_collblocks) {
        uint64_t num_blocks = 0;
        readExact(in, reinterpret_cast<char*>(&num_blocks), sizeof(num_blocks));
        if (num_blocks > 0) {
            // consume blocks but ignore; merge/extract do not use them
            std::vector<uint32_t> flat_blocks(num_blocks * 2);
            readExact(in, reinterpret_cast<char*>(flat_blocks.data()), static_cast<std::streamsize>(flat_blocks.size() * sizeof(uint32_t)));
        }
    }

    std::vector<Mum> mums;
    mums.reserve(n_mums);
    for (uint64_t r = 0; r < n_mums; ++r) {
        Mum m;
        m.length = lengths[r];
        m.offsets.resize(n_seqs);
        m.strands.resize(n_seqs);
        for (uint64_t c = 0; c < n_seqs; ++c) {
            size_t idx = static_cast<size_t>(r) * static_cast<size_t>(n_seqs) + static_cast<size_t>(c);
            int64_t start = flat_starts[idx];
            if (noPartials && start == -1) {
                throw std::runtime_error("Cannot parse partial bumbl: -1 start encountered");
            }
            m.offsets[c] = start;
            if (!packed.empty()) {
                uint8_t byte = packed[idx / 8];
                uint8_t bit  = (byte >> (7 - (idx % 8))) & 0x1; // MSB-first
                m.strands[c] = (bit != 0);
            } else {
                m.strands[c] = false;
            }
        }
        mums.emplace_back(std::move(m));
    }

    return mums;
}

// Streaming APIs: emit (length, firstOffset) per MUM without materializing all data
// Text .mums streamer (first offset only)
inline void stream_mums_first(
        const std::string& path,
        bool noPartials = true,
        const std::function<void(uint32_t, int64_t)>& consumer) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Failed to open MUMs file: " + path);
    }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        // fast parse: length until first tab
        size_t t1 = line.find('\t');
        if (t1 == std::string::npos) throw std::runtime_error("Malformed MUMs line");
        uint32_t length = static_cast<uint32_t>(std::stoul(line.substr(0, t1)));
        // offsets field between first tab and second tab
        size_t t2 = line.find('\t', t1 + 1);
        if (t2 == std::string::npos) throw std::runtime_error("Malformed MUMs line: missing strands field");
        const std::string offsets_str = line.substr(t1 + 1, t2 - (t1 + 1));
        // parse first offset only
        size_t comma = offsets_str.find(',');
        int64_t first = -1;
        if (comma ==  std::string::npos) {
            // single entry (could be empty)
            if (offsets_str.empty()) {
                if (noPartials) throw std::runtime_error("Cannot parse partial MUMs: empty offset encountered");
                first = -1;
            } else {
                first = static_cast<int64_t>(std::stoll(offsets_str));
            }
        } else {
            // take substring before first comma (could be empty)
            if (comma == 0) {
                if (noPartials) throw std::runtime_error("Cannot parse partial MUMs: empty offset encountered");
                first = -1;
            } else {
                first = static_cast<int64_t>(std::stoll(offsets_str.substr(0, comma)));
            }
        }
        consumer(length, first);
    }
}

// Binary .bumbl streamer (first offset only)
inline void stream_bumbl_first(
        const std::string& path,
        bool noPartials = true,
        const std::function<void(uint32_t, int64_t)>& consumer) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Failed to open bumbl file: " + path);
    }
    uint16_t flags = 0;
    readExact(in, reinterpret_cast<char*>(&flags), sizeof(flags));
    bool flag_partial    = (flags & static_cast<uint16_t>(1u << 13)) != 0;
    bool flag_length32   = (flags & static_cast<uint16_t>(1u << 15)) != 0;
    uint64_t n_seqs = 0, n_mums = 0;
    readExact(in, reinterpret_cast<char*>(&n_seqs), sizeof(n_seqs));
    readExact(in, reinterpret_cast<char*>(&n_mums), sizeof(n_mums));
    if (noPartials && flag_partial) {
        throw std::runtime_error("Cannot parse partial bumbl: header indicates partial MUMs");
    }
    std::vector<uint32_t> lengths(n_mums);
    if (flag_length32) {
        readExact(in, reinterpret_cast<char*>(lengths.data()), static_cast<std::streamsize>(n_mums * sizeof(uint32_t)));
    } else {
        std::vector<uint16_t> lens16(n_mums);
        readExact(in, reinterpret_cast<char*>(lens16.data()), static_cast<std::streamsize>(n_mums * sizeof(uint16_t)));
        for (size_t i = 0; i < n_mums; ++i) lengths[i] = lens16[i];
    }
    // current position is starts base
    std::streampos starts_base = in.tellg();
    const std::streamoff row_stride = static_cast<std::streamoff>(n_seqs * sizeof(int64_t));
    for (uint64_t i = 0; i < n_mums; ++i) {
        std::streampos pos = starts_base + static_cast<std::streamoff>(i) * row_stride; // first offset
        in.seekg(pos);
        int64_t start = -1;
        readExact(in, reinterpret_cast<char*>(&start), sizeof(start));
        if (noPartials && start == -1) {
            throw std::runtime_error("Cannot parse partial bumbl: empty start encountered");
        }
        consumer(lengths[static_cast<size_t>(i)], start);
    }
}


} // namespace mumsio


