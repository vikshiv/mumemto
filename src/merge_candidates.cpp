#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <tuple>
#include <filesystem>

using namespace std;

struct Mum {
    int length;
    vector<uint32_t> offsets;
    vector<bool> strands;
};

std::vector<uint16_t> readThresholds(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary | std::ios::ate);
    if (!file) {
        throw std::runtime_error("Failed to open file");
    }

    std::streamsize size = file.tellg();
    if (size % sizeof(uint16_t) != 0) {
        throw std::runtime_error("File size is not a multiple of uint16_t");
    }

    file.seekg(0, std::ios::beg);

    std::vector<uint16_t> buffer(size / sizeof(uint16_t));
    file.read(reinterpret_cast<char*>(buffer.data()), size);

    return buffer;
}

tuple<vector<Mum>, vector<bool>, vector<uint16_t>> parse_candidate(const string& path) {

    vector<uint16_t> next_best = readThresholds(path + ".thresh");
    vector<bool> mum_bv(next_best.size(), false);
    vector<Mum> mums;

    ifstream mums_file(path + ".mums");
    string line;
    while (getline(mums_file, line)) {
        istringstream iss(line);
        int length;
        string offsets_str, strands_str;
        iss >> length >> offsets_str >> strands_str;

        vector<uint32_t> offsets;
        stringstream offsets_stream(offsets_str);
        string offset;
        while (getline(offsets_stream, offset, ',')) {
            offsets.push_back(stoul(offset));
        }

        vector<bool> strands;
        stringstream strands_stream(strands_str);
        string strand;
        while (getline(strands_stream, strand, ',')) {
            strands.push_back(strand == "+");
        }

        mums.emplace_back(Mum{length, move(offsets), move(strands)});
        mum_bv[mums.back().offsets[0]] = true;
    }

    // Sort the mums by the first value in offsets
    std::sort(mums.begin(), mums.end(), [](const Mum& a, const Mum& b) {
        return a.offsets[0] < b.offsets[0];
    });

    return {move(mums), move(mum_bv), move(next_best)};
}

vector<uint32_t> fix_neg_strand(const Mum& cand, size_t cand_offset, int new_len, int old_len) {
    int len_diff = old_len - new_len;
    vector<uint32_t> offsets = cand.offsets;
    for (size_t i = 0; i < cand.offsets.size(); ++i) {
        offsets[i] += (!cand.strands[i]) ? len_diff : cand_offset;
    }
    return offsets;
}

tuple<vector<Mum>, vector<bool>, vector<uint16_t>> merge_partitions(
        const tuple<vector<Mum>, vector<bool>, vector<uint16_t>>& p1,
        const tuple<vector<Mum>, vector<bool>, vector<uint16_t>>& p2) {
    

    const auto& [p1_mums, p1_mum_bv, p1_nb] = p1;
    const auto& [p2_mums, p2_mum_bv, p2_nb] = p2;

    size_t mum_idx1 = 0, mum_idx2 = 0;
    const Mum* cur_mum1 = nullptr;
    const Mum* cur_mum2 = nullptr;
    size_t last_mum1 = 0, last_mum2 = 0;
    vector<Mum> new_mums;
    vector<bool> new_mum_bv(p1_nb.size(), false);
    vector<uint16_t> new_nb(p1_nb.size(), 0);

    for (size_t i = 0; i < p1_mum_bv.size(); ++i) {
        if (p1_nb[i] > 0 && p2_nb[i] > 0){
            new_nb[i] = max(p1_nb[i], p2_nb[i]);
        }
        if (p1_mum_bv[i]) {
            cur_mum1 = &p1_mums[mum_idx1++];
            last_mum1 = i;
        }
        if (p2_mum_bv[i]) {
            cur_mum2 = &p2_mums[mum_idx2++];
            last_mum2 = i;
        }
        if (cur_mum1 && cur_mum2 && (p1_mum_bv[i] || p2_mum_bv[i]) && (p1_nb[i] > 0 && p2_nb[i] > 0)) {
            // std::cout << "merging " << i << " " << last_mum1 << " " << last_mum2 << std::endl;
            // std::cout << "cur_mum1: " << cur_mum1->length << " " << cur_mum1->offsets[0] << " " << cur_mum1->strands[0] << std::endl;
            // std::cout << "cur_mum2: " << cur_mum2->length << " " << cur_mum2->offsets[0] << " " << cur_mum2->strands[0] << std::endl;   
            // std::cout << "p1_nb: " << p1_nb[i] << " " << p2_nb[i] << std::endl;
            // std::cout << "new_nb: " << new_nb[i] << std::endl;
            int s1_len = cur_mum1->length - (i - last_mum1);
            int s2_len = cur_mum2->length - (i - last_mum2);
            int new_len = min(s1_len, s2_len);
            // std::cout << "new_len: " << new_len << std::endl;
            if (new_len > new_nb[i] && new_len >= 20) {
                vector<uint32_t> new_offsets2 = fix_neg_strand(*cur_mum2, i - last_mum2, new_len, s2_len);
                vector<uint32_t> new_offsets1 = fix_neg_strand(*cur_mum1, i - last_mum1, new_len, s1_len);
                vector<uint32_t> combined_offsets = new_offsets1;
                combined_offsets.insert(combined_offsets.end(), new_offsets2.begin() + 1, new_offsets2.end());

                vector<bool> combined_strands = cur_mum1->strands;
                combined_strands.insert(combined_strands.end(), cur_mum2->strands.begin() + 1, cur_mum2->strands.end());

                new_mums.emplace_back(Mum{new_len, move(combined_offsets), move(combined_strands)});
                new_mum_bv[new_offsets1[0]] = true;
            }

            // std::cout << "new_mums: " << new_mums.back().length << " " << new_mums.back().offsets[0] << " " << new_mums.back().strands[0] << std::endl;
            // std::exit(0);
        }
    }

    return {move(new_mums), move(new_mum_bv), move(new_nb)};
}

string get_path(const string& path) {
    string base_path = path;
    if (path.size() >= 7 && path.compare(path.size() - 7, 7, ".thresh") == 0) {
        base_path = path.substr(0, path.size() - 7);
    } else if (path.size() >= 5 && path.compare(path.size() - 5, 5, ".mums") == 0) {
        base_path = path.substr(0, path.size() - 5);
    }
    // Check that both .thresh and .mums files exist
    string thresh_path = base_path + ".thresh";
    string mums_path = base_path + ".mums";
    
    if (!filesystem::exists(thresh_path)) {
        throw runtime_error("Could not find threshold file: " + thresh_path);
    }
    if (!filesystem::exists(mums_path)) {
        throw runtime_error("Could not find MUMs file: " + mums_path);
    }

    return base_path;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <input_paths>... -o <output_prefix> [-v]" << endl;
        return 1;
    }

    vector<string> paths;
    string output_prefix = "merged";
    bool verbose = false;

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-o" && i + 1 < argc) {
            output_prefix = argv[++i];
        } else if (arg == "-v") {
            verbose = true;
        } else {
            paths.push_back(arg);
        }
    }

    if (verbose) {
        cout << "Processing " << paths.size() << " input files" << endl;
    }

    assert(paths.size() >= 2 && "requires at least two input files");

    tuple<vector<Mum>, vector<bool>, vector<uint16_t>> left_mums;
    tuple<vector<Mum>, vector<bool>, vector<uint16_t>> right_mums;

    std::string base_left = get_path(paths[0]);
    std::string base_right = get_path(paths[1]);
    std::cerr << "merging " << filesystem::path(base_left).filename().string() << " and " << filesystem::path(base_right).filename().string() << endl;
        
    left_mums = parse_candidate(base_left);
    right_mums = parse_candidate(base_right);

    left_mums = merge_partitions(left_mums, right_mums);

    if (paths.size() > 2) {
        for (size_t i = 2; i < paths.size(); ++i) {
            base_right = get_path(paths[i]);
            std::cerr << "merging in " << filesystem::path(base_right).filename().string() << endl;
            right_mums = parse_candidate(base_right);
            left_mums = merge_partitions(left_mums, right_mums);
        }
    }

    string output_path = output_prefix + ".mums";
    if (verbose) {
        cout << "Writing results to " << output_path << endl;
    }

    ofstream output_file(output_path);
    for (const auto& mum : std::get<0>(left_mums)) {
        output_file << mum.length << "\t";
        for (size_t i = 0; i < mum.offsets.size(); ++i) {
            output_file << mum.offsets[i];
            if (i < mum.offsets.size() - 1) output_file << ",";
        }
        output_file << "\t";
        for (size_t i = 0; i < mum.strands.size(); ++i) {
            output_file << mum.strands[i];
            if (i < mum.strands.size() - 1) output_file << ",";
        }
        output_file << "\n";
    }

    return 0;
}
