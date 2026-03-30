#include <filesystem>
#include <string>
#include <vector>

// These helpers are referenced by `RefBuilder` constructors (CLI-style paths).
// The mumemto shared library includes `RefBuilder` for the in-memory API as well,
// so we provide minimal definitions here to avoid unresolved symbols when linking
// user code against `libmumemto.so`.

bool endsWith(const std::string& str, const std::string& suffix) {
    return str.size() >= suffix.size() &&
           0 == str.compare(str.size() - suffix.size(), suffix.size(), suffix);
}

int is_file(std::string path) {
    return std::filesystem::is_regular_file(path) ? 1 : 0;
}

std::vector<std::string> split(std::string input, char delim) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : input) {
        if (c == delim) {
            out.push_back(std::move(cur));
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    out.push_back(std::move(cur));
    return out;
}

