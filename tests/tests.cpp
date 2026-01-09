#include "../src/csr.cpp"
#include "../src/random_keys.cpp"
#include "../src/heuristics.cpp"
#include "../src/peripheral_vertices.cpp"
#include <filesystem>


std::vector<std::string> getFilesInDirectory(const std::string& directory_path) {
    std::vector<std::string> files;
    namespace fs = std::filesystem;
    
    try {
        for (const auto& entry : fs::recursive_directory_iterator(directory_path)) {
            if (entry.is_regular_file() && entry.path().extension() == ".mtx") {
                files.push_back(entry.path().string());
            }
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
    }
    
    return files;
}

bool failed = false;
usize n = 0;
// If parameter is not true, test fails
#define IS_TRUE(x) { if (!(x)) { std::cout << "\n - " << __FUNCTION__ << ":\n\t * Failed on line :" << __LINE__ << std::endl; failed = true;} else {n++;}}

void testCSRFromFile2() {
    // %   V         = [ 3 1 3 5 4 1 5 4 ]
    // %   COL_INDEX = [ 1 3 0 3 3 0 1 2 ]
    // %   ROW_INDEX = [ 0 2 4 5 8  ] 
    std::string file = "input/test1.mtx";
    CSR csr(file);
    IS_TRUE(csr.n_nz == 8);
    IS_TRUE(csr.m == 4);
    IS_TRUE(csr.row_index.size() == (csr.m + 1));
    IS_TRUE(csr.col_index.size() == csr.n_nz);
    // Check elements of csr.values
    // std::vector<real> expected_v = {3, 1, 3, 5, 4, 1, 5, 4};
    // IS_TRUE(csr.values == expected_v);
    // Check elements of csr.col_index
    std::vector<usize> expected_col_index = {1, 3, 0, 3, 3, 0, 1, 2};
    IS_TRUE(csr.col_index == expected_col_index);
    // Check elements of csr.row_index
    std::vector<usize> expected_row_index = {0, 2, 4, 5, 8};
    IS_TRUE(csr.row_index == expected_row_index);
}

void testCSRFromFile3() {
//  V         = [  ]
//  COL_INDEX = [ 2 5 4 5 0 3 4 2 1 2 5 0 1 4 ]
//  ROW_INDEX = [ 0 2 4 7 8 11 14 ] 
    std::string file = "input/test2.mtx";
    CSR csr(file);
    IS_TRUE(csr.n_nz == 14);
    IS_TRUE(csr.m == 6);
    IS_TRUE(csr.row_index.size() == (csr.m + 1));
    IS_TRUE(csr.col_index.size() == csr.n_nz);
    // Check elements of csr.col_index
    std::vector<usize> expected_col_index = {2, 5, 4, 5, 0, 3, 4, 2, 1, 2, 5, 0, 1, 4};
    IS_TRUE(csr.col_index == expected_col_index);
    // Check elements of csr.row_index
    std::vector<usize> expected_row_index = {0, 2, 4, 7, 8, 11, 14};
    IS_TRUE(csr.row_index == expected_row_index);
}


void testLastLevel() {
    std::string file = "input/test1.mtx";
    CSR csr(file);

    auto ll_e = csr.getLastLevelAndEccentricity(0);
    std::vector<usize> v_exp = {2};
    IS_TRUE(v_exp == ll_e.first);
    IS_TRUE(ll_e.second == 2);

    ll_e = csr.getLastLevelAndEccentricity(1);
    v_exp = {2};
    IS_TRUE(v_exp == ll_e.first);
    IS_TRUE(ll_e.second == 2);

    ll_e = csr.getLastLevelAndEccentricity(2);
    v_exp = {0, 1};
    IS_TRUE(v_exp == ll_e.first);
    IS_TRUE(ll_e.second == 2);

    ll_e = csr.getLastLevelAndEccentricity(3);
    v_exp = {0, 1, 2};
    IS_TRUE(v_exp == ll_e.first);
    IS_TRUE(ll_e.second == 1);

    /// test2
    file = "input/test2.mtx";
    csr = CSR(file);

    ll_e = csr.getLastLevelAndEccentricity(0);
    v_exp = {1, 3, 4};
    sort(ll_e.first.begin(), ll_e.first.end());
    IS_TRUE(v_exp == ll_e.first);
    IS_TRUE(ll_e.second == 2);

    ll_e = csr.getLastLevelAndEccentricity(1);
    v_exp = {3};
    sort(ll_e.first.begin(), ll_e.first.end());
    IS_TRUE(v_exp == ll_e.first);
    IS_TRUE(ll_e.second == 3);

    ll_e = csr.getLastLevelAndEccentricity(2);
    v_exp = {1, 5};
    sort(ll_e.first.begin(), ll_e.first.end());
    IS_TRUE(v_exp == ll_e.first);
    IS_TRUE(ll_e.second == 2);

    ll_e = csr.getLastLevelAndEccentricity(3);
    v_exp = {1, 5};
    sort(ll_e.first.begin(), ll_e.first.end());
    IS_TRUE(v_exp == ll_e.first);
    IS_TRUE(ll_e.second == 3);

    /// MST (CLRS)
    file = "input/mst.mtx";
    csr = CSR(file);

    ll_e = csr.getLastLevelAndEccentricity(0);
    v_exp = {4};
    sort(ll_e.first.begin(), ll_e.first.end());
    IS_TRUE(v_exp == ll_e.first);
    IS_TRUE(ll_e.second == 4);

    ll_e = csr.getLastLevelAndEccentricity(1);
    v_exp = {4};
    sort(ll_e.first.begin(), ll_e.first.end());
    IS_TRUE(v_exp == ll_e.first);
    IS_TRUE(ll_e.second == 3);

    ll_e = csr.getLastLevelAndEccentricity(2);
    v_exp = {0, 4, 6, 7};
    sort(ll_e.first.begin(), ll_e.first.end());
    IS_TRUE(v_exp == ll_e.first);
    IS_TRUE(ll_e.second == 2);

    ll_e = csr.getLastLevelAndEccentricity(3);
    v_exp = {0, 7};
    sort(ll_e.first.begin(), ll_e.first.end());
    IS_TRUE(v_exp == ll_e.first);
    IS_TRUE(ll_e.second == 3);

    ll_e = csr.getLastLevelAndEccentricity(4);
    v_exp = {0};
    sort(ll_e.first.begin(), ll_e.first.end());
    IS_TRUE(v_exp == ll_e.first);
    IS_TRUE(ll_e.second == 4);

    ll_e = csr.getLastLevelAndEccentricity(5);
    v_exp = {0};
    sort(ll_e.first.begin(), ll_e.first.end());
    IS_TRUE(v_exp == ll_e.first);
    IS_TRUE(ll_e.second == 3);

    ll_e = csr.getLastLevelAndEccentricity(6);
    v_exp = {0, 1, 2, 3, 4};
    sort(ll_e.first.begin(), ll_e.first.end());
    IS_TRUE(v_exp == ll_e.first);
    IS_TRUE(ll_e.second == 2);

    ll_e = csr.getLastLevelAndEccentricity(7);
    v_exp = {3, 4};
    sort(ll_e.first.begin(), ll_e.first.end());
    IS_TRUE(v_exp == ll_e.first);
    IS_TRUE(ll_e.second == 3);

    ll_e = csr.getLastLevelAndEccentricity(8);
    v_exp = {4};
    sort(ll_e.first.begin(), ll_e.first.end());
    IS_TRUE(v_exp == ll_e.first);
    IS_TRUE(ll_e.second == 3);
}

void testDiameter() {
    std::string file = "input/test1.mtx";
    CSR csr(file, true);

    usize diameter = csr.getDiameter();
    IS_TRUE(csr.min_degree == 1);
    IS_TRUE(csr.max_degree == 3);
    IS_TRUE(std::ceil(static_cast<realT>(csr.max_degree) / 2.0) == 2);
    IS_TRUE(std::ceil((static_cast<float>(csr.m) - 1.0) / diameter) == 2);

    csr = CSR("input/test2.mtx", true);
    diameter = csr.getDiameter();
    IS_TRUE(csr.min_degree == 1);
    IS_TRUE(csr.max_degree == 3);
    IS_TRUE(std::ceil(static_cast<realT>(csr.max_degree) / 2.0) == 2);
    IS_TRUE(std::ceil((static_cast<float>(csr.m) - 1.0) / diameter) == 2);

    csr = CSR("input/mst.mtx", true);
    diameter = csr.getDiameter();
    IS_TRUE(csr.min_degree == 2);
    IS_TRUE(csr.max_degree == 4);
    IS_TRUE(std::ceil(static_cast<realT>(csr.max_degree) / 2.0) == 2);
    IS_TRUE(std::ceil((static_cast<float>(csr.m) - 1.0) / diameter) == 2);
}


void testProfile2() {
    CSR csr("input/test1.mtx");
    csr.evaluateProfile();
    IS_TRUE(csr.profile == 4);

    csr.labels = {3, 2, 1, 0};
    csr.evaluateProfile(); 
    IS_TRUE(csr.profile == 6);

    csr.labels = {3, 1, 2, 0};
    csr.evaluateProfile();
    IS_TRUE(csr.profile == 6);

    csr.labels = {1, 0, 3, 2};
    csr.evaluateProfile();
    IS_TRUE(csr.profile == 4);

    csr.labels = {1, 3, 0, 2};
    csr.evaluateProfile();
    IS_TRUE(csr.profile == 4);
}

void testProfile3() {
    CSR csr("input/test2.mtx");
    csr.evaluateProfile();
    IS_TRUE(csr.profile == 11);

    csr.labels = {4, 0, 1, 5, 3, 2};
    csr.evaluateProfile(); 
    IS_TRUE(csr.profile == 12);

    csr.labels = {1, 3, 5, 4, 0, 2};
    csr.evaluateProfile();
    IS_TRUE(csr.profile == 10);

    csr.labels = {3, 5, 1, 4, 2, 0};
    csr.evaluateProfile();
    IS_TRUE(csr.profile == 13);

    csr.labels = {3, 1, 4, 5, 2, 0};
    csr.evaluateProfile();
    IS_TRUE(csr.profile == 9);
}

void testEncoderDecoder() {
    CSR csr(5, 10);
    csr.random_keys.resize(5);
    csr.labels.resize(5);
    csr.tmp_rk.resize(5);
    csr.indexed_rk.resize(5);

    /// Decoder
    std::vector<usize> l = {0, 1, 2, 3, 4};
    csr.labels = l;
    encoder(csr);
    decoder(csr);
    IS_TRUE(l == csr.labels);

    l = {4, 2, 0, 3, 1};
    csr.labels = l;
    encoder(csr);
    decoder(csr);
    IS_TRUE(l == csr.labels);

    std::random_device rd;
    std::mt19937 gen(rd());
    shuffle(l.begin(), l.end(), gen);
    csr.labels = l;
    encoder(csr);
    decoder(csr);
    IS_TRUE(l == csr.labels);

    /// Decoder
    l = {0, 1, 2, 3, 4};
    csr.labels = l;
    encoder(csr);
    decoder(csr);
    IS_TRUE(l == csr.labels);

    l = {4, 2, 0, 3, 1};
    csr.labels = l;
    encoder(csr);
    decoder(csr);
    IS_TRUE(l == csr.labels);

    shuffle(l.begin(), l.end(), gen);
    csr.labels = l;
    encoder(csr);
    decoder(csr);
    IS_TRUE(l == csr.labels);
}


void testEccentricityNWidth() {
    CSR csr("input/test1.mtx");

    auto [eccentricity, width] = csr.getEccentricityNWidth(0);
    IS_TRUE(eccentricity == 2);
    IS_TRUE(width == 2);
    std::tie(eccentricity, width) = csr.getEccentricityNWidth(1);
    IS_TRUE(eccentricity == 2);
    IS_TRUE(width == 2);
    std::tie(eccentricity, width) = csr.getEccentricityNWidth(2);
    IS_TRUE(eccentricity == 2);
    IS_TRUE(width == 2);
    std::tie(eccentricity, width) = csr.getEccentricityNWidth(3);
    IS_TRUE(eccentricity == 1);
    IS_TRUE(width == 3);

    csr = CSR("input/test2.mtx");
    std::tie(eccentricity, width) = csr.getEccentricityNWidth(0);
    IS_TRUE(eccentricity == 2);
    IS_TRUE(width == 3);
    std::tie(eccentricity, width) = csr.getEccentricityNWidth(1);
    IS_TRUE(eccentricity == 3);
    IS_TRUE(width == 2);
    std::tie(eccentricity, width) = csr.getEccentricityNWidth(2);
    IS_TRUE(eccentricity == 2);
    IS_TRUE(width == 3);
    std::tie(eccentricity, width) = csr.getEccentricityNWidth(3);
    IS_TRUE(eccentricity == 3);
    IS_TRUE(width == 2);
    std::tie(eccentricity, width) = csr.getEccentricityNWidth(4);
    IS_TRUE(eccentricity == 2);
    IS_TRUE(width == 3);
    std::tie(eccentricity, width) = csr.getEccentricityNWidth(5);
    IS_TRUE(eccentricity == 3);
    IS_TRUE(width == 3);
}

int main() {
    std::cout << "\n Executing unit tests...";

    // "realT" must be double (in types.hpp) for the BRKGA's encoder to work
    testEncoderDecoder();
    testCSRFromFile2();
    testCSRFromFile3();
    testLastLevel(); 
    testDiameter();
    testEccentricityNWidth();
    testProfile2();
    testProfile3();

    if (!failed)
        std::cout << "\n All " << n <<" tests passed." << std::endl;
}