#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include <omp.h>

#define COEF 1.234 // Need to be great then 1
#define MAX_POSSIBLE_NEIB 5

enum VertexType { Square, Lower, Upper };

enum LogLevel { NoLog, TimeLog, InfoLog, ArrayLog };

/* Get cell number and vertex type by vertex number
 * # Arguments:
 * * vertex - vertex number
 * * k1 - square cells sequence length
 * * k2 - triangular cells sequence length
 * # Return values:
 * cell number
 * vertex type
 */
inline std::pair<size_t, VertexType> get_cell_n_type(
    size_t vertex,
    size_t k1,
    size_t k2
) {
    size_t div = vertex / (k1 + k2 * 2);
    size_t mod = vertex % (k1 + k2 * 2);
    size_t cell = vertex - div * k2;
    VertexType type = Square;
    if (mod >= k1) {
        cell -= (mod - k1 + 1) / 2;
        type = (mod - k1 + 1) % 2 ? Lower : Upper;
    }
    return {cell, type};
}

/* Get vertex number by cell number
 * # Arguments:
 * * cell - cell number
 * * k1 - square cells sequence length
 * * k2 - triangular cells sequence length
 * * type - expected vertex type (valuable for triangular)
 */
inline size_t get_vertex(
    size_t cell,
    size_t k1,
    size_t k2,
    VertexType type
) {
    size_t div = cell / (k1 + k2);
    size_t mod = cell % (k1 + k2);
    size_t vertex = cell + div * k2;
    if (mod >= k1) { // Triangular
        vertex += mod - k1;
        if (type == Upper) {
            vertex++;
        }
    }
    return vertex;
}

/* Generate CSR portrait by grid params
 * # Arguments:
 * * nx - grid hieght
 * * ny - grid width 
 * * k1 - square cells sequence length
 * * k2 - triangular cells sequence length
 * * ll - log level
 * # Return values:
 * * IA
 * * JA
 */
auto gen(
    size_t nx,
    size_t ny,
    size_t k1,
    size_t k2,
    LogLevel ll
) {
    // Get ia/ja sizes
    size_t div = (nx * ny) / (k1 + k2);
    size_t mod = (nx * ny) % (k1 + k2);
    // Square vertices number
    size_t ns = k1 * div;
    // Triangular vertices number
    size_t nt = k2 * div * 2;
    if (mod > 0) {
        ns += k1 > mod ? mod : k1;
        mod = k1 > mod ? 0 : mod - k1;
        nt += 2 * mod;
    }
    // Vertices number
    size_t nv = ns + nt;
    // Edges number (doubled)
    size_t ne2 = 4 * ns + 3 * nt - 2 * (nx + ny);
    // IA size
    size_t size_ia = nv + 1;
    // JA size
    size_t size_ja = nv + ne2;
    if (ll >= InfoLog) {
        std::cout << "Square cells number: " << ns << "\n";
        std::cout << "Triangular cells number: " << nt << "\n";
        std::cout << "Cells number: " << nv << "\n";
        std::cout << "Edges number: " << ne2 << "\n";
        std::cout << "IA size: " << size_ia << "\n";
        std::cout << "JA size: " << size_ja << "\n";
    }
    std::vector<size_t> ia(size_ia);
    ia[0] = 0;
    // Memory: 4|V|, instead |E| for ja :(
    std::vector<std::vector<size_t>> dist_ja(nv);
    #pragma omp parallel
    {
        #pragma omp for
        for (auto& dist_val : dist_ja) {
            dist_val = std::vector<size_t>(MAX_POSSIBLE_NEIB);
        }
    }
    // Fill vecs
    // Go through vertices nums
    double start = omp_get_wtime();
    #pragma omp parallel
    {
        #pragma omp for
        for (size_t v = 0; v < nv; v++) {
            // NOTE: Define and count this here by parallel issues
            // NOTE: In start of iteration we are always on square vertex
            // Cur cell number and vertex type
            auto [c, t] = get_cell_n_type(v, k1, k2);
            size_t i = c / ny; // Cur cell row number
            size_t j = c % ny; // Cur cell col number
            unsigned char neibs = 0; // Vertex neighbor number
            switch (t) {
                case Square:
                    if (i > 0) { // Upper neighbor
                        size_t neighbor = get_vertex(c - ny, k1, k2, Lower);
                        dist_ja[v][neibs++] = neighbor;
                    }
                    if (j > 0) { // Left neighbor
                        dist_ja[v][neibs++] = v - 1;
                    }
                    dist_ja[v][neibs++] = v; // Self
                    if (j < ny - 1) { // Right neighbor
                        dist_ja[v][neibs++] = v + 1;
                    }
                    if (i < nx - 1) { // Lower neighbor
                        size_t neighbor = get_vertex(c + ny, k1, k2, Upper);
                        dist_ja[v][neibs++] = neighbor;
                    }
                    break;
                case Lower:
                    if (j > 0) { // Left neighbor
                        dist_ja[v][neibs++] = v - 1;
                    }
                    dist_ja[v][neibs++] = v; // Self
                    dist_ja[v][neibs++] = v + 1; // Pair upper triangular cell
                    if (i < nx - 1) { // Lower neighbor
                        size_t neighbor = get_vertex(c + ny, k1, k2, Upper);
                        dist_ja[v][neibs++] = neighbor;
                    }
                    break;
                case Upper:
                    if (i > 0) { // Upper neighbor
                        size_t neighbor = get_vertex(c - ny, k1, k2, Lower);
                        dist_ja[v][neibs++] = neighbor;
                    }
                    dist_ja[v][neibs++] = v - 1; // Pair lower triangular cell
                    dist_ja[v][neibs++] = v; // Self
                    if (j < ny - 1) { // Right neighbor
                        dist_ja[v][neibs++] = v + 1;
                    }
                    break;
            }
            ia[v + 1] = neibs;
        }
    }
    double end = omp_get_wtime();
    if (ll >= TimeLog) {
        std::cout << "Iter time:\t" << end - start << std::endl;
    }
    // Adjust ia vals
    int prev = 0;
    for (int i = 1; i < size_ia; i++) {
        ia[i] += prev;
        prev = ia[i];
    }
    // Concat dist_ja to ja
    std::vector<size_t> ja(size_ja);
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < nv; i++) {
            auto pos = ia[i];
            auto begin = dist_ja[i].begin();
            auto end = begin + (ia[i + 1] - ia[i]);
            std::copy(begin, end, ja.begin() + pos);
        }
    }
    return std::pair(std::move(ia), std::move(ja));
}

inline float filler(size_t i, size_t j) {
    return std::cos(i * j + i + j);
}

auto fill(
    std::vector<size_t>& ia,
    std::vector<size_t>& ja
) {
    std::vector<float> a(ja.size());
    std::vector<float> b(ia.size() - 1);
    #pragma omp parallel
    {
        #pragma omp for
        for (size_t i = 0; i < ia.size() - 1; i++) {
            size_t diag_ind;
            float sum = 0;
            for (size_t ind = ia[i]; ind < ia[i + 1]; ind++) {
                size_t j = ja[ind];
                if (j != i) {
                    a[ind] = filler(i, j);
                    sum += std::abs(a[ind]);
                } else {
                    diag_ind = ind;
                }
            }
            a[diag_ind] = COEF * sum;
            b[i] = std::sin(i);
        }
    }
    return std::pair(std::move(a), std::move(b));
}

int main(int argc, char** argv) {
    if (argc != 7) {
        // Help
        std::cout << "Usage: " << argv[0] << " Nx Ny K1 K2 Tn Dl\n";
        std::cout << "Where:\n";
        std::cout << "Nx is positive int that represents grid hieght\n";
        std::cout << "Ny is positive int that represents grid width\n";
        std::cout << "K1 is positive int that represents square cells sequence length\n";
        std::cout << "K2 is positive int that represents triangular cells sequence length\n";
        std::cout << "Tn is tread number";
        std::cout << "Ll is log level:\n";
        std::cout << "\t<=" << NoLog << " - no logs\n";
        std::cout << "\t>=" << TimeLog << " - show time\n";
        std::cout << "\t>=" << InfoLog << " - show info\n"; // TODO: change description
        std::cout << "\t>=" << ArrayLog << " - show arrays\n";
        return 1;
    }
    size_t nx = std::stoll(argv[1]);
    size_t ny = std::stoll(argv[2]);
    size_t k1 = std::stoll(argv[3]);
    size_t k2 = std::stoll(argv[4]);
    size_t tn = std::stoll(argv[5]);
    LogLevel ll = (LogLevel) std::stoi(argv[6]);

    omp_set_num_threads(tn);

    // Gen ia/ja
    if (ll >= InfoLog) {
        std::cout << "Generate IA/JA" << std::endl;
    }
    double start = omp_get_wtime();
    auto [ia, ja] = gen(nx, ny, k1, k2, ll);
    double end_gen = omp_get_wtime();
    if (ll >= ArrayLog) {
        std::cout << "IA:" << std::endl;
        for (const auto val: ia) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        std::cout << "JA:" << std::endl;
        for (const auto val: ja) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }

    // Fill a
    if (ll >= InfoLog) {
        std::cout << "Fill A:" << std::endl;
    }
    double start_fill = omp_get_wtime();
    auto [a, b] = fill(ia, ja);
    double end = omp_get_wtime();
    if (ll >= TimeLog) {
        std::cout << "Gen time:\t" << end_gen - start << std::endl;
        std::cout << "Fill time:\t" << end - start_fill << std::endl;
        std::cout << "Total time:\t" << end - start << std::endl;
    }
    if (ll >= ArrayLog) {
        std::cout << "A" << std::endl;
        for (const auto val: a) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        std::cout << "b" << std::endl;
        for (const auto val: b) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}
