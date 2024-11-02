#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include <omp.h>

#define COEF 1.234 // Need to be great then 1
#define MAX_POSSIBLE_NEIB 5

// Enum to distinguish between upper triangular and lower triangular vortex
enum TriangularType { Upper, Lower };

/* Get vertex number by cell number
 * # Arguments:
 * * cell - cell number
 * * k1 - square cells sequence length
 * * k2 - triangular cells sequence length
 * * type - expected triangular type
 */
size_t get_vertex(
    size_t cell,
    size_t k1,
    size_t k2,
    TriangularType type
) {
    size_t div = cell / (k1 + k2);
    size_t mod = cell % (k1 + k2);
    size_t vertex = cell + div * k2;
    if (mod >= k2) { // Triangular
        vertex += mod - k2;
        vertex += type == TriangularType::Upper ? 1 : 0;
    }
    return vertex;
}

/* Generate CSR portrait by grid params
 * # Arguments:
 * * nx - grid hieght
 * * ny - grid width 
 * * k1 - square cells sequence length
 * * k2 - triangular cells sequence length
 * * tn - tread number
 * * dl - debug level
 * # Return value
 * * IA
 * * JA
 */
auto gen(
    size_t nx,
    size_t ny,
    size_t k1,
    size_t k2,
    size_t tn,
    size_t dl
) {
    // Get ia/ja sizes
    size_t div = (nx * ny) / (k1 + k2);
    size_t mod = (nx * ny) % (k1 + k2);
    // Square vertices number
    size_t ns = k1 * div;
    // Triangular vertices number
    size_t nt = k2 * div * 2;
    size_t iters = div;
    if (mod > 0) {
        iters++;
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
    if (dl >= 1) {
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
    #pragma omp parallel
    {
        #pragma omp for
        for (size_t iter = 0; iter < iters; iter++) {
            // NOTE: Define and count this here by parallel issues
            // NOTE: In start of iteration we are always on square vertex
            size_t v = iter * (k1 + k2 * 2); // Cur vertex num
            size_t c = iter * (k1 + k2); // Cur cell number
            size_t i = c / ny; // Cur cell row number
            size_t j = c % ny; // Cur cell col number
            // Some crazy loop unrolling
            // Square cells
            for (int counter = 0; counter < k1 && v < nv; counter++) {
                unsigned char rel = 0;
                if (i > 0) { // Upper neighbor
                    size_t neighbor = get_vertex(c - ny, k1, k2, TriangularType::Lower);
                    dist_ja[v][rel++] = neighbor;
                }
                if (j > 0) { // Left neighbor
                    dist_ja[v][rel++] = v - 1;
                }
                dist_ja[v][rel++] = v; // Self
                if (j < ny - 1) { // Right neighbor
                    dist_ja[v][rel++] = v + 1;
                }
                if (i < nx - 1) { // Lower neighbor
                    size_t neighbor = get_vertex(c + ny, k1, k2, TriangularType::Upper);
                    dist_ja[v][rel++] = neighbor;
                }
                ia[v + 1] = rel;
                // Turn to next cell
                j = (j + 1) % ny;
                if (j == 0) {
                    i++;
                }
                v++;
                c++;
            }
            // Triangular cells
            for (int counter = 0; counter < k2 && v < nv; counter++) {
                // Lower
                unsigned char rel = 0;
                if (j > 0) { // Left neighbor
                    dist_ja[v][rel++] = v - 1;
                }
                dist_ja[v][rel++] = v; // Self
                dist_ja[v][rel++] = v + 1; // Pair upper triangular cell
                if (i < nx - 1) { // Lower neighbor
                    size_t neighbor = get_vertex(c + ny, k1, k2, TriangularType::Upper);
                    dist_ja[v][rel++] = neighbor;
                }
                ia[v + 1] = rel;
                // Upper
                rel = 0;
                if (i > 0) { // Upper neighbor
                    size_t neighbor = get_vertex(c - ny, k1, k2, TriangularType::Lower);
                    dist_ja[v + 1][rel++] = neighbor;
                }
                dist_ja[v + 1][rel++] = v; // Pair lower triangular cell
                dist_ja[v + 1][rel++] = v + 1; // Self
                if (j < ny - 1) { // Right neighbor
                    dist_ja[v + 1][rel++] = v + 2;
                }
                ia[v + 2] = rel;
                // Turn to next cell
                j = (j + 1) % ny;
                if (j == 0) {
                    i++;
                }
                v += 2;
                c++;
            }
        }
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
    std::vector<size_t>& ja,
    size_t tn
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
        std::cout << "Dl is debug level, may be 0, 1 and 2\n";
        return 1;
    }
    size_t nx = std::stoll(argv[1]);
    size_t ny = std::stoll(argv[2]);
    size_t k1 = std::stoll(argv[3]);
    size_t k2 = std::stoll(argv[4]);
    size_t tn = std::stoll(argv[5]);
    size_t dl = std::stoll(argv[6]);

    omp_set_num_threads(tn);

    // Gen ia/ja
    if (dl >= 1) {
        std::cout << "Generate IA/JA" << std::endl;
    }
    //double start = omp_get_wtime();
    auto [ia, ja] = gen(nx, ny, k1, k2, tn, dl);
    if (dl >= 2) {
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
    if (dl >= 1) {
        std::cout << "Fill A:" << std::endl;
    }
    auto [a, b] = fill(ia, ja, tn);
    //double end = omp_get_wtime();
    //std::cout << end - start << std::endl;
    if (dl >= 2) {
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
