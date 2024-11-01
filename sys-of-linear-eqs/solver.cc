#include <iostream>
#include <utility>
#include <vector>

// Enum to distinguish between upper triangular and lower triangular vortex
enum TriangularType { Upper, Lower };

/* Get vertex number my cell number
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
 * # Return value
 * * IA
 * * JA
 */
std::pair<std::vector<size_t>, std::vector<size_t>>
gen(
    size_t nx,
    size_t ny,
    size_t k1,
    size_t k2
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
    /*
    std::cout << "Square cells number: " << ns << "\n";
    std::cout << "Triangular cells number: " << nt << "\n";
    std::cout << "Cells number: " << nv << "\n";
    std::cout << "Edges number: " << ne2 << "\n";
    std::cout << "IA size: " << size_ia << "\n";
    std::cout << "JA size: " << size_ja << "\n";
    */
    std::vector<size_t> ia(size_ia);
    std::vector<size_t> ja(size_ja);
    // Fill vecs
    // Go through vertices nums
    size_t v = 0; // Cur vertex num
    size_t c = 0; // Cur cell number
    size_t i = 0; // Cur cell row number
    size_t j = 0; // Cur cell col number
    size_t ja_counter = 0;
    size_t ia_counter = 1;
    ia[0] = 0;
    while (v < nv) {
        // Some crazy loop unrolling
        // Square cells
        for (int counter = 0; counter < k1 && v < nv; counter++) {
            int rel = 1;
            if (i > 0) { // Upper neighbor
                size_t neighbor = get_vertex(c - ny, k1, k2, TriangularType::Lower);
                ja[ja_counter++] = neighbor;
                rel++;
            }
            if (j > 0) { // Left neighbor
                size_t neighbor = get_vertex(c - 1, k1, k2, TriangularType::Upper);
                ja[ja_counter++] = neighbor;
                rel++;
            }
            ja[ja_counter++] = v; // Self
            if (j < ny - 1) { // Right neighbor
                size_t neighbor = get_vertex(c + 1, k1, k2, TriangularType::Lower);
                ja[ja_counter++] = neighbor;
                rel++;
            }
            if (i < nx - 1) { // Lower neighbor
                size_t neighbor = get_vertex(c + ny, k1, k2, TriangularType::Upper);
                ja[ja_counter++] = neighbor;
                rel++;
            }
            ia[ia_counter++] = ia[ia_counter - 1] + rel;
            j = (j + 1) % ny;
            if (j == 0) {
                i += 1;
            }
            v += 1;
            c += 1;
        }
        // Triangular cells
        for (int counter = 0; counter < k2 && v < nv; counter++) {
            // Lower
            int rel = 2;
            if (j > 0) { // Left neighbor
                size_t neighbor = get_vertex(c - 1, k1, k2, TriangularType::Upper);
                ja[ja_counter++] = neighbor;
                rel++;
            }
            ja[ja_counter++] = v; // Self
            ja[ja_counter++] = v + 1; // Upper
            if (i < nx - 1) { // Lower neighbor
                size_t neighbor = get_vertex(c + ny, k1, k2, TriangularType::Upper);
                ja[ja_counter++] = neighbor;
                rel++;
            }
            ia[ia_counter++] = ia[ia_counter - 1] + rel;
            // Upper
            rel = 2;
            if (i > 0) { // Upper neighbor
                size_t neighbor = get_vertex(c - ny, k1, k2, TriangularType::Lower);
                ja[ja_counter++] = neighbor;
                rel++;
            }
            ja[ja_counter++] = v ; // Lower
            ja[ja_counter++] = v + 1; // Self
            if (j < ny - 1) { // Right neighbor
                size_t neighbor = get_vertex(c + 1, k1, k2, TriangularType::Lower);
                ja[ja_counter++] = neighbor;
                rel++;
            }
            ia[ia_counter++] = ia[ia_counter - 1] + rel;
            j = (j + 1) % ny;
            if (j == 0) {
                i += 1;
            }
            v += 2;
            c += 1;
        }
    }
    return {ia, ja};
}

int main(int argc, char** argv) {
    if (argc != 5) {
        // Help
        std::cout << "Usage: " << argv[0] << " Nx Ny K1 K2\n";
        std::cout << "Where:\n";
        std::cout << "Nx is positive int that represents grid hieght\n";
        std::cout << "Ny is positive int that represents grid width\n";
        std::cout << "K1 is positive int that represents square cells sequence length\n";
        std::cout << "K2 is positive int that represents triangular cells sequence length\n";
        return 1;
    }
    size_t nx = std::stoll(argv[1]);
    size_t ny = std::stoll(argv[2]);
    size_t k1 = std::stoll(argv[3]);
    size_t k2 = std::stoll(argv[4]);
    auto [ia, ja] = gen(nx, ny, k1, k2);
    /*
    std::cout << "IA: ";
    for (const auto a: ia) {
        std::cout << a << " ";
    }
    std::cout << "\n";
    std::cout << "JA: ";
    for (const auto a: ja) {
        std::cout << a << " ";
    }
    std::cout << "\n";
    */
    return 0;
}
