#include <iostream>
#include <utility>
#include <vector>

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
    std::cout << "div: " << div << "\n";
    std::cout << "mod: " << mod << "\n";
    // Square vertices number
    size_t ns = k1 * div;
    // Triangular vertices number
    size_t nt = k2 * div * 2;
    if (mod > 0) {
        ns += k1 > mod ? mod : k1;
        mod = k1 > mod ? 0 : mod - k1;
        nt += 2 * mod;
    }
    std::cout << "Square cells number: " << ns << "\n";
    std::cout << "Triangular cells number: " << nt << "\n";
    // Vertices number
    size_t nv = ns + nt;
    std::cout << "Cells number: " << nv << "\n";
    // Edges number (doubled)
    size_t ne2 = 4 * ns + 3 * nt - 2 * (nx + ny);
    std::cout << "Edges number: " << ne2 << "\n";
    // IA size
    size_t size_ia = nv + 1;
    std::cout << "IA size: " << size_ia << "\n";
    // JA size
    size_t size_ja = nv + ne2;
    std::cout << "JA size: " << size_ja << "\n";
    std::vector<int> ia(size_ia);
    std::vector<int> ja(size_ja);
    return {{}, {}};
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
    return 0;
}
