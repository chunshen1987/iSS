// Prototype: the FSSW yield loop (bulk_deltaf_kind = 20) on CPU and CUDA.
// See doc/GPU_feasibility.md for the results.
//
// It is standalone and not part of the iSS build. It times the per-cell,
// per-species yield of FSSW::calculate_dN_analytic (22-moment bulk delta f,
// division-free form) with OpenMP and with CUDA in FP64 and FP32, and the
// CUDA start-up and transfers, and compares the results with the CPU.
// The cells are random but span the T and mu of the benchmark surfaces.
// The species are the particle lines of the UrQMD table with mass above
// 0.1 GeV and their antiparticles (369, iSS samples 321 of them).
//
//   nvcc -O3 -std=c++17 -arch=sm_75 -Xcompiler -fopenmp,-O3 \
//       gpu_yield_bench.cu -o gpu_yield_bench
//   ./gpu_yield_bench [n_cells] [iSS_tables/pdg-urqmd_v3.3+.dat]
#include <cuda_runtime.h>
#include <omp.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#define CK(x) do { cudaError_t e = (x); if (e != cudaSuccess) { \
    printf("CUDA %s at %d\n", cudaGetErrorString(e), __LINE__); exit(1);} } while (0)

struct Species { double mass, prefactor; int sign, baryon, strange, charge; };
struct Cell { float T, muB, muS, muQ, dsigma0, bulkPi; double beta; double c[6]; };

static const double X_MIN = 0.5, X_MAX = 400, DX = 0.05;
static const int NTAB = static_cast<int>((X_MAX - X_MIN)/DX) + 1;
__constant__ Species d_species[400];

double now() {
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
}

template <typename R>
__host__ __device__ inline R yield(const Species &sp, const Cell &cell,
                                   const R *table) {
    const R T = cell.T, beta = static_cast<R>(cell.beta);
    const R mu = sp.baryon*cell.muB + sp.strange*cell.muS + sp.charge*cell.muQ;
    const R mass = sp.mass;
    const R lambda = exp(beta*mu);
    const int order = (mass < R(0.7) && T > R(0.05)) ? 10 : 1;
    R fug = lambda, theta = 1, Neq = 0, t1 = 0, t2 = 0, t3 = 0;
    const R mbeta = mass*beta;
    for (int n = 1; n <= order; ++n) {
        const R arg = n*mbeta, inv_n = R(1)/n;
        const R pos = (arg - R(X_MIN))*R(1/DX);
        const int i = static_cast<int>(pos);
        const R f = pos - i;
        const R *lo = table + 4*i, *hi = lo + 4;
        const R K1 = (1-f)*lo[0] + f*hi[0];
        const R K2 = (1-f)*lo[1] + f*hi[1];
        const R K3 = (1-f)*lo[2] + f*hi[2];
        Neq += theta*inv_n*fug*K2;
        t1 += theta*fug*K2;
        t2 += theta*fug*(mbeta*K1 + 3*K2*inv_n);
        t3 += theta*fug*(mbeta*K2 + 3*K3*inv_n);
        fug *= lambda; theta *= -sp.sign;
    }
    Neq *= mass*mass*T;
    t1 *= mass*mass*T; t2 *= mass*mass*T*T; t3 *= mass*mass*mass*T*T;
    const R unit = R(1/(0.19733*0.19733*0.19733));
    const R ds = cell.dsigma0, Pi = cell.bulkPi;
    const R dNbulk = unit*sp.prefactor*ds*Pi*(
        t1*mass*mass*R(cell.c[2]) + t2*(sp.baryon*R(cell.c[3])
        + sp.strange*R(cell.c[4]) + sp.charge*R(cell.c[5]))
        + t3*(R(cell.c[1]) - R(cell.c[2])));
    const R total = unit*sp.prefactor*ds*Neq + dNbulk;
    return total > 0 ? total : 0;
}

template <typename R>
__global__ void yield_kernel(const Cell *cells, long ncell, int nsp,
                             const R *table, R *out) {
    const long l = blockIdx.x*static_cast<long>(blockDim.x) + threadIdx.x;
    if (l >= ncell) return;
    const Cell cell = cells[l];
    for (int s = 0; s < nsp; ++s)
        out[s*ncell + l] = yield<R>(d_species[s], cell, table);
}

int main(int argc, char **argv) {
    const long ncell = argc > 1 ? atol(argv[1]) : 300204;
    // species: particle lines of the UrQMD table plus antiparticles
    std::vector<Species> sp;
    std::ifstream pdg(argc > 2 ? argv[2] : "iSS_tables/pdg-urqmd_v3.3+.dat");
    if (!pdg) { printf("cannot open the particle table\n"); return 1; }
    std::string line;
    while (std::getline(pdg, line)) {
        std::istringstream is(line);
        long id; std::string name; double m, w; int g, b, s, c, bt, iso, q, nd;
        if (!(is >> id >> name)) continue;
        if (name.find_first_not_of("0123456789-.") == std::string::npos) continue;
        if (!(is >> m >> w >> g >> b >> s >> c >> bt >> iso >> q >> nd)) continue;
        if (m < 0.1) continue;
        Species x{m, g/(2*M_PI*M_PI), b != 0 ? 1 : -1, b, s, q};
        sp.push_back(x);
        if (b != 0 || s != 0 || q != 0) { x.baryon = -b; x.strange = -s; x.charge = -q; sp.push_back(x); }
    }
    const int nsp = static_cast<int>(std::min<size_t>(sp.size(), 400));
    // cells: varying T and mu, as on the benchmark surfaces
    std::mt19937 gen(1);
    std::uniform_real_distribution<double> u(0., 1.);
    std::vector<Cell> cells(ncell);
    for (auto &c : cells) {
        c.T = 0.135 + 0.025*u(gen); c.muB = 0.3*u(gen); c.muS = c.muB/3;
        c.muQ = -0.02*u(gen); c.dsigma0 = 0.0156*(0.5 + u(gen));
        c.bulkPi = -0.015*(1 + 0.3*u(gen)); c.beta = 1./c.T;
        for (int k = 0; k < 6; ++k) c.c[k] = 1 + u(gen);
    }
    std::vector<double> table(4*(NTAB + 1));
    for (int i = 0; i < NTAB; ++i) {
        const double x = X_MIN + i*DX;
        table[4*i] = std::cyl_bessel_k(1., x);
        table[4*i+1] = std::cyl_bessel_k(2., x);
        table[4*i+2] = std::cyl_bessel_k(3., x);
    }
    std::vector<float> tablef(table.begin(), table.end());
    printf("%ld cells, %d species\n", ncell, nsp);

    // CPU: species outer, cells inner, as in FSSW
    std::vector<double> cpu(static_cast<size_t>(nsp)*ncell);
    for (int rep = 0; rep < 2; ++rep) {
        double t0 = now();
        for (int s = 0; s < nsp; ++s) {
            #pragma omp parallel for schedule(static)
            for (long l = 0; l < ncell; ++l)
                cpu[s*ncell + l] = yield<double>(sp[s], cells[l], table.data());
        }
        printf("CPU %d threads: %.3f s\n", omp_get_max_threads(), now() - t0);
    }

    double t0 = now();
    CK(cudaFree(0));
    printf("CUDA context: %.3f s\n", now() - t0);
    CK(cudaMemcpyToSymbol(d_species, sp.data(), nsp*sizeof(Species)));
    Cell *d_cells; double *d_tab, *d_out; float *d_tabf, *d_outf;
    CK(cudaMalloc(&d_cells, ncell*sizeof(Cell)));
    CK(cudaMalloc(&d_tab, table.size()*sizeof(double)));
    CK(cudaMalloc(&d_tabf, table.size()*sizeof(float)));
    CK(cudaMalloc(&d_out, cpu.size()*sizeof(double)));
    CK(cudaMalloc(&d_outf, cpu.size()*sizeof(float)));
    t0 = now();
    CK(cudaMemcpy(d_cells, cells.data(), ncell*sizeof(Cell), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_tab, table.data(), table.size()*8, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_tabf, tablef.data(), table.size()*4, cudaMemcpyHostToDevice));
    printf("upload cells + tables: %.3f s\n", now() - t0);

    const int block = 128;
    const int grid = static_cast<int>((ncell + block - 1)/block);
    std::vector<double> gpu(cpu.size());
    std::vector<float> gpuf(cpu.size());
    for (int rep = 0; rep < 2; ++rep) {
        t0 = now();
        yield_kernel<double><<<grid, block>>>(d_cells, ncell, nsp, d_tab, d_out);
        CK(cudaDeviceSynchronize());
        double t1 = now();
        CK(cudaMemcpy(gpu.data(), d_out, gpu.size()*8, cudaMemcpyDeviceToHost));
        double t2 = now();
        yield_kernel<float><<<grid, block>>>(d_cells, ncell, nsp, d_tabf, d_outf);
        CK(cudaDeviceSynchronize());
        double t3 = now();
        CK(cudaMemcpy(gpuf.data(), d_outf, gpuf.size()*4, cudaMemcpyDeviceToHost));
        double t4 = now();
        printf("GPU fp64 kernel %.3f s, copy back %.3f s (%.0f MB) | "
               "fp32 kernel %.3f s, copy back %.3f s\n", t1 - t0, t2 - t1,
               gpu.size()*8/1e6, t3 - t2, t4 - t3);
    }
    double e64 = 0, e32 = 0, s64 = 0, s32 = 0, scpu = 0;
    for (size_t i = 0; i < cpu.size(); ++i) {
        if (cpu[i] > 0) {
            e64 = std::max(e64, std::abs(gpu[i] - cpu[i])/cpu[i]);
            e32 = std::max(e32, std::abs(gpuf[i] - cpu[i])/cpu[i]);
        }
    }
    for (int s = 0; s < nsp; ++s) {   // per-species totals, as FSSW sums them
        double a = 0, b = 0, c = 0;
        for (long l = 0; l < ncell; ++l) { a += cpu[s*ncell+l]; b += gpu[s*ncell+l]; c += gpuf[s*ncell+l]; }
        if (a > 0) { s64 = std::max(s64, std::abs(b-a)/a); s32 = std::max(s32, std::abs(c-a)/a); }
        scpu += a;
    }
    printf("max rel. diff per cell: fp64 %.1e, fp32 %.1e; per species total: fp64 %.1e, fp32 %.1e\n",
           e64, e32, s64, s32);
    printf("total yield %.1f\n", scpu);
}
