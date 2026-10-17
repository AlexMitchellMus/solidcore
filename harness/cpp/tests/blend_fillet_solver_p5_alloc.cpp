// SPDX-License-Identifier: Apache-2.0
// Standalone regression for allocation exceptions during rollover fixture parsing.
// This overrides global new only in this test executable, never in the shim library.
// From the repository root (SOLIDCORE_CPP_DIR can point at a shared reference build):
// nice -n 10 c++ -std=c++17 -O2 -ffp-contract=off -fno-fast-math -DNDEBUG \
//   -DNOMINMAX -DSM_GRAPHICS_CALLBACKS=1 -Iharness/cpp/shim \
//   -Iupstream/source/SMLib/inc -Iupstream/source/SM_API/inc \
//   harness/cpp/tests/blend_fillet_solver_p5_alloc.cpp \
//   harness/cpp/shim/blend_fillet_solver_p5.cpp \
//   -L"${SOLIDCORE_CPP_DIR:-bin/cpp}/lib" -lSM_API -lSMLib -pthread \
//   -o bin/blend_fillet_solver_p5_alloc
// nice -n 10 bin/blend_fillet_solver_p5_alloc

#include "sc_shim.h"
#include <cstdio>
#include <cstdlib>
#include <new>
static bool fail_allocation = false;
void* operator new(std::size_t size) {
    if (fail_allocation) throw std::bad_alloc();
    if (void* p = std::malloc(size)) return p;
    throw std::bad_alloc();
}
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
extern "C" int32_t sc_blend_fillet_solver_p5_eval(const double*,int64_t,double*,int64_t,int64_t*);
int main() {
    double in[] = {0,0,-1,0,1,0,1,102,1,0,1,103};
    double out[128]; int64_t len=0;
    try {
        fail_allocation = true;
        int32_t status = sc_blend_fillet_solver_p5_eval(in,12,out,128,&len);
        fail_allocation = false;
        std::printf("allocation failure contained: status=%d\n",status);
        return status == SC_ERR_EXCEPTION ? 0 : 2;
    } catch (...) {
        fail_allocation = false;
        std::puts("allocation exception escaped the C boundary");
        return 1;
    }
}
