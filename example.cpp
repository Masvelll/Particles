#include <extern/pybind11/pybind11/pybind11.h>

namespace py = pybind11;

int64_t multiply(int64_t a, int64_t b){
    return a * b;
};

PYBIND11_MODULE(mol_dynamics, handle){
    handle.doc() = "This is the module docs.";
    handle.def("update", &multiply);
}
