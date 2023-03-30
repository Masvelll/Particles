#include <pybind11/pybind11.h>

namespace py = pybind11;

int64_t multiply(int64_t a, int64_t b){
    return a * b;
};

PYBIND11_MODULE(example, m){
    m.doc() = "This is the module docs.";
    m.def("multiply", &multiply, "A function that multiplies two numbers");
}

//    To make file
//c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3-config --includes) -I extern/pybind11/include example.cpp -o example$(python3-config --extension-suffix)
