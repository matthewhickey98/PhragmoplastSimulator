#include "mex.hpp"
#include "mexAdapter.hpp"

using namespace matlab::data;
using matlab::mex::ArgumentList;

class Simulate : public matlab::mex::Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs)
    {

        // Implement function
        TypedArray<double> doubleArray = std::move(inputs[0]);
        for (auto& elem : doubleArray)
        {
            elem *= 2;
        }

        // Assign outputs
        outputs[0] = doubleArray;
    }