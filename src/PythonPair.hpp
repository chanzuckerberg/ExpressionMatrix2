// Some magic needed to convert Python tuples with two elements
// to std::pair and vice versa.
// Adapted from http://cci.lbl.gov/cctbx_sources/boost_adaptbx/std_pair_conversion.h

// WE ARE ACTUALLY NOT USING THIS - IT IS NOT INCLUDED ANYWERE


#ifndef CZI_EXPRESSION_MATRIX2_PYTHON_PAIR_HPP
#define CZI_EXPRESSION_MATRIX2_PYTHON_PAIR_HPP

#include "utility.hpp"

#include <boost/python/tuple.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/to_python_converter.hpp>

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {



        template <typename T, typename U> class StdPairToPythonTupleLowLevel {
        public:
            static PyObject* convert(const pair<T,U>& p)
            {
                return boost::python::incref(boost::python::make_tuple(p.first, p.second).ptr());
            }

            static PyTypeObject const *get_pytype()
            {
                return &PyTuple_Type;
            }
        };



        template <typename T, typename U> class StdPairToPythonTuple {
        public:
            StdPairToPythonTuple()
            {
                boost::python::to_python_converter<pair<T,U>, StdPairToPythonTupleLowLevel<T,U>, true>();
            }
        };



        template <typename T, typename U> class StdPairFromPythonTuple {
        public:
            StdPairFromPythonTuple()
            {
                using namespace boost::python::converter;
                registry::push_back(
                    &convertible,
                    &construct,
                    boost::python::type_id< pair<T,U> >(),
                    get_pytype
                );
            }

            static const PyTypeObject* get_pytype()
            {
                return &PyTuple_Type;
            }

            static void *convertible(PyObject* o) {
                using namespace boost::python;
                if (!PyTuple_Check(o) || PyTuple_GET_SIZE(o) != 2) {
                    return 0;
                }
                return o;
            }

            static void construct(
                PyObject* o,
                boost::python::converter::rvalue_from_python_stage1_data* data)
            {
                using boost::python::extract;
                using namespace boost::python::converter;
                PyObject *first  = PyTuple_GET_ITEM(o, 0);
                PyObject* second = PyTuple_GET_ITEM(o, 1);
                void* storage = ((rvalue_from_python_storage<std::pair<T,U> >*) data)->storage.bytes;
                new (storage) pair<T,U>(extract<T>(first), extract<U>(second));
                data->convertible = storage;
            }
        };



        template <typename T, typename U> class StdPairToAndFromPythonTuple {
        public:
            StdPairToAndFromPythonTuple()
            {
                StdPairToPythonTuple<T,U>();
                StdPairFromPythonTuple<T,U>();
            }
        };


    }
}

#endif
