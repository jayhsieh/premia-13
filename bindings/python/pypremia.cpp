//  Copyright Joel de Guzman 2002-2004. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//  Hello World Example from the tutorial
//  [Joel de Guzman 10/9/2002]

#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/list.hpp>
#include <boost/python/exception_translator.hpp>

namespace py = boost::python;

#include <premia/import.h>
#include <premia/exception.h>
#include <premia/runtime/core.h>

void translate(premia::api::exception const & e)
{
	// Use the Python 'C' API to set up an exception object
	PyErr_SetString(PyExc_RuntimeError, e.what());
}

BOOST_PYTHON_MODULE(libpypremia)
{
  using namespace boost::python;
  using namespace premia::api;

  register_exception_translator<premia::api::exception>(&translate);

  def("init_premia", init_premia);
  def("setCurrentAsset", setCurrentAsset);
  def("setCurrentModel", setCurrentModel);
  def("setCurrentMethod", setCurrentMethod);
  def("setCurrentOption", setCurrentOption);
  def("write_double", write_double);
  def("write_long", write_long);
  def("write_int", write_int);
  def("write_enum", write_enum);
  def("write_array", write_array);
  def("write_array_double", write_array_double);
  def("write_array_end", write_array_end);
  def("write_filename", write_filename);
  def("stopWriteParameters", stopWriteParameters);

  def("ignore_double", ignore_double);
  def("ignore_long", ignore_long);
  def("ignore_int", ignore_int);
  def("ignore_enum", ignore_enum);
  def("ignore_array", ignore_array);
  def("ignore_array_double", ignore_array_double);
  def("ignore_array_end", ignore_array_end);
  def("ignore_filename", ignore_filename);

  def("readCurrentModel", readCurrentModel);
  def("readCurrentMethod", readCurrentMethod);
  def("readCurrentOption", readCurrentOption);
  def("read_double", read_double);
  def("read_long", read_long);
  def("read_int", read_int);
  def("read_enum", read_enum);
  def("read_array_size", read_array_size);
  def("read_array_double", read_array_double);
  def("read_array_end", read_array_end);
  def("read_filename", read_filename);
  def("stopReadParameters", stopWriteParameters);

  def("reset", reset);
  def("getResultDouble", get_result_double);
  def("getResultBool", get_result_bool);
  def("getResultArraySize", get_result_array_size);
  def("getResultArrayItem", get_result_array_item);
  def("compute_3", compute);
  def("init", init_premia);
} 
