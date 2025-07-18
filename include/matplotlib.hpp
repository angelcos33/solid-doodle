#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <numpy/arrayobject.h>

#include <matrix.hpp>

// Interpreter is implemented as a singleton (NOT thread-safe)
class Interpreter
{
public:
  PyObject* mPythonFuncContourf;
  PyObject* mPythonFuncFigure;
  PyObject* mPythonFuncImshow;
  PyObject* mPythonFuncPlot;
  PyObject* mPythonFuncShow;
  PyObject* mPythonFuncSave;

  Interpreter(Interpreter& other) = delete;
  void operator=(const Interpreter&) = delete;

  ~Interpreter() { Py_Finalize(); }

  static Interpreter* getInstance();

  PyObject* safe_import(PyObject* module, std::string fname)
  {
    PyObject* fn = PyObject_GetAttrString(module, fname.c_str());

    if (!fn) throw std::runtime_error(std::string("Couldn't find required function: ") + fname);

    if (!PyFunction_Check(fn)) throw std::runtime_error(fname + std::string(" is unexpectedly not a PyFunction."));

    return fn;
  }

private:
  static Interpreter* pInstance;

  Interpreter()
  {
    Py_Initialize();

    import_numpy(); // Initialize numpy C-API

    PyObject* matplotlib = PyImport_ImportModule("matplotlib");
    if (!matplotlib) {
      PyErr_Print();
      throw std::runtime_error("Error loading module matplotlib!");
    }
    PyObject* pyplot = PyImport_ImportModule("matplotlib.pyplot");
    if (!pyplot) {
      PyErr_Print();
      throw std::runtime_error("Error loading module matplotlib.pyplot!");
    }

    mPythonFuncPlot = safe_import(pyplot, "plot");
    mPythonFuncShow = safe_import(pyplot, "show");
    mPythonFuncSave = safe_import(pyplot, "savefig");
    mPythonFuncContourf = safe_import(pyplot, "contourf");
    mPythonFuncFigure = safe_import(pyplot, "figure");
    mPythonFuncImshow = safe_import(pyplot, "imshow");
  }

  void* import_numpy()
  {
    import_array(); // Initialize C-API of NumPy required for Python3
    return NULL;
  }
};

Interpreter* Interpreter::pInstance = nullptr;

Interpreter* Interpreter::getInstance()
{
  if (pInstance == nullptr) pInstance = new Interpreter();

  return pInstance;
}

// Type selector for Numpy types
template<typename T> struct select_np_type
{
  const static NPY_TYPES type = NPY_NOTYPE;
}; // Default
template<> struct select_np_type<double>
{
  const static NPY_TYPES type = NPY_DOUBLE;
};
template<> struct select_np_type<float>
{
  const static NPY_TYPES type = NPY_FLOAT;
};
template<> struct select_np_type<int>
{
  const static NPY_TYPES type = NPY_INT;
};

template<typename T> PyObject* to_array(const std::vector<T>& v)
{
  NPY_TYPES type = select_np_type<T>::type;
  const npy_intp dims = static_cast<npy_intp>(v.size());
  PyObject* array = PyArray_SimpleNewFromData(1, &dims, type, static_cast<void*>(const_cast<T*>(v.data())));
  return array;
}

template<typename T> PyObject* to_2d_array(const Matrix<T>& m)
{
  NPY_TYPES type = select_np_type<T>::type;
  const npy_intp dims[2] = {static_cast<npy_intp>(m.rows()), static_cast<npy_intp>(m.cols())};
  PyObject* array = PyArray_SimpleNewFromData(2, dims, type, static_cast<void*>(const_cast<T*>(m.data())));
  return array;
}

inline void figure()
{
  Interpreter::getInstance();

  PyObject* res = PyObject_CallNoArgs(Interpreter::getInstance()->mPythonFuncFigure);

  if (!res) {
    PyErr_Print();
    throw std::runtime_error("Call to figure() failed");
  }
  Py_DECREF(res);
}

template<typename T> bool plot(const std::vector<T>& x, const std::vector<T>& y)
{
  assert(x.size() == y.size());

  Interpreter::getInstance();

  PyObject* xList = to_array<T>(x);
  PyObject* yList = to_array<T>(y);

  // Construct positional arguments
  PyObject* args = PyTuple_New(2);
  PyTuple_SetItem(args, 0, xList);
  PyTuple_SetItem(args, 1, yList);

  PyObject* res = PyObject_Call(Interpreter::getInstance()->mPythonFuncPlot, args, NULL);

  Py_DECREF(args);
  if (res) Py_DECREF(res);

  return res;
}

template<typename T>
bool plot(const std::vector<T>& x, const std::vector<T>& y, const std::map<std::string, std::string>& keywords)
{
  assert(x.size() == y.size());

  Interpreter::getInstance();

  PyObject* xList = to_array<T>(x);
  PyObject* yList = to_array<T>(y);

  // Construct positional arguments
  PyObject* args = PyTuple_New(2);
  PyTuple_SetItem(args, 0, xList);
  PyTuple_SetItem(args, 1, yList);

  // Construct kwargs
  PyObject* kwargs = PyDict_New();
  for (const auto& [key, val] : keywords) {
    PyDict_SetItemString(kwargs, key.c_str(), PyUnicode_FromString(val.c_str()));
  }

  PyObject* res = PyObject_Call(Interpreter::getInstance()->mPythonFuncPlot, args, kwargs);

  Py_DECREF(args);
  Py_DECREF(kwargs);
  if (res) Py_DECREF(res);

  return res;
}

template<typename T> void contourf(const std::vector<T>& X, const std::vector<T>& Y, const Matrix<T>& Z)
{
  m_assert(X.size() == Z.cols() && Y.size() == Z.rows(),
    "Invalid shape for input data. X.size()==Z.cols() and Y.size()==Z.rows()");

  Interpreter::getInstance();

  PyObject* XList = to_array<T>(X);
  PyObject* YList = to_array<T>(Y);

  PyObject* ZNp = to_2d_array<T>(Z);

  // Construct positional arguments
  PyObject* args = PyTuple_New(3);
  PyTuple_SetItem(args, 0, XList);
  PyTuple_SetItem(args, 1, YList);
  PyTuple_SetItem(args, 2, ZNp);

  PyObject* res = PyObject_Call(Interpreter::getInstance()->mPythonFuncContourf, args, NULL);
  Py_DECREF(args);

  if (!res) {
    PyErr_Print();
    throw std::runtime_error("Call to contourf() failed");
  }
  Py_DECREF(res);
}

template<typename T> void imshow(const Matrix<T>& Z)
{
  Interpreter::getInstance();

  PyObject* ZNp = to_2d_array<T>(Z);

  PyObject* args = PyTuple_New(1);
  PyTuple_SetItem(args, 0, ZNp);

  PyObject* res = PyObject_Call(Interpreter::getInstance()->mPythonFuncImshow, args, NULL);
  Py_DECREF(args);

  if (!res) {
    PyErr_Print();
    throw std::runtime_error("Call to imshow() failed");
  }
  Py_DECREF(res);
}

inline void show()
{
  Interpreter::getInstance();

  PyObject* res = PyObject_CallNoArgs(Interpreter::getInstance()->mPythonFuncShow);

  if (!res) throw std::runtime_error("Call to show() failed");
  Py_DECREF(res);
}

inline void save(const std::string& filepath)
{
  Interpreter::getInstance();

  // Construct positional args
  PyObject* args = PyTuple_New(1);
  PyTuple_SetItem(args, 0, PyUnicode_FromString(filepath.c_str()));

  PyObject* res = PyObject_Call(Interpreter::getInstance()->mPythonFuncSave, args, NULL);

  Py_DECREF(args);
  if (!res) {
    PyErr_Print();
    throw std::runtime_error("Call to savefig() failed");
  }
  Py_DECREF(res);
}