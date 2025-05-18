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

// Interpreter is implemented as a singleton (NOT thread-safe)
class Interpreter
{
public:
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

template<typename Numeric> PyObject* to_list(const std::vector<Numeric>& v)
{
  PyObject* list = PyList_New(static_cast<Py_ssize_t>(v.size()));
  if (!list) throw std::runtime_error("Couldn't convert vector to list!");

  for (size_t i = 0; i < v.size(); ++i) {
    PyList_SetItem(list, static_cast<Py_ssize_t>(i), PyFloat_FromDouble(v.at(i)));
  }
  return list;
}

template<typename Numeric> bool plot(const std::vector<Numeric>& x, const std::vector<Numeric>& y)
{
  assert(x.size() == y.size());

  Interpreter::getInstance();

  PyObject* xList = to_list(x);
  PyObject* yList = to_list(y);

  // Construct positional arguments
  PyObject* args = PyTuple_New(2);
  PyTuple_SetItem(args, 0, xList);
  PyTuple_SetItem(args, 1, yList);

  PyObject* res = PyObject_Call(Interpreter::getInstance()->mPythonFuncPlot, args, NULL);

  Py_DECREF(args);
  if (res) Py_DECREF(res);

  return res;
}

template<typename Numeric>
bool plot(const std::vector<Numeric>& x,
  const std::vector<Numeric>& y,
  const std::map<std::string, std::string>& keywords)
{
  assert(x.size() == y.size());

  Interpreter::getInstance();

  PyObject* xList = to_list(x);
  PyObject* yList = to_list(y);

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