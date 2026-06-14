#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <numpy/arrayobject.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <vector>

// acceleration(dx, dy, dz, num, pos_x, pos_y, pos_z, i, masses, eps, G)

static PyObject* acceleration(PyObject* self, PyObject* args)
{
    double dx, dy, dz;
    int num;
    double eps, G;

    PyObject *pos_x_obj;
    PyObject *pos_y_obj;
    PyObject *pos_z_obj;
    PyObject *mass_obj;

    if (!PyArg_ParseTuple(args,
        "dddiOOOOdd",
        &dx, &dy, &dz,
        &num,
        &pos_x_obj,
        &pos_y_obj,
        &pos_z_obj,
        &mass_obj,
        &eps,
        &G))
    {
        return NULL;
    }

    Py_ssize_t n = PyList_Size(pos_x_obj);

    double ax = 0.0, ay = 0.0, az = 0.0;

    for (Py_ssize_t j = 0; j < n; j++)
    {
        if (j == num) continue;

        double rx = PyFloat_AsDouble(PyList_GetItem(pos_x_obj, j)) - dx;
        double ry = PyFloat_AsDouble(PyList_GetItem(pos_y_obj, j)) - dy;
        double rz = PyFloat_AsDouble(PyList_GetItem(pos_z_obj, j)) - dz;

        double mj = PyFloat_AsDouble(PyList_GetItem(mass_obj, j));

        double r2 = rx*rx + ry*ry + rz*rz + eps*eps;
        double r_soft = sqrt(r2);
        double inv_r3 = 1.0 / (r_soft * r_soft * r_soft);

        ax += G * mj * rx * inv_r3;
        ay += G * mj * ry * inv_r3;
        az += G * mj * rz * inv_r3;
    }

    return Py_BuildValue("ddd", ax, ay, az);
}

static PyObject* Yoshida4thIntegrator(PyObject* self, PyObject* args) 
{
    
    PyObject *mass_obj;
    double eps, G;

    double w0 = -1.702414289193;
    double w1 = 1.3512071919597;

    double d_coeff[3] = 
        {
            (w1) , 
            (w0) , 
            (w1)
        };

    double c_coeff[4] = 
        {
            (w1/2.0) , 
            ((w0+w1)/2.0) , 
            ((w0+w1)/2.0) , 
            (w1/2.0)
        };

    PyObject* pos_x; 
    PyObject* pos_y; 
    PyObject* pos_z; 

    PyObject* vel_x; 
    PyObject* vel_y; 
    PyObject* vel_z;

    double dt_seconds;

    if (!PyArg_ParseTuple(args,
        "OOOOOOOddd",
        &pos_x, &pos_y, &pos_z,
        &vel_x, &vel_y, &vel_z,
        &mass_obj,
        &eps,
        &G, 
        &dt_seconds))
    {
        return NULL;
    }

    PyObject* seq_px = PySequence_Fast(pos_x, "pos_x must be a sequence");
    PyObject* seq_py = PySequence_Fast(pos_y, "pos_y must be a sequence");
    PyObject* seq_pz = PySequence_Fast(pos_z, "pos_z must be a sequence");
    PyObject* seq_vx = PySequence_Fast(vel_x, "vel_x must be a sequence");
    PyObject* seq_vy = PySequence_Fast(vel_y, "vel_y must be a sequence");
    PyObject* seq_vz = PySequence_Fast(vel_z, "vel_z must be a sequence");
    PyObject* seq_masses = PySequence_Fast(mass_obj, "masses must be a sequence");

    if (!seq_px || !seq_py || !seq_pz || !seq_vx || !seq_vy || !seq_vz || !seq_masses)
    {
        Py_XDECREF(seq_px);
        Py_XDECREF(seq_py);
        Py_XDECREF(seq_pz);
        Py_XDECREF(seq_vx);
        Py_XDECREF(seq_vy);
        Py_XDECREF(seq_vz);
        Py_XDECREF(seq_masses);
        return NULL;
    }

    Py_ssize_t n = PySequence_Fast_GET_SIZE(seq_px);
    if (
        PySequence_Fast_GET_SIZE(seq_py) != n ||
        PySequence_Fast_GET_SIZE(seq_pz) != n ||
        PySequence_Fast_GET_SIZE(seq_vx) != n ||
        PySequence_Fast_GET_SIZE(seq_vy) != n ||
        PySequence_Fast_GET_SIZE(seq_vz) != n ||
        PySequence_Fast_GET_SIZE(seq_masses) != n)
    {
        Py_DECREF(seq_px);
        Py_DECREF(seq_py);
        Py_DECREF(seq_pz);
        Py_DECREF(seq_vx);
        Py_DECREF(seq_vy);
        Py_DECREF(seq_vz);
        Py_DECREF(seq_masses);
        PyErr_SetString(PyExc_ValueError, "all state arrays must have the same length");
        return NULL;
    }

    std::vector<double> px(n), py(n), pz(n), vx(n), vy(n), vz(n), masses(n);

    for (Py_ssize_t i = 0; i < n; i++)
    {
        px[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq_px, i));
        py[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq_py, i));
        pz[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq_pz, i));
        vx[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq_vx, i));
        vy[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq_vy, i));
        vz[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq_vz, i));
        masses[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq_masses, i));
        if (PyErr_Occurred())
        {
            Py_DECREF(seq_px);
            Py_DECREF(seq_py);
            Py_DECREF(seq_pz);
            Py_DECREF(seq_vx);
            Py_DECREF(seq_vy);
            Py_DECREF(seq_vz);
            Py_DECREF(seq_masses);
            return NULL;
        }
    }

    double* axs = (double*) malloc(n * sizeof(double));
    double* ays = (double*) malloc(n * sizeof(double));
    double* azs = (double*) malloc(n * sizeof(double));

    for (int num = 0; num < 3; num++)
    {
        double c = c_coeff[num];
        double d = d_coeff[num];

        for (Py_ssize_t j = 0; j < n; j++)
        {
            px[j] += c * vx[j] * dt_seconds;
            py[j] += c * vy[j] * dt_seconds;
            pz[j] += c * vz[j] * dt_seconds;
        }

        for (Py_ssize_t j = 0; j < n; j++)
        {
            axs[j] = 0;
            ays[j] = 0;
            azs[j] = 0;

            for (Py_ssize_t i = 0; i < n; i++)
            {
                if (i == j) continue;

                double rx = px[i] - px[j];
                double ry = py[i] - py[j];
                double rz = pz[i] - pz[j];

                double r2 = rx*rx + ry*ry + rz*rz + eps*eps;
                double r_soft = sqrt(r2);
                if (r2 == 0) continue;
                double inv_r3 = 1.0 / (r_soft * r_soft * r_soft);

                axs[j] += G * masses[i] * rx * inv_r3;
                ays[j] += G * masses[i] * ry * inv_r3;
                azs[j] += G * masses[i] * rz * inv_r3;
            }
        }

        for (Py_ssize_t j = 0; j < n; j++)
        {
            vx[j] += d * dt_seconds * axs[j];
            vy[j] += d * dt_seconds * ays[j];
            vz[j] += d * dt_seconds * azs[j];
        }
    }

    for (Py_ssize_t j = 0; j < n; j++)
    {

        px[j] += c_coeff[3] * vx[j] * dt_seconds;
        py[j] += c_coeff[3] * vy[j] * dt_seconds;
        pz[j] += c_coeff[3] * vz[j] * dt_seconds;
    }

    free(axs);
    free(ays);
    free(azs);

    PyObject* out_px = PyList_New(n);
    PyObject* out_py = PyList_New(n);
    PyObject* out_pz = PyList_New(n);
    PyObject* out_vx = PyList_New(n);
    PyObject* out_vy = PyList_New(n);
    PyObject* out_vz = PyList_New(n);

    if (!out_px || !out_py || !out_pz || !out_vx || !out_vy || !out_vz)
    {
        Py_XDECREF(out_px);
        Py_XDECREF(out_py);
        Py_XDECREF(out_pz);
        Py_XDECREF(out_vx);
        Py_XDECREF(out_vy);
        Py_XDECREF(out_vz);
        Py_DECREF(seq_px);
        Py_DECREF(seq_py);
        Py_DECREF(seq_pz);
        Py_DECREF(seq_vx);
        Py_DECREF(seq_vy);
        Py_DECREF(seq_vz);
        Py_DECREF(seq_masses);
        return NULL;
    }

    for (Py_ssize_t i = 0; i < n; i++)
    {
        PyList_SET_ITEM(out_px, i, PyFloat_FromDouble(px[i]));
        PyList_SET_ITEM(out_py, i, PyFloat_FromDouble(py[i]));
        PyList_SET_ITEM(out_pz, i, PyFloat_FromDouble(pz[i]));
        PyList_SET_ITEM(out_vx, i, PyFloat_FromDouble(vx[i]));
        PyList_SET_ITEM(out_vy, i, PyFloat_FromDouble(vy[i]));
        PyList_SET_ITEM(out_vz, i, PyFloat_FromDouble(vz[i]));
    }

    Py_DECREF(seq_px);
    Py_DECREF(seq_py);
    Py_DECREF(seq_pz);
    Py_DECREF(seq_vx);
    Py_DECREF(seq_vy);
    Py_DECREF(seq_vz);
    Py_DECREF(seq_masses);

    return Py_BuildValue("NNNNNN", out_px, out_py, out_pz, out_vx, out_vy, out_vz);
}


static PyMethodDef module_methods[] = {
    {"acceleration", acceleration, METH_VARARGS, "Compute acceleration"},
    {"Yoshida4thIntegrator", Yoshida4thIntegrator, METH_VARARGS, "Compute velocities and positions"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef core_module = {
    PyModuleDef_HEAD_INIT,
    "core_calculations",
    "Core calculations module",
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit_core_calculations(void)
{
    return PyModule_Create(&core_module);
}
