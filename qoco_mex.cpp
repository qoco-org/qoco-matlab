/* Copyright (c) 2025, Govind M. Chari <govindchari1@gmail.com> */
/* Copyright (c) 2017, Bartolomeo Stellato */
/* This source code is licensed under the BSD 3-Clause License  */

#include "qoco_mex.hpp"
#include <iostream>
const char *QOCO_SETTINGS_FIELDS[] = {"max_iters",
                                      "bisect_iters",
                                      "ruiz_iters",
                                      "iter_ref_iters",
                                      "kkt_static_reg",
                                      "kkt_dynamic_reg",
                                      "abstol",
                                      "reltol",
                                      "abstol_inacc",
                                      "reltol_inacc",
                                      "verbose"};

const char *QOCO_RESULT_FIELDS[] = {"x",
                                    "z",
                                    "y",
                                    "z",
                                    "iters",
                                    "setup_time_sec",
                                    "solve_time_sec",
                                    "obj",
                                    "pres",
                                    "dres",
                                    "gap",
                                    "status"};

const char *CSC_FIELDS[] = {"m",
                            "n",
                            "nnz",
                            "i",
                            "p",
                            "x"};

// wrapper class for all qoco data and settings
class QocoData
{
public:
  QocoData()
  {
    solver = NULL;
  }
  QOCOSolver *solver; // Solver
};

// internal utility functions
void castToDoubleArr(QOCOFloat *arr, double *arr_out, QOCOInt len);
void setToNaN(double *arr_out, QOCOInt len);
void copyMxStructToSettings(const mxArray *, QOCOSettings *);
void copyUpdatedSettingsToWork(const mxArray *, QOCOSettings *);
void castCintToDoubleArr(QOCOInt *arr, double *arr_out, QOCOInt len);
QOCOInt *copyToCintVector(mwIndex *vecData, QOCOInt numel);
QOCOInt *copyDoubleToCintVector(double *vecData, QOCOInt numel);
QOCOInt *roundDoubleIntoInt(double *vecData, QOCOInt numel);
QOCOFloat *copyToCfloatVector(double *vecData, QOCOInt numel);
// mxArray *copyInfoToMxStruct(OSQPInfo *info);
mxArray *copySettingsToMxStruct(QOCOSettings *settings);
// mxArray *copyCscMatrixToMxStruct(csc *M);
// mxArray *copyDataToMxStruct(OSQPWorkspace *work);
// mxArray *copyLinsysSolverToMxStruct(OSQPWorkspace *work);
// mxArray *copyScalingToMxStruct(OSQPWorkspace *work);
// mxArray *copyWorkToMxStruct(OSQPWorkspace *work);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  QocoData *qocoData; // QOCO data identifier

  // Exitflag
  QOCOInt exitflag = 1;

  // Static string for static methods
  char stat_string[64];

  // Get the command string
  char cmd[64];
  if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
    mexErrMsgTxt("First input should be a command string less than 64 characters long.");

  // new object
  if (!strcmp("new", cmd))
  {
    // Check parameters
    if (nlhs != 1)
    {
      mexErrMsgTxt("New: One output expected.");
    }
    // Return a handle to a new C++ wrapper instance
    plhs[0] = convertPtr2Mat<QocoData>(new QocoData);
    return;
  }

  // Check for a second input, which should be the class instance handle or string 'static'
  if (nrhs < 2)
    mexErrMsgTxt("Second input should be a class instance handle or the string 'static'.");

  if (mxGetString(prhs[1], stat_string, sizeof(stat_string)))
  {
    // If we are dealing with non-static methods, get the class instance pointer from the second input
    qocoData = convertMat2Ptr<QocoData>(prhs[1]);
  }
  else
  {
    if (strcmp("static", stat_string))
    {
      mexErrMsgTxt("Second argument for static functions is string 'static'");
    }
  }

  // delete the object and its data
  if (!strcmp("delete", cmd))
  {

    // clean up the problem workspace
    if (qocoData->solver)
    {
      qoco_cleanup(qocoData->solver);
    }

    // clean up the handle object
    destroyObject<QocoData>(prhs[1]);
    // Warn if other commands were ignored
    if (nlhs != 0 || nrhs != 2)
      mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
    return;
  }

  // report the current settings
  if (!strcmp("current_settings", cmd))
  {

    // throw an error if this is called before solver is configured
    if (!qocoData->solver)
    {
      mexErrMsgTxt("Solver is uninitialized.  No settings have been configured.");
    }
    // report the current settings
    plhs[0] = copySettingsToMxStruct(qocoData->solver->settings);
    return;
  }

  // write_settings
  if (!strcmp("update_settings", cmd))
  {

    // overwrite the current settings.  Mex function is responsible
    // for disallowing overwrite of selected settings after initialization,
    // and for all error checking
    // throw an error if this is called before solver is configured
    if (!qocoData->solver)
    {
      mexErrMsgTxt("Solver is uninitialized.  No settings have been configured.");
    }

    copyUpdatedSettingsToWork(prhs[2], qocoData->solver->settings);
    return;
  }

  // report the default settings
  if (!strcmp("default_settings", cmd))
  {

    // Warn if other commands were ignored
    if (nrhs > 2)
      mexWarnMsgTxt("Default settings: unexpected number of arguments.");

    // Create a Settings structure in default form and report the results
    // Useful for external solver packages (e.g. Yalmip) that want to
    // know which solver settings are supported
    QOCOSettings *defaults = (QOCOSettings *)mxCalloc(1, sizeof(QOCOSettings));
    set_default_settings(defaults);
    plhs[0] = copySettingsToMxStruct(defaults);
    mxFree(defaults);
    return;
  }

  // setup
  if (!strcmp("setup", cmd))
  {

    // throw an error if this is called more than once
    if (qocoData->solver)
    {
      mexErrMsgTxt("Solver is already initialized with problem data.");
    }

    // Create data and settings containers
    QOCOSettings *settings = new QOCOSettings;

    // handle the problem data first.  Matlab-side
    // class wrapper is responsible for ensuring that
    // P and A are sparse matrices,  everything
    // else is a dense vector and all inputs are
    // of compatible dimension

    const QOCOInt n = (QOCOInt)mxGetScalar(prhs[2]);
    const QOCOInt m = (QOCOInt)mxGetScalar(prhs[3]);
    const QOCOInt p = (QOCOInt)mxGetScalar(prhs[4]);
    const mxArray *P = prhs[5];
    const mxArray *c = prhs[6];
    const mxArray *A = prhs[7];
    const mxArray *b = prhs[8];
    const mxArray *G = prhs[9];
    const mxArray *h = prhs[10];
    const QOCOInt l = (QOCOInt)mxGetScalar(prhs[11]);
    const QOCOInt nsoc = (QOCOInt)mxGetScalar(prhs[12]);
    const mxArray *q = prhs[13];

    // Create Data Structure
    QOCOFloat *cvec = copyToCfloatVector(mxGetPr(c), n);
    QOCOFloat *bvec = copyToCfloatVector(mxGetPr(b), p);
    QOCOFloat *hvec = copyToCfloatVector(mxGetPr(h), m);
    QOCOInt *qvec = roundDoubleIntoInt(mxGetPr(q), nsoc); // Warning: This only works with QOCOInt is int32

    // Matrix P:  nnz = P->p[n]
    QOCOCscMatrix *Pcsc = new QOCOCscMatrix;
    QOCOInt *Pp = copyToCintVector(mxGetJc(P), n + 1);
    QOCOInt *Pi = copyToCintVector(mxGetIr(P), Pp[n]);
    QOCOFloat *Px = copyToCfloatVector(mxGetPr(P), Pp[n]);
    qoco_set_csc(Pcsc, n, n, Pp[n], Px, Pp, Pi);

    // // Matrix A: nnz = A->p[n]
    QOCOCscMatrix *Acsc = new QOCOCscMatrix;
    QOCOInt *Ap = copyToCintVector(mxGetJc(A), n + 1);
    QOCOInt *Ai = copyToCintVector(mxGetIr(A), Ap[n]);
    QOCOFloat *Ax = copyToCfloatVector(mxGetPr(A), Ap[n]);
    qoco_set_csc(Acsc, p, n, Ap[n], Ax, Ap, Ai);

    // Matrix G: nnz = G->p[n]
    QOCOCscMatrix *Gcsc = new QOCOCscMatrix;
    QOCOInt *Gp = copyToCintVector(mxGetJc(G), n + 1);
    QOCOInt *Gi = copyToCintVector(mxGetIr(G), Gp[n]);
    QOCOFloat *Gx = copyToCfloatVector(mxGetPr(G), Gp[n]);
    qoco_set_csc(Gcsc, m, n, Gp[n], Gx, Gp, Gi);

    // Create Settings
    const mxArray *mxSettings = prhs[14];
    if (mxIsEmpty(mxSettings))
    {
      // use defaults
      set_default_settings(settings);
    }
    else
    {
      // populate settings structure from mxArray input
      copyMxStructToSettings(mxSettings, settings);
    }

    // Setup solver.
    qocoData->solver = new QOCOSolver;
    exitflag = qoco_setup(qocoData->solver, n, m, p, Pcsc, cvec, Acsc, bvec, Gcsc, hvec, l, nsoc, qvec, settings);

    // cleanup temporary structures
    if (cvec)
      qoco_free(cvec);
    if (bvec)
      qoco_free(bvec);
    if (hvec)
      qoco_free(hvec);
    // if (qvec)
    //   qoco_free(qvec);
    if (Px)
      qoco_free(Px);
    if (Pi)
      qoco_free(Pi);
    if (Pp)
      qoco_free(Pp);
    if (Pcsc)
      delete Pcsc;
    if (Ax)
      qoco_free(Ax);
    if (Ai)
      qoco_free(Ai);
    if (Ap)
      qoco_free(Ap);
    if (Acsc)
      delete Acsc;
    if (Gx)
      qoco_free(Gx);
    if (Gi)
      qoco_free(Gi);
    if (Gp)
      qoco_free(Gp);
    if (Gcsc)
      delete Gcsc;

    // Settings
    if (settings)
      delete settings;

    // Report error (if any)
    if (exitflag != 0)
    {
      mexErrMsgTxt("Invalid problem setup");
    }
    return;
  }

  if (!strcmp("get_dimensions", cmd))
  {

    // throw an error if this is called before solver is configured
    if (!qocoData->solver)
    {
      mexErrMsgTxt("Solver has not been initialized.");
    }

    plhs[0] = mxCreateDoubleScalar(qocoData->solver->work->data->m);
    plhs[1] = mxCreateDoubleScalar(qocoData->solver->work->data->n);
    plhs[2] = mxCreateDoubleScalar(qocoData->solver->work->data->p);

    return;
  }

  // solve
  if (!strcmp("solve", cmd))
  {
    if (nlhs != 12 || nrhs != 2)
    {
      mexErrMsgTxt("Solve : wrong number of inputs / outputs");
    }
    if (!qocoData->solver)
    {
      mexErrMsgTxt("No problem data has been given.");
    }

    // solve the problem
    qoco_solve(qocoData->solver);

    // Allocate space for solution
    // primal variables
    plhs[0] = mxCreateDoubleMatrix(qocoData->solver->work->data->n, 1, mxREAL);
    // slack variables
    plhs[1] = mxCreateDoubleMatrix(qocoData->solver->work->data->m, 1, mxREAL);
    // dual variables for equality constraints
    plhs[2] = mxCreateDoubleMatrix(qocoData->solver->work->data->p, 1, mxREAL);
    // dual variables for conic constraints
    plhs[3] = mxCreateDoubleMatrix(qocoData->solver->work->data->m, 1, mxREAL);

    // copy results to mxArray outputs
    castToDoubleArr(qocoData->solver->sol->x, mxGetPr(plhs[0]), qocoData->solver->work->data->n);
    castToDoubleArr(qocoData->solver->sol->s, mxGetPr(plhs[1]), qocoData->solver->work->data->m);
    castToDoubleArr(qocoData->solver->sol->y, mxGetPr(plhs[2]), qocoData->solver->work->data->p);
    castToDoubleArr(qocoData->solver->sol->z, mxGetPr(plhs[3]), qocoData->solver->work->data->m);

    // copy the rest of the information.
    plhs[4] = mxCreateDoubleScalar(qocoData->solver->sol->iters);
    plhs[5] = mxCreateDoubleScalar(qocoData->solver->sol->setup_time_sec);
    plhs[6] = mxCreateDoubleScalar(qocoData->solver->sol->solve_time_sec);
    plhs[7] = mxCreateDoubleScalar(qocoData->solver->sol->obj);
    plhs[8] = mxCreateDoubleScalar(qocoData->solver->sol->pres);
    plhs[9] = mxCreateDoubleScalar(qocoData->solver->sol->dres);
    plhs[10] = mxCreateDoubleScalar(qocoData->solver->sol->gap);
    plhs[11] = mxCreateDoubleScalar(qocoData->solver->sol->status);

    return;
  }

  if (!strcmp("constant", cmd))
  { // Return solver constants

    char constant[32];
    mxGetString(prhs[2], constant, sizeof(constant));

    if (!strcmp("QOCO_UNSOLVED", constant))
    {
      plhs[0] = mxCreateDoubleScalar(QOCO_UNSOLVED);
      return;
    }

    if (!strcmp("QOCO_SOLVED", constant))
    {
      plhs[0] = mxCreateDoubleScalar(QOCO_SOLVED);
      return;
    }

    if (!strcmp("QOCO_SOLVED_INACCURATE", constant))
    {
      plhs[0] = mxCreateDoubleScalar(QOCO_SOLVED_INACCURATE);
      return;
    }

    if (!strcmp("QOCO_NUMERICAL_ERROR", constant))
    {
      plhs[0] = mxCreateDoubleScalar(QOCO_NUMERICAL_ERROR);
      return;
    }

    if (!strcmp("QOCO_MAX_ITER", constant))
    {
      plhs[0] = mxCreateDoubleScalar(QOCO_MAX_ITER);
      return;
    }

    mexErrMsgTxt("Constant not recognized.");

    return;
  }

  // Got here, so command not recognized
  mexErrMsgTxt("Command not recognized.");
}

QOCOFloat *copyToCfloatVector(double *vecData, QOCOInt numel)
{
  // This memory needs to be freed!
  QOCOFloat *out = (QOCOFloat *)qoco_malloc(numel * sizeof(QOCOFloat));

  // copy data
  for (QOCOInt i = 0; i < numel; i++)
  {
    out[i] = (QOCOFloat)vecData[i];
  }
  return out;
}

QOCOInt *roundDoubleIntoInt(double *vecData, QOCOInt numel)
{
  // This memory needs to be freed!
  QOCOInt *out = (QOCOInt *)qoco_malloc(numel * sizeof(QOCOInt));

  // round data
  for (QOCOInt i = 0; i < numel; i++)
  {
    out[i] = round(vecData[i]);
  }
  return out;
}

QOCOInt *copyToCintVector(mwIndex *vecData, QOCOInt numel)
{
  // This memory needs to be freed!
  QOCOInt *out = (QOCOInt *)qoco_malloc(numel * sizeof(QOCOInt));

  // copy data
  for (QOCOInt i = 0; i < numel; i++)
  {
    out[i] = (QOCOInt)vecData[i];
  }
  return out;
}

QOCOInt *copyDoubleToCintVector(double *vecData, QOCOInt numel)
{
  // This memory needs to be freed!
  QOCOInt *out = (QOCOInt *)qoco_malloc(numel * sizeof(QOCOInt));

  // copy data
  for (QOCOInt i = 0; i < numel; i++)
  {
    out[i] = (QOCOInt)vecData[i];
  }
  return out;
}

void castCintToDoubleArr(QOCOInt *arr, double *arr_out, QOCOInt len)
{
  for (QOCOInt i = 0; i < len; i++)
  {
    arr_out[i] = (double)arr[i];
  }
}

void castToDoubleArr(QOCOFloat *arr, double *arr_out, QOCOInt len)
{
  for (QOCOInt i = 0; i < len; i++)
  {
    arr_out[i] = (double)arr[i];
  }
}

void setToNaN(double *arr_out, QOCOInt len)
{
  QOCOInt i;
  for (i = 0; i < len; i++)
  {
    arr_out[i] = mxGetNaN();
  }
}

mxArray *copySettingsToMxStruct(QOCOSettings *settings)
{

  int nfields = sizeof(QOCO_SETTINGS_FIELDS) / sizeof(QOCO_SETTINGS_FIELDS[0]);
  mxArray *mxPtr = mxCreateStructMatrix(1, 1, nfields, QOCO_SETTINGS_FIELDS);

  // map the QOCO_SETTINGS fields one at a time into mxArrays
  // matlab handles everything as a double
  mxSetField(mxPtr, 0, "max_iters", mxCreateDoubleScalar(settings->max_iters));
  mxSetField(mxPtr, 0, "bisect_iters", mxCreateDoubleScalar(settings->bisect_iters));
  mxSetField(mxPtr, 0, "ruiz_iters", mxCreateDoubleScalar(settings->ruiz_iters));
  mxSetField(mxPtr, 0, "iter_ref_iters", mxCreateDoubleScalar(settings->iter_ref_iters));
  mxSetField(mxPtr, 0, "kkt_static_reg", mxCreateDoubleScalar(settings->kkt_static_reg));
  mxSetField(mxPtr, 0, "kkt_dynamic_reg", mxCreateDoubleScalar(settings->kkt_dynamic_reg));
  mxSetField(mxPtr, 0, "abstol", mxCreateDoubleScalar(settings->abstol));
  mxSetField(mxPtr, 0, "reltol", mxCreateDoubleScalar(settings->reltol));
  mxSetField(mxPtr, 0, "abstol_inacc", mxCreateDoubleScalar(settings->abstol_inacc));
  mxSetField(mxPtr, 0, "reltol_inacc", mxCreateDoubleScalar(settings->reltol_inacc));
  mxSetField(mxPtr, 0, "verbose", mxCreateDoubleScalar(settings->verbose));

  return mxPtr;
}

// ======================================================================

void copyMxStructToSettings(const mxArray *mxPtr, QOCOSettings *settings)
{

  // this function assumes that only a complete and validated structure
  // will be passed.  matlab mex-side function is responsible for checking
  // structure validity

  // map the QOCO_SETTINGS fields one at a time into mxArrays
  // matlab handles everything as a double
  settings->max_iters = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "max_iters"));
  settings->bisect_iters = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "bisect_iters"));
  settings->ruiz_iters = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "ruiz_iters"));
  settings->iter_ref_iters = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "iter_ref_iters"));
  settings->kkt_static_reg = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "kkt_static_reg"));
  settings->kkt_dynamic_reg = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "kkt_dynamic_reg"));
  settings->abstol = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "abstol"));
  settings->reltol = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "reltol"));
  settings->abstol_inacc = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "abstol_inacc"));
  settings->reltol_inacc = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "reltol_inacc"));
  settings->verbose = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "verbose"));
}

void copyUpdatedSettingsToWork(const mxArray *mxPtr, QOCOSettings *settings)
{
  settings->max_iters = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "max_iters"));
  settings->bisect_iters = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "sigma"));
  settings->ruiz_iters = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "scaling"));
  settings->iter_ref_iters = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "adaptive_rho"));
  settings->kkt_static_reg = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "adaptive_rho_interval"));
  settings->kkt_dynamic_reg = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "adaptive_rho_tolerance"));
  settings->abstol = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "adaptive_rho_fraction"));
  settings->reltol = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "max_iter"));
  settings->abstol_inacc = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "eps_abs"));
  settings->reltol_inacc = (QOCOFloat)mxGetScalar(mxGetField(mxPtr, 0, "eps_rel"));
  settings->verbose = (QOCOInt)mxGetScalar(mxGetField(mxPtr, 0, "eps_dual_inf"));
}