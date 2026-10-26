//-----------------------------------------------------------------------------
//   swmm5.h
//
//   Project: EPA SWMM5
//   Version: 5.1
//   Date:    03/24/14  (Build 5.1.001)
//            08/01/16  (Build 5.1.011)
//   Author:  L. Rossman
//
//   Prototypes for SWMM5 functions exported to swmm5.dll.
//
//-----------------------------------------------------------------------------

#ifndef SWMM5_H
#define SWMM5_H


// --- define WINDOWS

#undef WINDOWS
#ifdef _WIN32
  #define WINDOWS
#endif
#ifdef __WIN32__
  #define WINDOWS
#endif

// --- define DLLEXPORT

#ifdef WINDOWS
    #define DLLEXPORT __declspec(dllexport) __stdcall
#else
    #define DLLEXPORT
#endif

// --- use "C" linkage for C++ programs

#ifdef __cplusplus
extern "C" {
#endif

int  DLLEXPORT   swmm_run(char* f1, char* f2, char* f3);
int  DLLEXPORT   swmm_open(char* f1, char* f2, char* f3);
int  DLLEXPORT   swmm_start(int saveFlag);
int  DLLEXPORT   swmm_step(double* elapsedTime);
int  DLLEXPORT   swmm_end(void);
int  DLLEXPORT   swmm_report(void);
int  DLLEXPORT   swmm_getMassBalErr(float* runoffErr, float* flowErr,
                 float* qualErr);
int  DLLEXPORT   swmm_close(void);
int  DLLEXPORT   swmm_getVersion(void);
int  DLLEXPORT   swmm_getError(char* errMsg, int msgLen);
int  DLLEXPORT   swmm_getWarnings(void);

int  DLLEXPORT   swmm_getBasicJunctions(double **node_x, double **node_y, double **node_elev, double **node_depth);
int  DLLEXPORT   swmm_inputinflow(double *node_q);
int  DLLEXPORT   swmm_resetTimeStep(double timeStep);
int  DLLEXPORT   swmm_getAdvancedJunctions(double *node_z, double *node_outfall);
int  DLLEXPORT   swmm_setOutfallExternalStage(int nodeIndex, double stageMetres);
int  DLLEXPORT   swmm_beginCoupledExchange(void);
int  DLLEXPORT   swmm_getCoupledAudit(double* values);
int  DLLEXPORT   swmm_setOutfallReceivingCell(int nodeIndex, int cellIndex,
                 double bedMetres, double depthMetres, double areaSquareMetres);
int  DLLEXPORT   swmm_advanceCoupled(double seconds, double* elapsedTime);
int  DLLEXPORT   swmm_getNumThreads(void);   /* Dynamic-wave thread count. */
int  DLLEXPORT   swmm_hydraulicRestartSupported(void);
int  DLLEXPORT   swmm_hydraulicRestart(void* stream, int reading);
double DLLEXPORT swmm_getOutfallExchangeVolume(int nodeIndex);
double DLLEXPORT swmm_getNodeAvailableCoupledVolume(int nodeIndex);
double DLLEXPORT swmm_getNodeFloodExchangeVolume(int nodeIndex);

#ifdef __cplusplus
}   // matches the linkage specification from above */
#endif


#endif //SWMM5_H
