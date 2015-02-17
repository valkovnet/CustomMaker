#include "AEConfig.h"
#include "entry.h"

#ifdef AE_OS_WIN
	#include <windows.h>
	#include <stdio.h>
	#include <string.h>
#elif defined AE_OS_MAC
	#include <wchar.h>
#endif

#include "AE_GeneralPlug.h"
#include "AE_Effect.h"
#include "A.h"
#include "AE_EffectUI.h"
#include "SPSuites.h"
#include "AE_AdvEffectSuites.h"
#include "AE_EffectCBSuites.h"
#include "AEGP_SuiteHandler.h"
#include "AE_Macros.h"

#include "MoveMaker_vector.h"
#include "MoveMaker_types.h"

#ifdef AE_OS_WIN
	#define PROXY_FOOTAGE_PATH	L"C:\\proxy.jpg"
	#define LAYERED_PATH		L"C:\\noel_clown_nose.psd"
#elif defined AE_OS_MAC
	#define PROXY_FOOTAGE_PATH	L"proxy.jpg"
	#define LAYERED_PATH		L"noel_clown_nose.psd"
#endif

#define ERROR_MISSING_FOOTAGE	"Footage not found! Make sure you've copied proxy.jpg from the MoveMaker folder and noel_clown_nose.psd from the Projector folder to the root of your primary hard disk."
#define DUMP_SUCCEEDED_WIN		"MoveMaker successfully completed work"
#define DUMP_SUCCEEDED_MAC		"MoveMaker to admin of primary hard drive."
#define DUMP_FAILED				"MoveMaker failed."

// This entry point is exported through the PiPL (.r file)
extern "C" DllExport AEGP_PluginInitFuncPrototype EntryPointFunc;

A_FpLong Dot4D(VectorType4D v1, VectorType4D v2);

#define Artie_MAX_POLYGONS 10

typedef struct {
	AEGP_RenderLayerContextH			layer_contexts[Artie_MAX_POLYGONS];
	A_long								count;
} Artie_LayerContexts, *Artie_LayerContextsP, **Artie_LayerContextsH;


#define	Artie_TRANSFORM_POINT4(_IN, _MATRIX, _OUT)													\
		do {																							\
																										\
			(_OUT).P[0] =	(_IN).P[0] * (_MATRIX).mat[0][0] + (_IN).P[1] * (_MATRIX).mat[1][0] +			\
							(_IN).P[2] * (_MATRIX).mat[2][0] + (_IN).P[3] * (_MATRIX).mat[3][0];			\
																										\
			(_OUT).P[1] =	(_IN).P[0] * (_MATRIX).mat[0][1] + (_IN).P[1] * (_MATRIX).mat[1][1] +			\
							(_IN).P[2] * (_MATRIX).mat[2][1] + (_IN).P[3] * (_MATRIX).mat[3][1];			\
																										\
			(_OUT).P[2] =	(_IN).P[0] * (_MATRIX).mat[0][2] + (_IN).P[1] * (_MATRIX).mat[1][2] +			\
							(_IN).P[2] * (_MATRIX).mat[2][2] + (_IN).P[3] * (_MATRIX).mat[3][2];			\
																										\
			(_OUT).P[3] =	(_IN).P[0] * (_MATRIX).mat[0][3] + (_IN).P[1] * (_MATRIX).mat[1][3] +			\
							(_IN).P[2] * (_MATRIX).mat[2][3] + (_IN).P[3] * (_MATRIX).mat[3][3];			\
																										\
			if ((_OUT).P[3] != 0) {																	\
				A_FpLong	inv_homoF = 1.0 / (_OUT).P[3];												\
				(_OUT).P[0] *=  inv_homoF;																\
				(_OUT).P[1] *=  inv_homoF;																\
				(_OUT).P[2] *=  inv_homoF;																\
				(_OUT).P[3] =   1.0;																	\
			}																							\
		} while (0)

#define	Artie_TRANSFORM_VECTOR4(_IN, _MATRIX, _OUT)													\
		do {																							\
																										\
			(_OUT).V[0] =	(_IN).V[0] * (_MATRIX).mat[0][0] + (_IN).V[1] * (_MATRIX).mat[1][0] +		\
							(_IN).V[2] * (_MATRIX).mat[2][0] + (_IN).V[3] * (_MATRIX).mat[3][0];		\
																										\
			(_OUT).V[1] =	(_IN).V[0] * (_MATRIX).mat[0][1] + (_IN).V[1] * (_MATRIX).mat[1][1] +		\
							(_IN).V[2] * (_MATRIX).mat[2][1] + (_IN).V[3] * (_MATRIX).mat[3][1];		\
																										\
			(_OUT).V[2] =	(_IN).V[0] * (_MATRIX).mat[0][2] + (_IN).V[1] * (_MATRIX).mat[1][2] +		\
							(_IN).V[2] * (_MATRIX).mat[2][2] + (_IN).V[3] * (_MATRIX).mat[3][2];		\
																										\
			(_OUT).V[3] =	(_IN).V[0] * (_MATRIX).mat[0][3] + (_IN).V[1] * (_MATRIX).mat[1][3] +		\
							(_IN).V[2] * (_MATRIX).mat[2][3] + (_IN).V[3] * (_MATRIX).mat[3][3];		\
																										\
			if ((_OUT).V[3] != 0) {																		\
				A_FpLong	inv_homoF = 1.0 / (_OUT).V[3];												\
				(_OUT).V[0] *=  inv_homoF;																\
				(_OUT).V[1] *=  inv_homoF;																\
				(_OUT).V[2] *=  inv_homoF;																\
				(_OUT).V[3] =   1.0;																	\
			}																							\
		} while (0)



#define			AEFX_COPY_STRUCT(FROM, TO)	\
	do {									\
		long _t = sizeof(FROM);				\
		char *_p = (char*)&(FROM);			\
		char *_q = (char*)&(TO);			\
		while (_t--) {						\
			*_q++ = *_p++;					\
		}									\
	} while (0);										