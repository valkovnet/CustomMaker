#include "MoveMaker.h"
#include "AE_GeneralPlug.h"

static AEGP_Command	S_dump_proj_cmd = 0,
						S_other_cmd = 0;

static AEGP_PluginID	S_my_id				= 0;

static SPBasicSuite		*sP		=	0;

static AEGP_PluginID		S_aegp_plugin_id;

static SPBasicSuite			*S_sp_basic_suiteP;

#define Artie_MAX_LIGHTS				10
#define Artie_DEFAULT_FIELD_OF_VIEW		1.0
#define Artie_RAYTRACE_THRESHOLD		0.1
#define Artie_MAX_RAY_DEPTH				3

#ifndef ABS
#define ABS(_A) ((_A) < 0 ? -(_A) : (_A))
#endif

#define THE_PHOTO_PATH L"C:\\Temp\\af\\logo.png"

#define FOLDER_NAME	L"Custom folder container"

#define RENDER_PATH L"C:\\Temp\\af\\"

static A_char *workingPath = "C:\\Temp\\af\\";

const A_char *pathToRenderVideo = "c:\\Temp\\af\\Working:output.avi";

const A_UTF16Char *pathToTemplateProject = reinterpret_cast<const A_UTF16Char *>(L"C:\\Temp\\af\\Online Shop Promo\\Online Shop Promo.aep");

const A_UTF16Char *pathToWorkingProject = reinterpret_cast<const A_UTF16Char *>(L"C:\\Temp\\af\\Online Shop Promo\\WorkingProject.aep");

const char* layerFirstText = "01 Intro text";

const char* layerSecondText = "02 Intro text";

const char* layerWithImage = "03 Import Your logo";

const char* layerLastText = "04 Intro text";

const char* layerToRender = "Render";

const char* layerToChangeBackground = "Red Solid 1";

VectorType4D 
Vneg(
VectorType4D 	V)
{
	V.V[X] = -V.V[X];
	V.V[Y] = -V.V[Y];
	V.V[Z] = -V.V[Z];
	
	/* since V is aFP vector, V.V[W] is always 0. */
	return V;
}

A_FpLong 
Dot4D(
	VectorType4D 	v1F, 
	VectorType4D 	v2F)
{
  return v1F.V[0] * v2F.V[0] + v1F.V[1] * v2F.V[1] + v1F.V[2] * v2F.V[2];
}

static VectorType4D 
Pdiff(
	PointType4D P1,
	PointType4D P0)
{
  VectorType4D V;

  NormalizePoint(&P0);
  NormalizePoint(&P1);

  V.V[0] = P1.P[0] - P0.P[0];
  V.V[1] = P1.P[1] - P0.P[1];
  V.V[2] = P1.P[2] - P0.P[2];
  V.V[3] = 0.0;

  return V;
}

A_Err
IdentityMatrix4(
	A_Matrix4 	*matrixP)
{
	AEFX_CLR_STRUCT(*matrixP);

	matrixP->mat[0][0] = 1.0;
	matrixP->mat[1][1] = 1.0;
	matrixP->mat[2][2] = 1.0;
	matrixP->mat[3][3] = 1.0;

  return A_Err_NONE; 
}

A_Err
TranslateMatrix4(
	A_FpLong 	x,
	A_FpLong 	y,
	A_FpLong 	z,
	A_Matrix4 	*resultP)
{
	IdentityMatrix4(resultP);

	resultP->mat[3][0] = x;
	resultP->mat[3][1] = y;
	resultP->mat[3][2] = z;
	resultP->mat[3][3] = 1.0;

	return A_Err_NONE;
}

A_Err
ScaleMatrix4(
	A_FpLong x,
	A_FpLong y,
	A_FpLong z,
	A_Matrix4 *resultP)
{
	IdentityMatrix4(resultP);

	resultP->mat[0][0] = x;
	resultP->mat[1][1] = y;
	resultP->mat[2][2] = z;
	resultP->mat[3][3] = 1.0;

	return A_Err_NONE;
}


A_Err 
InverseMatrix4(
	const	A_Matrix4	*m,
			A_Matrix4 *resultP)
{
	A_Err err = A_Err_NONE;

	A_FpLong d00, d01, d02, d03,
			 d10, d11, d12, d13,
			 d20, d21, d22, d23,
			 d30, d31, d32, d33,
			 m00, m01, m02, m03,
			 m10, m11, m12, m13,
			 m20, m21, m22, m23,
			 m30, m31, m32, m33,
			 D;

	m00 = m->mat[0][0];  m01 = m->mat[0][1];  m02 = m->mat[0][2];  m03 = m->mat[0][3];
	m10 = m->mat[1][0];  m11 = m->mat[1][1];  m12 = m->mat[1][2];  m13 = m->mat[1][3];
	m20 = m->mat[2][0];  m21 = m->mat[2][1];  m22 = m->mat[2][2];  m23 = m->mat[2][3];
	m30 = m->mat[3][0];  m31 = m->mat[3][1];  m32 = m->mat[3][2];  m33 = m->mat[3][3];

	d00 = m11*m22*m33 + m12*m23*m31 + m13*m21*m32 - m31*m22*m13 - m32*m23*m11 - m33*m21*m12;
	d01 = m10*m22*m33 + m12*m23*m30 + m13*m20*m32 - m30*m22*m13 - m32*m23*m10 - m33*m20*m12;
	d02 = m10*m21*m33 + m11*m23*m30 + m13*m20*m31 - m30*m21*m13 - m31*m23*m10 - m33*m20*m11;
	d03 = m10*m21*m32 + m11*m22*m30 + m12*m20*m31 - m30*m21*m12 - m31*m22*m10 - m32*m20*m11;

	d10 = m01*m22*m33 + m02*m23*m31 + m03*m21*m32 - m31*m22*m03 - m32*m23*m01 - m33*m21*m02;
	d11 = m00*m22*m33 + m02*m23*m30 + m03*m20*m32 - m30*m22*m03 - m32*m23*m00 - m33*m20*m02;
	d12 = m00*m21*m33 + m01*m23*m30 + m03*m20*m31 - m30*m21*m03 - m31*m23*m00 - m33*m20*m01;
	d13 = m00*m21*m32 + m01*m22*m30 + m02*m20*m31 - m30*m21*m02 - m31*m22*m00 - m32*m20*m01;

	d20 = m01*m12*m33 + m02*m13*m31 + m03*m11*m32 - m31*m12*m03 - m32*m13*m01 - m33*m11*m02;
	d21 = m00*m12*m33 + m02*m13*m30 + m03*m10*m32 - m30*m12*m03 - m32*m13*m00 - m33*m10*m02;
	d22 = m00*m11*m33 + m01*m13*m30 + m03*m10*m31 - m30*m11*m03 - m31*m13*m00 - m33*m10*m01;
	d23 = m00*m11*m32 + m01*m12*m30 + m02*m10*m31 - m30*m11*m02 - m31*m12*m00 - m32*m10*m01;

	d30 = m01*m12*m23 + m02*m13*m21 + m03*m11*m22 - m21*m12*m03 - m22*m13*m01 - m23*m11*m02;
	d31 = m00*m12*m23 + m02*m13*m20 + m03*m10*m22 - m20*m12*m03 - m22*m13*m00 - m23*m10*m02;
	d32 = m00*m11*m23 + m01*m13*m20 + m03*m10*m21 - m20*m11*m03 - m21*m13*m00 - m23*m10*m01;
	d33 = m00*m11*m22 + m01*m12*m20 + m02*m10*m21 - m20*m11*m02 - m21*m12*m00 - m22*m10*m01;

	D = m00*d00 - m01*d01 + m02*d02 - m03*d03;

	if (D) {
		resultP->mat[0][0] =  d00/D; resultP->mat[0][1] = -d10/D;  resultP->mat[0][2] =  d20/D; resultP->mat[0][3] = -d30/D;
		resultP->mat[1][0] = -d01/D; resultP->mat[1][1] =  d11/D;  resultP->mat[1][2] = -d21/D; resultP->mat[1][3] =  d31/D;
		resultP->mat[2][0] =  d02/D; resultP->mat[2][1] = -d12/D;  resultP->mat[2][2] =  d22/D; resultP->mat[2][3] = -d32/D;
		resultP->mat[3][0] = -d03/D; resultP->mat[3][1] =  d13/D;  resultP->mat[3][2] = -d23/D; resultP->mat[3][3] =  d33/D;
	} else {
		resultP->mat[0][0] =  d00; resultP->mat[0][1] = -d10;  resultP->mat[0][2] =  d20; resultP->mat[0][3] = -d30;
		resultP->mat[1][0] = -d01; resultP->mat[1][1] =  d11;  resultP->mat[1][2] = -d21; resultP->mat[1][3] =  d31;
		resultP->mat[2][0] =  d02; resultP->mat[2][1] = -d12;  resultP->mat[2][2] =  d22; resultP->mat[2][3] = -d32;
		resultP->mat[3][0] = -d03; resultP->mat[3][1] =  d13;  resultP->mat[3][2] = -d23; resultP->mat[3][3] =  d33;
	}
	return err;
}

A_Err
MultiplyMatrix4(
	const A_Matrix4	*A,
	const A_Matrix4	*B,
		  A_Matrix4	*resultP)
{
	A_Err err = A_Err_NONE;

	A_Matrix4	tmp;

	for (register A_u_long iLu = 0; iLu < 4; iLu++){
		for (register int jLu = 0; jLu < 4; jLu++) {
			tmp.mat[iLu][jLu] = A->mat[iLu][0] * B->mat[0][jLu] + 
								A->mat[iLu][1] * B->mat[1][jLu] +
								A->mat[iLu][2] * B->mat[2][jLu] +
								A->mat[iLu][3] * B->mat[3][jLu];
		}
	}
	AEFX_COPY_STRUCT(tmp, *resultP);
	return err;
}

A_Err
TransformPoint4(
	const PointType4D	*pointP,
	const A_Matrix4		*transformP,
		  PointType4D	*resultP)
{
	A_Err err = A_Err_NONE;
	PointType4D	tmp;

	Artie_TRANSFORM_POINT4(*pointP, *transformP, tmp);

	*resultP  = tmp;

	return err;
}

A_Err
TransformVector4(
	const	VectorType4D *vectorP,
	const	A_Matrix4	 *matrixP,
			VectorType4D *resultP)
{
	A_Err err = A_Err_NONE;
	VectorType4D tmp;

	Artie_TRANSFORM_VECTOR4(*vectorP, *matrixP, tmp);

	*resultP = tmp;

	return err;
}

A_Err 
TransposeMatrix4(
	const A_Matrix4 *matrix1,		
	 	  A_Matrix4 *resultP) 		

{
	register A_u_long iLu, jLu;
	A_Matrix4	tmp;

	for (iLu = 0 ; iLu < 4 ; iLu++){
		for (jLu = 0 ; jLu < 4 ; jLu++) {
			tmp.mat[iLu][jLu] = matrix1->mat[jLu][iLu];
		}
	}

	AEFX_COPY_STRUCT(tmp, *resultP);
	
	return A_Err_NONE;
}

Ray CreateRay(PointType4D P0, PointType4D P1)
{
  Ray R;

  R.P0 = P0;
  R.P1 = P1;
  return R;
}

PointType4D 
PVadd(
	PointType4D 	P, 
	VectorType4D 	V)
{
  PointType4D p;

  NormalizePoint(&P);
  p.P[X] = P.P[X] + V.V[X];
  p.P[Y] = P.P[Y] + V.V[Y];
  p.P[Z] = P.P[Z] + V.V[Z];
  p.P[W] = 1.0;

  return p;
}

VectorType4D 
Vadd(
	VectorType4D 	V1, 
	VectorType4D 	V2)
{
  VectorType4D V;

  V.V[X] = V1.V[X] + V2.V[X];
  V.V[Y] = V1.V[Y] + V2.V[Y];
  V.V[Z] = V1.V[Z] + V2.V[Z];
  V.V[W] = 0.0;

  return V;
}

Ray
TransformRay(
	const Ray		*rayP,
	const A_Matrix4	*xformP)
{
	Ray	R;

	TransformPoint4(&rayP->P0, xformP, &R.P0);
	TransformPoint4(&rayP->P1, xformP, &R.P1);

	return R;
}

static void 
NormalizePoint(
	PointType4D *P1)
{
  register double *p = P1->P;
  register double  w = p[W];

  if (w != 1.0){
    /* We are assuming that the order in P is:  X, Y, Z, W */
    *p++ /= w;  
    *p++ /= w; 
    *p++ /= w;
    *p = 1.0;
  }
}

static VectorType4D 
Normalize(
	A_FpLong	(*sqrt)(A_FpLong),  
	VectorType4D v)
{
  VectorType4D vout;
  A_FpLong l = sqrt( v.V[0]*v.V[0] + v.V[1]*v.V[1] + v.V[2]*v.V[2]);
 
  if (l < Artie_EPSILON){
    vout = v;
  } else {
    vout.V[0] = v.V[0]/l;
    vout.V[1] = v.V[1]/l;
    vout.V[2] = v.V[2]/l;
    vout.V[3] = 0.0;
  }
  return vout;
}

/*	Return the vector s*V. */

static VectorType4D 
Vscale(
	PF_FpLong sF,
	VectorType4D V)
{
  V.V[X] *= sF;  
  V.V[Y] *= sF;  
  V.V[Z] *= sF; 
  V.V[W] = 0.0;

  return V;
}

static A_Err
PrintAndDisposeStream(
	AEGP_StreamRefH		streamH,
	A_char				*format_stringZ,
	A_char				*indent_stringZ,
	A_char				*stream_nameZ,
	FILE 				*out)
{
	A_Err 				err = A_Err_NONE, err2;
	AEGP_StreamType		stream_type;
	A_long				num_kfsL;
	A_Time				sampleT = { 0, 100 };
	AEGP_StreamValue	val;
	AEGP_StreamValue	*sample_valP = &val;
	A_char				unitsAC[AEGP_MAX_STREAM_UNITS_SIZE + 1];
	
	AEGP_SuiteHandler	suites(sP);
	
	AEFX_CLR_STRUCT(val);

	ERR(suites.StreamSuite2()->AEGP_GetStreamName(streamH, FALSE, stream_nameZ));
	
	if (!err && stream_nameZ) {		// don't bother to print blank name...
		fprintf(out, format_stringZ, indent_stringZ, stream_nameZ);

		ERR(suites.StreamSuite2()->AEGP_GetStreamType(streamH, &stream_type));

		ERR(suites.StreamSuite2()->AEGP_GetStreamUnitsText(streamH, FALSE, unitsAC));
		ERR(suites.KeyframeSuite3()->AEGP_GetStreamNumKFs(streamH, &num_kfsL));
		
		if (!err && num_kfsL) {
			A_Boolean pre_expressionB	=	TRUE;
			
			if (AEGP_StreamType_NO_DATA != stream_type)	{
				ERR(suites.StreamSuite2()->AEGP_GetNewStreamValue(S_my_id, 
																	streamH, 
																	AEGP_LTimeMode_LayerTime, 
																	&sampleT, 
																	pre_expressionB, 
																	sample_valP));
			}
			if (!err) {
				switch (stream_type) 
				{
					case AEGP_StreamType_TwoD_SPATIAL:
					case AEGP_StreamType_TwoD:
						fprintf(out, "\t\tX %1.2f %s Y %1.2f %s", sample_valP->val.two_d.x, unitsAC, sample_valP->val.two_d.y, unitsAC);
						break;		
					case AEGP_StreamType_OneD:
						fprintf(out, "\t\t%1.2f %s", sample_valP->val.one_d, unitsAC);
						break;		
					case AEGP_StreamType_COLOR:				
						fprintf(out, "\t\tR %1.2f G %1.2f B %1.2f", sample_valP->val.color.redF, sample_valP->val.color.greenF, sample_valP->val.color.blueF);
						break;		
					case AEGP_StreamType_ARB:
						fprintf(out, "arb data");
						break;		
					case AEGP_StreamType_LAYER_ID:
						fprintf(out, "\t\tID%d", sample_valP->val.layer_id);
						break;		
					case AEGP_StreamType_MASK_ID:
						fprintf(out, "\t\tID%d", sample_valP->val.mask_id);
						break;							
				}
				if (AEGP_StreamType_NO_DATA != stream_type)	{
					ERR(suites.StreamSuite2()->AEGP_DisposeStreamValue(sample_valP));
				}
			}
		}		

		if (!err) {
			fprintf(out, "\n");	
		}
	}

	ERR2(suites.StreamSuite2()->AEGP_DisposeStream(streamH));

	return err;
}

static A_Err
Artie_GetPlaneEquation(
	Artie_Poly		*polyP,
	A_FpLong		*aFP,
	A_FpLong		*bFP,
	A_FpLong		*cFP,
	A_FpLong		*dFP)
{
	A_Err			err = A_Err_NONE;
	A_FpLong		sum_x = 0, 
					sum_y = 0, 
					sum_z = 0;

	for (A_long iL = 0; iL < 4; iL++) {
		sum_x += polyP->vertsA[iL].coord.P[0];
		sum_y += polyP->vertsA[iL].coord.P[1];
		sum_z += polyP->vertsA[iL].coord.P[2];
	}
	sum_x *= 0.25;
	sum_y *= 0.25;
	sum_z *= 0.25;

	*aFP = polyP->normal.V[0];
	*bFP = polyP->normal.V[1];
	*cFP = polyP->normal.V[2];

	*dFP = -(sum_x * (*aFP) + sum_y * (*bFP) + sum_z * (*cFP));
	
	return err;
}

static A_Err
Artie_RayIntersectPlane(
	Artie_Poly		*polyP,
	Ray				*rayP,
	A_FpLong		*tPF)
{
	A_Err			err 	= A_Err_NONE;
	A_FpLong		topF 	= 0.0, 
					bottomF = 0.0,
					aFP 	= 0.0, 
					bFP 	= 0.0, 
					cFP		= 0.0,
					dFP 	= 0.0;
					
	VectorType4D	dir;

	ERR(Artie_GetPlaneEquation(polyP, &aFP, &bFP, &cFP, &dFP));
	if (!err) {
		dir		= Pdiff(rayP->P1, rayP->P0);
		bottomF = aFP * dir.V[0] +bFP * dir.V[1] + cFP * dir.V[2];

		if (bottomF){
			topF = -(dFP + aFP*rayP->P0.P[0] + bFP*rayP->P0.P[1] + cFP*rayP->P0.P[2]);
			*tPF = topF / bottomF;
		} else {
			*tPF = Artie_MISS;
		}
	}
	return err;
}

static A_Err
Artie_RayIntersectLayer(
	Artie_Poly  	*polyP,
	Ray				*rayP,
	A_FpLong		*tPF,
	A_FpLong		*uPF,
	A_FpLong		*vPF)
{
	A_Err		err 		= A_Err_NONE;
	A_FpLong	u0F 		= 0.0, 
				v0F			= 0.0,
				u1F			= 0.0,
				v1F			= 0.0, 
				u2F			= 0.0, 
				v2F			= 0.0,
				alphaF		= 0.0,
				betaF		= 0.0,
				gammaF		= 0.0;	
	A_Boolean	insideB 	= FALSE;
	A_long		i1L			= 0,
			 	i2			= 0,
			 	kL			= 0;
	PointType4D		P;
	VectorType4D	dir;

	*uPF = *vPF = -1;

	ERR(Artie_RayIntersectPlane(polyP, rayP, tPF));

	if ( *tPF != Artie_MISS) {
		// find uv-coords
		dir = Pdiff(rayP->P1, rayP->P0);
		P.P[0] = rayP->P0.P[0] + dir.V[0] * *tPF;
		P.P[1] = rayP->P0.P[1] + dir.V[1] * *tPF;
		P.P[2] = rayP->P0.P[2] + dir.V[2] * *tPF;

		if (0 == polyP->dominant) {
			i1L = 1;
			i2 = 2;
		} else if ( 1 == polyP->dominant) {
			i1L = 0;
			i2 = 2;
		} else {
			i1L = 0;
			i2 = 1;
		}

		u0F = P.P[i1L] - polyP->vertsA[0].coord.P[i1L];
		v0F = P.P[i2] - polyP->vertsA[0].coord.P[i2];
		
		kL = 2;
		do{
			u1F = polyP->vertsA[kL-1].coord.P[i1L] - polyP->vertsA[0].coord.P[i1L];
			v1F = polyP->vertsA[kL-1].coord.P[i2] - polyP->vertsA[0].coord.P[i2];
			u2F = polyP->vertsA[kL  ].coord.P[i1L] - polyP->vertsA[0].coord.P[i1L];
			v2F = polyP->vertsA[kL  ].coord.P[i2] - polyP->vertsA[0].coord.P[i2];
			
			if (u1F == 0){
				betaF = u0F/u2F;
				if (betaF >= 0.0 && betaF <= 1.0){
					alphaF 	= (v0F - betaF*v2F)/v1F;
					insideB = ((alphaF>=0.) && (alphaF + betaF) <= 1);
				}
			} else {
				betaF = (v0F * u1F - u0F * v1F)/(v2F * u1F - u2F * v1F);
				if (betaF >= 0. && betaF <= 1.) {
					alphaF 	= (u0F - betaF*u2F) / u1F;
					insideB = ((alphaF >= 0.) && (alphaF + betaF) <= 1);
				}
			}
		} while (!insideB && (++kL < 4));

		if (insideB){
			gammaF = 1.0 - alphaF - betaF;
			*uPF = gammaF * polyP->vertsA[0].txtur[0] + alphaF * polyP->vertsA[kL-1].txtur[0] + betaF * polyP->vertsA[kL].txtur[0];
			*vPF = gammaF * polyP->vertsA[0].txtur[1] + alphaF * polyP->vertsA[kL-1].txtur[1] + betaF * polyP->vertsA[kL].txtur[1];
		} else {
			*tPF = Artie_MISS;
		}
	}
	return err;
}

static A_Err
Artie_FillIntersection(
	Ray				*rayP,
	Artie_Poly		*polyP,
	A_FpLong		t,
	A_FpLong		u,
	A_FpLong		v,
	Intersection	*intersectionP)
{
	A_Err err = A_Err_NONE;
	
	AEGP_SuiteHandler suites(S_sp_basic_suiteP);	
	AEGP_WorldSuite2 *wsP = suites.WorldSuite2();
	
	A_long	widthL	=	0,
			heightL	=	0;
			
	A_u_long	rowbytesLu = 0;
	PF_Pixel8*	baseP = NULL;
	
	ERR(wsP->AEGP_GetSize(polyP->texture, &widthL, &heightL));
	ERR(wsP->AEGP_GetBaseAddr8(polyP->texture, &baseP));
	ERR(wsP->AEGP_GetRowBytes(polyP->texture, &rowbytesLu));

	if ( t != Artie_MISS) {
		intersectionP->hitB 						= TRUE;
		intersectionP->intersect_point_coord_on_ray = t;

		// build reflected ray
		VectorType4D	dir, reflected_dir;
		A_FpLong		dotF;

		dir								= Pdiff(rayP->P1, rayP->P0);
		dotF							= Dot4D(dir, polyP->normal);
		reflected_dir					= Vadd(polyP->normal, Vscale(-2.0, dir));
		intersectionP->ray.P0			= PVadd(rayP->P0, Vscale(t, dir));
		intersectionP->ray.P0			= PVadd(intersectionP->ray.P0, Vscale(0.0001, reflected_dir));
		intersectionP->ray.P1			= PVadd(intersectionP->ray.P0, reflected_dir);
		intersectionP->ray.depth		= rayP->depth + 1;
		intersectionP->ray.contribution	= rayP->contribution * polyP->material.ksF;
		
		intersectionP->reflectance_percentF = polyP->material.ksF;
		
		if ( u >= 0 && v >= 0 && u <= 1 && v <= 1) {
			A_long iL = static_cast<A_long>(u * (widthL-1));
			A_long j = static_cast<A_long>(v * (heightL-1));
			PF_Pixel8 *pP = baseP + iL;
			intersectionP->color_struck_surface = *(PF_Pixel *)((char *)pP + j * rowbytesLu);
		}
	}
	return err;
}

static 
Intersection
SceneIntersection(
	Artie_Scene		*sceneP, 
	Ray				*rayP)
{
	A_Err			err			= A_Err_NONE;
	Intersection	returned_I, 
					I[Artie_MAX_POLYGONS];
	A_FpLong		u = 0.0, 
					v = 0.0;

	returned_I.hitB = FALSE;
	returned_I.intersect_point_coord_on_ray = Artie_MISS;
	returned_I.reflectance_percentF = 0.0;

	for(A_long iL = 0; iL < sceneP->num_polysL; ++iL){
		I[iL].hitB = 0;
		err = Artie_RayIntersectLayer(&sceneP->polygons[iL], rayP, &I[iL].intersect_point_coord_on_ray, &u, &v);
		if (!err && I[iL].intersect_point_coord_on_ray != Artie_MISS){
			I[iL].hitB = 1;
		}
	}
	
	if (!err) {
		for(A_long iL = 0; iL < sceneP->num_polysL; iL++) {
			if (I[iL].hitB) {
				if (I[iL].intersect_point_coord_on_ray < returned_I.intersect_point_coord_on_ray && I[iL].intersect_point_coord_on_ray > 0){
					err = Artie_FillIntersection(	rayP, 
													&sceneP->polygons[iL], 
													I[iL].intersect_point_coord_on_ray, 
													u, 
													v, 
													&returned_I);
				}
			}
		}
	}
	return returned_I;
}

static PF_Pixel 
Raytrace(
	const PR_InData				*in_dataP,	
		  Artie_Scene			*sceneP,
		  Ray					*rayP)
{
	Intersection    I;
	VectorType4D    Nhat,Rvec, Rrefl;
	Ray             RefRay;
	PF_Pixel        cFP, crefl;
	A_FpLong        cosThetaI;
	AEGP_SuiteHandler suites(in_dataP->pica_basicP);	
 
	if (rayP->contribution < Artie_RAYTRACE_THRESHOLD){
		cFP.red 	= 0;
		cFP.green	= 0;
		cFP.blue 	= 0;
		cFP.alpha	= 0;
		return cFP;
	}

	if (rayP->depth > Artie_MAX_RAY_DEPTH) {
		cFP.red 	= 0;
		cFP.green	= 0;
		cFP.blue 	= 0;
		cFP.alpha	= 0;
		return cFP;
	}

	I = SceneIntersection(sceneP, rayP);

	if (I.hitB){
		cFP.red	= I.color_struck_surface.red;
		cFP.green	= I.color_struck_surface.green;
		cFP.blue	= I.color_struck_surface.blue;
		cFP.alpha	= I.color_struck_surface.alpha;

		Nhat = Normalize(suites.ANSICallbacksSuite1()->sqrt, Pdiff(I.ray.P1, I.ray.P0));
		Rvec = Normalize(suites.ANSICallbacksSuite1()->sqrt, Pdiff(rayP->P1, rayP->P0));

		/* Compute the direction of the reflected ray */
		cosThetaI = -Dot4D(Rvec, Nhat);
		
		Rrefl = Vadd(Rvec, Vscale(2.0 * cosThetaI, Nhat));

		if (I.reflectance_percentF){

			/* Cast aFP reflected ray */

			RefRay = CreateRay(I.ray.P0, PVadd(I.ray.P0, Rrefl));
			RefRay.contribution = rayP->contribution * I.reflectance_percentF;
			RefRay.depth = rayP->depth + 1;
			crefl = Raytrace(in_dataP,  sceneP, &RefRay);

			cFP.red 	+= (A_u_char)(I.reflectance_percentF * crefl.red);
			cFP.green	+= (A_u_char)(I.reflectance_percentF * crefl.green);
			cFP.blue 	+= (A_u_char)(I.reflectance_percentF * crefl.blue);
		}
	} else {
		cFP.red		= 0;
		cFP.green	= 0;
		cFP.blue	= 0;
		cFP.alpha	= 0;
	}
	return cFP;
}

static 
PF_Pixel ThrowRay(
	const	PR_InData				*in_dataP,	
			Artie_Scene				*sceneP,
			A_FpLong				xF, 
			A_FpLong				yF,
			A_FpLong				zF)
{
  Ray		R;
  PointType4D P0, P1;

  P0.P[X] = 0;
  P0.P[Y] = 0;
  P0.P[Z] = 0; 
  P0.P[W] = 1;

  P1.P[X] = xF;
  P1.P[Y] = yF; 
  P1.P[Z] = zF; 
  P1.P[W] = 1;

  R = CreateRay(P0, P1);
  R.contribution 	= 1.0;
  R.depth 			= 0;
  return Raytrace(in_dataP,  sceneP, &R);
}

static A_Err
Artie_SampleImage(
	const PR_InData				*in_dataP,	
	PR_RenderContextH			render_contextH, 
	Artie_Scene					*sceneP,
	Artie_Camera				*cameraP,
	AEGP_WorldH					*rwP)
{
	A_Err				err		= A_Err_NONE;
	A_FpLong			x		= 0, 
						y		= 0, 
						z		= cameraP->focal_lengthF;
	PF_Pixel8			*pixelP	= NULL,
						*baseP	= NULL;
	
	A_long				widthL	= 0,
						heightL	= 0;
						
	A_u_long			rowbytesLu = 0;
						
	
	AEGP_SuiteHandler	suites(in_dataP->pica_basicP);
	
	ERR(suites.WorldSuite2()->AEGP_GetSize(*rwP, &widthL, &heightL));
	ERR(suites.WorldSuite2()->AEGP_GetBaseAddr8(*rwP, &baseP));
	ERR(suites.WorldSuite2()->AEGP_GetRowBytes(*rwP, &rowbytesLu));
	
	if (z < 0){
		z *= -1;
	}
	
	for (A_long iL = 0; iL < heightL; iL++) {
		y	= -heightL/2.0 + iL + 0.5;

		pixelP = (PF_Pixel8*)( (char *)baseP + iL * rowbytesLu);	
			
		for (A_long jL = 0; jL < widthL; jL++){
			x = -widthL / 2.0 + jL + 0.5;
			*pixelP++ = ThrowRay(	in_dataP, sceneP, x, y, z);
		}
	}
	return err;
}

static A_Err
Artie_CreateDefaultCamera(
	const PR_InData				*in_dataP,
	AEGP_CompH					compH,
	AEGP_DownsampleFactor		*dsfP,
	AEGP_MemHandle				*cameraPH)
{
	A_Err					err				= A_Err_NONE;
	Artie_Camera			*cameraP		= NULL;
	AEGP_ItemH				itemH			= NULL;
	A_long					widthL 			= 0, 
							heightL			= 0;
	A_FpLong				comp_origin_xF 	= 0.0, 
							comp_origin_yF 	= 0.0, 
							min_dimensionF 	= 0.0;

	A_Matrix4				matrix;

	AEGP_SuiteHandler 	suites(in_dataP->pica_basicP);	

	*cameraPH = NULL;

	ERR(suites.MemorySuite1()->AEGP_NewMemHandle(	in_dataP->aegp_plug_id, 
													"Camera Data", 
													sizeof(Artie_Camera), 
													AEGP_MemFlag_CLEAR, 
													cameraPH));

	ERR(suites.MemorySuite1()->AEGP_LockMemHandle(*cameraPH, reinterpret_cast<void**>(&cameraP)));
	ERR(suites.CompSuite4()->AEGP_GetItemFromComp(compH, &itemH));
	ERR(suites.ItemSuite6()->AEGP_GetItemDimensions(itemH, &widthL, &heightL));

	//  comp origin is the middle of the comp in x and y, and z = 0.

	if (!err){
		comp_origin_xF = widthL / 2.0;
		comp_origin_yF = heightL / 2.0;
		min_dimensionF = MIN(comp_origin_xF, comp_origin_yF);
		A_FpLong tanF			= suites.ANSICallbacksSuite1()->tan(0.5 * Artie_DEFAULT_FIELD_OF_VIEW);
		cameraP->type			= AEGP_CameraType_PERSPECTIVE;
		cameraP->res_xLu		= static_cast<A_u_long>(widthL);
		cameraP->res_yLu		= static_cast<A_u_long>(heightL);
		cameraP->focal_lengthF	= min_dimensionF / tanF;			
		cameraP->dsf			= *dsfP;

		ERR(IdentityMatrix4(&matrix));
		ERR(TranslateMatrix4(comp_origin_xF, 
							comp_origin_yF, 
							-cameraP->focal_lengthF, 
							&matrix));
		ERR(InverseMatrix4(&matrix, &cameraP->view_matrix));
		ERR(suites.MemorySuite1()->AEGP_UnlockMemHandle(*cameraPH));
	}
	return err;
}

static A_Err
Artie_DisposeCamera(
	const PR_InData			*in_dataP,	
	AEGP_MemHandle			cameraH)
{
	A_Err				err		 = A_Err_NONE;
	Artie_Camera		*cameraP = NULL;
	AEGP_SuiteHandler 	suites(in_dataP->pica_basicP);	

	if (!cameraH) {
		err = PF_Err_UNRECOGNIZED_PARAM_TYPE;
		in_dataP->msg_func(err, "Trying to dispose NULL camera");
	}
	ERR(suites.MemorySuite1()->AEGP_LockMemHandle(cameraH, reinterpret_cast<void**>(&cameraP)));
	ERR(suites.MemorySuite1()->AEGP_FreeMemHandle(cameraH));

	return err;
}

static A_Err
Artie_DisposePolygonTexture(
	const PR_InData			*in_dataP,
	PR_RenderContextH		render_contextH,
	Artie_Poly				*polygonP)
{
	A_Err		err	= A_Err_NONE;
	AEGP_SuiteHandler suites(in_dataP->pica_basicP);	
	AEGP_CanvasSuite5	*canP	= suites.CanvasSuite5();

	if (polygonP->texture){
		err = canP->AEGP_DisposeTexture(render_contextH, polygonP->layer_contextH, polygonP->texture);
	}
	return err;
}

static A_Err
Artie_DisposeScene(
	const PR_InData				*in_dataP,	
	PR_RenderContextH			render_contextH,
	AEGP_MemHandle				sceneH)
{
	A_Err			err 	= A_Err_NONE, 
					err2 	= A_Err_NONE;
	Artie_Scene		*sceneP = NULL;	
	
	AEGP_SuiteHandler suites(in_dataP->pica_basicP);	

	if (sceneH) {
		ERR(suites.MemorySuite1()->AEGP_LockMemHandle(sceneH, reinterpret_cast<void**>(&sceneP)));
		if (!err) {
			for(A_long iL = 0; iL< sceneP->num_polysL; iL++){
				ERR(Artie_DisposePolygonTexture(in_dataP, render_contextH, &sceneP->polygons[iL]));
			}
		}
		ERR2(suites.MemorySuite1()->AEGP_FreeMemHandle(sceneH));
	}
	return err;
}

static A_Err
Artie_DisposeLights(
	const PR_InData				*in_dataP,	
	AEGP_MemHandle				lightsH)
{
	AEGP_SuiteHandler suites(in_dataP->pica_basicP);	
	return suites.MemorySuite1()->AEGP_FreeMemHandle(lightsH);
}

static A_Err
Artie_CreateLayerCamera(
	const PR_InData				*in_dataP,
	AEGP_CompH					compH,
	AEGP_DownsampleFactor		*dsfP,
	A_Time						comp_time,
	A_Rect						*roiRP0,
	AEGP_LayerH					camera_layerH,
	AEGP_MemHandle				*cameraPH)
{
	A_Err				err				= A_Err_NONE,
						err2			= A_Err_NONE;
	Artie_Camera		*cameraP		= NULL;
	AEGP_ItemH			itemH			= NULL;
	A_long				widthL			= 0, 
						heightL			= 0;	
	A_Ratio				pix_aspectR		= {1,1};
	A_FpLong			comp_origin_xF	= 0.0, 
						comp_origin_yF	= 0.0;
	A_Matrix4			matrix;
	
	AEGP_SuiteHandler suites(in_dataP->pica_basicP);	

	*cameraPH = NULL;

	ERR(suites.MemorySuite1()->AEGP_NewMemHandle(	in_dataP->aegp_plug_id, 
													"Camera Data", 
													sizeof(Artie_Camera), 
													AEGP_MemFlag_CLEAR,  
													cameraPH));	

	ERR(suites.MemorySuite1()->AEGP_LockMemHandle(*cameraPH, reinterpret_cast<void**>(&cameraP)));
	ERR(suites.CompSuite4()->AEGP_GetItemFromComp(compH, &itemH));
	ERR(suites.ItemSuite6()->AEGP_GetItemDimensions(itemH, &widthL, &heightL));
	ERR(suites.ItemSuite6()->AEGP_GetItemPixelAspectRatio(itemH, &pix_aspectR));

	//  comp origin is the middle of the comp in x and y, and z = 0.
	if (!err) {
		comp_origin_xF					=	widthL / 2.0;
		comp_origin_yF					=	heightL / 2.0;
		cameraP->view_aspect_ratioF		=	widthL / heightL;
		cameraP->pixel_aspect_ratioF	=	pix_aspectR.num / pix_aspectR.den;
		cameraP->dsf					=	*dsfP;
	}

	ERR(suites.LayerSuite5()->AEGP_GetLayerToWorldXform(camera_layerH, &comp_time, &matrix));
	ERR(InverseMatrix4(&matrix, &cameraP->view_matrix));
	ERR(suites.CameraSuite2()->AEGP_GetCameraType(camera_layerH, &cameraP->type));

	if (!err) {
		cameraP->res_xLu	= (unsigned long) widthL;				
		cameraP->res_yLu	= (unsigned long) heightL;
		cameraP->dsf		= *dsfP;
		
		if (roiRP0){
			cameraP->roiR		=	*roiRP0;
		}
	}
	if (cameraP){
		ERR2(suites.MemorySuite1()->AEGP_UnlockMemHandle(*cameraPH));
	}
	if (*cameraPH){
		ERR2(Artie_DisposeCamera(in_dataP, *cameraPH));
		*cameraPH = NULL;
	}
	return err;
}


static A_Err
Artie_BuildCamera(
	const PR_InData				*in_dataP,	
	PR_RenderContextH			render_contextH,		
	AEGP_MemHandle				*cameraPH,
	AEGP_SuiteHandler suites)
{
	A_Err					err				= A_Err_NONE;
	AEGP_LayerH				camera_layerH	= NULL;
	A_Time					render_time		= {0,1}, 
							time_step		= {0,1};
	AEGP_DownsampleFactor	dsf				= {1,1};
	AEGP_CompH				compH 			= NULL;
	A_LegacyRect			roiLR			= {0,0,0,0};
	A_Rect					roiR			= {0,0,0,0};

	*cameraPH = NULL;

	ERR(suites.CanvasSuite5()->AEGP_GetCompToRender(render_contextH, &compH));
	ERR(suites.CanvasSuite5()->AEGP_GetRenderDownsampleFactor(render_contextH, &dsf));
	ERR(suites.CanvasSuite5()->AEGP_GetCompRenderTime(render_contextH, &render_time, &time_step));
	ERR(suites.CameraSuite2()->AEGP_GetCamera(render_contextH, &render_time, &camera_layerH));
	ERR(suites.CanvasSuite5()->AEGP_GetROI(render_contextH, &roiLR));
	if (!err && !camera_layerH){
		ERR(Artie_CreateDefaultCamera(in_dataP,  compH, &dsf, cameraPH));
	} else 	{ 
		roiR.top	= roiLR.top;
		roiR.left	= roiLR.left;
		roiR.right	= roiLR.right;
		roiR.bottom	= roiLR.bottom;
		ERR(Artie_CreateLayerCamera(	in_dataP,  
										compH, 
										&dsf, 
										render_time, 
										&roiR,
										camera_layerH, 
										cameraPH));
	}
	return err;
}

static A_Err
Artie_CreateListOfLayerContexts(
	const PR_InData				*in_dataP,
	PR_RenderContextH			render_contextH,
	A_long						*num_layersPL,
	Artie_LayerContexts			*layer_contexts)	
{
	A_Err						err 					= A_Err_NONE;
	AEGP_RenderLayerContextH	layer_contextH 			= NULL;
	A_long						local_layers_to_renderL	= 0;

	AEGP_SuiteHandler			suites(in_dataP->pica_basicP);
	AEGP_CanvasSuite5			*csP 	=	 suites.CanvasSuite5();	

	*num_layersPL =	0;
	
	err = csP->AEGP_GetNumLayersToRender(render_contextH, &local_layers_to_renderL);
	
	if (!err) {
		for(A_long iL = 0; (!err && (iL < local_layers_to_renderL) && (iL < Artie_MAX_POLYGONS)); iL++){
			ERR(csP->AEGP_GetNthLayerContextToRender(render_contextH, iL, &layer_contextH));
			if (!err){
				layer_contexts->layer_contexts[iL] = layer_contextH;
				layer_contexts->count++;
			}
		}
		*num_layersPL = local_layers_to_renderL;
	}
	return err;
}

static A_Err
Artie_BuildLight(
	const PR_InData			*in_dataP,			/* >> */
	PR_RenderContextH		render_contextH,	/* >> */
	AEGP_LayerH				layerH,				/* >> */
	Artie_Light				*lightP)			/* <> */
{
	A_Err				err	= A_Err_NONE;
	A_Time				comp_time, comp_time_step;
	AEGP_StreamVal		stream_val;
	Artie_Light			light;

	AEGP_SuiteHandler	suites(in_dataP->pica_basicP);

	ERR(suites.CanvasSuite5()->AEGP_GetCompRenderTime(render_contextH, &comp_time, &comp_time_step));
	ERR(suites.LightSuite2()->AEGP_GetLightType(layerH, &lightP->type));
	ERR(suites.StreamSuite2()->AEGP_GetLayerStreamValue(layerH, 
														AEGP_LayerStream_COLOR,
														AEGP_LTimeMode_CompTime,
														&comp_time,
														FALSE,
														&stream_val,
														NULL));

	light.color.alpha 	= (A_u_char)(PF_MAX_CHAN8 * stream_val.color.alphaF 	+ 0.5);
	light.color.red 	= (A_u_char)(PF_MAX_CHAN8 * stream_val.color.redF 		+ 0.5);
	light.color.green 	= (A_u_char)(PF_MAX_CHAN8 * stream_val.color.greenF 	+ 0.5);
	light.color.blue 	= (A_u_char)(PF_MAX_CHAN8 * stream_val.color.blueF 		+ 0.5);

	ERR(suites.StreamSuite2()->AEGP_GetLayerStreamValue(	layerH, 
															AEGP_LayerStream_INTENSITY,
															AEGP_LTimeMode_CompTime,
															&comp_time,
															FALSE,
															&stream_val,
															NULL));
	if (!err){
		light.intensityF = stream_val.one_d;
		err = suites.StreamSuite2()->AEGP_GetLayerStreamValue( 	layerH, 
																AEGP_LayerStream_CASTS_SHADOWS,
																AEGP_LTimeMode_CompTime,
																&comp_time,
																FALSE,
																&stream_val,
																NULL);
		light.casts_shadowsB = (stream_val.one_d != 0.0);
	}
	
	// Normalize intensity

	if (!err){
		lightP->intensityF *= 0.01;

		light.direction.x 	= 0;
		light.direction.y 	= 0;
		light.direction.z 	= 1;
		light.position.x	= 0;
		light.position.y 	= 0;
		light.position.z 	= 0;

		// GetLightAttenuation was also unimplemented and only returned this constant.
		light.attenuation.x = 1;
		light.attenuation.y = 0;
		light.attenuation.z = 0;
	
		err = suites.StreamSuite2()->AEGP_GetLayerStreamValue( layerH, 
																AEGP_LayerStream_CONE_ANGLE,
																AEGP_LTimeMode_CompTime,
																&comp_time,
																FALSE,
																&stream_val,
																NULL);
		light.cone_angleF = stream_val.one_d;

		ERR(suites.StreamSuite2()->AEGP_GetLayerStreamValue( layerH, 
											AEGP_LayerStream_CONE_FEATHER,
											AEGP_LTimeMode_CompTime,
											&comp_time,
											FALSE,
											&stream_val,
											NULL));
		light.spot_exponentF = stream_val.one_d;
		
		ERR(suites.LayerSuite5()->AEGP_GetLayerToWorldXform(	layerH, 
											&comp_time, 
											&light.light_to_world_matrix));
		ERR(suites.StreamSuite2()->AEGP_GetLayerStreamValue(	layerH, 
															AEGP_LayerStream_SHADOW_DARKNESS,
															AEGP_LTimeMode_CompTime,
															&comp_time,
															FALSE,
															&stream_val,
															NULL));
		light.shadow_densityF = stream_val.one_d;
		if (!err) {	
			light.shadow_densityF *= 0.01;
		}

		if (!err && light.type != AEGP_LightType_PARALLEL) 	{
			err = suites.StreamSuite2()->AEGP_GetLayerStreamValue(	layerH, 
																	AEGP_LayerStream_SHADOW_DIFFUSION,
																	AEGP_LTimeMode_CompTime,
																	&comp_time,
																	FALSE,
																	&stream_val,
																	NULL);
			light.radiusF = stream_val.one_d;
		}
	}
	return err;
}	

static A_Err
Artie_GetSourceLayerSize(
	const PR_InData				*in_dataP,				/* >> */
	PR_RenderContextH			render_contextH,		/* >> */
	AEGP_LayerH					layerH,					/* >> */
	A_FpLong					*widthPF,				/* << */
	A_FpLong					*heightPF,				/* << */
	AEGP_DownsampleFactor		*dsfP)
{
	A_Err 					err 			= 	A_Err_NONE;
	AEGP_ItemH				source_itemH 	= 	NULL;
	AEGP_DownsampleFactor	dsf				=	{1, 1};
	A_long					widthL			=	0, 
							heightL			=	0;

	AEGP_SuiteHandler suites(in_dataP->pica_basicP);	

	if (widthPF		==	NULL	||
		heightPF	==	NULL	||
		dsfP		==	NULL){
		err = A_Err_PARAMETER;
	}
	if (!err){
		if (widthPF){
			*widthPF	= 0.0;
		}
		if (heightPF){
			*heightPF	= 0.0;
		}
		if (dsfP) {
			dsfP->xS	=
			dsfP->yS	= 1;
		}

		// This doesn't return a valid item if the layer is a text layer
		// Use AEGP_GetCompRenderTime, AEGP_GetRenderLayerBounds, AEGP_GetRenderDownsampleFactor instead?
		ERR(suites.LayerSuite5()->AEGP_GetLayerSourceItem(layerH, &source_itemH));
		
		if (!err && source_itemH){
			err = suites.ItemSuite6()->AEGP_GetItemDimensions(source_itemH, &widthL, &heightL);	
			
			if (!err){
				*widthPF	= widthL 	/ (A_FpLong)dsf.xS;
				*heightPF	= heightL 	/ (A_FpLong)dsf.yS;
				
				ERR(suites.CanvasSuite5()->AEGP_GetRenderDownsampleFactor(render_contextH, &dsf));
			
				if (!err && dsfP){
					*dsfP = dsf;
				}
			}
		}
	}
	return err;
}

static A_Err
Artie_BuildLights(
	const PR_InData				*in_dataP,	
	PR_RenderContextH			render_contextH,		
	AEGP_MemHandle				*lightsPH)
{
	A_Err				err			=	A_Err_NONE, 
						err2		= A_Err_NONE;
	A_long				layers_to_renderL	=	0;
	Artie_Light			*lightP 	= 	NULL;
	A_Boolean			is_activeB	=	FALSE;
	AEGP_LayerH			layerH		=	NULL;;
	A_Time				comp_time		=	{0, 1}, 
						comp_time_step	=	{0, 1};
	AEGP_CompH			compH 		= 	NULL;
	AEGP_ObjectType		layer_type	=	AEGP_ObjectType_NONE;
	
	AEGP_SuiteHandler	suites(in_dataP->pica_basicP);

	ERR(suites.MemorySuite1()->AEGP_NewMemHandle(	S_aegp_plugin_id, 
													"Light Data", 
													sizeof(Artie_Light), 
													AEGP_MemFlag_CLEAR, 
													lightsPH));

	ERR(suites.MemorySuite1()->AEGP_LockMemHandle(*lightsPH, reinterpret_cast<void**>(&lightP)));
	ERR(suites.CanvasSuite5()->AEGP_GetCompToRender(render_contextH, &compH));
	ERR(suites.LayerSuite5()->AEGP_GetCompNumLayers(compH, &layers_to_renderL));
	
	for (A_long iL = 0; !err && iL < layers_to_renderL && iL < Artie_MAX_LIGHTS; iL++){		
		ERR(suites.LayerSuite5()->AEGP_GetCompLayerByIndex(compH, iL, &layerH));
		ERR(suites.CanvasSuite5()->AEGP_GetCompRenderTime(render_contextH, &comp_time, &comp_time_step));
		ERR(suites.LayerSuite5()->AEGP_IsVideoActive(layerH, 
													AEGP_LTimeMode_CompTime,
													&comp_time, 
													&is_activeB));
		ERR(suites.LayerSuite5()->AEGP_GetLayerObjectType(layerH, &layer_type));
		
		if (!err && (AEGP_ObjectType_LIGHT == layer_type)  && is_activeB){
			err = Artie_BuildLight(in_dataP,  render_contextH, layerH,  lightP);
		}
		if (!err){
			lightP->num_lightsL++;
		}
	}
	if (lightP){
		ERR(suites.MemorySuite1()->AEGP_UnlockMemHandle(*lightsPH));
	}	
	if (lightsPH){	
		ERR2(Artie_DisposeLights(in_dataP,  *lightsPH));
	}
	return err;
}

static A_Err
Artie_GetPolygonTexture(
	const PR_InData			*in_dataP,
	PR_RenderContextH		render_contextH,
	AEGP_RenderHints		render_hints,
	A_Boolean				*is_visiblePB,
	Artie_Poly				*polygonP)
{
	A_Err					err	= A_Err_NONE;
	A_FpLong				widthF = 0.0, heightF = 0.0;
	AEGP_DownsampleFactor	dsf = {1, 1};
	AEGP_SuiteHandler suites(in_dataP->pica_basicP);	

	*is_visiblePB = FALSE;
	
	if (!polygonP->texture){
		// no texture map yet

		// if we're in wireframe, there is no texture.
		// we still need the dimensions, so we'll get the layers source dimensions and use that.
		// we'll also need to correct for the comp down sampling as all textures have this correction

		if ( polygonP->aegp_quality == AEGP_LayerQual_WIREFRAME){
		
			ERR(Artie_GetSourceLayerSize(	in_dataP, 
											render_contextH, 
											polygonP->layerH, 
											&widthF, 
											&heightF, 
											&dsf));
			if (!err){
				polygonP->origin.x			= 0;
				polygonP->origin.y			= 0;
			}
				
		} else {
			ERR(suites.CanvasSuite5()->AEGP_RenderTexture(	render_contextH,
															polygonP->layer_contextH,
															AEGP_RenderHints_NONE,
															NULL, 
															NULL,
															NULL,
															&polygonP->texture));

		}
		if (!err && polygonP->texture) {
			*is_visiblePB = TRUE;
		}
		if (!err && *is_visiblePB){
			A_FpLong	widthL 	= 0, 
						heightL = 0;
			
			err = Artie_GetSourceLayerSize(	in_dataP, 
											render_contextH, 
											polygonP->layerH, 
											&widthL, 
											&heightL, 
											&dsf);
			if (!err){
				polygonP->origin.x /= (A_FpLong)dsf.xS;	
				polygonP->origin.y /= (A_FpLong)dsf.yS;
			}
		}

		// construct the polygon vertices in local ( pixel) space.
		if (!err && *is_visiblePB) 	{
		
			A_long	widthL = 0, heightL = 0;
			
			ERR(suites.WorldSuite2()->AEGP_GetSize(polygonP->texture, &widthL, &heightL));
			
			widthF 	= static_cast<A_FpLong>(widthL);
			heightF = static_cast<A_FpLong>(heightL);

			// counterclockwise specification -- for texture map
			//vertex 0
			polygonP->vertsA[0].coord.P[0] = -polygonP->origin.x;				
			polygonP->vertsA[0].coord.P[1] = -polygonP->origin.y;   	
			polygonP->vertsA[0].coord.P[2] = 0;
			polygonP->vertsA[0].coord.P[3] = 1;

			polygonP->vertsA[0].txtur[0] = 0.0;
			polygonP->vertsA[0].txtur[1] = 0.0;
			polygonP->vertsA[0].txtur[2] = 1;


			// vertex 1
			polygonP->vertsA[1].coord.P[0] = -polygonP->origin.x;				
			polygonP->vertsA[1].coord.P[1] = heightF - polygonP->origin.y;				
			polygonP->vertsA[1].coord.P[2] = 0.0;
			polygonP->vertsA[1].coord.P[3] = 1;

			polygonP->vertsA[1].txtur[0] = 0.0;
			polygonP->vertsA[1].txtur[1] = 1.0;
			polygonP->vertsA[1].txtur[2] = 1;


			//vertex 2
			polygonP->vertsA[2].coord.P[0] = widthF  - polygonP->origin.x ;	
			polygonP->vertsA[2].coord.P[1] = heightF - polygonP->origin.y ;				
			polygonP->vertsA[2].coord.P[2] = 0;
			polygonP->vertsA[2].coord.P[3] = 1;

			polygonP->vertsA[2].txtur[0] = 1.0;
			polygonP->vertsA[2].txtur[1] = 1.0;
			polygonP->vertsA[2].txtur[2] = 1;

		
			// vertex 3
			polygonP->vertsA[3].coord.P[0] = (widthF  - polygonP->origin.x);	
			polygonP->vertsA[3].coord.P[1] = -polygonP->origin.y;
			polygonP->vertsA[3].coord.P[2] = 0;
			polygonP->vertsA[3].coord.P[3] = 1;


			polygonP->vertsA[3].txtur[0] = 1.0;
			polygonP->vertsA[3].txtur[1] = 0.0;
			polygonP->vertsA[3].txtur[2] = 1;

		}
	}
	return err;
}

static A_Err
Artie_BuildPolygon(
	const PR_InData				*in_dataP,
	PR_RenderContextH			render_contextH,
	AEGP_RenderLayerContextH	layer_contextH,
	AEGP_LayerH					layerH,
	A_Matrix4					*xform,
	A_FpLong					widthF,
	A_FpLong					heightF,
	PF_CompositeMode			*composite_mode,
	Artie_Material				*material,
	AEGP_TrackMatte				*track_matte,
	AEGP_LayerQuality			*aegp_quality,
	Artie_Poly					*polygonP)			/* <> */
{
	A_Err			err	   = A_Err_NONE;
	AEGP_SuiteHandler suites(in_dataP->pica_basicP);	
	
	AEFX_CLR_STRUCT(*polygonP);

	polygonP->aegp_quality		= *aegp_quality;
	polygonP->layerH			= layerH;
	polygonP->layer_contextH	= layer_contextH;
	polygonP->normal.V[X]		= 0;
	polygonP->normal.V[Y]		= 0;
	polygonP->normal.V[Z]		= -1;
	polygonP->normal.V[W]		= 0;
	polygonP->material			= *material;	
	polygonP->world_matrix		= *xform;

	
	// fill in vertices and texture map
	err = Artie_GetPolygonTexture(	in_dataP,  
									render_contextH,
									AEGP_RenderHints_NONE,
									&polygonP->is_visibleB, 
									polygonP);
	return err;
}

static A_Err
Artie_AddPolygonToScene(
	const PR_InData				*in_dataP,
	PR_RenderContextH			render_contextH,
	A_long						indexL,
	AEGP_RenderLayerContextH	layer_contextH,
	AEGP_LayerH					layerH, 
	AEGP_MemHandle				sceneH)
{
	
	A_Err				err				= A_Err_NONE, 
						err2			= A_Err_NONE;
	Artie_Scene			*sceneP			= NULL;
	PF_CompositeMode	composite_mode;
	AEGP_TrackMatte		track_matte;
	A_Matrix4			xform;
	AEGP_LayerQuality	aegp_quality;
	AEGP_LayerFlags		layer_flags;
	AEGP_LayerTransferMode layer_transfer_mode;
	Artie_Poly			poly;
	A_FpLong			opacityF		= 0.0,
						widthF 			= 0.0, 
						heightF 		= 0.0;
	A_long				seedL 			= 0;
	Artie_Material		material;
	A_Time				comp_time		= {0,1},
						comp_time_step	= {0,1};
						
	AEGP_StreamVal		stream_val;
	AEGP_SuiteHandler 	suites(in_dataP->pica_basicP);	

	ERR(suites.MemorySuite1()->AEGP_LockMemHandle(sceneH, reinterpret_cast<void**>(&sceneP)));
	ERR(suites.LayerSuite5()->AEGP_GetLayerFlags(layerH, &layer_flags));
	ERR(suites.LayerSuite5()->AEGP_GetLayerTransferMode(layerH, &layer_transfer_mode));
	ERR(suites.LayerSuite5()->AEGP_GetLayerQuality(layerH, &aegp_quality));
	ERR(suites.CanvasSuite5()->AEGP_GetCompRenderTime(render_contextH, &comp_time, &comp_time_step));
	ERR(suites.LayerSuite5()->AEGP_GetLayerToWorldXform(layerH, &comp_time, &xform));
	ERR(suites.CanvasSuite5()->AEGP_GetCompRenderTime(render_contextH, &comp_time, &comp_time_step));
	ERR(suites.CanvasSuite5()->AEGP_GetRenderOpacity(render_contextH, layer_contextH, &comp_time, &opacityF));
	ERR(suites.LayerSuite5()->AEGP_GetLayerDancingRandValue(layerH, &comp_time, &seedL));
	if (!err) {
		composite_mode.xfer			= layer_transfer_mode.mode;
		composite_mode.rand_seed	= seedL;
		composite_mode.opacity		= (unsigned char) (255.0 * opacityF / 100.0 + 0.5);
		composite_mode.rgb_only		= (layer_transfer_mode.flags & AEGP_TransferFlag_PRESERVE_ALPHA) != 0;
		track_matte					= layer_transfer_mode.track_matte;
	}
	// get polygon material
	ERR(suites.StreamSuite2()->AEGP_GetLayerStreamValue(layerH, 
														AEGP_LayerStream_AMBIENT_COEFF, 
														AEGP_LTimeMode_CompTime,
														&comp_time,	
														FALSE,	
														&stream_val, 
														NULL));
	if (!err){
		material.kaF = stream_val.one_d;
		err = suites.StreamSuite2()->AEGP_GetLayerStreamValue(	layerH, 
																AEGP_LayerStream_DIFFUSE_COEFF, 
																AEGP_LTimeMode_CompTime,
																&comp_time,	
																FALSE,	
																&stream_val, 
																NULL);
		material.kdF = stream_val.one_d;

	}
	ERR(suites.StreamSuite2()->AEGP_GetLayerStreamValue(	layerH, 
															AEGP_LayerStream_SPECULAR_INTENSITY,	
															AEGP_LTimeMode_CompTime,
															&comp_time,	
															FALSE,	
															&stream_val, 
															NULL));

	// Normalize coeffs
	if (!err) {
		material.ksF = stream_val.one_d;
		material.kaF *= 0.01;	
		material.kdF *= 0.01;	
		material.ksF *= 0.01;	
	}
	ERR(suites.StreamSuite2()->AEGP_GetLayerStreamValue(layerH, 
														AEGP_LayerStream_SPECULAR_SHININESS,
														AEGP_LTimeMode_CompTime,
														&comp_time,	
														FALSE,	
														&stream_val, 
														NULL));

	if (!err) {
		material.specularF = stream_val.one_d;

		AEGP_DownsampleFactor	dsf	= {1,1};
		
		err = Artie_GetSourceLayerSize(	in_dataP, 
										render_contextH, 
										layerH, 
										&widthF, 
										&heightF, 
										&dsf);
	}

	ERR(Artie_BuildPolygon(	in_dataP, 
							render_contextH,
							layer_contextH,
							layerH,  
							&xform, 
							widthF, 
							heightF, 
							&composite_mode, 
							&material, 
							&track_matte, 
							&aegp_quality,  
							&poly));
	if (!err) {
		sceneP->polygons[indexL] = poly;
		sceneP->num_polysL++;
	}
	if (sceneP) {
		ERR2(suites.MemorySuite1()->AEGP_UnlockMemHandle(sceneH));
	}
	return err;
}

static A_Err
Artie_BuildScene(
	const PR_InData				*in_dataP,	
	Artie_LayerContexts			layer_contexts,
	PR_RenderContextH			render_contextH,
	AEGP_MemHandle				*scenePH)
{
	A_Err			err					= A_Err_NONE, 
					err2				= A_Err_NONE;

	Artie_Scene		*sceneP				= NULL;
	AEGP_LayerH		layerH				= NULL;
	AEGP_CompH		compH				= NULL;
	AEGP_ItemH		source_itemH		= NULL;
	AEGP_ItemFlags	item_flags			= 0;
	A_long			layers_to_renderL	= 0;

	AEGP_MemHandle	lightsH				= NULL;
	
	AEGP_SuiteHandler	suites(in_dataP->pica_basicP);

	ERR(suites.MemorySuite1()->AEGP_NewMemHandle(	in_dataP->aegp_plug_id, 
													"Scene Data", 
													sizeof(Artie_Scene), 
													AEGP_MemFlag_CLEAR,  
													scenePH)); 

	ERR(suites.CanvasSuite5()->AEGP_GetCompToRender(render_contextH, &compH));

	ERR(Artie_CreateListOfLayerContexts(in_dataP, 
										render_contextH,
										&layers_to_renderL, 
										&layer_contexts));

	ERR(Artie_BuildLights(in_dataP, render_contextH, &lightsH));
	ERR(suites.MemorySuite1()->AEGP_LockMemHandle(*scenePH, reinterpret_cast<void**>(&sceneP)));

	for(A_long iL = 0; iL < layers_to_renderL; iL++){			
		ERR(suites.CanvasSuite5()->AEGP_GetLayerFromLayerContext(	render_contextH, 
																	layer_contexts.layer_contexts[iL], 
																	&layerH));
			
		ERR(suites.LayerSuite5()->AEGP_GetLayerSourceItem(layerH, &source_itemH));
		ERR(suites.ItemSuite6()->AEGP_GetItemFlags(source_itemH, &item_flags));
	
		if (item_flags & AEGP_ItemFlag_HAS_VIDEO){ 
			ERR(Artie_AddPolygonToScene(	in_dataP, 
											render_contextH,
											iL,
											layer_contexts.layer_contexts[iL], 
											layerH, 
											*scenePH));
		}
	}
	ERR2(suites.MemorySuite1()->AEGP_UnlockMemHandle(*scenePH));
	return err;
}

static void
copyConvertStringLiteralIntoUTF16(
	const wchar_t* inputString,
	A_UTF16Char* destination)
{
#ifdef AE_OS_MAC
	int length = wcslen(inputString);
	CFRange	range = {0, 256};
	range.length = length;
	CFStringRef inputStringCFSR = CFStringCreateWithBytes(kCFAllocatorDefault,
														  reinterpret_cast<const UInt8 *>(inputString),
														  length * sizeof(wchar_t),
														  kCFStringEncodingUTF32LE,
														  false);
	CFStringGetBytes(inputStringCFSR,
					 range,
					 kCFStringEncodingUTF16,
					 0,
					 false,
					 reinterpret_cast<UInt8 *>(destination),
					 length * (sizeof (A_UTF16Char)),
					 NULL);
	destination[length] = 0; // Set NULL-terminator, since CFString calls don't set it
	CFRelease(inputStringCFSR);
#elif defined AE_OS_WIN
	size_t length = wcslen(inputString);
	wcscpy_s(reinterpret_cast<wchar_t*>(destination), length + 1, inputString);
#endif
}

static void
copyConvertStringLiteralInto(
	const wchar_t* inputString,
	A_char* destination)
{
	size_t length = wcslen(inputString);
	wcscpy_s(reinterpret_cast<wchar_t*>(destination), length + 1, inputString);
}

static A_Err 
	AddNewImageToProject(
	AEGP_SuiteHandler	suites,
	AEGP_ProjectH		projH)
{
	A_Err 				err		= A_Err_NONE;

	AEGP_ItemH root_itemH = NULL;
	ERR(suites.ProjSuite5()->AEGP_GetProjectRootFolder(projH, &root_itemH));

	AEGP_ItemH new_folderH = NULL;
	A_UTF16Char	folder_name[256];
	copyConvertStringLiteralIntoUTF16(FOLDER_NAME, folder_name);
	ERR(suites.ItemSuite8()->AEGP_CreateNewFolder(folder_name,
						root_itemH,
						&new_folderH));

	if (!err && *new_folderH) {

		AEGP_FootageLayerKey	key1, key2;
	
		AEFX_CLR_STRUCT(key1);
		AEFX_CLR_STRUCT(key2);
	
		key1.layer_idL =
		key2.layer_idL = AEGP_LayerID_UNKNOWN;
	
		key1.layer_indexL = 0;
		key2.layer_indexL = 1;
	
		A_UTF16Char	psd_path[256];
		copyConvertStringLiteralIntoUTF16(THE_PHOTO_PATH, psd_path);

		AEGP_FootageH footageH = NULL;
		ERR(suites.FootageSuite5()->AEGP_NewFootage(S_my_id, 
														psd_path,
														&key1,
														NULL,
														FALSE, 
														NULL, 
														&footageH));																	

		AEGP_ItemH layer1_itemH = NULL;
		ERR(suites.FootageSuite5()->AEGP_AddFootageToProject(footageH, new_folderH, &layer1_itemH));									

		A_char scriptAC[AEGP_MAX_MARKER_URL_SIZE] = {"var newSource = null; \
														for (var i = 1; i < app.project.items.length; i++) { \
															var theComposition = app.project.item(i); \
															if (theComposition.name == 'logo.png') { \
																newSource = theComposition; \
															} \
														} \
														if (newSource) { \
															for (var i = 1; i < app.project.items.length; i++) { \
																var theComposition = app.project.item(i); \
																if (theComposition.name == '03 Import Your logo') { \
																	var layer = theComposition.layer(1); \
																	layer.replaceSource(newSource,  true) \
																} \
															} \
														}\0"};

		A_Boolean encodingB = TRUE; 
		AEGP_MemHandle resultMemH = NULL, errorMemH = NULL;									

		ERR(suites.UtilitySuite5()->AEGP_ExecuteScript(S_my_id, scriptAC, encodingB, &resultMemH, &errorMemH));
	}

	return err;
}

static A_Err 
	SetTextToLayer(
		AEGP_SuiteHandler	suites,
		AEGP_PluginID		S_my_id,
		AEGP_LayerH			layerH,		
		AEGP_StreamRefH		streamH,
		const A_u_short		unicode[],
		A_long lengthL
	) 
{

	A_Err 				err		= A_Err_NONE;

	AEGP_StreamValue	value;
	AEFX_CLR_STRUCT(value);
								
	AEGP_StreamType		stream_type				= AEGP_StreamType_NO_DATA;
								
	A_Time				timeT					= {0,1};
	AEGP_LTimeMode		time_mode				= AEGP_LTimeMode_LayerTime;

	ERR(suites.DynamicStreamSuite2()->AEGP_GetNewStreamRefForLayer(S_my_id, layerH, &streamH));
	ERR(suites.StreamSuite2()->AEGP_GetNewLayerStream(S_my_id, 
										layerH,
										AEGP_LayerStream_SOURCE_TEXT,
										&streamH));

	ERR(suites.StreamSuite2()->AEGP_GetStreamType(streamH, &stream_type));

	if (!err && (stream_type == AEGP_StreamType_TEXT_DOCUMENT)){
		ERR(suites.LayerSuite5()->AEGP_GetLayerCurrentTime(layerH, time_mode, &timeT));
		ERR(suites.StreamSuite2()->AEGP_GetNewStreamValue(	S_my_id,
															streamH,
															AEGP_LTimeMode_LayerTime,
															&timeT,
															TRUE,
															&value));																									
				
		//const A_u_short unicode[] = {0x0042, 0x0042, 0x0042};  // "BBB"
		//const A_u_short unicode[] = {0x0048, 0x0065, 0x006c, 0x006c, 0x006f, 0x0020, 0x0021 }; // "Hello !"

		//const A_u_short unicode[] = {0x0042, 0x0065, 0x0020, 0x0042, 0x0065, 0x0020, 0x0042, 0x0065, 0x0020, 0x0021 }; // "Be be be !"
		//const A_u_short unicode[] = { 0x0053, 0x0061, 0x006d, 0x0070, 0x006c, 0x0065, 0x0020, 0x0074, 0x0065, 0x0078, 0x0074 };
		//A_long	lengthL = sizeof(unicode) / sizeof(A_u_short);				
		//A_long	lengthL = sizeof(unicode) + 1;		

		ERR(suites.TextDocumentSuite1()->AEGP_SetText(value.val.text_documentH, unicode, lengthL));
		ERR(suites.StreamSuite2()->AEGP_SetStreamValue(S_my_id, streamH, &value));									
	}

	return err;
}

static A_Err
RecursiveDump(
	FILE 				*out,
	A_long				indent_levelsL,
	AEGP_ItemH			itemH,
	AEGP_ItemH			parent_itemH0,
	AEGP_ItemH			*prev_itemPH)
{
	A_Err 				err		= A_Err_NONE, 
						err2	= A_Err_NONE;

	A_char				item_nameAC[AEGP_MAX_ITEM_NAME_SIZE]	= {'\0'},
						type_nameAC[AEGP_MAX_ITEM_NAME_SIZE]	= {'\0'},
						layer_nameAC[AEGP_MAX_ITEM_NAME_SIZE]	= {'\0'},
						effect_nameAC[AEGP_MAX_ITEM_NAME_SIZE]	= {'\0'},
						stream_nameAC[AEGP_MAX_ITEM_NAME_SIZE]	= {'\0'},
						indent_stringAC[128];
	
	AEGP_ItemH			curr_parent_itemH = parent_itemH0;	
	
	A_long				iL	=	0,
						jL	=	0, 
						kL	=	0;
	
	AEGP_EffectRefH		effectH = NULL;

	AEGP_SuiteHandler	suites(sP);	
	
	indent_stringAC[0] = '\0';
	
	for (iL = 0; !err && iL < indent_levelsL; iL++) {
		suites.ANSICallbacksSuite1()->sprintf(indent_stringAC, "%s\t", indent_stringAC);
	}

	while (	!err	&& 
			(itemH) &&
			((parent_itemH0 == 0) || (parent_itemH0 == curr_parent_itemH))) {
			
		*prev_itemPH = itemH;
		
		AEGP_ProjectH	projH	=	NULL;
		
		ERR(suites.ProjSuite5()->AEGP_GetProjectByIndex(0, &projH));
		ERR(suites.ItemSuite6()->AEGP_GetNextProjItem(projH, itemH, &itemH));

		if (itemH) {
			AEGP_ItemType		item_type;
			
			ERR(suites.ItemSuite6()->AEGP_GetItemParentFolder(itemH, &curr_parent_itemH));

			if ((parent_itemH0 == 0) || (parent_itemH0 == curr_parent_itemH)) {
				
				ERR(suites.ItemSuite6()->AEGP_GetItemType(itemH, &item_type));
				ERR(suites.ItemSuite6()->AEGP_GetTypeName(item_type, type_nameAC));
				ERR(suites.ItemSuite6()->AEGP_GetItemName(itemH, item_nameAC));

				if (!err && (item_type != AEGP_ItemType_FOLDER)) {
					A_Time 				durationT, currT;
					A_long				widthL, heightL;
					AEGP_MemHandle		pathH;

					ERR(suites.ItemSuite6()->AEGP_GetItemDuration(itemH, &durationT));
					ERR(suites.ItemSuite6()->AEGP_GetItemCurrentTime(itemH, &currT));
					ERR(suites.ItemSuite6()->AEGP_GetItemDimensions(itemH, &widthL, &heightL));

					fprintf(out, "%s%s: %s Duration: %1.2f Current: %1.2f Width: %d Height: %d ", indent_stringAC,
							type_nameAC, item_nameAC, (A_FpLong)(durationT.value) / (durationT.scale), (A_FpLong)(currT.value) / (currT.scale),
							widthL, heightL);
							
					if (!err && (item_type == AEGP_ItemType_FOOTAGE)) {
						AEGP_FootageH			footageH;
						AEGP_FootageInterp		interp;
						A_long					num_footage_filesL;
						
						ERR(suites.FootageSuite5()->AEGP_GetMainFootageFromItem(itemH, &footageH));
						ERR(suites.FootageSuite5()->AEGP_GetFootageInterpretation(itemH, false, &interp));
						ERR(suites.FootageSuite5()->AEGP_GetFootageNumFiles(footageH, &num_footage_filesL, 0));
						ERR(suites.FootageSuite5()->AEGP_GetFootagePath(footageH, 0, AEGP_FOOTAGE_MAIN_FILE_INDEX, &pathH));

						// Set background
						if (strcmp(item_nameAC, layerToChangeBackground) == 0) {													
							AEGP_ColorVal color = {0.0, 1.0, 1.0, 0.0};
							ERR(suites.FootageSuite5()->AEGP_SetSolidFootageColor(itemH, false, &color));							
						}

						if (!err) {
							fprintf(out, "Num Files: %d\n", num_footage_filesL);
							ERR(suites.MemorySuite1()->AEGP_FreeMemHandle(pathH));
						}
					} else if (item_type == AEGP_ItemType_COMP) {
						A_Time 			lay_inT, lay_durT, comp_inT, comp_durT;
						AEGP_CompH		compH;
						AEGP_LayerH		layerH;
						A_long			num_layersL, num_effectsL, num_paramsL;
						
						ERR(suites.CompSuite4()->AEGP_GetCompFromItem(itemH, &compH));
						ERR(suites.LayerSuite5()->AEGP_GetCompNumLayers(compH, &num_layersL));						

						fprintf(out, "Num Layers: %d\n", num_layersL);

						for (iL = 0; !err && iL < num_layersL; iL++) {
						
							AEGP_StreamRefH		streamH = NULL;

							ERR(suites.LayerSuite5()->AEGP_GetCompLayerByIndex(compH, iL, &layerH));
							ERR(suites.LayerSuite5()->AEGP_GetLayerName(layerH, layer_nameAC, 0));
							ERR(suites.LayerSuite5()->AEGP_GetLayerInPoint(layerH, AEGP_LTimeMode_LayerTime, &lay_inT));
							ERR(suites.LayerSuite5()->AEGP_GetLayerDuration(layerH, AEGP_LTimeMode_LayerTime, &lay_durT));
							ERR(suites.LayerSuite5()->AEGP_GetLayerInPoint(layerH, AEGP_LTimeMode_CompTime, &comp_inT));
							ERR(suites.LayerSuite5()->AEGP_GetLayerDuration(layerH, AEGP_LTimeMode_CompTime, &comp_durT));
							ERR(suites.LayerSuite5()->AEGP_GetLayerCurrentTime(layerH, AEGP_LTimeMode_LayerTime, &currT));
							ERR(suites.EffectSuite2()->AEGP_GetLayerNumEffects(layerH, &num_effectsL));														

							fprintf(out, "%s\t%d: %s LayIn: %1.2f Layer Duration: %1.2f CompIn: %1.2f CompDur: %1.2f Current: %1.2f Num Effects: %d\n", indent_stringAC, iL + 1, layer_nameAC,
								(A_FpLong)(lay_inT.value) / (lay_inT.scale), (A_FpLong)(lay_durT.value) / (lay_durT.scale),
								(A_FpLong)(comp_inT.value) / (comp_inT.scale), (A_FpLong)(comp_durT.value) / (comp_durT.scale),
								(A_FpLong)(currT.value) / (currT.scale), num_effectsL);

							ERR(suites.StreamSuite2()->AEGP_GetNewLayerStream(S_my_id, layerH, AEGP_LayerStream_ANCHORPOINT, &streamH));

							ERR(PrintAndDisposeStream(streamH, "%s\t\t%s", indent_stringAC, stream_nameAC, out));
							
							ERR(suites.StreamSuite2()->AEGP_GetNewLayerStream(S_my_id, layerH, AEGP_LayerStream_POSITION, &streamH));
							
							ERR(PrintAndDisposeStream(streamH, "%s\t\t%s", indent_stringAC, stream_nameAC, out));
							
							ERR(suites.StreamSuite2()->AEGP_GetNewLayerStream(S_my_id, layerH, AEGP_LayerStream_SCALE, &streamH));
							
							ERR(PrintAndDisposeStream(streamH, "%s\t\t%s", indent_stringAC, stream_nameAC, out));
							
							ERR(suites.StreamSuite2()->AEGP_GetNewLayerStream(S_my_id, layerH, AEGP_LayerStream_ROTATION, &streamH));

							ERR(PrintAndDisposeStream(streamH, "%s\t\t%s", indent_stringAC, stream_nameAC, out));
							
							ERR(suites.StreamSuite2()->AEGP_GetNewLayerStream(S_my_id, layerH, AEGP_LayerStream_OPACITY, &streamH));
							
							ERR(PrintAndDisposeStream(streamH, "%s\t\t%s", indent_stringAC, stream_nameAC, out));
							
							ERR(suites.StreamSuite2()->AEGP_GetNewLayerStream(S_my_id, layerH, AEGP_LayerStream_AUDIO, &streamH));

							ERR(PrintAndDisposeStream(streamH, "%s\t\t%s", indent_stringAC, stream_nameAC, out));
							
							ERR(suites.StreamSuite2()->AEGP_GetNewLayerStream(S_my_id, layerH, AEGP_LayerStream_MARKER, &streamH));
							
							ERR(PrintAndDisposeStream(streamH, "%s\t\t%s", indent_stringAC, stream_nameAC, out));
							
							// Add layer to render
							if (strcmp(item_nameAC, layerToRender) == 0) {							
								AEGP_RenderQueueState current_state  = AEGP_RenderQueueState_RENDERING;
																	
								ERR(suites.RenderQueueSuite1()->AEGP_AddCompToRenderQueue(compH, pathToRenderVideo));																			
								ERR(suites.RenderQueueSuite1()->AEGP_GetRenderQueueState(&current_state));
							}

							//TODO:
							if (strcmp(item_nameAC, layerWithImage) == 0) {
								ERR(AddNewImageToProject(suites, projH));
							}
	
							// Set text to first layer
							if (strcmp(item_nameAC, layerFirstText) == 0) {
								const A_u_short unicode[] = { 0x0048, 0x0065, 0x006c, 0x006c, 0x006f, 0x0020, 0x0021 };
								A_long	lengthL = sizeof(unicode) / sizeof(A_u_short);

								ERR(SetTextToLayer(suites, S_my_id, layerH, streamH, unicode, lengthL)); 

								// Set background
								AEGP_ColorVal color = {1.0, 1.0, 0.0, 0.0};								
								ERR(suites.CompSuite9()->AEGP_SetCompBGColor(compH, &color));
							}

							if (strcmp(item_nameAC, layerSecondText) == 0) {	

								const A_u_short unicode[] = { 0x0053, 0x0061, 0x006d, 0x0070, 0x006c, 0x0065 };
								A_long	lengthL = sizeof(unicode) / sizeof(A_u_short);
								ERR(SetTextToLayer(suites, S_my_id, layerH, streamH, unicode, lengthL)); 
							}
							
							if (strcmp(item_nameAC, layerLastText) == 0) {		
								
								const A_u_short unicode[] = { 0x0077, 0x0077, 0x0077, 0x002e, 0x0079, 0x0061, 0x002e, 0x0072, 0x0075 };
								A_long	lengthL = sizeof(unicode) / sizeof(A_u_short);
								ERR(SetTextToLayer(suites, S_my_id, layerH, streamH, unicode, lengthL)); 								
							}							
							
							for (jL = 0; !err && jL < num_effectsL; jL++) {									
								ERR(suites.EffectSuite2()->AEGP_GetLayerEffectByIndex(S_my_id, layerH, jL, &effectH));

								if (!err && effectH) {
									AEGP_InstalledEffectKey	key;
									
									ERR(suites.EffectSuite2()->AEGP_GetInstalledKeyFromLayerEffect(effectH, &key));
									ERR(suites.EffectSuite2()->AEGP_GetEffectName(key, effect_nameAC));									
									
									if (!strcmp(effect_nameAC, "Shifter")) {
										A_Time					fake_timeT = {0, 100};

										err = suites.EffectSuite2()->AEGP_EffectCallGeneric(S_my_id, effectH, &fake_timeT, (void*)"MoveMaker Shifter");
									}
								}
								fprintf(out, "%s\t\t%d: %s\n", indent_stringAC, jL + 1, effect_nameAC);
								
								ERR(suites.StreamSuite2()->AEGP_GetEffectNumParamStreams(effectH, &num_paramsL));
								
								for (kL = 1; !err && kL < num_paramsL; ++kL) {		// start at 1 to skip initial layer param
									
									ERR(suites.StreamSuite2()->AEGP_GetNewEffectStreamByIndex(S_my_id, effectH, kL, &streamH));
									ERR(PrintAndDisposeStream(streamH, "%s\t\t\t%s", indent_stringAC, stream_nameAC, out));
								}
								ERR2(suites.EffectSuite2()->AEGP_DisposeEffect(effectH));
							}
						}
					}

				} else {
					fprintf(out, "%s%s: %s\n", indent_stringAC, type_nameAC, item_nameAC);
					
					err = RecursiveDump(out, 
										indent_levelsL + 1, 
										itemH, 
										itemH, 
										prev_itemPH);

					itemH = *prev_itemPH;
				}
			}
		}
	}
	return err;
}


static A_Err
DumpProj(void)
{
	A_Err 				err 	= A_Err_NONE, 
						err2 	= A_Err_NONE;
	FILE 				*out 	= NULL;
	AEGP_ProjectH		projH	= NULL;
	AEGP_ItemH			itemH 	= NULL;
	A_char				proj_nameAC[AEGP_MAX_PROJ_NAME_SIZE];
	A_char				path_nameAC[AEGP_MAX_PROJ_NAME_SIZE + 16];		// so we can add "_Dump.txt"
	AEGP_SuiteHandler	suites(sP);

	suites.ANSICallbacksSuite1()->strcpy(path_nameAC, workingPath);
	
	// Open template project from path pathToTemplateProejct
	ERR(suites.ProjSuite5()->AEGP_OpenProjectFromPath(pathToTemplateProject, &projH));

	// Save new project from template
	ERR(suites.ProjSuite5()->AEGP_SaveProjectToPath(projH, pathToWorkingProject));

	// Reopen new project
	ERR(suites.ProjSuite5()->AEGP_OpenProjectFromPath(pathToWorkingProject, &projH));

	ERR(suites.UtilitySuite3()->AEGP_StartUndoGroup("Keepin' On"));

	if (!err) {
		ERR(suites.ProjSuite5()->AEGP_GetProjectByIndex(0, &projH));
		ERR(suites.ProjSuite5()->AEGP_GetProjectRootFolder(projH, &itemH));

		if (!err && itemH) {
			err = suites.ProjSuite5()->AEGP_GetProjectName(projH, proj_nameAC);
			strcat(path_nameAC, proj_nameAC);
			strcat(path_nameAC, "_Dump.txt");
			if (!err) {
				out = fopen(path_nameAC, "w");
				if (out) {
					AEGP_ItemH			prev_itemH = 0;
			
					err = RecursiveDump(out, 0, itemH, 0, &prev_itemH);

					fprintf(out, "\n\n");
					fclose(out);
				}
			}

			// Save project			
			ERR(suites.ProjSuite5()->AEGP_SaveProjectAs(projH, pathToWorkingProject));

			// Go to render
			ERR(suites.RenderQueueSuite1()->AEGP_SetRenderQueueState(AEGP_RenderQueueState_RENDERING));
		}
	}
	ERR2(suites.UtilitySuite3()->AEGP_EndUndoGroup());

	if (!err)
	{
//#ifdef AE_OS_WIN
//		suites.UtilitySuite5()->AEGP_ReportInfo(S_my_id, DUMP_SUCCEEDED_WIN);
//#elif defined AE_OS_MAC
//		suites.UtilitySuite5()->AEGP_ReportInfo(S_my_id, DUMP_SUCCEEDED_MAC);
//#endif
	}
	else
	{
		suites.UtilitySuite5()->AEGP_ReportInfo(S_my_id, DUMP_FAILED);
	}

	return err;
}

static A_Err
MakeCuter(void)
{
	A_Err 				err 	= A_Err_NONE;
	AEGP_SuiteHandler	suites(sP);
	AEGP_ItemH			active_itemH;
	ERR(suites.ItemSuite6()->AEGP_GetActiveItem(&active_itemH));
	
	AEGP_ItemFlags	item_flags;
	char			info_string_1AC[256];
	char			info_string_2AC[256];
	char			active_item_nameAC[AEGP_MAX_ITEM_NAME_SIZE] = "ain't got nothin";
	
	if (active_itemH){					
		ERR(suites.ItemSuite6()->AEGP_GetItemName(active_itemH, active_item_nameAC));
		if (!err) {
			suites.ANSICallbacksSuite1()->sprintf(info_string_1AC, "Active Item: %s", active_item_nameAC);
		}
		ERR(suites.ItemSuite6()->AEGP_GetItemFlags(active_itemH, &item_flags));
		if (!err) {
			suites.ANSICallbacksSuite1()->sprintf(	info_string_2AC, "m: %s hp: %s up: %s mp: %s v: %s a: %s s: %s",
										(item_flags & AEGP_ItemFlag_MISSING) ? "Y" : "N",
										(item_flags & AEGP_ItemFlag_HAS_PROXY) ? "Y" : "N",
										(item_flags & AEGP_ItemFlag_USING_PROXY) ? "Y" : "N",
										(item_flags & AEGP_ItemFlag_MISSING_PROXY) ? "Y" : "N",
										(item_flags & AEGP_ItemFlag_HAS_VIDEO) ? "Y" : "N",
										(item_flags & AEGP_ItemFlag_HAS_AUDIO) ? "Y" : "N",
										(item_flags & AEGP_ItemFlag_STILL) ? "Y" : "N");						
		}
	}
			
	/*	
		the code below creates footage, adds the footage to project (returns item) and then
		adds the item to the active item if it's a composition.
	*/
	if (!err) {
		AEGP_ProjectH		projH;
		AEGP_FootageH		footageH, proxy_footageH;
		AEGP_ItemH			footage_itemH, root_itemH;
		AEGP_LayerH			added_layerH;
		A_UTF16Char			layer_pathZ[AEGP_MAX_PATH_SIZE],
							proxy_pathZ[AEGP_MAX_PATH_SIZE];

		copyConvertStringLiteralIntoUTF16(LAYERED_PATH, layer_pathZ);
		copyConvertStringLiteralIntoUTF16(PROXY_FOOTAGE_PATH, proxy_pathZ);

		ERR(suites.FootageSuite5()->AEGP_NewFootage(	S_my_id, 
														layer_pathZ,
														NULL,
														NULL,
														FALSE,
														NULL,  
														&footageH));

		if (err)
		{
			suites.UtilitySuite5()->AEGP_ReportInfo(S_my_id, ERROR_MISSING_FOOTAGE);
		}

		ERR(suites.ProjSuite5()->AEGP_GetProjectByIndex(0, &projH));
		ERR(suites.ProjSuite5()->AEGP_GetProjectRootFolder(projH, &root_itemH));
		if (!err && root_itemH) {
			err = suites.FootageSuite5()->AEGP_AddFootageToProject(footageH, root_itemH, &footage_itemH);
		}
		if (!err && footage_itemH) {
			err = suites.FootageSuite5()->AEGP_NewFootage(	S_my_id, 
															proxy_pathZ,
															NULL,
															NULL,
															FALSE,  
															NULL,
															&proxy_footageH);

			ERR(suites.FootageSuite5()->AEGP_SetItemProxyFootage(proxy_footageH, footage_itemH));
			
			ERR(suites.ItemSuite6()->AEGP_SetItemUseProxy(footage_itemH, TRUE));

			if (!err) {
				AEGP_CompH		compH;
				A_Boolean		is_add_validB = FALSE;
				AEGP_ItemType	item_type;

				ERR(suites.ItemSuite6()->AEGP_GetItemType(active_itemH, &item_type));	
				
				if (item_type == AEGP_ItemType_COMP) {
					ERR(suites.CompSuite4()->AEGP_GetCompFromItem(active_itemH, &compH));
					ERR(suites.LayerSuite5()->AEGP_IsAddLayerValid(footage_itemH, compH, &is_add_validB));
				}
				
				if (is_add_validB) {
					ERR(suites.LayerSuite5()->AEGP_AddLayer(footage_itemH, compH, &added_layerH));
					ERR(suites.LayerSuite5()->AEGP_ReorderLayer(added_layerH, AEGP_REORDER_LAYER_TO_END));
					ERR(suites.LayerSuite5()->AEGP_SetLayerQuality(added_layerH, AEGP_LayerQual_BEST));
				}
			}
		}
		ERR(suites.CommandSuite1()->AEGP_SetMenuCommandName(S_other_cmd, "i changed my name"));
		ERR(suites.CommandSuite1()->AEGP_CheckMarkMenuCommand(S_other_cmd, TRUE));
	}
	return err;
}


static A_Err
UpdateMenuHook(
	AEGP_GlobalRefcon		plugin_refconPV,		/* >> */
	AEGP_UpdateMenuRefcon	refconPV,				/* >> */
	AEGP_WindowType			active_window)			/* >> */
{
	A_Err 				err = A_Err_NONE;
	AEGP_SuiteHandler	suites(sP);

	if (S_dump_proj_cmd) {
		err = suites.CommandSuite1()->AEGP_EnableCommand(S_dump_proj_cmd);
	}

	if (!err && S_other_cmd) {
		AEGP_ItemH		active_itemH;
		
		ERR(suites.ItemSuite6()->AEGP_GetActiveItem(&active_itemH));
		if (!err && active_itemH) {
			ERR(suites.CommandSuite1()->AEGP_EnableCommand(S_other_cmd));
			ERR(suites.CommandSuite1()->AEGP_CheckMarkMenuCommand(S_other_cmd, TRUE));
		}
	} else {
		ERR(suites.CommandSuite1()->AEGP_CheckMarkMenuCommand(S_other_cmd, FALSE));
	}
	return err;
}
						
static A_Err
CommandHook(
	AEGP_GlobalRefcon	plugin_refconPV,		/* >> */
	AEGP_CommandRefcon	refconPV,				/* >> */
	AEGP_Command		command,				/* >> */
	AEGP_HookPriority	hook_priority,			/* >> */
	A_Boolean			already_handledB,		/* >> */
	A_Boolean			*handledPB)				/* << */
{
	A_Err 				err = A_Err_NONE;
	AEGP_SuiteHandler	suites(sP);

	*handledPB = FALSE;

	if ((command == S_dump_proj_cmd) || (command == S_other_cmd)) {

		if (command == S_dump_proj_cmd) {
			err = DumpProj();
			*handledPB = TRUE;
		} else if (command == S_other_cmd) {
			err = MakeCuter();
			*handledPB = TRUE;
		}
	}
	return err;
}

A_Err
EntryPointFunc(
	struct SPBasicSuite		*pica_basicP,			/* >> */
	A_long				 	major_versionL,			/* >> */		
	A_long					minor_versionL,			/* >> */		
	AEGP_PluginID			aegp_plugin_id,			/* >> */
	AEGP_GlobalRefcon		*global_refconP)		/* << */
{
	A_Err 				err = A_Err_NONE;

	S_my_id				= aegp_plugin_id;

	sP 	= pica_basicP;
	
	AEGP_SuiteHandler	suites(pica_basicP);
	
	ERR(suites.CommandSuite1()->AEGP_GetUniqueCommand(&S_dump_proj_cmd));
	ERR(suites.CommandSuite1()->AEGP_InsertMenuCommand(S_dump_proj_cmd, "Custom Move Maker", AEGP_Menu_EDIT, AEGP_MENU_INSERT_SORTED));

	ERR(suites.CommandSuite1()->AEGP_GetUniqueCommand(&S_other_cmd));
	ERR(suites.CommandSuite1()->AEGP_InsertMenuCommand(S_other_cmd, "Custom Move Maker: Make Cuter", AEGP_Menu_EDIT, AEGP_MENU_INSERT_AT_BOTTOM));

	ERR(suites.RegisterSuite5()->AEGP_RegisterCommandHook(	S_my_id, 
															AEGP_HP_BeforeAE, 
															AEGP_Command_ALL, 
															CommandHook, 
															NULL));

	ERR(suites.RegisterSuite5()->AEGP_RegisterUpdateMenuHook(S_my_id, UpdateMenuHook, NULL));
	return err;
}
