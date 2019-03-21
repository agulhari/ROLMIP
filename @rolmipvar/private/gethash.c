#define NLHS 1
#define NRHS 3
#define pEXPONENT 0
#define pEXPTABLE 1
#define pJUMP 2
#define pINDEX 0 
#define pargEXPONENT prhs[pEXPONENT]
#define pargEXPTABLE prhs[pEXPTABLE]
#define pargJUMP prhs[pJUMP]
#define pargINDEX plhs[pINDEX]

#include "mex.h"

#define GETINDEX(type) for (ii = 0; ii < dimensionsEXPTABLEElement[0]; ++ii) {\
				ismatch = true;\
				for (jj = 0; jj < dimensionsEXPTABLEElement[1]; ++jj) {\
					if (exptableContsimplexNumeric##type[jj*(dimensionsEXPTABLEElement[0]) + ii] != exponentContsimplexNumeric##type[jj]) {\
						ismatch = false;\
						break;\
					}\
				}\
				if (ismatch) {\
					index = ((double)ii)*jumpContsimplexNumeric[countsimplex];\
					break;\
				}\
			}

static int csr_tocsc(const mwIndex n_row, const mwIndex n_col, const mwIndex *Ap, const mwIndex *Aj, const double *Ax, mwIndex *Bp, mwIndex *Bi, double *Bx);
static int csc_tocsr(const mwIndex n_row, const mwIndex n_col, const mwIndex *Ap, const mwIndex *Aj, const double *Ax, mwIndex *Bp, mwIndex *Bi, double *Bx);
   
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/*function index = gethash(exponent,exptable,jump)
		index = 0;
		for contsimplex=1:(length(exponent))
			if (length(exptable{contsimplex}) > 0)
				for ii = 1:size(exptable{contsimplex}, 1)
					if all(exptable{contsimplex}(ii, :) == exponent{contsimplex})
						index = index + (ii - 1)*jump(contsimplex);
						break;
					end
				end
			end
		end
		index = index + 1;
		return;
	end*/
	unsigned int index = 0;
	boolean_T ismatch = false, ismatchcolumns = false, ismatchvalue = false, isdouble = false, isfloat = false, isint8 = false, isint16 = false, isint32 = false, isint64 = false, isuint8 = false, isuint16 = false, isuint32 = false, isuint64 = false;
	double *indexptr = NULL;
	unsigned int countsimplex = 0, ii = 0, jj = 0, kk = 0, ll = 0;
	const mwSize *dimensionsEXPONENT = NULL, *dimensionsEXPTABLE = NULL, *dimensionsJUMP = NULL, *dimensionsEXPONENTElement = NULL, *dimensionsEXPTABLEElement = NULL;
	//mwSize nzmaxEXPONENT = 0, nzmaxEXPTABLE = 0;
	size_t lengthExponent = 0;
	const mxArray *exptableContsimplex = NULL, *exponentContsimplex = NULL;
	const double *exponentContsimplexNumericDouble = NULL, *exptableContsimplexNumericDouble = NULL, *jumpContsimplexNumeric = NULL;
	const float *exponentContsimplexNumericFloat = NULL, *exptableContsimplexNumericFloat = NULL;
	const char *exponentContsimplexNumericInt8 = NULL, *exptableContsimplexNumericInt8 = NULL;
	const short *exponentContsimplexNumericInt16 = NULL, *exptableContsimplexNumericInt16 = NULL;
	const int *exponentContsimplexNumericInt32 = NULL, *exptableContsimplexNumericInt32 = NULL;
	const long long *exponentContsimplexNumericInt64 = NULL, *exptableContsimplexNumericInt64 = NULL;
	const unsigned char *exponentContsimplexNumericUint8 = NULL, *exptableContsimplexNumericUint8 = NULL;
	const unsigned short *exponentContsimplexNumericUint16 = NULL, *exptableContsimplexNumericUint16 = NULL;
	const unsigned int *exponentContsimplexNumericUint32 = NULL, *exptableContsimplexNumericUint32 = NULL;
	const unsigned long long *exponentContsimplexNumericUint64 = NULL, *exptableContsimplexNumericUint64 = NULL;
	const double *exponentContsimplexNumericSparseDouble = NULL, *exptableContsimplexNumericSparseDouble = NULL;
	const mwIndex *exponentContsimplexNumericSparseIR = NULL, *exptableContsimplexNumericSparseIR = NULL;
	const mwIndex *exponentContsimplexNumericSparseJC = NULL, *exptableContsimplexNumericSparseJC = NULL;
	//uint8_T *isequalcolumn = NULL;
    
	/* Check for proper number of input and output arguments */
	if (nrhs != NRHS) { 
		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Not enough input arguments.");
		return;
	} 
	if (nlhs > NLHS) {
		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Too many output arguments.");
		return;
	}
	if (pargEXPONENT == NULL || pargEXPTABLE == NULL || pargJUMP == NULL) {
		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Invalid input arguments.");
		return;
	}
	if (!mxIsCell(pargEXPONENT)) {
		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Exponents must be of type 'cell'.");
		return;
	}
	if (!mxIsCell(pargEXPTABLE)) {
		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Exptable must be of type 'cell'.");
		return;
	}
	if (!mxIsNumeric(pargJUMP) || !mxIsDouble(pargJUMP) || mxIsComplex(pargJUMP)) {
		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Jump must be of type 'double'.");
		return;
	}
	if (mxIsEmpty(pargEXPONENT) || mxIsEmpty(pargEXPTABLE) || mxIsEmpty(pargJUMP)) {
		pargINDEX = mxCreateDoubleMatrix(1, 1, mxREAL);
		if (pargINDEX == NULL) {
			mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Out of Memory.");
			return;
		}
		indexptr = mxGetPr(pargINDEX);
		indexptr[0] = 1;
		return;
	}
	if (mxGetNumberOfDimensions(pargEXPONENT) > 2) {
		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Exponent must be a tow dimensional cell array.");
		return;
	}
	if (mxGetNumberOfDimensions(pargEXPTABLE) > 2) {
		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Exptable must be a tow dimensional cell array.");
		return;
	}
	if (mxGetNumberOfDimensions(pargJUMP) > 2) {
		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Jump must be a tow dimensional cell array.");
		return;
	}
	dimensionsEXPONENT = mxGetDimensions(pargEXPONENT);
	dimensionsEXPTABLE = mxGetDimensions(pargEXPTABLE);
	dimensionsJUMP = mxGetDimensions(pargJUMP);
	if (dimensionsEXPONENT[0] != dimensionsEXPTABLE[0] || dimensionsEXPONENT[0] != dimensionsJUMP[0]) {
		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "All input arguments must have the same number of rows.");
		return;
	}
	if (dimensionsEXPONENT[1] != dimensionsEXPTABLE[1] || dimensionsEXPONENT[1] != dimensionsJUMP[1]) {
		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "All input arguments must have the same number of columns.");
		return;
	}
	lengthExponent = dimensionsEXPONENT[0]*dimensionsEXPONENT[1];
	jumpContsimplexNumeric = mxGetPr(pargJUMP);
	for (countsimplex = 0; countsimplex < lengthExponent; ++countsimplex) {
		exponentContsimplex = mxGetCell(pargEXPONENT, countsimplex);
		if (exponentContsimplex == NULL) {
			mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Element %d of exponent is invalid.", countsimplex + 1);
			return;
		}
		exptableContsimplex = mxGetCell(pargEXPTABLE, countsimplex);
		if (exptableContsimplex == NULL) {
			mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Element %d of exponent is invalid.", countsimplex + 1);
			return;
		}
		if (mxIsComplex(exponentContsimplex)) {
			mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Element %d of exponent must not be complex.", countsimplex + 1);
			return;
		}
		if (mxIsComplex(exptableContsimplex)) {
			mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Element %d of exptable must not be complex.", countsimplex + 1);
			return;
		}
		if (mxIsEmpty(exptableContsimplex) && mxIsEmpty(exptableContsimplex)) {
			index = 0;
			break;
		}
		if (mxGetNumberOfDimensions(exponentContsimplex) > 2) {
			mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Elements of exponent must be two dimensional matrices.");
			return;
		}
		dimensionsEXPONENTElement = mxGetDimensions(exponentContsimplex);
		if (mxGetNumberOfDimensions(exptableContsimplex) > 2) {
			mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Elements of exptable must be two dimensional matrices.");
			return;
		}
		dimensionsEXPTABLEElement = mxGetDimensions(exptableContsimplex);
		if (dimensionsEXPONENTElement[1] != dimensionsEXPTABLEElement[1]) {
			mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Number of columns of exponent and exptable must match for element %d.", countsimplex + 1);
			return;
		}
		if (1 != dimensionsEXPONENTElement[0]) {
			mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Exponent must have one column for element %d.", countsimplex + 1);
			return;
		}
		if (dimensionsEXPTABLEElement[0] <= 0) {
			mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Exptable must have at least one column for element %d.", countsimplex + 1);
			return;
		}
		if (mxIsSparse(exponentContsimplex) && mxIsSparse(exptableContsimplex)) {
			int success = 0;
			mwSize nzmaxEXPTABLE = 0;
			double *exptableContsimplexNumericSparseDoubleCSR = NULL;
			mwIndex *exptableContsimplexNumericSparseJCCSR = NULL, *exptableContsimplexNumericSparseIRCSR = NULL;
			//double *exptableContsimplexNumericSparseDoubleCSC = NULL;
			//mwIndex *exptableContsimplexNumericSparseJCCSC = NULL, *exptableContsimplexNumericSparseIRCSC = NULL;
			exponentContsimplexNumericSparseDouble = mxGetPr(exponentContsimplex);
			exptableContsimplexNumericSparseDouble = mxGetPr(exptableContsimplex);
			if (exponentContsimplexNumericSparseDouble == NULL || exptableContsimplexNumericSparseDouble == NULL) {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
				return;
			}
			exponentContsimplexNumericSparseIR = mxGetIr(exponentContsimplex);
			exptableContsimplexNumericSparseIR = mxGetIr(exptableContsimplex);
			if (exponentContsimplexNumericSparseIR == NULL || exptableContsimplexNumericSparseIR == NULL) {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
				return;
			}
			exponentContsimplexNumericSparseJC = mxGetJc(exponentContsimplex);
			exptableContsimplexNumericSparseJC = mxGetJc(exptableContsimplex);
			if (exponentContsimplexNumericSparseJC == NULL || exptableContsimplexNumericSparseJC == NULL) {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
				return;
			}
			nzmaxEXPTABLE = mxGetNzmax(exptableContsimplex);
			exptableContsimplexNumericSparseDoubleCSR = (double*) mxCalloc(nzmaxEXPTABLE, sizeof(double));
			if (exptableContsimplexNumericSparseDoubleCSR == NULL) {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Out of Memory.");
			}
			exptableContsimplexNumericSparseIRCSR = (mwIndex*) mxCalloc(nzmaxEXPTABLE, sizeof(mwIndex));
			if (exptableContsimplexNumericSparseIRCSR == NULL) {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Out of Memory.");
			}
			exptableContsimplexNumericSparseJCCSR = (mwIndex*) mxCalloc(dimensionsEXPTABLEElement[0] + 1, sizeof(mwIndex));
			if (exptableContsimplexNumericSparseJCCSR == NULL) {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Out of Memory.");
			}
			success = csc_tocsr(dimensionsEXPTABLEElement[0], dimensionsEXPTABLEElement[1], exptableContsimplexNumericSparseJC, exptableContsimplexNumericSparseIR, exptableContsimplexNumericSparseDouble, exptableContsimplexNumericSparseJCCSR, exptableContsimplexNumericSparseIRCSR, exptableContsimplexNumericSparseDoubleCSR);
			if (success != 0) {
				if (exptableContsimplexNumericSparseDoubleCSR != NULL) {
					mxFree(exptableContsimplexNumericSparseDoubleCSR);
					exptableContsimplexNumericSparseDoubleCSR = NULL;
				}
				if (exptableContsimplexNumericSparseJCCSR != NULL) {
					mxFree(exptableContsimplexNumericSparseJCCSR);
					exptableContsimplexNumericSparseJCCSR = NULL;
				}
				if (exptableContsimplexNumericSparseIRCSR != NULL) {
					mxFree(exptableContsimplexNumericSparseIRCSR);
					exptableContsimplexNumericSparseIRCSR = NULL;
				}
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Could not convert sparse to CSR.");
			}
			/*exptableContsimplexNumericSparseDoubleCSC = (double*) mxCalloc(nzmaxEXPTABLE, sizeof(double));
			if (exptableContsimplexNumericSparseDoubleCSC == NULL) {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Out of Memory.");
			}
			exptableContsimplexNumericSparseIRCSC = (mwIndex*) mxCalloc(nzmaxEXPTABLE, sizeof(mwIndex));
			if (exptableContsimplexNumericSparseIRCSC == NULL) {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Out of Memory.");
			}
			exptableContsimplexNumericSparseJCCSC = (mwIndex*) mxCalloc(dimensionsEXPTABLEElement[1] + 1, sizeof(mwIndex));
			if (exptableContsimplexNumericSparseJCCSC == NULL) {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Out of Memory.");
			}
			success = csr_tocsc(dimensionsEXPTABLEElement[0], dimensionsEXPTABLEElement[1], exptableContsimplexNumericSparseJCCSR, exptableContsimplexNumericSparseIRCSR, exptableContsimplexNumericSparseDoubleCSR, exptableContsimplexNumericSparseJCCSC, exptableContsimplexNumericSparseIRCSC, exptableContsimplexNumericSparseDoubleCSC);
			mexPrintf("sparse original\n");
			for (ii = 0; ii < dimensionsEXPTABLEElement[1]; ++ii) {
				for (jj = exptableContsimplexNumericSparseJC[ii]; jj < exptableContsimplexNumericSparseJC[ii + 1]; ++jj) {
					mexPrintf("i:  %d,   j:  %d,   v:  %f\n", exptableContsimplexNumericSparseIR[jj], jj, exptableContsimplexNumericSparseDouble[jj]);
				}
			}
			mexPrintf("\n\nsparse CSR\n");
			for (ii = 0; ii < dimensionsEXPTABLEElement[0]; ++ii) {
				for (jj = exptableContsimplexNumericSparseJCCSR[ii]; jj < exptableContsimplexNumericSparseJCCSR[ii + 1]; ++jj) {
					mexPrintf("i:  %d,   j:  %d,   v:  %f\n", exptableContsimplexNumericSparseIRCSR[jj], ii, exptableContsimplexNumericSparseDoubleCSR[jj]);
				}
			}
			mexPrintf("sparse double converted\n");
			for (ii = 0; ii < dimensionsEXPTABLEElement[1]; ++ii) {
				for (jj = exptableContsimplexNumericSparseJCCSC[ii]; jj < exptableContsimplexNumericSparseJCCSC[ii + 1]; ++jj) {
					mexPrintf("i:  %d,   j:  %d,   v:  %f\n", exptableContsimplexNumericSparseIRCSC[jj], jj, exptableContsimplexNumericSparseDoubleCSC[jj]);
				}
			}*/
			for (ii = 0; ii < dimensionsEXPTABLEElement[0]; ++ii) {
				ismatch = true;
				ismatchvalue = true;
				ismatchcolumns = true;
				jj = 0;
				for (kk = 0; kk < dimensionsEXPTABLEElement[1]; ++kk) {
					// table contains value in current column
					//if (kk <= exptableContsimplexNumericSparseJCCSR[ii] && kk < exptableContsimplexNumericSparseJCCSR[ii + 1]) {
					if (kk >= exptableContsimplexNumericSparseIRCSR[exptableContsimplexNumericSparseJCCSR[ii]] && exptableContsimplexNumericSparseJCCSR[ii + 1] > 0 && kk <= exptableContsimplexNumericSparseIRCSR[exptableContsimplexNumericSparseJCCSR[ii + 1] - 1]) {
						if (exptableContsimplexNumericSparseIRCSR[jj + exptableContsimplexNumericSparseJCCSR[ii]] == kk) {
							// row in exponent is empty but not empty in exptable
							//if (exponentContsimplexNumericSparseJC[exptableContsimplexNumericSparseIRCSR[kk] + 1] - exponentContsimplexNumericSparseJC[exptableContsimplexNumericSparseIRCSR[kk]] <= 0) {
							if (exponentContsimplexNumericSparseJC[kk + 1] - exponentContsimplexNumericSparseJC[kk] <= 0) {
								ismatchvalue = false;
								ismatchcolumns = false;
								break;
							}
							// value is equal
							//if (exponentContsimplexNumericSparseDouble[exponentContsimplexNumericSparseJC[exptableContsimplexNumericSparseIRCSR[kk]]] != exptableContsimplexNumericSparseDoubleCSR[kk]) {
							if (exponentContsimplexNumericSparseDouble[exponentContsimplexNumericSparseJC[kk]] != exptableContsimplexNumericSparseDoubleCSR[jj + exptableContsimplexNumericSparseJCCSR[ii]]) {
								ismatchvalue = false;
								ismatchcolumns = false;
								break;
							}
							++jj;
						}
					}
					else {
						// exponent has value not in table
						if (exponentContsimplexNumericSparseJC[kk + 1] - exponentContsimplexNumericSparseJC[kk] > 0) {
							ismatchcolumns = false;
							break;
						}
					}
					if (!ismatchvalue) {
						ismatch = false;
						break;
					}
					if (!ismatchcolumns) {
						ismatchvalue = false;
						ismatch = false;
						break;
					}
				}
				if (!ismatchcolumns) {
					ismatch = false;
				}
				if (!ismatchvalue) {
					ismatch = false;
				}
				if (ismatch) {
					index = ((double)ii)*jumpContsimplexNumeric[countsimplex];
					break;
				}
			}
			if (exptableContsimplexNumericSparseDoubleCSR != NULL) {
				mxFree(exptableContsimplexNumericSparseDoubleCSR);
				exptableContsimplexNumericSparseDoubleCSR = NULL;
			}
			if (exptableContsimplexNumericSparseJCCSR != NULL) {
				mxFree(exptableContsimplexNumericSparseJCCSR);
				exptableContsimplexNumericSparseJCCSR = NULL;
			}
			if (exptableContsimplexNumericSparseIRCSR != NULL) {
				mxFree(exptableContsimplexNumericSparseIRCSR);
				exptableContsimplexNumericSparseIRCSR = NULL;
			}
		}
		else if (!mxIsSparse(exponentContsimplex) && mxIsSparse(exptableContsimplex)) {
			/*Sparse table and double exponent*/
			int success = 0;
			mwSize nzmaxEXPTABLE = 0;
			double *exptableContsimplexNumericSparseDoubleCSR = NULL;
			mwIndex *exptableContsimplexNumericSparseJCCSR = NULL, *exptableContsimplexNumericSparseIRCSR = NULL;
			exponentContsimplexNumericDouble = mxGetPr(exponentContsimplex);
			exptableContsimplexNumericSparseDouble = mxGetPr(exptableContsimplex);
			if (exponentContsimplexNumericDouble == NULL || exptableContsimplexNumericSparseDouble == NULL) {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
				return;
			}
			exptableContsimplexNumericSparseIR = mxGetIr(exptableContsimplex);
			if (exptableContsimplexNumericSparseIR == NULL) {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
				return;
			}
			exptableContsimplexNumericSparseJC = mxGetJc(exptableContsimplex);
			if (exptableContsimplexNumericSparseJC == NULL) {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
				return;
			}
			nzmaxEXPTABLE = mxGetNzmax(exptableContsimplex);
			exptableContsimplexNumericSparseDoubleCSR = (double*) mxCalloc(nzmaxEXPTABLE, sizeof(double));
			if (exptableContsimplexNumericSparseDoubleCSR == NULL) {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Out of Memory.");
			}
			exptableContsimplexNumericSparseIRCSR = (mwIndex*) mxCalloc(nzmaxEXPTABLE, sizeof(mwIndex));
			if (exptableContsimplexNumericSparseIRCSR == NULL) {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Out of Memory.");
			}
			exptableContsimplexNumericSparseJCCSR = (mwIndex*) mxCalloc(dimensionsEXPTABLEElement[0] + 1, sizeof(mwIndex));
			if (exptableContsimplexNumericSparseJCCSR == NULL) {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Out of Memory.");
			}
			success = csc_tocsr(dimensionsEXPTABLEElement[0], dimensionsEXPTABLEElement[1], exptableContsimplexNumericSparseJC, exptableContsimplexNumericSparseIR, exptableContsimplexNumericSparseDouble, exptableContsimplexNumericSparseJCCSR, exptableContsimplexNumericSparseIRCSR, exptableContsimplexNumericSparseDoubleCSR);
			if (success != 0) {
				if (exptableContsimplexNumericSparseDoubleCSR != NULL) {
					mxFree(exptableContsimplexNumericSparseDoubleCSR);
					exptableContsimplexNumericSparseDoubleCSR = NULL;
				}
				if (exptableContsimplexNumericSparseJCCSR != NULL) {
					mxFree(exptableContsimplexNumericSparseJCCSR);
					exptableContsimplexNumericSparseJCCSR = NULL;
				}
				if (exptableContsimplexNumericSparseIRCSR != NULL) {
					mxFree(exptableContsimplexNumericSparseIRCSR);
					exptableContsimplexNumericSparseIRCSR = NULL;
				}
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Could not convert sparse to CSR.");
			}
			for (ii = 0; ii < dimensionsEXPTABLEElement[0]; ++ii) {
				ismatch = true;
				ismatchvalue = true;
				ismatchcolumns = true;
				kk = 0;
				for (jj = 0; jj < dimensionsEXPTABLEElement[1]; ++jj) {
					// table contains value in current column
					if (jj >= exptableContsimplexNumericSparseIRCSR[exptableContsimplexNumericSparseJCCSR[ii]] && exptableContsimplexNumericSparseJCCSR[ii + 1] > 0 && jj <= exptableContsimplexNumericSparseIRCSR[exptableContsimplexNumericSparseJCCSR[ii + 1] - 1]) {
					//if (jj <= exptableContsimplexNumericSparseJCCSR[ii] && jj < exptableContsimplexNumericSparseJCCSR[ii + 1]) {
						if (exptableContsimplexNumericSparseIRCSR[kk + exptableContsimplexNumericSparseJCCSR[ii]] == jj) {
							// value is equal
							if (exponentContsimplexNumericDouble[jj] != exptableContsimplexNumericSparseDoubleCSR[kk + exptableContsimplexNumericSparseJCCSR[ii]]) {
								ismatchvalue = false;
								ismatchcolumns = false;
								break;
							}
							++kk;
						}
					}
					else {
						// value not in table and exponent must be 0
						if (exponentContsimplexNumericDouble[jj] != 0.0) {
							ismatchvalue = false;
							break;
						}
					}
					/*for (ll = exptableContsimplexNumericSparseJC[jj]; ll < exptableContsimplexNumericSparseJC[jj + 1]; ++ll) {
						if (exptableContsimplexNumericSparseIR[ll] == ii) {
							ismatchcolumns = true;
							//if (exptableContsimplexNumericSparseJC[jj + 1] - exptableContsimplexNumericSparseJC[jj] != exponentContsimplexNumericSparseJC[jj + 1] - exponentContsimplexNumericSparseJC[jj]) {
							//	ismatchvalue = false;
							//	break;
							//}
							if (exptableContsimplexNumericSparseDouble[ll] != exponentContsimplexNumericDouble[jj]) {
								ismatchvalue = false;
								break;
							}
						}
						else if (exptableContsimplexNumericSparseIR[ll] > ii) {
							break;
						}
					}*/
					//if (!ismatchcolumns && exponentContsimplexNumericDouble[jj] != 0.0) {
					//	ismatchvalue = false;
					//	ismatch = false;
					//	break;
					//}
					if (!ismatchvalue) {
						ismatch = false;
						break;
					}
					if (!ismatchcolumns) {
						ismatch = false;
						break;
					}
					//if (!ismatch) {
					//	break;
					//}
				}
				if (!ismatchcolumns) {
					ismatch = false;
				}
				if (!ismatchvalue) {
					ismatch = false;
				}
				if (ismatch) {
					index = ((double)ii)*jumpContsimplexNumeric[countsimplex];
					break;
				}
			}
			if (exptableContsimplexNumericSparseDoubleCSR != NULL) {
				mxFree(exptableContsimplexNumericSparseDoubleCSR);
				exptableContsimplexNumericSparseDoubleCSR = NULL;
			}
			if (exptableContsimplexNumericSparseJCCSR != NULL) {
				mxFree(exptableContsimplexNumericSparseJCCSR);
				exptableContsimplexNumericSparseJCCSR = NULL;
			}
			if (exptableContsimplexNumericSparseIRCSR != NULL) {
				mxFree(exptableContsimplexNumericSparseIRCSR);
				exptableContsimplexNumericSparseIRCSR = NULL;
			}
		}
		else {
			if (mxIsSparse(exponentContsimplex) || mxIsSparse(exptableContsimplex)) {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Element %d of exptable and exponent must be full and not sparse.", countsimplex + 1);
				return;
			}
			isdouble = mxIsDouble(exponentContsimplex) && mxIsDouble(exptableContsimplex);
			isfloat = mxIsSingle(exponentContsimplex) && mxIsSingle(exptableContsimplex);
			isint8 = mxIsInt8(exponentContsimplex) && mxIsInt8(exptableContsimplex);
			isint16 = mxIsInt16(exponentContsimplex) && mxIsInt16(exptableContsimplex);
			isint32 = mxIsInt32(exponentContsimplex) && mxIsInt32(exptableContsimplex);
			isint64 = mxIsInt64(exponentContsimplex) && mxIsInt64(exptableContsimplex);
			isuint8 = mxIsUint8(exponentContsimplex) && mxIsUint8(exptableContsimplex);
			isuint16 = mxIsUint16(exponentContsimplex) && mxIsUint16(exptableContsimplex);
			isuint32 = mxIsUint32(exponentContsimplex) && mxIsUint32(exptableContsimplex);
			isuint64 = mxIsUint64(exponentContsimplex) && mxIsUint64(exptableContsimplex);
			/*if (!mxIsDouble(exponentContsimplex)) {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Element %d of exponent must be of type 'double'.", countsimplex + 1);
				return;
			}
			if (!mxIsDouble(exptableContsimplex)) {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Element %d of exptable must be of type 'double'.", countsimplex + 1);
				return;
			}*/
			if (!isdouble && !isfloat && !isint8 && !isint16 && !isint32 && !isint64 && !isuint8 && !isuint16 && !isuint32 && !isuint64) {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Element %d of exptable must be numeric and of the same type as exponent.", countsimplex + 1);
				return;
			}
			if (isdouble) {
				exponentContsimplexNumericDouble = mxGetPr(exponentContsimplex);
				exptableContsimplexNumericDouble = mxGetPr(exptableContsimplex);
				if (exponentContsimplexNumericDouble == NULL || exptableContsimplexNumericDouble == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEX(Double)
			}
			else if (isfloat) {
				exponentContsimplexNumericFloat = (float*) mxGetData(exponentContsimplex);
				exptableContsimplexNumericFloat = (float*) mxGetData(exptableContsimplex);
				if (exponentContsimplexNumericFloat == NULL || exptableContsimplexNumericFloat == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEX(Float)
			}
			else if (isint8) {
				exponentContsimplexNumericInt8 = (char*) mxGetData(exponentContsimplex);
				exptableContsimplexNumericInt8 = (char*) mxGetData(exptableContsimplex);
				if (exponentContsimplexNumericInt8 == NULL || exptableContsimplexNumericInt8 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEX(Int8)
			}
			else if (isint16) {
				exponentContsimplexNumericInt16 = (short*) mxGetData(exponentContsimplex);
				exptableContsimplexNumericInt16 = (short*) mxGetData(exptableContsimplex);
				if (exponentContsimplexNumericInt16 == NULL || exptableContsimplexNumericInt16 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEX(Int16)
			}
			else if (isint32) {
				exponentContsimplexNumericInt32 = (int*) mxGetData(exponentContsimplex);
				exptableContsimplexNumericInt32 = (int*) mxGetData(exptableContsimplex);
				if (exponentContsimplexNumericInt32 == NULL || exptableContsimplexNumericInt32 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEX(Int32)
			}
			else if (isint64) {
				exponentContsimplexNumericInt64 = (long long*) mxGetData(exponentContsimplex);
				exptableContsimplexNumericInt64 = (long long*) mxGetData(exptableContsimplex);
				if (exponentContsimplexNumericInt64 == NULL || exptableContsimplexNumericInt64 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEX(Int64)
			}
			else if (isuint8) {
				exponentContsimplexNumericUint8 = (unsigned char*) mxGetData(exponentContsimplex);
				exptableContsimplexNumericUint8 = (unsigned char*) mxGetData(exptableContsimplex);
				if (exponentContsimplexNumericUint8 == NULL || exptableContsimplexNumericUint8 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEX(Uint8)
			}
			else if (isuint16) {
				exponentContsimplexNumericUint16 = (unsigned short*) mxGetData(exponentContsimplex);
				exptableContsimplexNumericUint16 = (unsigned short*) mxGetData(exptableContsimplex);
				if (exponentContsimplexNumericUint16 == NULL || exptableContsimplexNumericUint16 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEX(Uint16)
			}
			else if (isuint32) {
				exponentContsimplexNumericUint32 = (unsigned int*) mxGetData(exponentContsimplex);
				exptableContsimplexNumericUint32 = (unsigned int*) mxGetData(exptableContsimplex);
				if (exponentContsimplexNumericUint32 == NULL || exptableContsimplexNumericUint32 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEX(Uint32)
			}
			else if (isuint64) {
				exponentContsimplexNumericUint64 = (unsigned long long*) mxGetData(exponentContsimplex);
				exptableContsimplexNumericUint64 = (unsigned long long*) mxGetData(exptableContsimplex);
				if (exponentContsimplexNumericUint64 == NULL || exptableContsimplexNumericUint64 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEX(Uint64)
			}
			else {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Unknown data for element %d.", countsimplex + 1);
				return;
			}
		}
	}
	index = index + 1;
	pargINDEX = mxCreateDoubleMatrix(1, 1, mxREAL);
	if (pargINDEX == NULL) {
		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Out of Memory.");
		return;
	}
	indexptr = mxGetPr(pargINDEX);
	indexptr[0] = index;
}


/*https://github.com/scipy/scipy/blob/3b36a574dc657d1ca116f6e230be694f3de31afc/scipy/sparse/sparsetools/csr.h#L376
 * Compute B = A for CSR matrix A, CSC matrix B
 *
 * Also, with the appropriate arguments can also be used to:
 *   - compute B = A^t for CSR matrix A, CSR matrix B
 *   - compute B = A^t for CSC matrix A, CSC matrix B
 *   - convert CSC->CSR
 *
 * Input Arguments:
 *   I  n_row         - number of rows in A
 *   I  n_col         - number of columns in A
 *   I  Ap[n_row+1]   - row pointer
 *   I  Aj[nnz(A)]    - column indices
 *   T  Ax[nnz(A)]    - nonzeros
 *
 * Output Arguments:
 *   I  Bp[n_col+1] - column pointer
 *   I  Bj[nnz(A)]  - row indices
 *   T  Bx[nnz(A)]  - nonzeros
 *
 * Note:
 *   Output arrays Bp, Bj, Bx must be preallocated
 *
 * Note: 
 *   Input:  column indices *are not* assumed to be in sorted order
 *   Output: row indices *will be* in sorted order
 *
 *   Complexity: Linear.  Specifically O(nnz(A) + max(n_row,n_col))
 * 
 */
int csr_tocsc(const mwIndex n_row,
				const mwIndex n_col, 
				const mwIndex *Ap, 
				const mwIndex *Aj, 
				const double *Ax,
						mwIndex *Bp,
						mwIndex *Bi,
						double *Bx)
{  
	mwIndex ii = 0, n = 0, col = 0, cumsum = 0, row = 0, jj = 0, dest = 0, last;
	const mwIndex nnz = Ap[n_row];
	mwIndex temp = 0;

	if (n_row <= 0) {
		return 1;
	}
	if (n_col <= 0) {
		return 1;
	}
	if (Ap == NULL) {
		return 1;
	}
	if (Aj == NULL) {
		return 1;
	}
	if (Ax == NULL) {
		return 1;
	}
	if (Bp == NULL) {
		return 1;
	}
	if (Bi == NULL) {
		return 1;
	}
	if (Bx == NULL) {
		return 1;
	}

	//compute number of non-zero entries per column of A 
	for (ii = 0; ii < n_col; ++ii) {
		//std::fill(Bp, Bp + n_col, 0);
		Bp[ii] = 0;
	}

	for (n = 0; n < nnz; n++){            
		Bp[Aj[n]]++;
	}

	//cumsum the nnz per column to get Bp[]
	for (col = 0, cumsum = 0; col < n_col; col++){     
		temp  = Bp[col];
		Bp[col] = cumsum;
		cumsum += temp;
	}
	Bp[n_col] = nnz; 

	for (row = 0; row < n_row; row++){
		for (jj = Ap[row]; jj < Ap[row+1]; jj++){
			col  = Aj[jj];
			dest = Bp[col];

			Bi[dest] = row;
			Bx[dest] = Ax[jj];

			Bp[col]++;
		}
	}  

	for (col = 0, last = 0; col <= n_col; col++){
		temp  = Bp[col];
		Bp[col] = last;
		last    = temp;
	}
	return 0;
}
int csc_tocsr(const mwIndex n_row,
               const mwIndex n_col, 
               const mwIndex *Ap, 
               const mwIndex *Ai, 
               const double *Ax,
                     mwIndex *Bp,
                     mwIndex *Bj,
                     double *Bx)
{ return csr_tocsc(n_col, n_row, Ap, Ai, Ax, Bp, Bj, Bx); }


////https://github.com/fourier/libspmatrix/blob/master/inc/sp_matrix.h
//#define CRS 0
//#define CCS 1
//#define MTX(m,i,j,v) sp_matrix_element_add((m),(i),(j),(v));
//typedef struct
//{
//  int width;                    /* size of an array */
//  int last_index;               /* last stored index, i.e. if width = 20
//                                 * it will be 9 if only 10 nonzero elements
//                                 * stored */
//  int  *indexes;                /* array of column/row indexes */
//  double *values;               /* array of values */
//} indexed_array;
//typedef indexed_array* indexed_array_ptr;
//typedef struct
//{
//  int rows_count;
//  int cols_count;
//  indexed_array_ptr storage;
//  int ordered;                               /* if matrix was finalized */
//  int storage_type;          /* Storage type */
//} sp_matrix;
//typedef sp_matrix* sp_matrix_ptr;
//double sp_matrix_element_add(sp_matrix_ptr self,int i, int j, double value)
//{
//  int index,new_width,I,J;
//  int* indexes = (int*)0;
//  double* values = (double*)0;
//  /* check for matrix and if i,j are proper indicies */
//  if (self == NULL) {
//  
//  }
//  if (i >= 0 && i < self->rows_count) {
//  
//  }
//  if (j >= 0 && j < self->cols_count) {
//  
//  }
//  /* set I and J to be i and j in case of CRS or j and i otherwise */
//  I = self->storage_type == CRS ? i : j;
//  J = self->storage_type == CRS ? j : i;
//  /* loop by nonzero columns in row/col i */
//  for (index = 0; index <= self->storage[I].last_index; ++ index)
//    if (self->storage[I].indexes[index] == J)
//    {
//      /* nonzerod element found, add to it */
//      self->storage[I].values[index] += value;
//      return self->storage[I].values[index];
//    }
//  /* needed to add a new element to the row/col */
//    
//  /*
//   * check if bandwidth is not exceed and reallocate memory
//   * if necessary
//   */
//  if (self->storage[I].last_index == self->storage[I].width - 1)
//  {
//    new_width = self->storage[I].width*2;
//    if (new_width <= 0)             /* avoid crashes on bad bandwidth */
//      new_width = 1;
//    indexes = (int*)realloc(self->storage[I].indexes,new_width*sizeof(int));
//	if (indexes == NULL) {
//	
//	}
//    self->storage[I].indexes = indexes;
//    values = (double*)realloc(self->storage[I].values,new_width*sizeof(double));
//	if (values == NULL) {
//	
//	}
//    self->storage[I].values = values;
//    self->storage[I].width = new_width;
//    //self->ordered = NOT_ORDERED;
//  }
//  /* add an element to the row/col */
//  self->storage[I].last_index++;
//  self->storage[I].values[self->storage[I].last_index] = value;
//  self->storage[I].indexes[self->storage[I].last_index] = J;
//  return value;
//}
//
//void sp_matrix_init(sp_matrix_ptr mtx,
//                    int rows,
//                    int cols,
//                    int bandwidth,
//                    int type)
//{
//  int i,n;
//  if (mtx)
//  {
//    mtx->rows_count = rows;
//    mtx->cols_count = cols;
//    //mtx->ordered = NOT_ORDERED;
//    mtx->storage_type = type;
//    n = type == CRS ? rows : cols;
//    mtx->storage = (indexed_array*)malloc(sizeof(indexed_array)*n);
//    /* create rows or cols with fixed bandwidth */
//    for (i = 0; i < n; ++ i)
//    {
//      mtx->storage[i].width = bandwidth;
//      mtx->storage[i].last_index = -1;
//      mtx->storage[i].indexes = (int*)malloc(sizeof(int)*bandwidth);
//      mtx->storage[i].values = (double*)malloc(sizeof(double)*bandwidth);
//      memset(mtx->storage[i].indexes,0,sizeof(int)*bandwidth);
//      memset(mtx->storage[i].values,0,sizeof(double)*bandwidth);
//    }
//  }
//}
//int sp_matrix_convert(sp_matrix_ptr mtx_from, sp_matrix_ptr mtx_to, int type) {
//	int i,j;
//	if (type == mtx_from->storage_type) {
//		return 0;
//	}
//	sp_matrix_init(mtx_to,
//					mtx_from->rows_count,
//					mtx_from->cols_count,
//					mtx_from->storage[0].width,
//					type);
//	if (type == CCS) {/* CRS -> CCS */
//		for (i = 0; i < mtx_from->rows_count; ++ i) {
//			for (j = 0; j <= mtx_from->storage[i].last_index; ++ j)
//			MTX(mtx_to,i,mtx_from->storage[i].indexes[j], mtx_from->storage[i].values[j]);
//		}
//	}
//	else {/* CCS -> CRS */
//		for (i = 0; i < mtx_from->cols_count; ++ i) {
//			for (j = 0; j <= mtx_from->storage[i].last_index; ++ j)
//				MTX(mtx_to,mtx_from->storage[i].indexes[j],i, mtx_from->storage[i].values[j]);
//		}
//	}
//	return 1;
//}
/*______________________________________________________________________________________________________________________________________________________________*/

//#include <string.h>
///* Type Definitions */
//#ifndef struct_emxArray_int32_T
//#define struct_emxArray_int32_T
//
//struct emxArray_int32_T
//{
//  int32_T *data;
//  int32_T *size;
//  int32_T allocatedSize;
//  int32_T numDimensions;
//  boolean_T canFreeData;
//};
//
//#endif                                 /*struct_emxArray_int32_T*/
//
//#ifndef typedef_emxArray_int32_T
//#define typedef_emxArray_int32_T
//
//typedef struct emxArray_int32_T emxArray_int32_T;
//
//#endif                                 /*typedef_emxArray_int32_T*/
//
//#ifndef struct_emxArray_real_T
//#define struct_emxArray_real_T
//
//struct emxArray_real_T
//{
//  real_T *data;
//  int32_T *size;
//  int32_T allocatedSize;
//  int32_T numDimensions;
//  boolean_T canFreeData;
//};
//
//#endif                                 /*struct_emxArray_real_T*/
//
//#ifndef typedef_emxArray_real_T
//#define typedef_emxArray_real_T
//
//typedef struct emxArray_real_T emxArray_real_T;
//
//#endif                                 /*typedef_emxArray_real_T*/
//
//typedef struct {
//  emxArray_real_T *d;
//  emxArray_int32_T *colidx;
//  emxArray_int32_T *rowidx;
//  int32_T m;
//  int32_T n;
//  int32_T maxnz;
//} coder_internal_sparse;
//
//#ifndef typedef_coder_internal_sparse_1
//#define typedef_coder_internal_sparse_1
//
//typedef struct {
//  emxArray_real_T *d;
//  emxArray_int32_T *colidx;
//  emxArray_int32_T *rowidx;
//  int32_T n;
//  int32_T maxnz;
//} coder_internal_sparse_1;
//
//#endif                                 /*typedef_coder_internal_sparse_1*/
//
//#ifndef struct_emxArray_boolean_T
//#define struct_emxArray_boolean_T
//
//struct emxArray_boolean_T
//{
//  boolean_T *data;
//  int32_T *size;
//  int32_T allocatedSize;
//  int32_T numDimensions;
//  boolean_T canFreeData;
//};
//
//#endif                                 /*struct_emxArray_boolean_T*/
//
//#ifndef typedef_emxArray_boolean_T
//#define typedef_emxArray_boolean_T
//
//typedef struct emxArray_boolean_T emxArray_boolean_T;
//
//#endif                                 /*typedef_emxArray_boolean_T*/
//void emxEnsureCapacity_int32_T(emxArray_int32_T *emxArray,
//  int32_T oldNumel)
//{
//  int32_T newNumel;
//  int32_T i;
//  void *newData;
//  if (oldNumel < 0) {
//    oldNumel = 0;
//  }
//
//  newNumel = 1;
//  for (i = 0; i < emxArray->numDimensions; i++) {
//    newNumel = (int32_T)newNumel*(uint32_T)emxArray->size[i];
//  }
//
//  if (newNumel > emxArray->allocatedSize) {
//    i = emxArray->allocatedSize;
//    if (i < 16) {
//      i = 16;
//    }
//
//    while (i < newNumel) {
//      if (i > 1073741823) {
//        i = MAX_int32_T;
//      } else {
//        i <<= 1;
//      }
//    }
//
//    newData = mxCalloc((uint32_T)i, sizeof(int32_T));
//    if (newData == NULL) {
//		mexErrMsgIdAndTxt("", "Out of Memory.");
//    }
//
//    if (emxArray->data != NULL) {
//      memcpy(newData, emxArray->data, sizeof(int32_T) * oldNumel);
//      if (emxArray->canFreeData) {
//        mxFree(emxArray->data);
//      }
//    }
//
//    emxArray->data = (int32_T *)newData;
//    emxArray->allocatedSize = i;
//    emxArray->canFreeData = true;
//  }
//}
//
//void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int32_T oldNumel) {
//  int32_T newNumel;
//  int32_T i;
//  void *newData;
//  if (oldNumel < 0) {
//    oldNumel = 0;
//  }
//
//  newNumel = 1;
//  for (i = 0; i < emxArray->numDimensions; i++) {
//    newNumel = (uint32_T)newNumel*(uint32_T)emxArray->size[i];
//  }
//
//  if (newNumel > emxArray->allocatedSize) {
//    i = emxArray->allocatedSize;
//    if (i < 16) {
//      i = 16;
//    }
//
//    while (i < newNumel) {
//      if (i > 1073741823) {
//        i = MAX_int32_T;
//      } else {
//        i <<= 1;
//      }
//    }
//
//    newData = mxCalloc((uint32_T)i, sizeof(real_T));
//    if (newData == NULL) {
//		mexErrMsgIdAndTxt("", "Out of Memory.");
//    }
//
//    if (emxArray->data != NULL) {
//      memcpy(newData, emxArray->data, sizeof(real_T) * oldNumel);
//      if (emxArray->canFreeData) {
//        mxFree(emxArray->data);
//      }
//    }
//
//    emxArray->data = (real_T *)newData;
//    emxArray->allocatedSize = i;
//    emxArray->canFreeData = true;
//  }
//}
//
//int ismembersparse(const coder_internal_sparse table, const coder_internal_sparse_1 search, emxArray_boolean_T *b_ismemb) {
//	int32_T i0;
//	int32_T low_ip1;
//	emxArray_int32_T *eq_rowidx;
//	emxArray_real_T *obj_d;
//	emxArray_int32_T *obj_colidx;
//	emxArray_boolean_T *S;
//	emxArray_int32_T *y_colidx;
//	emxArray_boolean_T *expl_temp;
//	emxArray_int32_T *b_expl_temp;
//	int32_T ii;
//	int32_T sn;
//	int32_T i1;
//	int32_T colNnz;
//	int32_T k;
//	int32_T cidx;
//	int32_T nInt;
//	int32_T high_i;
//	boolean_T moreAToDo;
//	int32_T low_i;
//	int32_T mid_i;
//	boolean_T moreBToDo;
//
//	/* ismemb = ismember(table, search, 'rows'); */
//	/* 'ismemb:3' ismemb = false(size(table, 1), 1); */
//	i0 = b_ismemb->size[0];
//	if (table.m < 0) {
//		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Exptable must have positive number of rows.");
//		return;
//	}
//
//	b_ismemb->size[0] = table.m;
//	emxEnsureCapacity_boolean_T(b_ismemb, i0);
//	low_ip1 = table.m;
//	if (table.m < 0) {
//		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Exptable must have positive number of rows.");
//		return;
//	}
//
//	for (i0 = 0; i0 < low_ip1; i0++) {
//		b_ismemb->data[i0] = false;
//	}
//
//	/* 'ismemb:4' for ii = 1:size(table, 1) */
//	i0 = table.m;
//	//emxInit_int32_T(sp, &eq_rowidx, 1, &e_emlrtRTEI, true);
//	eq_rowidx->data = (int32_T*) NULL;
//	eq_rowidx->numDimensions = 1;
//	eq_rowidx->size = (int32_T*) mxCalloc(1, sizeof(int32_T));
//	if (eq_rowidx->size == NULL) {
//		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Out of Memory.");
//		return;
//	}
//	eq_rowidx->allocatedSize = 0;
//	eq_rowidx->canFreeData = true;
//
//	//emxInit_real_T(sp, &obj_d, 1, &f_emlrtRTEI, true);
//	obj_d->data = (real_T*) NULL;
//	obj_d->numDimensions = 1;
//	obj_d->size = (int32_T*) mxCalloc(1, sizeof(int32_T));
//	if (obj_d->size == NULL) {
//		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Out of Memory.");
//		return;
//	}
//	obj_d->allocatedSize = 0;
//	obj_d->canFreeData = true;
//
//	//emxInit_int32_T(sp, &obj_colidx, 1, &f_emlrtRTEI, true);
//	obj_colidx->data = (int32_T*) NULL;
//	obj_colidx->numDimensions = 1;
//	obj_colidx->size = (int32_T*) mxCalloc(1, sizeof(int32_T));
//	if (obj_colidx->size == NULL) {
//		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Out of Memory.");
//		return;
//	}
//	obj_colidx->allocatedSize = 0;
//	obj_colidx->canFreeData = true;
//
//	//emxInit_boolean_T(sp, &S, 2, &d_emlrtRTEI, true);
//	S->data = (boolean_T*) NULL;
//	S->numDimensions = 1;
//	S->size = (int32_T*) mxCalloc(2, sizeof(int32_T));
//	if (S->size == NULL) {
//		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Out of Memory.");
//		return;
//	}
//	S->allocatedSize = 0;
//	S->canFreeData = true;
//
//	//emxInit_int32_T(sp, &y_colidx, 1, &g_emlrtRTEI, true);
//	y_colidx->data = (int32_T*) NULL;
//	y_colidx->numDimensions = 1;
//	y_colidx->size = (int32_T*) mxCalloc(1, sizeof(int32_T));
//	if (y_colidx->size == NULL) {
//		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Out of Memory.");
//		return;
//	}
//	y_colidx->allocatedSize = 0;
//	y_colidx->canFreeData = true;
//
//	//emxInit_boolean_T(sp, &expl_temp, 1, &h_emlrtRTEI, true);
//	expl_temp->data = (boolean_T*) NULL;
//	expl_temp->numDimensions = 1;
//	expl_temp->size = (int32_T*) mxCalloc(2, sizeof(int32_T));
//	if (expl_temp->size == NULL) {
//		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Out of Memory.");
//		return;
//	}
//	expl_temp->allocatedSize = 0;
//	expl_temp->canFreeData = true;
//
//	//emxInit_int32_T(sp, &b_expl_temp, 1, &h_emlrtRTEI, true);
//	b_expl_temp->data = (int32_T*) NULL;
//	b_expl_temp->numDimensions = 1;
//	b_expl_temp->size = (int32_T*) mxCalloc(1, sizeof(int32_T));
//	if (b_expl_temp->size == NULL) {
//		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Out of Memory.");
//		return;
//	}
//	b_expl_temp->allocatedSize = 0;
//	b_expl_temp->canFreeData = true;
//
//	for (ii = 0; ii < i0; ii++) {
//		/* 'ismemb:5' eq = table(ii, :) == search; */
//		if (!(1.0 + (real_T)ii <= table.m)) {
//			mexErrMsgIdAndTxt("Coder:builtins:IndexOutOfBounds", "Coder:builtins:IndexOutOfBounds");
//			return;
//		}
//
//		sn = table.n;
//		if (0 > table.n) {
//			mexErrMsgIdAndTxt("Coder:builtins:IndexOutOfBounds", "Coder:builtins:IndexOutOfBounds");
//			return;
//		}
//
//		obj_d->size[0] = 0;
//		i1 = obj_colidx->size[0];
//		low_ip1 = table.n + 1;
//		if (low_ip1 < 0) {
//			mexErrMsgIdAndTxt("ROLMIP:gethash:input", "low_ip1 < 0");
//		}
//
//		obj_colidx->size[0] = low_ip1;
//		emxEnsureCapacity_int32_T(obj_colidx, i1);
//		if (low_ip1 < 0) {
//			mexErrMsgIdAndTxt("ROLMIP:gethash:input", "low_ip1 < 0");
//		}
//
//		for (i1 = 0; i1 < low_ip1; i1++) {
//			obj_colidx->data[i1] = 0;
//		}
//
//		obj_colidx->data[0] = 1;
//		colNnz = 1;
//		k = 0;
//		if ((1 <= table.n) && (table.n > 2147483646)) {
//			check_forloop_overflow_error(&f_st);
//		}
//
//		for (cidx = 0; cidx < sn; cidx++) {
//			if (table.colidx->data[cidx] < table.colidx->data[cidx + 1]) {
//				if (ii + 1 < table.rowidx->data[table.colidx->data[cidx] - 1]) {
//					nInt = table.colidx->data[cidx];
//					moreAToDo = false;
//				} else {
//					high_i = table.colidx->data[cidx + 1];
//					low_i = table.colidx->data[cidx];
//					low_ip1 = table.colidx->data[cidx];
//					while (high_i > low_ip1 + 1) {
//					mid_i = (low_i >> 1) + (high_i >> 1);
//					if (((low_i & 1) == 1) && ((high_i & 1) == 1)) {
//						mid_i++;
//					}
//
//					if (ii + 1 >= table.rowidx->data[mid_i - 1]) {
//						low_i = mid_i;
//						low_ip1 = mid_i;
//					} else {
//						high_i = mid_i;
//					}
//					}
//
//					nInt = low_i + 1;
//					moreAToDo = (table.rowidx->data[low_i - 1] == ii + 1);
//				}
//			} else if (table.colidx->data[cidx] == table.colidx->data[cidx + 1]) {
//				nInt = table.colidx->data[cidx];
//				moreAToDo = false;
//			} else {
//				nInt = 1;
//				moreAToDo = false;
//			}
//
//			if (moreAToDo) {
//				i1 = obj_d->size[0];
//				mid_i = obj_d->size[0];
//				obj_d->size[0] = i1 + 1;
//				emxEnsureCapacity_real_T(obj_d, mid_i);
//				obj_d->data[i1] = table.d->data[nInt - 2];
//				obj_d->data[k] = table.d->data[nInt - 2];
//				k++;
//				colNnz++;
//			}
//
//			obj_colidx->data[cidx + 1] = colNnz;
//		}
//
//		if (obj_colidx->data[obj_colidx->size[0] - 1] - 1 == 0) {
//			i1 = obj_d->size[0];
//			obj_d->size[0] = 1;
//			emxEnsureCapacity_real_T(obj_d, i1);
//			obj_d->data[0] = 0.0;
//		}
//
//		if (table.n != search.n) {
//			mexErrMsgIdAndTxt("MATLAB:dimagree", "MATLAB:dimagree", 0);
//		}
//
//		i1 = S->size[0] * S->size[1];
//		S->size[0] = 1;
//		S->size[1] = search.n;
//		emxEnsureCapacity_boolean_T(S, i1);
//		low_ip1 = search.n;
//		for (i1 = 0; i1 < low_ip1; i1++) {
//			S->data[i1] = true;
//		}
//
//		i1 = search.n;
//		for (mid_i = 0; mid_i < i1; mid_i++) {
//			high_i = obj_colidx->data[mid_i];
//			low_ip1 = search.colidx->data[mid_i] - 1;
//			moreAToDo = (obj_colidx->data[mid_i] < obj_colidx->data[1 + mid_i]);
//			moreBToDo = (search.colidx->data[mid_i] < search.colidx->data[1 + mid_i]);
//			while (moreAToDo || moreBToDo) {
//			while ((high_i < obj_colidx->data[1 + mid_i]) && ((!moreBToDo) || (1 <
//						search.rowidx->data[low_ip1]))) {
//				S->data[mid_i] = (obj_d->data[high_i - 1] == 0.0);
//				high_i++;
//			}
//
//			moreAToDo = (high_i < obj_colidx->data[1 + mid_i]);
//			while ((low_ip1 + 1 < search.colidx->data[1 + mid_i]) && ((!moreAToDo) ||
//					(search.rowidx->data[low_ip1] < 1))) {
//				S->data[mid_i] = (0.0 == search.d->data[low_ip1]);
//				low_ip1++;
//			}
//
//			while ((high_i < obj_colidx->data[1 + mid_i]) && (low_ip1 + 1 <
//					search.colidx->data[1 + mid_i]) && (1 == search.rowidx->
//					data[low_ip1])) {
//				S->data[mid_i] = (obj_d->data[high_i - 1] == search.d->data[low_ip1]);
//				low_ip1++;
//				high_i++;
//			}
//
//			moreAToDo = (high_i < obj_colidx->data[1 + mid_i]);
//			moreBToDo = (low_ip1 + 1 < search.colidx->data[1 + mid_i]);
//			}
//		}
//
//		nInt = S->size[1];
//		if (S->size[1] >= MAX_int32_T) {
//			mexErrMsgIdAndTxt("Coder:toolbox:SparseMaxSize", "Coder:toolbox:SparseMaxSize", 2, 12, MAX_int32_T);
//		}
//
//		high_i = 0;
//		i1 = S->size[1];
//		for (k = 0; k < i1; k++) {
//			if (S->data[k]) {
//			high_i++;
//			}
//		}
//
//		low_i = S->size[1];
//		if (high_i >= 1) {
//		} else {
//			high_i = 1;
//		}
//
//		i1 = obj_colidx->size[0];
//		obj_colidx->size[0] = S->size[1] + 1;
//		emxEnsureCapacity_int32_T(obj_colidx, i1);
//		low_ip1 = S->size[1];
//		for (i1 = 0; i1 <= low_ip1; i1++) {
//			obj_colidx->data[i1] = 0;
//		}
//
//		obj_colidx->data[0] = 1;
//		i1 = eq_rowidx->size[0];
//		eq_rowidx->size[0] = high_i;
//		emxEnsureCapacity_int32_T(eq_rowidx, i1);
//		for (i1 = 0; i1 < high_i; i1++) {
//			eq_rowidx->data[i1] = 0;
//		}
//
//		eq_rowidx->data[0] = 1;
//		mid_i = 1;
//		for (high_i = 0; high_i < nInt; high_i++) {
//			if (S->data[high_i]) {
//			eq_rowidx->data[mid_i - 1] = 1;
//			mid_i++;
//			}
//
//			obj_colidx->data[high_i + 1] = mid_i;
//		}
//
//		/* 'ismemb:6' ismemb(ii, 1) = full(all(eq)); */
//		moreAToDo = (S->size[1] == 1);
//		if (moreAToDo || (S->size[1] != 1)) {
//		} else {
//			mexErrMsgIdAndTxt("Coder:toolbox:eml_all_or_any_autoDimIncompatibility", "Coder:toolbox:eml_all_or_any_autoDimIncompatibility", 0);
//		}
//
//		sparse_allOrAny(&c_st, obj_colidx, low_i, expl_temp, y_colidx, b_expl_temp);
//		i1 = b_ismemb->size[0];
//		mid_i = (int32_T)(1U + ii);
//		if ((mid_i < 1) || (mid_i > i1)) {
//			emlrtDynamicBoundsCheckR2012b(mid_i, 1, i1, &emlrtBCI, sp);
//		}
//
//		b_ismemb->data[mid_i - 1] = false;
//		high_i = y_colidx->data[1] - 1;
//		i1 = y_colidx->data[0];
//		for (nInt = i1; nInt <= high_i; nInt++) {
//			low_ip1 = b_ismemb->size[0];
//			if ((mid_i < 1) || (mid_i > low_ip1)) {
//			emlrtDynamicBoundsCheckR2012b(mid_i, 1, low_ip1, &emlrtBCI, sp);
//			}
//
//			b_ismemb->data[mid_i - 1] = true;
//		}
//	}
//
//	emxFree_int32_T(&b_expl_temp);
//	emxFree_boolean_T(&expl_temp);
//	emxFree_int32_T(&y_colidx);
//	emxFree_boolean_T(&S);
//	emxFree_int32_T(&obj_colidx);
//	emxFree_real_T(&obj_d);
//	emxFree_int32_T(&eq_rowidx);
//}
