/*
MEX function for calculation of gethash functionality.
Does basically, what sum(find(ismember(exptable, exponent, 'rows')).*jump) does, if both arguments are numeric and dimensions are appropriate.
All types of numeric exptable and exponent arguments are allowed (including sparse) inside the cell elements.

function index = gethash(exponent,exptable,jump)
	index = 0;
	for contsimplex = 1:(length(exponent))
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

#define GETINDEX(type) /*for ii = 1:size(exptable{contsimplex}, 1)*/\
		for (ii = 0; ii < dimensionsEXPTABLEElement[0]; ++ii) {\
			ismatch = true;\
			/*if all(exptable{contsimplex}(ii, :) == exponent{contsimplex})*/\
			for (jj = 0; jj < dimensionsEXPTABLEElement[1]; ++jj) {\
				if (exptableContsimplexNumeric##type[jj*(dimensionsEXPTABLEElement[0]) + ii] != exponentContsimplexNumeric##type[jj]) {\
					ismatch = false;\
					break;\
				}\
			}\
			if (ismatch) {\
				/*index = index + (ii - 1)*jump(contsimplex);*/\
				index = index + ((double)ii)*jumpContsimplexNumeric[countsimplex];\
				break;\
			}\
		}
#define GETINDEXSPARSE(type) /*for ii = 1:size(exptable{contsimplex}, 1)*/\
		for (ii = 0; ii < dimensionsEXPTABLEElement[0]; ++ii) {\
			ismatch = true;\
			ismatchvalue = true;\
			ismatchcolumns = true;\
			kk = 0;\
			/*if all(exptable{contsimplex}(ii, :) == exponent{contsimplex})*/\
			for (jj = 0; jj < dimensionsEXPTABLEElement[1]; ++jj) {\
				/*table contains value in current column*/\
				if (jj >= exptableContsimplexNumericSparseIRCSR[exptableContsimplexNumericSparseJCCSR[ii]] && exptableContsimplexNumericSparseJCCSR[ii + 1] > 0 && jj <= exptableContsimplexNumericSparseIRCSR[exptableContsimplexNumericSparseJCCSR[ii + 1] - 1]) {\
					if (exptableContsimplexNumericSparseIRCSR[kk + exptableContsimplexNumericSparseJCCSR[ii]] == jj) {\
						/*value is equal*/\
						if (((double) exponentContsimplexNumeric##type[jj]) != exptableContsimplexNumericSparseDoubleCSR[kk + exptableContsimplexNumericSparseJCCSR[ii]]) {\
							ismatchvalue = false;\
							ismatchcolumns = false;\
							break;\
						}\
						++kk;\
					}\
				}\
				else {\
					/*value not in table and exponent must be 0*/\
					if ((double) exponentContsimplexNumeric##type[jj] != 0.0) {\
						ismatchvalue = false;\
						break;\
					}\
				}\
				if (!ismatchvalue) {\
					ismatch = false;\
					break;\
				}\
				if (!ismatchcolumns) {\
					ismatch = false;\
					break;\
				}\
			}\
			if (!ismatchcolumns) {\
				ismatch = false;\
			}\
			if (!ismatchvalue) {\
				ismatch = false;\
			}\
			if (ismatch) {\
				/*index = index + (ii - 1)*jump(contsimplex);*/\
				index = index + ((double)ii)*jumpContsimplexNumeric[countsimplex];\
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
	mwSize countsimplex = 0, ii = 0, jj = 0, kk = 0, ll = 0;
	const mwSize *dimensionsEXPONENT = NULL, *dimensionsEXPTABLE = NULL, *dimensionsJUMP = NULL, *dimensionsEXPONENTElement = NULL, *dimensionsEXPTABLEElement = NULL;
	//mwSize nzmaxEXPONENT = 0, nzmaxEXPTABLE = 0;
	size_t lengthExponent = 0;
	const mxArray *exptableContsimplex = NULL, *exponentContsimplex = NULL;
	const double *exponentContsimplexNumericDouble = NULL, *exptableContsimplexNumericDouble = NULL, *jumpContsimplexNumeric = NULL;
	const float *exponentContsimplexNumericFloat = NULL, *exptableContsimplexNumericFloat = NULL;
	const INT8_T *exponentContsimplexNumericInt8 = NULL, *exptableContsimplexNumericInt8 = NULL;
	const INT16_T *exponentContsimplexNumericInt16 = NULL, *exptableContsimplexNumericInt16 = NULL;
	const INT32_T *exponentContsimplexNumericInt32 = NULL, *exptableContsimplexNumericInt32 = NULL;
	const INT64_T *exponentContsimplexNumericInt64 = NULL, *exptableContsimplexNumericInt64 = NULL;
	const UINT8_T *exponentContsimplexNumericUint8 = NULL, *exptableContsimplexNumericUint8 = NULL;
	const UINT16_T *exponentContsimplexNumericUint16 = NULL, *exptableContsimplexNumericUint16 = NULL;
	const UINT32_T *exponentContsimplexNumericUint32 = NULL, *exptableContsimplexNumericUint32 = NULL;
	const UINT64_T *exponentContsimplexNumericUint64 = NULL, *exptableContsimplexNumericUint64 = NULL;
	const double *exponentContsimplexNumericSparseDouble = NULL, *exptableContsimplexNumericSparseDouble = NULL;
	const mwIndex *exponentContsimplexNumericSparseIR = NULL, *exptableContsimplexNumericSparseIR = NULL;
	const mwIndex *exponentContsimplexNumericSparseJC = NULL, *exptableContsimplexNumericSparseJC = NULL;
    
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
	/*for contsimplex=1:(length(exponent))*/
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
		/*if (length(exptable{contsimplex}) > 0)*/
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
			/*for ii = 1:size(exptable{contsimplex}, 1)*/
			for (ii = 0; ii < dimensionsEXPTABLEElement[0]; ++ii) {
				ismatch = true;
				ismatchvalue = true;
				ismatchcolumns = true;
				jj = 0;
				/*if all(exptable{contsimplex}(ii, :) == exponent{contsimplex})*/
				for (kk = 0; kk < dimensionsEXPTABLEElement[1]; ++kk) {
					// table contains value in current column
					if (kk >= exptableContsimplexNumericSparseIRCSR[exptableContsimplexNumericSparseJCCSR[ii]] && exptableContsimplexNumericSparseJCCSR[ii + 1] > 0 && kk <= exptableContsimplexNumericSparseIRCSR[exptableContsimplexNumericSparseJCCSR[ii + 1] - 1]) {
						if (exptableContsimplexNumericSparseIRCSR[jj + exptableContsimplexNumericSparseJCCSR[ii]] == kk) {
							// row in exponent is empty but not empty in exptable
							if (exponentContsimplexNumericSparseJC[kk + 1] - exponentContsimplexNumericSparseJC[kk] <= 0) {
								ismatchvalue = false;
								ismatchcolumns = false;
								break;
							}
							// value is equal
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
					/*index = index + (ii - 1)*jump(contsimplex);*/
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
			exptableContsimplexNumericSparseDouble = mxGetPr(exptableContsimplex);
			if (exptableContsimplexNumericSparseDouble == NULL) {
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
			isdouble = mxIsDouble(exponentContsimplex);
			isfloat = mxIsSingle(exponentContsimplex);
			isint8 = mxIsInt8(exponentContsimplex);
			isint16 = mxIsInt16(exponentContsimplex);
			isint32 = mxIsInt32(exponentContsimplex);
			isint64 = mxIsInt64(exponentContsimplex);
			isuint8 = mxIsUint8(exponentContsimplex);
			isuint16 = mxIsUint16(exponentContsimplex);
			isuint32 = mxIsUint32(exponentContsimplex);
			isuint64 = mxIsUint64(exponentContsimplex);
			if (!isdouble && !isfloat && !isint8 && !isint16 && !isint32 && !isint64 && !isuint8 && !isuint16 && !isuint32 && !isuint64) {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Element %d of exptable must be numeric and of the same type as exponent.", countsimplex + 1);
				return;
			}
			if (isdouble) {
				exponentContsimplexNumericDouble = mxGetPr(exponentContsimplex);
				if (exponentContsimplexNumericDouble == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEXSPARSE(Double)
			}
			else if (isfloat) {
				exponentContsimplexNumericFloat = (float*) mxGetData(exponentContsimplex);
				if (exponentContsimplexNumericFloat == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEXSPARSE(Float)
			}
			else if (isint8) {
				exponentContsimplexNumericInt8 = (INT8_T*) mxGetData(exponentContsimplex);
				if (exponentContsimplexNumericInt8 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEXSPARSE(Int8)
			}
			else if (isint16) {
				exponentContsimplexNumericInt16 = (INT16_T*) mxGetData(exponentContsimplex);
				if (exponentContsimplexNumericInt16 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEXSPARSE(Int16)
			}
			else if (isint32) {
				exponentContsimplexNumericInt32 = (INT32_T*) mxGetData(exponentContsimplex);
				if (exponentContsimplexNumericInt32 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEXSPARSE(Int32)
			}
			else if (isint64) {
				exponentContsimplexNumericInt64 = (INT64_T*) mxGetData(exponentContsimplex);
				if (exponentContsimplexNumericInt64 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEXSPARSE(Int64)
			}
			else if (isuint8) {
				exponentContsimplexNumericUint8 = (UINT8_T*) mxGetData(exponentContsimplex);
				if (exponentContsimplexNumericUint8 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEXSPARSE(Uint8)
			}
			else if (isuint16) {
				exponentContsimplexNumericUint16 = (UINT16_T*) mxGetData(exponentContsimplex);
				if (exponentContsimplexNumericUint16 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEXSPARSE(Uint16)
			}
			else if (isuint32) {
				exponentContsimplexNumericUint32 = (UINT32_T*) mxGetData(exponentContsimplex);
				if (exponentContsimplexNumericUint32 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEXSPARSE(Uint32)
			}
			else if (isuint64) {
				exponentContsimplexNumericUint64 = (unsigned long long*) mxGetData(exponentContsimplex);
				if (exponentContsimplexNumericUint64 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEXSPARSE(Uint64)
			}
			else {
				mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Unknown data for element %d.", countsimplex + 1);
				return;
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
				exponentContsimplexNumericInt8 = (INT8_T*) mxGetData(exponentContsimplex);
				exptableContsimplexNumericInt8 = (INT8_T*) mxGetData(exptableContsimplex);
				if (exponentContsimplexNumericInt8 == NULL || exptableContsimplexNumericInt8 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEX(Int8)
			}
			else if (isint16) {
				exponentContsimplexNumericInt16 = (INT16_T*) mxGetData(exponentContsimplex);
				exptableContsimplexNumericInt16 = (INT16_T*) mxGetData(exptableContsimplex);
				if (exponentContsimplexNumericInt16 == NULL || exptableContsimplexNumericInt16 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEX(Int16)
			}
			else if (isint32) {
				exponentContsimplexNumericInt32 = (INT32_T*) mxGetData(exponentContsimplex);
				exptableContsimplexNumericInt32 = (INT32_T*) mxGetData(exptableContsimplex);
				if (exponentContsimplexNumericInt32 == NULL || exptableContsimplexNumericInt32 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEX(Int32)
			}
			else if (isint64) {
				exponentContsimplexNumericInt64 = (INT64_T*) mxGetData(exponentContsimplex);
				exptableContsimplexNumericInt64 = (INT64_T*) mxGetData(exptableContsimplex);
				if (exponentContsimplexNumericInt64 == NULL || exptableContsimplexNumericInt64 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEX(Int64)
			}
			else if (isuint8) {
				exponentContsimplexNumericUint8 = (UINT8_T*) mxGetData(exponentContsimplex);
				exptableContsimplexNumericUint8 = (UINT8_T*) mxGetData(exptableContsimplex);
				if (exponentContsimplexNumericUint8 == NULL || exptableContsimplexNumericUint8 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEX(Uint8)
			}
			else if (isuint16) {
				exponentContsimplexNumericUint16 = (UINT16_T*) mxGetData(exponentContsimplex);
				exptableContsimplexNumericUint16 = (UINT16_T*) mxGetData(exptableContsimplex);
				if (exponentContsimplexNumericUint16 == NULL || exptableContsimplexNumericUint16 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEX(Uint16)
			}
			else if (isuint32) {
				exponentContsimplexNumericUint32 = (UINT32_T*) mxGetData(exponentContsimplex);
				exptableContsimplexNumericUint32 = (UINT32_T*) mxGetData(exptableContsimplex);
				if (exponentContsimplexNumericUint32 == NULL || exptableContsimplexNumericUint32 == NULL) {
					mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
					return;
				}
				GETINDEX(Uint32)
			}
			else if (isuint64) {
				exponentContsimplexNumericUint64 = (UINT64_T*) mxGetData(exponentContsimplex);
				exptableContsimplexNumericUint64 = (UINT64_T*) mxGetData(exptableContsimplex);
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
	/*index = index + 1;*/
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
