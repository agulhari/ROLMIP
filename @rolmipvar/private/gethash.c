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
			//nzmaxEXPONENT = mxGetNzmax(exponentContsimplex);
			//nzmaxEXPTABLE = mxGetNzmax(exptableContsimplex);
			//if (nzmaxEXPONENT != nzmaxEXPTABLE) {
			//	mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
			//	return;
			//}
			for (ii = 0; ii < dimensionsEXPTABLEElement[0]; ++ii) {
				ismatch = true;
				for (jj = 0; jj < dimensionsEXPTABLEElement[1]; ++jj) {
					ismatchvalue = true;
					for (kk = exponentContsimplexNumericSparseJC[jj]; kk < exponentContsimplexNumericSparseJC[jj + 1]; ++kk) {
						ismatchcolumns = false;
						for (ll = exptableContsimplexNumericSparseJC[jj]; ll < exptableContsimplexNumericSparseJC[jj + 1]; ++ll) {
							if (exptableContsimplexNumericSparseIR[ll] == ii) {
								ismatchcolumns = true;
								//if (exptableContsimplexNumericSparseJC[jj + 1] - exptableContsimplexNumericSparseJC[jj] != exponentContsimplexNumericSparseJC[jj + 1] - exponentContsimplexNumericSparseJC[jj]) {
								//	ismatchvalue = false;
								//	break;
								//}
								if (exptableContsimplexNumericSparseDouble[ll] != exponentContsimplexNumericSparseDouble[kk]) {
									ismatchvalue = false;
									break;
								}
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
					if (!ismatch) {
						ismatchvalue = false;
						break;
					}
					if (!ismatchvalue) {
						ismatch = false;
						break;
					}
				}
				if (ismatch) {
					index = ((double)ii)*jumpContsimplexNumeric[countsimplex];
					break;
				}
			}
		}
		else if (!mxIsSparse(exponentContsimplex) && mxIsSparse(exptableContsimplex)) {
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
			//nzmaxEXPTABLE = mxGetNzmax(exptableContsimplex);
			//if (nzmaxEXPONENT != nzmaxEXPTABLE) {
			//	mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
			//	return;
			//}
			for (ii = 0; ii < dimensionsEXPTABLEElement[0]; ++ii) {
				ismatch = true;
				for (jj = 0; jj < dimensionsEXPTABLEElement[1]; ++jj) {
					ismatchvalue = true;
					ismatchcolumns = false;
					for (ll = exptableContsimplexNumericSparseJC[jj]; ll < exptableContsimplexNumericSparseJC[jj + 1]; ++ll) {
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
					}
					if (!ismatchcolumns && exponentContsimplexNumericDouble[jj] != 0.0) {
						ismatchvalue = false;
						ismatch = false;
						break;
					}
					if (!ismatchvalue) {
						ismatch = false;
						break;
					}
					if (!ismatch) {
						break;
					}
				}
				if (ismatch) {
					index = ((double)ii)*jumpContsimplexNumeric[countsimplex];
					break;
				}
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
