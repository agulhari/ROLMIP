/*=================================================================
 * mxmalloc.c
 * 
 * This function takes a MATLAB string as an argument and copies it in 
 * NULL terminated ANSI C string. Use this function only with 
 * MATLAB strings represented in single-byte encoding schemes.
 *
 * This is a MEX-file for MATLAB.  
 * Copyright 1984-2014 The MathWorks, Inc.
 *=================================================================*/
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
	boolean_T ismatch = false;
	double *indexptr = NULL;
	unsigned int countsimplex = 0, ii = 0, jj = 0;
	const mwSize *dimensionsEXPONENT = NULL, *dimensionsEXPTABLE = NULL, *dimensionsJUMP = NULL, *dimensionsEXPONENTElement = NULL, *dimensionsEXPTABLEElement = NULL;
	size_t lengthExponent = 0;
	const mxArray *exptableContsimplex = NULL, *exponentContsimplex = NULL;
	const double *exponentContsimplexNumeric = NULL, *exptableContsimplexNumeric = NULL, *jumpContsimplexNumeric = NULL;
    
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
		if (!mxIsDouble(exponentContsimplex) || mxIsComplex(exponentContsimplex)) {
			mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Element %d of exponent must be of type 'double'.", countsimplex + 1);
			return;
		}
		if (!mxIsDouble(exptableContsimplex) || mxIsComplex(exptableContsimplex)) {
			mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Element %d of exptable must be of type 'double'.", countsimplex + 1);
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
		exponentContsimplexNumeric = mxGetPr(exponentContsimplex);
		exptableContsimplexNumeric = mxGetPr(exptableContsimplex);
		if (exponentContsimplexNumeric == NULL || exptableContsimplexNumeric == NULL) {
			mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Numeric data for element %d is invalid.", countsimplex + 1);
			return;
		}
		for (ii = 0; ii < dimensionsEXPTABLEElement[0]; ++ii) {
			ismatch = true;
			for (jj = 0; jj < dimensionsEXPTABLEElement[1]; ++jj) {
				if (exptableContsimplexNumeric[jj*(dimensionsEXPTABLEElement[0]) + ii] != exponentContsimplexNumeric[jj]) {
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
	index = index + 1;
	pargINDEX = mxCreateDoubleMatrix(1, 1, mxREAL);
	if (pargINDEX == NULL) {
		mexErrMsgIdAndTxt("ROLMIP:gethash:input", "Out of Memory.");
		return;
	}
	indexptr = mxGetPr(pargINDEX);
	indexptr[0] = index;
}
