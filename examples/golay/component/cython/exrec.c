/*
 * Golay exRec counting functions.  This file contains the function countExRecForConfigs()
 * which is designed as a Python extension.
 */

#include <stdio.h>
#include <Python.h>
#include <gmp.h>

PyObject * PyLong_FromMPZ(mpz_t z);
char * mpz2str(mpz_t z);
void mpz_set_PyLong(mpz_t z, PyObject * pyLong);
PyObject * PyList_FromMpzList(mpz_t zList[], int len);
mpz_t * pyList2mpz(PyObject * pyList);
int * pyList2int(PyObject * pyList);
int ** pyTable2Int(PyObject * pyTable);
mpz_t ** pyTable2Mpz(PyObject * pyTable);
void freeIntTable(int ** table, int len);

void freeMpzTable(mpz_t ** table, PyObject* pyTable);

mpz_t ** countExRecForConfig(int * sLECaList,
						 int len_sLECaList,
						 int * sLECbList,
						 int len_sLECbList,
						 int * sCList,
						 int len_sCList,
						 mpz_t * countsLECa,
						 mpz_t * countsLECb,
						 mpz_t * countsC,
						 mpz_t ** lcountsTECa,
						 mpz_t ** lcountsTECb);

void checkPyError(void);

/* Hard-coded table indexed by (logical) syndrome that contains the result of decoding that syndrome.
 * A 1 represents a logical error, a 0 represents no logical error.
 */
static const char decodedSyndromes[] = {
0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,1,1,0,0,1,0,1,1,1,0,1,1,1,1,1,1,1,1,0,0,1,0,1,1,1,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,
1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,0,1,1,1,1,1,0,0,1,0,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,0,1,1,1,1,1,1,0,
1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,0,1,0,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,
1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,
1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,1,1,1,1,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,1,1,0,1,1,0,1,1,1,0,
1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,1,0,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,1,1,1,1,
1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,1,1,
1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,0,1,1,0,1,1,0,
1,1,1,1,1,1,0,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,0,1,0,1,1,0,1,1,1,1,1,1,0,1,1,0,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,
1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,0,1,1,0,1,1,0,1,1,1,1,1,1,1,0,1,
1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,0,1,0,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,
1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,0,1,1,0,1,1,0,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,0,1,0,0,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,
1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,
1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,0,0,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,0,1,1,
1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,
1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,0,1,1,1,1,1,0,1,1,0,1,1,0,1,
1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,
1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,0,1,1,1,0,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,0,1,1,1,1,0,1,1,1,0,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
1,1,1,1,1,1,0,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,0,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,
1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,0,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,
1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1,1,1,1,1,1,1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,0,0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,
0,0,0,0,0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,
0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,
0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,1,0,0,
1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,1,1,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,
0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,
0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,
0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0
};

/* Get entry i of Python list l as a long type */
#define PYGETLONG(l, i) PyLong_AsUnsignedLong(PyList_GetItem(l, i))

static const int MASK12 = (1<<12) - 1;



/*
 * Counts the ways in which logical X errors may occur in a Golay encoded CNOT 1-Rec for the given
 * set of failure configurations.
 * Returns a list of size 3.  The first element of the list contains counts for the logical XI error,
 * the second contains logical IX error counts, and the third contains logical XX error counts.  Each
 * element is itself a list of counts, indexed by the number of faults k.
 *
 * This function first converts all of the Python arguments into C data structures and then computes
 * each configuration in pure C.  For large configuration sets, this implementation is much faster
 * than an equivalent Python implementation.
 *
 * The resulting error counts are computed as arbitrary precision integers and returned as Python longs
 * (which support arbitrary precision natively).
 *
 * Input Arguments
 * ---------------
 * sLEClists - A list of syndromes, indexed by k, for which countsLEC[k][s] !=0
 * sClists - A list of syndromes, indexed by k, for which countsC[k][s] !=0
 * countsLEC - A list of counts for each output syndrome of the leading EC, indexed by k.
 * countsC - A list of counts for each output syndrome of the transversal CNOT, indexed by k.
 * countsECLogical - A list of pairs counts for each input syndrome to the trailing EC, indexed by k.
 * 					 each pair represents the counts for a correctable (FALSE) or uncorrectable (TRUE)
 * 					 error at the output of the trailing EC.
 * configs - A list of failure configurations to count.
 * 			 Each configuration specifies [kLECa, kLECb, kC, kTECa, kTECb]
 * kMax - max(sum(configs))
 */
PyObject * countExRecForConfigs(PyObject * sLEClists,
								PyObject * sClists,
								PyObject * countsLEC,
								PyObject * countsC,
								PyObject * countsECLogical,
								PyObject * configs,
								PyObject * kMax)
{

	//printf("0\n");
	/* Convert all of the Python objects into native C data types */
	int kMax_int = PyInt_AsLong(kMax);

	mpz_t maligCounts[4][3][kMax_int + 1];

	for(int ec=0; ec < 4; ec++) {
		for (int error=0; error < 3; error++) {
			for (int k=0; k < kMax_int + 1; k++) {
				mpz_init(maligCounts[ec][error][k]);
				mpz_set_ui(maligCounts[ec][error][k], 0);
			}
		}
	}

	int ** sLEClists_int = pyTable2Int(sLEClists);
	int ** sClists_int = pyTable2Int(sClists);
	mpz_t ** countsLEC_mpz = pyTable2Mpz(countsLEC);
	mpz_t ** countsC_mpz = pyTable2Mpz(countsC);

	int len = PyList_Size(countsECLogical);
	mpz_t *** countsECLogical_mpz = malloc(sizeof(mpz_t *) * len);
	for(int k=0; k < len; k++) {
		PyObject * countsLogicalK = PyList_GetItem(countsECLogical, k);
		countsECLogical_mpz[k] = pyTable2Mpz(countsLogicalK);
	}

	//printf("1\n");

	/* Now that everything is converted, do the actual computation on the assigned configurations */
	for(int iConfig=0; iConfig < PyList_Size(configs); iConfig++) {
		PyObject * config = PyList_GetItem(configs,iConfig);
		int kLECa = PyInt_AsLong(PyList_GetItem(config, 0));
		int kLECb = PyInt_AsLong(PyList_GetItem(config, 1));
		int kC = PyInt_AsLong(PyList_GetItem(config, 2));
		int kTECa = PyInt_AsLong(PyList_GetItem(config, 3));
		int kTECb = PyInt_AsLong(PyList_GetItem(config, 4));
		int k = kLECa + kLECb + kC + kTECa + kTECb;

		//printf("1.7 - (%d, %d, %d, %d, %d) - %d\n", kLECa, kLECb, kC, kTECa, kTECb, k);

		int len_a = PyList_Size(PyList_GetItem(sLEClists, kLECa));
		int len_b = PyList_Size(PyList_GetItem(sLEClists, kLECb));
		int len_c = PyList_Size(PyList_GetItem(sClists, kC));

		/* result is a 4x3 array indexed by [ec][error] */
		mpz_t ** result = countExRecForConfig(sLEClists_int[kLECa], len_a,
											 sLEClists_int[kLECb], len_b,
											 sClists_int[kC], len_c,
											 countsLEC_mpz[kLECa],
											 countsLEC_mpz[kLECb],
											 countsC_mpz[kC],
											 countsECLogical_mpz[kTECa],
											 countsECLogical_mpz[kTECb]
											 );

		for (int ec = 0; ec < 4; ec++) {
			for (int error = 0; error < 3; error++) {
				//printf("%d, %d, %d\n", ec, error, k);
				mpz_add(maligCounts[ec][error][k], maligCounts[ec][error][k], result[ec][error]);
				mpz_clear(result[ec][error]);
			}
			free(result[ec]);
		}

		free(result);
	}

	//printf("2\n");

	/* Package the result for Python */
	PyObject * resultList = PyList_New(4);
	for (int ec = 0; ec < 4; ec++) {
		PyObject * errorList = PyList_New(3);
		PyList_SetItem(resultList, ec, errorList);
		for (int error = 0; error < 3; error++) {
			PyList_SetItem(errorList, error, PyList_FromMpzList(maligCounts[ec][error], kMax_int + 1));
		}
	}


	/* Release all allocated memory */
	freeIntTable(sLEClists_int, PyList_Size(sLEClists));
	freeIntTable(sClists_int, PyList_Size(sClists));
	freeMpzTable(countsLEC_mpz, countsLEC);
	freeMpzTable(countsC_mpz, countsC);

	for(int k=0; k < len; k++) {
		PyObject * countsLogicalK = PyList_GetItem(countsECLogical, k);
		freeMpzTable(countsECLogical_mpz[k], countsLogicalK);
	}

	free(countsECLogical_mpz);

	for(int ec=0; ec < 4; ec++) {
		for (int error=0; error < 3; error++) {
			for (int k=0; k < kMax_int + 1; k++) {
				mpz_clear(maligCounts[ec][error][k]);
			}
		}
	}



	return resultList;

}

void checkPyError(void) {
	if(NULL != PyErr_Occurred()) {
		printf("Python Exception occurred\n");
	}
}

/* Converts a Python list of PyLong items into a corresponding array of mpz_t items. */
mpz_t * pyList2mpz(PyObject * pyList) {
	int len = PyList_Size(pyList);
	mpz_t * mpzList = malloc(sizeof(mpz_t) * len);
	for(int i=0; i < len; i++) {
		mpz_init(mpzList[i]);
		mpz_set_PyLong(mpzList[i], PyList_GetItem(pyList, i));
	}

	return mpzList;
}

/*
 * Converts a list of Python integers into a corresponding array of C integers.
 */
int * pyList2int(PyObject * pyList) {
	int len = PyList_Size(pyList);
	int * intList = malloc(sizeof(int) * len);
	for(int i=0; i < len; i++) {
		intList[i] = PYGETLONG(pyList, i);
	}

	return intList;
}

/*
 * Converts a 2-D Python list of integers into a corresponding 2-D array of C integers.
 */
int ** pyTable2Int(PyObject * pyTable) {
	int len = PyList_Size(pyTable);
	int ** intTable = malloc(sizeof(int *) * len);

	for(int i=0; i < len; i++) {
		intTable[i] = pyList2int(PyList_GetItem(pyTable, i));
	}

	return intTable;
}
/*
 * Converts a 2-D Python list of PyLong items into a corresponding 2-D array of C mpz_t items.
 */
mpz_t ** pyTable2Mpz(PyObject * pyTable) {
	int len = PyList_Size(pyTable);
	mpz_t ** mpzTable = malloc(sizeof(mpz_t *) * len);

	for(int i=0; i < len; i++) {
		mpzTable[i] = pyList2mpz(PyList_GetItem(pyTable, i));
	}

	return mpzTable;
}

/*
 * Free memory associated with a 2-D array of integers where the first dimension is len.
 */
void freeIntTable(int ** table, int len) {
	for(int i=0; i < len; i++) {
		free(table[i]);
	}

	free(table);
}

/*
 * Free memory associated with a 2-D array of mpz_t that corresponds to the 2-D Python list pyTable.
 */
void freeMpzTable(mpz_t ** table, PyObject * pyTable) {
	int len = PyList_Size(pyTable);
	for(int i=0; i < len; i++) {
		PyObject * pyList = PyList_GetItem(pyTable, i);
		int len2 = PyList_Size(pyList);
		for(int j=0; j < len2; j++) {
			mpz_clear(table[i][j]);
		}
		free(table[i]);
	}
	free(table);
}






/*
 * Counts the logical XI, IX, and XX errors of a 1-Rec for a particular failure configuration.
 * This function is designed as a Python extension.  However, it turns out to be much faster to
 * use countExRecForConfigs, instead.
 */
PyObject * countExRecForConfigPy(PyObject * sLECaList,
						 PyObject * sLECbList,
						 PyObject * sCList,
						 PyObject * countsLECa,
						 PyObject * countsLECb,
						 PyObject * countsC,
						 PyObject * lcountsTECa,
						 PyObject * lcountsTECb)
{
	mpz_t errorCountIX, errorCountXI, errorCountXX;
	mpz_t errorCountCIX, errorCountCXI, errorCountCXX;
	mpz_t errorCountAIX, errorCountAXI, errorCountAXX;
	mpz_t countLECa, countLECb, countC;
	mpz_t tECaBenign, tECaMalig, tECbBenign, tECbMalig;

	mpz_init(errorCountIX);
	mpz_init(errorCountXI);
	mpz_init(errorCountXX);

	mpz_init(errorCountCIX);
	mpz_init(errorCountCXI);
	mpz_init(errorCountCXX);

	mpz_init(errorCountAIX);
	mpz_init(errorCountAXI);
	mpz_init(errorCountAXX);

	mpz_init(countLECa);
	mpz_init(countLECb);
	mpz_init(countC);
	mpz_init(tECaBenign);
	mpz_init(tECaMalig);
	mpz_init(tECbBenign);
	mpz_init(tECbMalig);

	PyObject * tECaCount, * tECbCount;
	int sC, sCa, sCb, sLECa, sA, sB1, sLECb;
	int decodedLECa, expectedTECb;
	int iC, iLECa, iLECb;

	int len_sLECaList = PyList_Size(sLECaList);
	int len_sLECbList = PyList_Size(sLECbList);
	int len_sCList = PyList_Size(sCList);
	//printf("Computing %dx%dx%d\n", len_sCList, len_sLECaList, len_sLECbList);


	/* Iterate over all possible syndromes that can be produced at the input to the two trailing error
		corrections. Then, for each syndrome, determine the weighted count of malignant sets. */

	/* Transversal CNOT of the exRec */
	for(iC=0; iC < len_sCList; iC++)
	{

		mpz_set_ui(errorCountCIX, 0);
		mpz_set_ui(errorCountCXI, 0);
		mpz_set_ui(errorCountCXX, 0);

		sC = PYGETLONG(sCList, iC);
		sCa = sC >> 12;
		sCb = sC & MASK12;

		mpz_set_PyLong(countC, PyList_GET_ITEM(countsC, sC));
		//printf("At 1 - %d, %d, count=%s\n", iC, sC, mpz2str(countC));

		//printf("At 2.1\n");
		for(iLECa=0; iLECa < len_sLECaList; iLECa++)
		{
			//printf("At 2.2 - %d\n", iLECa);
			sLECa = PYGETLONG(sLECaList,iLECa);
			sA = sLECa ^ sCa;
			sB1 = sLECa ^ sCb;

			decodedLECa = decodedSyndromes[sLECa];

			mpz_set_PyLong(countLECa, PyList_GET_ITEM(countsLECa, sLECa));

			/* Get the counts for each logical outcome of the trailing EC */
			tECaCount = PyList_GET_ITEM(lcountsTECa, sA);
			mpz_set_PyLong(tECaBenign, PyList_GET_ITEM(tECaCount, decodedLECa));
			mpz_set_PyLong(tECaMalig, PyList_GET_ITEM(tECaCount, !decodedLECa));

			mpz_mul(tECaBenign, tECaBenign, countLECa);
			mpz_mul(tECaMalig, tECaMalig, countLECa);

			mpz_set_ui(errorCountAIX, 0);
			mpz_set_ui(errorCountAXI, 0);
			mpz_set_ui(errorCountAXX, 0);

			for(iLECb=0; iLECb < len_sLECbList; iLECb++)
			{
				//printf("At 2.3 - %d\n", iLECb);
				sLECb = PYGETLONG(sLECbList,iLECb);
				//sLECb = sLECbList_int[iLECb];
				expectedTECb = decodedSyndromes[sLECb] ^ decodedLECa;

				mpz_set_PyLong(countLECb, PyList_GET_ITEM(countsLECb, sLECb));
				//countLECb = countLECbList_mpz[iLECb];

				tECbCount = PyList_GET_ITEM(lcountsTECb, sB1 ^ sLECb);
				mpz_set_PyLong(tECbBenign, PyList_GET_ITEM(tECbCount, expectedTECb));
				mpz_set_PyLong(tECbMalig, PyList_GET_ITEM(tECbCount, !expectedTECb));

				mpz_mul(tECbBenign, tECbBenign, countLECb);
				mpz_mul(tECbMalig, tECbMalig, countLECb);

				mpz_add(errorCountAIX, errorCountAIX, tECbMalig);
				mpz_add(errorCountAXI, errorCountAXI, tECbBenign);
				mpz_add(errorCountAXX, errorCountAXX, tECbMalig);

			}


			mpz_mul(errorCountAIX, errorCountAIX, tECaBenign);
			mpz_mul(errorCountAXI, errorCountAXI, tECaMalig);
			mpz_mul(errorCountAXX, errorCountAXX, tECaMalig);

			//printf("ecAIX=%s\n", mpz2str(errorCountAIX));
			//printf("ecAXI=%s\n", mpz2str(errorCountAXI));
			//printf("ecAXX=%s\n", mpz2str(errorCountAXX));


			mpz_add(errorCountCIX, errorCountCIX, errorCountAIX);
			mpz_add(errorCountCXI, errorCountCXI, errorCountAXI);
			mpz_add(errorCountCXX, errorCountCXX, errorCountAXX);

		}

		mpz_mul(errorCountCIX, errorCountCIX, countC);
		mpz_mul(errorCountCXI, errorCountCXI, countC);
		mpz_mul(errorCountCXX, errorCountCXX, countC);

		mpz_add(errorCountIX, errorCountIX, errorCountCIX);
		mpz_add(errorCountXI, errorCountXI, errorCountCXI);
		mpz_add(errorCountXX, errorCountXX, errorCountCXX);
	}

	//printf("c-result is %s, %s, %s\n", mpz2str(errorCountIX), mpz2str(errorCountXI), mpz2str(errorCountXX));
	PyObject * resultList = PyList_New(3);
	PyList_SetItem(resultList, 0, PyLong_FromMPZ(errorCountIX));
	PyList_SetItem(resultList, 1, PyLong_FromMPZ(errorCountXI));
	PyList_SetItem(resultList, 2, PyLong_FromMPZ(errorCountXX));

	mpz_clear(errorCountIX);
	mpz_clear(errorCountXI);
	mpz_clear(errorCountXX);

	mpz_clear(errorCountCIX);
	mpz_clear(errorCountCXI);
	mpz_clear(errorCountCXX);

	mpz_clear(errorCountAIX);
	mpz_clear(errorCountAXI);
	mpz_clear(errorCountAXX);

	mpz_clear(countLECa);
	mpz_clear(countLECb);
	mpz_clear(tECaBenign);
	mpz_clear(tECaMalig);
	mpz_clear(tECbBenign);
	mpz_clear(tECbMalig);

	return resultList;
}


#define BENIGN 0
#define MALIG  1
#define EXISTS 1
#define NOT_EXIST 0


/*
 * Counts the logical XI, IX, and XX errors of a 1-Rec for a particular failure configuration.
 */
mpz_t ** countExRecForConfig(int * sLECaList,
						 int len_sLECaList,
						 int * sLECbList,
						 int len_sLECbList,
						 int * sCList,
						 int len_sCList,
						 mpz_t * countsLECa,
						 mpz_t * countsLECb,
						 mpz_t * countsC,
						 mpz_t ** lcountsTECa,
						 mpz_t ** lcountsTECb)
{
	int i, ec, error;
	mpz_t maligCountsC[4][3];
	mpz_t countsBlockA[2][2];
	mpz_t sumsBlockB[2][2];
	mpz_t tECaBenign, tECaMalig, tECbBenign, tECbMalig;
	mpz_t tempMPZ;


	/* Indexed as [ec][error], where:
	 * ec is a value 0-3.  A 1 in the MSB means the TEC on block A is intact.  A 1 in the LSB means
	 * the TEC on block B is intact.
	 * error is a value 0-2. 0 = IX, 1 = XI, 2 = XX
	 */
	mpz_t ** maligCounts = malloc(sizeof(mpz_t *) * 4);
	for (ec=0; ec < 4; ec++) {
		maligCounts[ec] = malloc(sizeof(mpz_t) * 3);
		for (error=0; error < 3; error++) {
			mpz_init(maligCounts[ec][error]);
			mpz_init(maligCountsC[ec][error]);
		}
	}

	for (ec = 0; ec < 2; ec++) {
		for (error = 0; error < 2; error++) {
			mpz_init(countsBlockA[ec][error]);
			mpz_init(sumsBlockB[ec][error]);
		}
	}

	mpz_init(tECaBenign);
	mpz_init(tECaMalig);
	mpz_init(tECbBenign);
	mpz_init(tECbMalig);
	mpz_init(tempMPZ);

	int sC, sCa, sCb, sLECa, sA, sB1, sLECb;
	int sAisMalig, sB, sBisMalig;
	int ecA, ecB, eA, eB;
	int decodedLECa, expectedTECb;
	int iC, iLECa, iLECb;

	//printf("Computing %dx%dx%d\n", len_sCList, len_sLECaList, len_sLECbList);


	/* Iterate over all possible syndromes that can be produced at the input to the two trailing error
		corrections. Then, for each syndrome, determine the weighted count of malignant sets. */


	/* Transversal CNOT of the exRec */
	for(iC=0; iC < len_sCList; iC++)
	{
		//printf("At 0 - %d", iC);
		for (ec=0; ec < 4; ec++) {
			for (error = 0; error < 3; error++) {
				mpz_set_ui(maligCountsC[ec][error], 0);
			}
		}
		//printf("At 0.1 - %d\n", iC);

		sC = sCList[iC];
		sCa = sC >> 12;
		sCb = sC & MASK12;

		//printf("At 0.2 - %d, %d\n", iC, sC);

		//printf("At 1 - %d, %d, count=%s\n", iC, sC, mpz2str(countsC[sC]));

		//printf("At 2.1\n");
		for(iLECa=0; iLECa < len_sLECaList; iLECa++)
		{
			//printf("At 2.2 - %d\n", iLECa);
			sLECa = sLECaList[iLECa];
			sA = sLECa ^ sCa;
			sB1 = sLECa ^ sCb;

			decodedLECa = decodedSyndromes[sLECa];
			sAisMalig = (decodedSyndromes[sA] ^ decodedLECa);

			mpz_set(countsBlockA[NOT_EXIST][sAisMalig], countsLECa[sLECa]);
			mpz_set_ui(countsBlockA[NOT_EXIST][!sAisMalig], 0);

			/* Get the counts for each logical outcome of the trailing EC */
			mpz_mul(countsBlockA[EXISTS][BENIGN], lcountsTECa[sA][decodedLECa], countsLECa[sLECa]);
			mpz_mul(countsBlockA[EXISTS][MALIG], lcountsTECa[sA][!decodedLECa], countsLECa[sLECa]);

			//printf("Block A: sLECa=%d, sC=%d, sA=%d, decLECa=%d, decSA=%d, sAisMalig=%d\n", sLECa, sC, sA, decodedLECa, decodedSyndromes[sA], sAisMalig);

			for (i = 0; i < 4; i++) {
				mpz_set_ui(((mpz_t *)sumsBlockB)[i], 0);
			}

			for(iLECb=0; iLECb < len_sLECbList; iLECb++)
			{
				//printf("At 2.3 - %d\n", iLECb);
				sLECb = sLECbList[iLECb];
				//sLECb = sLECbList_int[iLECb];
				expectedTECb = decodedSyndromes[sLECb] ^ decodedLECa;

				sB = sB1 ^ sLECb;
				sBisMalig = (decodedSyndromes[sB] ^ expectedTECb);

				/* First, the case in which the trailing EC has been removed. */
				mpz_add(sumsBlockB[NOT_EXIST][sBisMalig], sumsBlockB[NOT_EXIST][sBisMalig], countsLECb[sLECb]);

				/* Now the case in which the trailing EC is left intact. */
				mpz_mul(tECbBenign, lcountsTECb[sB][expectedTECb], countsLECb[sLECb]);
				mpz_mul(tECbMalig, lcountsTECb[sB][!expectedTECb], countsLECb[sLECb]);

				mpz_add(sumsBlockB[EXISTS][MALIG], sumsBlockB[EXISTS][MALIG], tECbMalig);
				mpz_add(sumsBlockB[EXISTS][BENIGN], sumsBlockB[EXISTS][BENIGN], tECbBenign);

			}

			// Use blockA=bit1, blockB=bit0
			for (ec=0; ec < 4; ec++) {
				ecA = ec >> 1;
				ecB = ec & 1;
				for (error=0; error < 3; error++) {
					eA = (error + 1) >> 1;
					eB = (error + 1) & 1;
					mpz_mul(tempMPZ, countsBlockA[ecA][eA], sumsBlockB[ecB][eB]);
					mpz_add(maligCountsC[ec][error], maligCountsC[ec][error], tempMPZ);
				}
			}
		}

		for (ec=0; ec < 4; ec++) {
			for (error=0; error < 3; error++) {
				mpz_mul(tempMPZ, maligCountsC[ec][error], countsC[sC]);
				mpz_add(maligCounts[ec][error], maligCounts[ec][error], tempMPZ);
			}
		}
	}

	//printf("c-result is %s, %s, %s\n", mpz2str(errorCountIX), mpz2str(errorCountXI), mpz2str(errorCountXX));

	for (ec = 0; ec < 4; ec++) {
		for (error = 0; error < 3; error++) {
			mpz_clear(maligCountsC[ec][error]);
		}
	}

	for (ec = 0; ec < 2; ec++) {
		for (error = 0; error < 2; error++) {
			mpz_clear(countsBlockA[ec][error]);
			mpz_clear(sumsBlockB[ec][error]);
		}
	}

	mpz_clear(tECaBenign);
	mpz_clear(tECaMalig);
	mpz_clear(tECbBenign);
	mpz_clear(tECbMalig);
	mpz_clear(tempMPZ);

	return maligCounts;
}






/*
 * Convert an array of mpz_t into a corresponding Python list.
 */
PyObject * PyList_FromMpzList(mpz_t zList[], int len) {
	PyObject * pyList = PyList_New(len);
	for(int i=0; i < len; i++) {
		PyList_SetItem(pyList, i, PyLong_FromMPZ(zList[i]));
	}

	return pyList;
}

/*
 * Convert an mpz_t into a corresponding Python long.
 */
PyObject * PyLong_FromMPZ(mpz_t z)
{
    PyObject* pi_result = NULL;
    char * buffer = malloc(mpz_sizeinbase(z, 10) + 2);
    mpz_get_str(buffer, 10, z);
    //printf("mpz=%s\n", buffer);
    pi_result = PyLong_FromString(buffer, NULL, 10);
    free(buffer);

    return pi_result;
}

/*
 * Set the value of an mpz_t based on the given Python long.
 */
void mpz_set_PyLong(mpz_t z, PyObject * pyLong)
{
	PyObject * pystr = PyObject_Str(pyLong);
	const char * str = PyString_AsString(pystr);
	mpz_set_str(z, str, 10);
	Py_DECREF(pystr);
}

char g_buffer[1024];

char * mpz2str(mpz_t z) {
	//printf("At mpz2str()\n");
	mpz_get_str(g_buffer, 10, z);
	return g_buffer;
}
