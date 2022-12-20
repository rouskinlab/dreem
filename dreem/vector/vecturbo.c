/*

>>>>==========================>>>>
>>>>======== VecTurbo ========>>>>
>>>>==========================>>>>

C extension module for faster vectoring than possible with pure Python.
Part of the DREEM Pipeline by Rouskin Lab
Author: Matty Allan
Date: 2023-??-??

For documentation on writing C extension modules for Python, see link:
https://docs.python.org/3/extending/extending.html

*/

// Include the Python API so that Python can call from this module.
#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Other included modules.
#include <stdio.h>
#include <regex.h>


/*
Define constants.
Note that these constants MUST match the constants defined in the Python code.
*/

// CIGAR operations
const char CIG_ALN = 'M';  // match or substitution
const char CIG_MAT = '=';  // match
const char CIG_SUB = 'X';  // substitution
const char CIG_DEL = 'D';  // deletion
const char CIG_INS = 'I';  // insertion
const char CIG_SCL = 'S';  // soft clipping
const char CIG_HCL = 'H';  // hard clipping


int main() {
    printf("VecTurbo module of the DREEM pipeline.\n");
    return 0;
}
