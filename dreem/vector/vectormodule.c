/*

DREEM Vector Module
===================
    C extension

Auth: Matty
Date: 2023-03-20

*/


#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


// Exceptions

static PyObject *VectorError;

/*
Return a description of each error code, to be used by VectorError.

Error code categories:
100s: syntax errors in SAM lines (excluding CIGAR operations)
200s: incompatible fields between two paired reads
300s: syntax errors in CIGAR operations
400s: incompatible sequences and CIGAR operations
500s: failure to find ambiguous insertions and deletions
*/
static const char *error_text(int error_code)
{
    switch (error_code)
    {
    case 0:
        return "Success";
    default:
        return "Error";
    }
}



// Base encodings

static const char BASEA = 'A';
static const char BASEC = 'C';
static const char BASEG = 'G';
static const char BASET = 'T';


// Mutation vector byte encodings

static const char MATCH = '\x01';
static const char DELET = '\x02';
static const char INS_5 = '\x04';
static const char INS_3 = '\x08';
static const char SUB_A = '\x10';
static const char SUB_C = '\x20';
static const char SUB_G = '\x40';
static const char SUB_T = '\x80';
static const char SUB_N = SUB_A | SUB_C | SUB_G | SUB_T;
static const char ANY_N = SUB_N | MATCH;


// Numeric constants

static const int NUMERIC_BASE = 10;  // SAM files use only base 10
static const int INT_BUFFER_SIZE = 21;  // up to 20 digits for 2^64 - 1


// SAM file parsing

static const char *SAM_SEP = "\t\n";
static const char CIG_ALIGN = 'M';
static const char CIG_DELET = 'D';
static const char CIG_INSRT = 'I';
static const char CIG_MATCH = '=';
static const char CIG_SCLIP = 'S';
static const char CIG_SUBST = 'X';
static const uint16_t FLAG_PAIRED = 1;
static const uint16_t FLAG_REV = 16;
static const uint16_t FLAG_1ST = 64;
static const uint16_t FLAG_2ND = 128;


typedef struct CigarOp
{
    char *op;      // type of operation
    uint32_t len;  // length of operation
} CigarOp;


/*
Parse a string (str) and store the value in an integer (number).
A basic wrapper around string-to-unsigned-long (strtoul) that also
determines whether a return value of 0 is correct (e.g. if str == "0")
or if it is because parsing failed (strtoul returns 0 upon failure).

Parameters
----------
ulong
    Non-nullable pointer to number in which to store the result
str
    Nullable pointer to string from which to parse the number

Returns
-------
int
    0 if successful, otherwise an error code (> 0)
*/
static int parse_ulong(unsigned long *ulong, char *str)
{
    // Return an error if strtok returned a NULL pointer.
    if (str == NULL) {return 1;}
    // Parse str and store numeric value in ulong.
    *ulong = strtoul(str, NULL, NUMERIC_BASE);
    // The parse succeeded if and only if the parsed number, formatted
    // as a string, matches the original string.
    char nstr[INT_BUFFER_SIZE];
    sprintf(nstr, "%lu", *ulong);
    // Return whether the strings match.
    return strcmp(str, nstr) != 0;
}


/*
Find the next operation in a CIGAR string. Store the kind of operation
in cigar->op and the length of the operation in cigar->len. If parsing
does not return a valid operation (which must have cigar->len > 0 and
cigar->op != NULL), then set cigar->op to NULL. This event happens when
the end of the CIGAR string is reached and does not signify an error.

Parameters
----------
op
    Non-nullable pointer to a CIGAR operation struct in which to store
    the parsed information from the CIGAR string
*/
static void get_next_cigar_op(CigarOp *cigar)
{
    // Parse as many characters of text as possible to a number, cast
    // to uint32, and store in cigar->len (length of CIGAR operation).
    // Parsing must start one character ahead of the last character to
    // be read from the CIGAR string (which is cigar->op, the type of
    // the previous operation); thus, parsing starts at cigar->op + 1.
    // The character after the last digit to be parsed is the type of
    // the new operation; thus, after cigar->len is parsed, cigar->op
    // is finally pointed at the character after the last parsed digit.
    cigar->len = (uint32_t)strtoul(cigar->op + 1, &cigar->op, NUMERIC_BASE);
    if (cigar->len == 0)
    {
        // If the parsed operation is invalid because the length is not
        // positive, then set cigar->op to NULL. This happens normally
        // upon reaching the end of the CIGAR string.
        cigar->op = NULL;
    }
}


typedef struct SamRead
{
    // Parsed attributes
    char *qname;     // query name
    uint16_t flag;   // bitwise flag
    char *rname;     // reference name
    uint32_t pos;    // mapping position
    uint8_t mapq;    // mapping quality
    char *cigar;     // CIGAR string
    char *seq;       // read sequence
    char *qual;      // read quality
    // Calculated attributes
    uint32_t len;    // read length
    uint8_t paired;  // whether read is paired
    uint8_t rev;     // whether read is reverse complemented
    uint8_t is1st;   // whether read is the 1st read
    uint8_t is2nd;   // whether read is the 2nd read
} SamRead;


/*
static void print_read(SamRead *read)
{
    printf("========================\n");
    printf("QNAM: %s\n", read->qname);
    printf("FLAG: %u\n", read->flag);
    printf("RNAM: %s\n", read->rname);
    printf("POSI: %u\n", read->pos);
    printf("MAPQ: %u\n", read->mapq);
    printf("CIGR: %s\n", read->cigar);
    printf("SEQU: %s\n", read->seq);
    printf("QUAL: %s\n", read->qual);
    printf("LENG: %u\n", read->len);
    printf("========================\n");
}
*/


/*
Parse one line from a SAM file and store the information in each field
in one of the field arguments passed to this function.

Parameters
----------
read
    Non-nullable pointer to the SAM read struct in which to store the
    parsed information.
line
    Non-nullable pointer to the text of a SAM-formatted line to parse.
    It must be tab-delimited and contain at least 11 fields (below).

Returns
-------
int
    0 if successful, otherwise an error code (> 0)
*/
static int parse_sam_line(SamRead *read, char *line)
{
    char *end;  // Point to the end of the current field.
    unsigned long temp_ulong;  // Hold parsed numbers before casting.
    // Query name
    if ((read->qname = strtok_r(line, SAM_SEP, &end)) == NULL) {return 101;}
    // Bitwise flag
    if (parse_ulong(&temp_ulong, strtok_r(NULL, SAM_SEP, &end))) {return 102;}
    read->flag = (uint16_t)temp_ulong;
    // Compute individual flag bits.
    read->paired = (read->flag & FLAG_PAIRED) > 0;
    read->rev = (read->flag & FLAG_REV) > 0;
    read->is1st = (read->flag & FLAG_1ST) > 0;
    read->is2nd = (read->flag & FLAG_2ND) > 0;
    // Reference name
    if ((read->rname = strtok_r(NULL, SAM_SEP, &end)) == NULL) {return 103;}
    // Mapping position
    if (parse_ulong(&temp_ulong, strtok_r(NULL, SAM_SEP, &end))) {return 104;}
    read->pos = (uint32_t)temp_ulong;
    // Mapping quality
    if (parse_ulong(&temp_ulong, strtok_r(NULL, SAM_SEP, &end))) {return 105;}
    read->mapq = (uint8_t)temp_ulong;
    // CIGAR string
    if ((read->cigar = strtok_r(NULL, SAM_SEP, &end)) == NULL) {return 106;}
    // Next reference (ignored)
    if (strtok_r(NULL, SAM_SEP, &end) == NULL) {return 107;}
    // Next position (ignored)
    if (strtok_r(NULL, SAM_SEP, &end) == NULL) {return 108;}
    // Template length (ignored)
    if (strtok_r(NULL, SAM_SEP, &end) == NULL) {return 109;}
    // Read sequence
    if ((read->seq = strtok_r(NULL, SAM_SEP, &end)) == NULL) {return 110;}
    // Read length (computed with pointer arithmetic)
    read->len = end - read->seq - 1;
    // Read quality
    if ((read->qual = strtok_r(NULL, SAM_SEP, &end)) == NULL) {return 111;}
    // Lengths of read and quality strings must match.
    if (read->len != end - read->qual - 1) {return 112;}
    // Parsing the line succeeded.
    return 0;
}


/* Encode a base character as a substitution. */
static char encode_subs(char base)
{
    switch (base)
    {
    case BASET:
        return SUB_T;
    case BASEG:
        return SUB_G;
    case BASEC:
        return SUB_C;
    case BASEA:
        return SUB_A;
    default:
        return SUB_N;
    }
}


/*
Encode a match as a match byte (if the read quality is sufficient),
otherwise as an ambiguous byte (otherwise).
*/
static char encode_match(char read_base,
                         char read_qual,
                         char min_qual)
{
    return read_qual >= min_qual ? MATCH : ANY_N ^ encode_subs(read_base);
}


static char encode_comp(char ref_base,
                        char read_base,
                        char read_qual,
                        char min_qual)
{
    if (read_qual < min_qual) {return ANY_N ^ encode_subs(ref_base);}
    return ref_base == read_base ? MATCH : encode_subs(read_base);
}


// Insertions and deletions

typedef struct Indel
{
    char insert;
} Indel;


/*
Truncate any part of the current operation that lies before (5' of) the
beginning (5' end) of the section of interest.
*/
static void truncate_before_start(uint32_t *op_len,
                                  char **op_start_read,
                                  char **op_start_qual,
                                  char **op_start_muts,
                                  char *muts)
{
    // Check if any portion of the operation lies before (5' of) the
    // beginning (5' end) of the section of interest.
    if (*op_start_muts < muts)
    {
        // Compute the length to truncate.
        uint32_t trunc = muts - *op_start_muts;
        // Decrease the remaining length of the operation.
        *op_len -= trunc;
        // Advance the start (5' end) of the operation in the read.
        *op_start_read += trunc;
        *op_start_qual += trunc;
        // Advance the start (5' end) of the operation in the section.
        *op_start_muts = muts;
    }
}


/*
Compute the mutation vector of a SamRead.

Parameters
----------
muts
    Non-nullable pointer to empty buffer to which to write the mutation
    calls. The length of muts must equal sect_len, otherwise the
    behavior is undefined.
sect_seq
    Non-nullable pointer to the sequence of the reference within the
    region of interest. The sequence may contain only the characters
    'A', 'C', 'G', and 'T' (lowercase not allowed), and its length must
    equal sect_len, otherwise the behavior is undefined.
sect_len
    Length of the section. Must equal lengths of muts and sect_seq,
    otherwise the behavior is undefined and may cause memory violations.
sect_end5
    Position of the 5' end of the section with respect to the beginning
    of the entire reference sequence (1-indexed). Must be positive.
read
    Read from a SAM file.
min_qual
    Minimum ASCII-encoded quality score to accept a base call.
ambid
    Whether to compute and label ambiguous insertions and deletions.

Returns
-------
On success: 0
On failure: >0
*/
static int vectorize_read(SamRead *read,
                          char *muts,
                          char *sect_seq,
                          uint32_t sect_len,
                          uint32_t sect_end5,
                          char min_qual,
                          uint8_t ambid)
{
    // Point to the position in the read sequence immediately after
    // (3' of) the end of the current operation.
    char *op_end_read = read->seq;
    // Point to the position in the read quality immediately after
    // (3' of) the end of the current operation.
    char *op_end_qual = read->qual;
    // Count how many positions of the read have been consumed.
    uint32_t cons_read = 0;
    // Point to the position in the mutation vector at the beginning
    // (5' side of) the current operation. It will lie outside the
    // bounds of the mutation vector if the current operation does not
    // start within the section of interest, so dereference only after
    // verifying that the position is in bounds.
    char *muts_pos = muts + (read->pos - sect_end5);
    // Point to the position in the mutation vector immediately after
    // (3' of) the end of the current operation. It will lie outside the
    // bounds of the mutation vector if the current operation does not
    // start within the section of interest, so dereference only after
    // verifying that the position is in bounds.
    char *op_end_muts = muts_pos;
    // Point to the position immediately after the end (3' side) of the
    // mutation vector.
    char *muts_end = muts + sect_len;
    // Point to recorded insertions and deletions.
    Indel *indels = NULL;
    // Initialize the CIGAR operation.
    CigarOp cigar;
    cigar.op = read->cigar - 1;  // point just before the CIGAR string
    int cigar_count = 0;  // count the number of operations parsed
    // Read the first operation from the CIGAR string.
    get_next_cigar_op(&cigar);
    // Read the entire CIGAR string one operation at a time.
    while (cigar.op != NULL)
    {
        cigar_count++;
        switch (cigar.op[0])
        {
        // Base(s) in the operation match the reference.
        case CIG_MATCH:
            // Check for an overshoot of the read sequence.
            if (cigar.len > read->len) {return 200;}
            // Advance the end (3' side) of the operation.
            op_end_muts += cigar.len;
            op_end_read += cigar.len;
            op_end_qual += cigar.len;
            cons_read += cigar.len;
            // Check if the operation overlaps the mutation vector.
            if (op_end_muts > muts && muts_pos < muts_end)
            {
                // Truncate any part of the operation before (5' of) the
                // mutation vector.
                truncate_before_start(&cigar.len, &read->seq, &read->qual,
                                      &muts_pos, muts);
                // Compute the mutational status of each position until
                // the 3' end of the operation or the 3' end of the
                // section, whichever comes first.
                char *op_limit_muts = (op_end_muts < muts_end ?
                                       op_end_muts : muts_end);
                while (muts_pos < op_limit_muts)
                {
                    // Mark the mutation vector as a match or ambiguous.
                    *muts_pos |= encode_match(read->seq[0],
                                              read->qual[0],
                                              min_qual);
                    // Advance the position in the vector and read.
                    muts_pos++;
                    read->seq++;
                    read->qual++;
                }
            }
            break;
        // Base(s) in the operation match or have a substitution.
        case CIG_ALIGN:
        case CIG_SUBST:
            // Check for an overshoot of the read sequence.
            if (cigar.len > read->len) {return 200;}
            // Advance the end (3' side) of the operation.
            op_end_muts += cigar.len;
            op_end_read += cigar.len;
            op_end_qual += cigar.len;
            cons_read += cigar.len;
            // Check if the operation overlaps the mutation vector.
            if (op_end_muts > muts && muts_pos < muts_end)
            {
                // Truncate any part of the operation before (5' of) the
                // mutation vector.
                truncate_before_start(&cigar.len, &read->seq, &read->qual,
                                      &muts_pos, muts);
                // Compute the mutational status of each position until
                // the 3' end of the operation or the 3' end of the
                // section, whichever comes first.
                char *op_limit_muts = (op_end_muts < muts_end ?
                                       op_end_muts : muts_end);
                while (muts_pos < op_limit_muts)
                {
                    // Mark the mutation vector as a match, substitution
                    // or ambiguous.
                    *muts_pos |= encode_comp(sect_seq[muts_pos - muts],
                                             read->seq[0],
                                             read->qual[0],
                                             min_qual);
                    // Advance the position in the vector and read.
                    muts_pos++;
                    read->seq++;
                    read->qual++;
                }
            }
            break;
        // Base(s) in the operation are deleted from the read.
        case CIG_DELET:
            // Advance the end (3' side) of the operation.
            op_end_muts += cigar.len;
            // Check if the operation overlaps the mutation vector.
            if (op_end_muts > muts && muts_pos < muts_end)
            {
                // Truncate any part of the operation before (5' of) the
                // mutation vector.
                truncate_before_start(&cigar.len, &read->seq, &read->qual,
                                      &muts_pos, muts);
                // Compute the mutational status of each position until
                // the 3' end of the operation or the 3' end of the
                // section, whichever comes first.
                char *op_limit_muts = (op_end_muts < muts_end ?
                                       op_end_muts : muts_end);
                while (muts_pos < op_limit_muts)
                {
                    // TODO: add deletion to list of indels
                    // Mark a deletion in the mutation vector.
                    *muts_pos |= DELET;
                    // Advance the position in the vector.
                    muts_pos++;
                }
            }
            break;
        case CIG_INSRT:
            // Check for an overshoot of the read sequence.
            if (cigar.len > read->len) {return 200;}
            // Advance the end (3' side) of the operation.
            op_end_read += cigar.len;
            op_end_qual += cigar.len;
            cons_read += cigar.len;
            // Check if the operation overlaps the mutation vector.
            // Note that for insertions only, test muts_pos <= muts_end
            // instead of muts_pos < muts_end because an insertion right
            // after the section is possible (albeit very unlikely).
            if (op_end_muts > muts && muts_pos <= muts_end)
            {
                // Mark the vector positions 5' and 3' of the insertion.
                if (muts_pos > muts) {*(muts_pos - 1) |= INS_5;}
                if (muts_pos < muts_end) {*muts_pos |= INS_3;}
                // Compute the mutational status of each position until
                // the 3' end of the operation.
                while (read->seq < op_end_read)
                {
                    // TODO: add insertion to list of indels
                    // Advance the position in the read.
                    read->seq++;
                    read->qual++;
                }
            }
            break;
        case CIG_SCLIP:
            // Check for an overshoot of the read sequence.
            if (cigar.len > read->len) {return 200;}
            // Skip over this operation.
            op_end_read += cigar.len;
            op_end_qual += cigar.len;
            read->seq = op_end_read;
            read->qual = op_end_qual;
            break;
        // The CIGAR operation was invalid.
        default:
            return 202;
        }
        // Read the next operation from the CIGAR string.
        get_next_cigar_op(&cigar);
    }
    // Return an error if no CIGAR operations were read.
    if (!cigar_count) {return 203;}
    // Return an error if the entire read has not been consumed.
    if (cons_read < read->len) {return 204;}
    if (cons_read > read->len) {return 200;}
    // Deallocate all recorded insertions and deletions (if any).
    free(indels);
    // The read was parsed to a mutation vector successfully.
    return 0;
}


static int vectorize_line(char *line,
                          char *muts,
                          char *ref_name,
                          char *sect_seq,
                          uint32_t sect_len,
                          uint32_t sect_end5,
                          char min_qual,
                          uint8_t ambid)
{
    // Declare a SamRead (not a pointer to a SamRead).
    SamRead read;
    // Parse the line and store its information in the read.
    int error = parse_sam_line(&read, line);
    // Return if an error occurred during parsing.
    if (error) {return error;}
    // Compute the mutation vector of the read and store it in muts.
    return vectorize_read(&read, muts, sect_seq, sect_len, sect_end5,
                          min_qual, ambid);
}


static int vectorize_pair(char *line1,
                          char *line2,
                          char *muts,
                          char *ref_name,
                          char *sect_seq,
                          uint32_t sect_len,
                          uint32_t sect_end5,
                          char min_qual,
                          uint8_t ambid)
{
    // Declare a pair of SamReads (not pointers to SamReads).
    SamRead read1, read2;
    // Parse line 1 and store its information in read 1.
    int error = parse_sam_line(&read1, line1);
    // Return if an error occurred during parsing.
    if (error) {return error;}
    // Parse line 2 and store its information in read 2.
    error = parse_sam_line(&read2, line2);
    // Return if an error occurred during parsing.
    if (error) {return error;}
    // Verify that the reads form a compatible pair.
    if (!((&read1)->paired && (&read2)->paired)) {return 501;}
    if (!((&read1)->is1st && (&read2)->is2nd)) {return 501;}
    if ((&read1)->rev == (&read2)->rev) {return 501;}
    if (strcmp((&read1)->qname, (&read2)->qname)) {return 502;}
    if (strcmp((&read1)->rname, (&read2)->rname)) {return 503;}
    if (strcmp((&read1)->rname, ref_name)) {return 504;}
    // Vectorize read 1; if it fails, return an error code. Otherwise,
    // vectorize read 2 and return 0 for success and >0 for failure.
    return (vectorize_read(&read1, muts, sect_seq, sect_len, sect_end5,
                           min_qual, ambid) ||
            vectorize_read(&read2, muts, sect_seq, sect_len, sect_end5,
                           min_qual, ambid));
}


/* Python interface to function vectorize_line */
static PyObject *py_vecline(PyObject *self, PyObject *args)
{
    // Arguments from Python function call
    Py_buffer line1;
    Py_buffer muts;
    Py_buffer sect_seq;
    uint32_t sect_len;
    uint32_t sect_end5;
    char *ref_name;
    char min_qual;
    uint8_t ambid;

    // Attempt to convert the arguments from the Python function call
    // into C data types. Return NULL upon failure.
    if (!PyArg_ParseTuple(args, "y*y*y*iisbp", &line1, &muts,
                          &sect_seq, &sect_len, &sect_end5, &ref_name,
                          &min_qual, &ambid)) {return NULL;}
    // Attempt to vectorize the line.
    int error = vectorize_line(line1.buf, muts.buf, ref_name,
                               sect_seq.buf, sect_len, sect_end5,
                               min_qual, ambid);
    if (error)
    {
        PyErr_SetString(VectorError, error_text(error));
        return NULL;
    }
    // Vectoring completed successfully.
    return PyLong_FromLong((long)42);
}


/* Python interface to function vectorize_line */
static PyObject *py_vecpair(PyObject *self, PyObject *args)
{
    // Arguments from Python function call
    Py_buffer line1;
    Py_buffer line2;
    Py_buffer muts;
    Py_buffer sect_seq;
    uint32_t sect_len;
    uint32_t sect_end5;
    char *ref_name;
    char min_qual;
    uint8_t ambid;

    // Attempt to convert the arguments from the Python function call
    // into C data types. Return NULL upon failure.
    if (!PyArg_ParseTuple(args, "y*y*y*y*iisbp", &line1, &line2, &muts,
                          &sect_seq, &sect_len, &sect_end5, &ref_name,
                          &min_qual, &ambid)) {return NULL;}
    // Attempt to vectorize the line.
    int error = vectorize_pair(line1.buf, line2.buf, muts.buf, ref_name,
                               sect_seq.buf, sect_len, sect_end5,
                               min_qual, ambid);
    if (error)
    {
        PyErr_SetString(VectorError, error_text(error));
        return NULL;
    }
    // Vectoring completed successfully.
    return PyLong_FromLong((long)42);
}


/* Python method table */
static PyMethodDef VectorMethods[] = {
    {"vectorize_line", py_vecline, METH_VARARGS,
     "Generate a mutation vector from a SAM line."},
    {"vectorize_pair", py_vecpair, METH_VARARGS,
     "Generate a mutation vector from a pair of SAM lines."},
    {NULL, NULL, 0, NULL}  // sentinel
};


/* Python module definition */
static struct PyModuleDef vectormodule = {
    PyModuleDef_HEAD_INIT,
    "vector",      // module name
    NULL,          // documentation (TODO)
    -1,            // module state
    VectorMethods  // method table
};


/* Python module initialization function */
PyMODINIT_FUNC PyInit_vector(void)
{
    // Initialize the module and ensure it exists.
    PyObject *module;
    module = PyModule_Create(&vectormodule);
    if (module == NULL) {return NULL;}

    // Define a new type of Python exception for vectoring.
    VectorError = PyErr_NewException("vector.error", NULL, NULL);
    // Add the exception type to the module.
    Py_XINCREF(VectorError);
    if (PyModule_AddObject(module, "error", VectorError) < 0)
    {
        // Adding the exception type failed. Stop and return NULL.
        Py_XDECREF(VectorError);
        Py_CLEAR(VectorError);
        Py_DECREF(module);
        return NULL;
    }

    // Initializing the module succeeded.
    return module;
}


int main()
{
    // This module must be called via its Python API, not run directly.
    return 0;
}
