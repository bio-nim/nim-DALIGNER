# vim: sw=2 ts=2 sts=2 tw=80 et:
{.compile: "QV.c".}
## ******************************************************************************************
## 
##   Compressor/decompressor for .quiv files: customized Huffman codes for each stream based on
##     the histogram of values occuring in a given file.  The two low complexity streams
##     (deletionQV and substitutionQV) use a Huffman coding of the run length of the prevelant
##     character.
## 
##   Author:   Gene Myers
##   Date:     Jan 18, 2014
##   Modified: July 25, 2014
## 
## ******************************************************************************************

##   The defined constant INTERACTIVE (set in DB.h) determines whether an interactive or
##     batch version of the routines in this library are compiled.  In batch mode, routines
##     print an error message and exit.  In interactive mode, the routines place the error
##     message in EPLACE (also defined in DB.h) and return an error value, typically NULL
##     if the routine returns a pointer, and an unusual integer value if the routine returns
##     an integer.
##   Below when an error return is described, one should understand that this value is returned
##     only if the routine was compiled in INTERACTIVE mode.
##   A PacBio compression scheme

type
  QVcoding* {.importc: "QVcoding", header: "QV.h".} = object
    delScheme* {.importc: "delScheme".}: pointer ##   Huffman scheme for deletion QVs
    insScheme* {.importc: "insScheme".}: pointer ##   Huffman scheme for insertion QVs
    mrgScheme* {.importc: "mrgScheme".}: pointer ##   Huffman scheme for merge QVs
    subScheme* {.importc: "subScheme".}: pointer ##   Huffman scheme for substitution QVs
    dRunScheme* {.importc: "dRunScheme".}: pointer ##   Huffman scheme for deletion run lengths (if delChar > 0)
    sRunScheme* {.importc: "sRunScheme".}: pointer ##   Huffman scheme for substitution run lengths (if subChar > 0)
    delChar* {.importc: "delChar".}: cint ##   If > 0, run-encoded deletion value
    subChar* {.importc: "subChar".}: cint ##   If > 0, run-encoded substitution value
    flip* {.importc: "flip".}: cint ##   Need to flip multi-byte integers
    prefix* {.importc: "prefix".}: cstring ##   Header line prefix
  

##  Read the next nlines of input, and QVentry returns a pointer to the first line if needed.
##    If end-of-input is encountered before any further input, -1 is returned.  If there is
##    an error than -2 is returned.  Otherwise the length of the line(s) read is returned.

proc Read_Lines*(input: ptr FILE; nlines: cint): cint {.cdecl, importc: "Read_Lines",
    header: "QV.h".}
proc QVentry*(): cstring {.cdecl, importc: "QVentry", header: "QV.h".}
##  Get and set the line counter for error reporting

proc Set_QV_Line*(line: cint) {.cdecl, importc: "Set_QV_Line", header: "QV.h".}
proc Get_QV_Line*(): cint {.cdecl, importc: "Get_QV_Line", header: "QV.h".}
##  Read up to the next num entries or until eof from the .quiva file on input and record
##    frequency statistics.  Copy these entries to the temporary file temp if != NULL.
##    If there is an error then -1 is returned, otherwise the number of entries read.

proc QVcoding_Scan*(input: ptr FILE; num: cint; temp: ptr FILE): cint {.cdecl,
    importc: "QVcoding_Scan", header: "QV.h".}
proc QVcoding_Scan1*(rlen: cint; del: cstring; tag: cstring; ins: cstring; mrg: cstring;
                    sub: cstring) {.cdecl, importc: "QVcoding_Scan1", header: "QV.h".}
##  Given QVcoding_Scan has been called at least once, create an encoding scheme based on
##    the accumulated statistics and return a pointer to it.  The returned encoding object
##    is *statically allocated within the routine.  If lossy is set then use a lossy scaling
##    for the insertion and merge streams.  If there is an error, then NULL is returned.

proc Create_QVcoding*(lossy: cint): ptr QVcoding {.cdecl, importc: "Create_QVcoding",
    header: "QV.h".}
##   Read/write a coding scheme to input/output.  The encoding object returned by the reader
##     is *statically* allocated within the routine.  If an error occurs while reading then
##     NULL is returned.

proc Read_QVcoding*(input: ptr FILE): ptr QVcoding {.cdecl, importc: "Read_QVcoding",
    header: "QV.h".}
proc Write_QVcoding*(output: ptr FILE; coding: ptr QVcoding) {.cdecl,
    importc: "Write_QVcoding", header: "QV.h".}
##   Free all the auxiliary storage associated with coding (but not the object itself!)

proc Free_QVcoding*(coding: ptr QVcoding) {.cdecl, importc: "Free_QVcoding",
                                        header: "QV.h".}
##   Assuming the file pointer is positioned just beyond an entry header line, read the
##     next set of 5 QV lines, compress them according to 'coding', and output.  If lossy
##     is set then the scheme is a lossy one.  A negative value is returned if an error
##     occurred, and the sequence length otherwise.

proc Compress_Next_QVentry*(input: ptr FILE; output: ptr FILE; coding: ptr QVcoding;
                           lossy: cint): cint {.cdecl,
    importc: "Compress_Next_QVentry", header: "QV.h".}
proc Compress_Next_QVentry1*(rlen: cint; del: cstring; tag: cstring; ins: cstring;
                            mrg: cstring; sub: cstring; output: ptr FILE;
                            coding: ptr QVcoding; lossy: cint) {.cdecl,
    importc: "Compress_Next_QVentry1", header: "QV.h".}
##   Assuming the input is position just beyond the compressed encoding of an entry header,
##     read the set of compressed encodings for the ensuing 5 QV vectors, decompress them,
##     and place their decompressed values into entry which is a 5 element array of character
##     pointers.  The parameter rlen computed from the preceeding header line, critically
##     provides the length of each of the 5 vectors.  A non-zero value is return only if an
##     error occured.

proc Uncompress_Next_QVentry*(input: ptr FILE; entry: cstringArray;
                             coding: ptr QVcoding; rlen: cint): cint {.cdecl,
    importc: "Uncompress_Next_QVentry", header: "QV.h".}
