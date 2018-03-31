# vim: sw=2 ts=2 sts=2 tw=80 et:
{.compile: "DB.c".}
## ******************************************************************************************
## 
##   Compressed data base module.  Auxiliary routines to open and manipulate a data base for
##     which the sequence and read information are separated into two separate files, and the
##     sequence is compressed into 2-bits for each base.  Support for tracks of additional
##     information, and trimming according to the current partition.  Eventually will also
##     support compressed quality information. 
## 
##   Author :  Gene Myers
##   Date   :  July 2013
##   Revised:  April 2014
## 
## ******************************************************************************************

from qv import QVcoding
## ******************************************************************************************
## 
##   UTILITIES
## 
## ******************************************************************************************
##   The following general utilities return NULL if any of their input pointers are NULL, or if they
##     could not perform their function (in which case they also print an error to stderr).

proc Malloc*(size: int64; mesg: cstring): pointer {.cdecl, importc: "Malloc",
    header: "DB.h".}
##   Guarded versions of malloc, realloc

proc Realloc*(`object`: pointer; size: int64; mesg: cstring): pointer {.cdecl,
    importc: "Realloc", header: "DB.h".}
##   and strdup, that output "mesg" to

proc Strdup*(string: cstring; mesg: cstring): cstring {.cdecl, importc: "Strdup",
    header: "DB.h".}
##   stderr if out of memory

proc Fopen*(path: cstring; mode: cstring): ptr FILE {.cdecl, importc: "Fopen",
    header: "DB.h".}
##  Open file path for "mode"

proc pPathTo(path: cstring): cstring {.cdecl, importc: "PathTo", header: "DB.h".}
##  Return path portion of file name "path"
proc PathTo*(path: string): string =
  return $pPathTo(path.cstring)

proc pRoot(path: cstring; suffix: cstring): cstring {.cdecl, importc: "Root", header: "DB.h".}
##  Return the root name, excluding suffix, of "path"
##  Catenate returns concatenation of path.sep.root.suffix in a *temporary* buffer
##  Numbered_Suffix returns concatenation of left.<num>.right in a *temporary* buffer
proc Root*(path: string; suffix: string): string =
  return $pRoot(path.cstring, suffix.cstring)

proc pCatenate(path: cstring; sep: cstring; root: cstring; suffix: cstring): cstring {.  cdecl, importc: "Catenate", header: "DB.h".}
proc Catenate*(path: string; sep: string; root: string; suffix: string): string =
  return $pCatenate(path.cstring, sep.cstring, root.cstring, suffix.cstring)

proc Numbered_Suffix*(left: cstring; num: cint; right: cstring): cstring {.cdecl,
    importc: "Numbered_Suffix", header: "DB.h".}
##  DB-related utilities

proc Print_Number*(num: int64; width: cint; `out`: FILE) {.cdecl,
    importc: "Print_Number", header: "DB.h".}
##   Print readable big integer

proc Number_Digits*(num: int64): cint {.cdecl, importc: "Number_Digits", header: "DB.h".}
##   Return # of digits in printed number

template COMPRESSED_LEN*(len: untyped): untyped =
  (((len) + 3) shr 2)

proc Compress_Read*(len: cint; s: cstring) {.cdecl, importc: "Compress_Read",
                                        header: "DB.h".}
##   Compress read in-place into 2-bit form

proc Uncompress_Read*(len: cint; s: cstring) {.cdecl, importc: "Uncompress_Read",
    header: "DB.h".}
##   Uncompress read in-place into numeric form

proc Print_Read*(s: cstring; width: cint) {.cdecl, importc: "Print_Read", header: "DB.h".}
proc Lower_Read*(s: cstring) {.cdecl, importc: "Lower_Read", header: "DB.h".}
##   Convert read from numbers to lowercase letters (0-3 to acgt)

proc Upper_Read*(s: cstring) {.cdecl, importc: "Upper_Read", header: "DB.h".}
##   Convert read from numbers to uppercase letters (0-3 to ACGT)

proc Number_Read*(s: cstring) {.cdecl, importc: "Number_Read", header: "DB.h".}
##   Convert read from letters to numbers

proc Letter_Arrow*(s: cstring) {.cdecl, importc: "Letter_Arrow", header: "DB.h".}
##   Convert arrow pw's from numbers to uppercase letters (0-3 to 1234)

proc Number_Arrow*(s: cstring) {.cdecl, importc: "Number_Arrow", header: "DB.h".}
##   Convert arrow pw string from letters to numbers
## 
## ******************************************************************************************
## 
##   DB IN-CORE DATA STRUCTURES
## 
## 
## ******************************************************************************************

const
  DB_QV* = 0x000003FF
  DB_CSS* = 0x00000400
  DB_BEST* = 0x00000800
  DB_ARROW* = 0x00000002
  DB_ALL* = 0x00000001

##   Fields have different interpretations if a .db versus a .dam

type
  DAZZ_READ* {.importc: "DAZZ_READ", header: "DB.h".} = object
    origin* {.importc: "origin".}: cint ##   Well # (DB), Contig # (DAM)
    rlen* {.importc: "rlen".}: cint ##   Length of the sequence (Last pulse = fpulse + rlen)
    fpulse* {.importc: "fpulse".}: cint ##   First pulse (DB), left index of contig in scaffold (DAM)
    boff* {.importc: "boff".}: int64 ##   Offset (in bytes) of compressed read in 'bases' file, or offset of
                                 ##     uncompressed bases in memory block
    coff* {.importc: "coff".}: int64 ##   Offset (in bytes) of compressed quiva streams in '.qvs' file (DB),
                                 ##   Offset (in bytes) of scaffold header string in '.hdr' file (DAM)
                                 ##   4 compressed shorts containing snr info if an arrow DB.
    flags* {.importc: "flags".}: cint ##   QV of read + flags above (DB only)
  

##   A track can be of 3 types:
##     data == NULL: there are nreads 'anno' records of size 'size'.
##     data != NULL && size == 4: anno is an array of nreads+1 int's and data[anno[i]..anno[i+1])
##                                     contains the variable length data
##     data != NULL && size == 8: anno is an array of nreads+1 int64's and data[anno[i]..anno[i+1])
##                                     contains the variable length data

type
  DAZZ_TRACK* {.importc: "DAZZ_TRACK", header: "DB.h".} = object
    next* {.importc: "next".}: ptr DAZZ_TRACK ##   Link to next track
    name* {.importc: "name".}: cstring ##   Symbolic name of track
    size* {.importc: "size".}: cint ##   Size in bytes of anno records
    anno* {.importc: "anno".}: pointer ##   over [0,nreads]: read i annotation: int, int64, or 'size' records
    data* {.importc: "data".}: pointer ##      data[anno[i] .. anno[i+1]-1] is data if data != NULL
  

##   The information for accessing QV streams is in a DAZZ_QV record that is a "pseudo-track"
##     named ".@qvs" and is always the first track record in the list (if present).  Since normal
##     track names cannot begin with a . (this is enforced), this pseudo-track is never confused
##     with a normal track.

type
  DAZZ_QV* {.importc: "DAZZ_QV", header: "DB.h".} = object
    next* {.importc: "next".}: ptr DAZZ_TRACK
    name* {.importc: "name".}: cstring
    ncodes* {.importc: "ncodes".}: cint ##   # of coding tables
    coding* {.importc: "coding".}: ptr QVcoding ##   array [0..ncodes-1] of coding schemes (see QV.h)
    table* {.importc: "table".}: ptr uint16 ##   for i in [0,db->nreads-1]: read i should be decompressed with
                                       ##     scheme coding[table[i]]
    quiva* {.importc: "quiva".}: ptr FILE ##   the open file pointer to the .qvs file
  

##   The DB record holds all information about the current state of an active DB including an
##     array of DAZZ_READS, one per read, and a linked list of DAZZ_TRACKs the first of which
##     is always a DAZZ_QV pseudo-track (if the QVs have been loaded).

type
  DAZZ_DB* {.importc: "DAZZ_DB", header: "DB.h".} = object
    ureads* {.importc: "ureads".}: cint ##   Total number of reads in untrimmed DB
    treads* {.importc: "treads".}: cint ##   Total number of reads in trimmed DB
    cutoff* {.importc: "cutoff".}: cint ##   Minimum read length in block (-1 if not yet set)
    allarr* {.importc: "allarr".}: cint ##   DB_ALL | DB_ARROW
    freq* {.importc: "freq".}: array[4, cfloat] ##   frequency of A, C, G, T, respectively
                                           ##   Set with respect to "active" part of DB (all vs block, untrimmed vs trimmed)
    maxlen* {.importc: "maxlen".}: cint ##   length of maximum read (initially over all DB)
    totlen* {.importc: "totlen".}: int64 ##   total # of bases (initially over all DB)
    nreads* {.importc: "nreads".}: cint ##   # of reads in actively loaded portion of DB
    trimmed* {.importc: "trimmed".}: cint ##   DB has been trimmed by cutoff/all
    part* {.importc: "part".}: cint ##   DB block (if > 0), total DB (if == 0)
    ufirst* {.importc: "ufirst".}: cint ##   Index of first read in block (without trimming)
    tfirst* {.importc: "tfirst".}: cint ##   Index of first read in block (with trimming)
                                    ##   In order to avoid forcing users to have to rebuild all thier DBs to accommodate
                                    ##     the addition of fields for the size of the actively loaded trimmed and untrimmed
                                    ##     blocks, an additional read record is allocated in "reads" when a DB is loaded into
                                    ##     memory (reads[-1]) and the two desired fields are crammed into the first two
                                    ##     integer spaces of the record.
    path* {.importc: "path".}: cstring ##   Root name of DB for .bps, .qvs, and tracks
    loaded* {.importc: "loaded".}: cint ##   Are reads loaded in memory?
    bases* {.importc: "bases".}: pointer ##   file pointer for bases file (to fetch reads from),
                                     ##     or memory pointer to uncompressed block of all sequences.
    reads* {.importc: "reads".}: ptr DAZZ_READ ##   Array [-1..nreads] of DAZZ_READ
    tracks* {.importc: "tracks".}: ptr DAZZ_TRACK ##   Linked list of loaded tracks
  

## ******************************************************************************************
## 
##   DB STUB FILE FORMAT = NFILE FDATA^nfile NBLOCK PARAMS BDATA^nblock
## 
## ******************************************************************************************

const
  MAX_NAME* = 10000
  DB_NFILE* = "files = %9d\x0A"
  DB_FDATA* = "  %9d %s %s\x0A"
  DB_NBLOCK* = "blocks = %9d\x0A"
  DB_PARAMS* = "size = %10lld cutoff = %9d all = %1d\x0A"
  DB_BDATA* = " %9d %9d\x0A"

## ******************************************************************************************
## 
##   DB ROUTINES
## 
## ******************************************************************************************
##  Suppose DB is the name of an original database.  Then there will be files .DB.idx, .DB.bps,
##     .DB.qvs, and files .DB.<track>.anno and DB.<track>.data where <track> is a track name
##     (not containing a . !).
##  A DAM is basically a DB except that:
##     1. there are no QV's, instead .coff points the '\0' terminated fasta header of the read
##           in the file .<dam>.hdr file
##     2. .origin contains the contig # of the read within a fasta entry (assembly sequences
##           contain N-separated contigs), and .fpulse the first base of the contig in the
##           fasta entry
##  Open the given database or dam, "path" into the supplied DAZZ_DB record "db". If the name has
##    a part # in it then just the part is opened.  The index array is allocated (for all or
##    just the part) and read in.
##  Return status of routine:
##     -1: The DB could not be opened for a reason reported by the routine to EPLACE
##      0: Open of DB proceeded without mishap
##      1: Open of DAM proceeded without mishap
# Note: path must be a writable string!
proc Open_DB*(path: cstring; db: ptr DAZZ_DB): cint {.cdecl, importc: "Open_DB",
    header: "DB.h".}
##  Trim the DB or part thereof and all loaded tracks according to the cutoff and all settings
##    of the current DB partition.  Reallocate smaller memory blocks for the information kept
##    for the retained reads.

proc Trim_DB*(db: ptr DAZZ_DB) {.cdecl, importc: "Trim_DB", header: "DB.h".}
##  Shut down an open 'db' by freeing all associated space, including tracks and QV structures,
##    and any open file pointers.  The record pointed at by db however remains (the user
##    supplied it and so should free it).

proc Close_DB*(db: ptr DAZZ_DB) {.cdecl, importc: "Close_DB", header: "DB.h".}
##  Return the size in bytes of the given DB

proc sizeof_DB*(db: ptr DAZZ_DB): int64 {.cdecl, importc: "sizeof_DB", header: "DB.h".}
##  If QV pseudo track is not already in db's track list, then load it and set it up.
##    The database must not have been trimmed yet.  -1 is returned if a .qvs file is not
##    present, and 1 is returned if an error (reported to EPLACE) occured and INTERACTIVE
##    is defined.  Otherwise a 0 is returned.

proc Load_QVs*(db: ptr DAZZ_DB): cint {.cdecl, importc: "Load_QVs", header: "DB.h".}
##  Remove the QV pseudo track, all space associated with it, and close the .qvs file.

proc Close_QVs*(db: ptr DAZZ_DB) {.cdecl, importc: "Close_QVs", header: "DB.h".}
##  Look up the file and header in the file of the indicated track.  Return:
##      1: Track is for trimmed DB
##      0: Track is for untrimmed DB
##     -1: Track is not the right size of DB either trimmed or untrimmed
##     -2: Could not find the track
##  In addition, if opened (0 or 1 returned), then kind points at an integer indicating
##    the type of track as follows:
##       CUSTOM  0 => a custom track
##       MASK    1 => a mask track

const
  CUSTOM_TRACK* = 0
  MASK_TRACK* = 1

proc Check_Track*(db: ptr DAZZ_DB; track: cstring; kind: ptr cint): cint {.cdecl,
    importc: "Check_Track", header: "DB.h".}
##  If track is not already in the db's track list, then allocate all the storage for it,
##    read it in from the appropriate file, add it to the track list, and return a pointer
##    to the newly created DAZZ_TRACK record.  If the track does not exist or cannot be
##    opened for some reason, then NULL is returned if INTERACTIVE is defined.  Otherwise
##    the routine prints an error message to stderr and exits if an error occurs, and returns
##    with NULL only if the track does not exist.

proc Load_Track*(db: ptr DAZZ_DB; track: cstring): ptr DAZZ_TRACK {.cdecl,
    importc: "Load_Track", header: "DB.h".}
##  If track is on the db's track list, then it is removed and all storage associated with it
##    is freed.

proc Close_Track*(db: ptr DAZZ_DB; track: cstring) {.cdecl, importc: "Close_Track",
    header: "DB.h".}
##  Allocate and return a buffer big enough for the largest read in 'db'.
##  **NB** free(x-1) if x is the value returned as *prefix* and suffix '\0'(4)-byte
##  are needed by the alignment algorithms.  If cannot allocate memory then return NULL
##  if INTERACTIVE is defined, or print error to stderr and exit otherwise.

proc New_Read_Buffer*(db: ptr DAZZ_DB): ptr char {.cdecl, importc: "New_Read_Buffer",
    header: "DB.h".}
##  Load into 'read' the i'th read in 'db'.  As a lower case ascii string if ascii is 1, an
##    upper case ascii string if ascii is 2, and a numeric string over 0(A), 1(C), 2(G), and 3(T)
##    otherwise.  A '\0' (or 4) is prepended and appended to the string so it has a delimeter
##    for traversals in either direction.  A non-zero value is returned if an error occured
##    and INTERACTIVE is defined.

proc Load_Read*(db: ptr DAZZ_DB; i: cint; read: cstring; ascii: cint): cint {.cdecl,
    importc: "Load_Read", header: "DB.h", discardable.}
##  Exactly the same as Load_Read, save the arrow information is loaded, not the DNA sequence,
##    and there is only a choice between numeric (0) or ascii (1);

proc Load_Arrow*(db: ptr DAZZ_DB; i: cint; read: cstring; ascii: cint): cint {.cdecl,
    importc: "Load_Arrow", header: "DB.h".}
##  Load into 'read' the subread [beg,end] of the i'th read in 'db' and return a pointer to the
##    the start of the subinterval (not necessarily = to read !!! ).  As a lower case ascii
##    string if ascii is 1, an upper case ascii string if ascii is 2, and a numeric string
##    over 0(A), 1(C), 2(G), and 3(T) otherwise.  A '\0' (or 4) is prepended and appended to
##    the string holding the substring so it has a delimeter for traversals in either direction.
##    A NULL pointer is returned if an error occured and INTERACTIVE is defined.

proc Load_Subread*(db: ptr DAZZ_DB; i: cint; beg: cint; `end`: cint; read: cstring;
                  ascii: cint): ptr char {.cdecl, importc: "Load_Subread",
                                       header: "DB.h".}
##  Allocate a set of 5 vectors large enough to hold the longest QV stream that will occur
##    in the database.  If cannot allocate memory then return NULL if INTERACTIVE is defined,
##    or print error to stderr and exit otherwise.

const
  DEL_QV* = 0
  DEL_TAG* = 1
  INS_QV* = 2
  SUB_QV* = 3
  MRG_QV* = 4

proc New_QV_Buffer*(db: ptr DAZZ_DB): cstringArray {.cdecl, importc: "New_QV_Buffer",
    header: "DB.h".}
##  Load into 'entry' the 5 QV vectors for i'th read in 'db'.  The deletion tag or characters
##    are converted to a numeric or upper/lower case ascii string as per ascii.  Return with
##    a zero, except when an error occurs and INTERACTIVE is defined in which case return wtih 1.

proc Load_QVentry*(db: ptr DAZZ_DB; i: cint; entry: cstringArray; ascii: cint): cint {.
    cdecl, importc: "Load_QVentry", header: "DB.h".}
##  Allocate a block big enough for all the uncompressed sequences, read them into it,
##    reset the 'off' in each read record to be its in-memory offset, and set the
##    bases pointer to point at the block after closing the bases file.  If ascii is
##    1 then the reads are converted to lowercase ascii, if 2 then uppercase ascii, and
##    otherwise the reads are left as numeric strings over 0(A), 1(C), 2(G), and 3(T).
##    Return with a zero, except when an error occurs and INTERACTIVE is defined in which
##    case return wtih 1.

proc Read_All_Sequences*(db: ptr DAZZ_DB; ascii: cint): cint {.cdecl,
    importc: "Read_All_Sequences", header: "DB.h".}
##  For the DB or DAM "path" = "prefix/root.[db|dam]", find all the files for that DB, i.e. all
##    those of the form "prefix/[.]root.part" and call actor with the complete path to each file
##    pointed at by path, and the suffix of the path by extension.  The . proceeds the root
##    name if the defined constant HIDE_FILES is set.  Always the first call is with the
##    path "prefix/root.[db|dam]" and extension "db" or "dam".  There will always be calls for
##    "prefix/[.]root.idx" and "prefix/[.]root.bps".  All other calls are for *tracks* and
##    so this routine gives one a way to know all the tracks associated with a given DB.
##    -1 is returned if the path could not be found, and 1 is returned if an error (reported
##    to EPLACE) occured and INTERACTIVE is defined.  Otherwise a 0 is returned.

proc List_DB_Files*(path: cstring;
                   actor: proc (path: cstring; extension: cstring) {.cdecl.}): cint {.
    cdecl, importc: "List_DB_Files", header: "DB.h".}
