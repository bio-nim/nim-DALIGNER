# vim: sw=2 ts=2 sts=2 tw=80 et:
{.compile: "DBX.c".}

##  Wrappers to extend HITS_DB.
##
##  Note that none of the extra fields are ever stored on-disk.

from os import nil
import
  DB

type
  HITS_DBX* {.importc: "HITS_DBX", header: "DBX.h".} = object
    db* {.importc: "db".}: HITS_DB

    ##  When "data" is non-null, it stores the entire DB
    ##  in memory, so we can avoid random-access disk operations.
    ##  But if null, then wrappers simply delegate.
    ##
    data* {.importc: "data".}: cstring


proc pOpen_DBX*(path: cstring; dbx: ptr HITS_DBX; preload: bool): cint {.cdecl,
    importc: "Open_DBX", header: "DBX.h".}

proc Open_DBX*(path: string, dbx: ptr HITS_DBX; preload: bool) =
  let ret = pOpen_DBX(path, dbx, preload)
  if ret != 0:
    os.raiseOSError(os.osLastError(), "Could not open DAZZ_DB '" & path & "'")

proc Load_ReadX*(dbx: ptr HITS_DBX; i: cint; read: cstring; ascii: cint): cint {.cdecl,
    importc: "Load_ReadX", header: "DBX.h", discardable.}
proc Close_DBX*(dbx: ptr HITS_DBX) {.cdecl, importc: "Close_DBX", header: "DBX.h".}
