# vim: sw=2 ts=2 sts=2 tw=80 et:
from os import DirSep
from strutils import rsplit
const thisdir = system.currentSourcePath.rsplit(DirSep, 1)[0]
{.passC: "-I" & thisdir.}
{.compile: "DBX.c".}

##  Wrappers to extend DAZZ_DB.
##
##  Note that none of the extra fields are ever stored on-disk.

from os import nil
import
  db 

type
  DAZZ_DBX* {.importc: "DAZZ_DBX", header: "DBX.h".} = object
    db* {.importc: "db".}: DAZZ_DB

    ##  When "data" is non-null, it stores the entire DB
    ##  in memory, so we can avoid random-access disk operations.
    ##  But if null, then wrappers simply delegate.
    ##
    data* {.importc: "data".}: cstring


proc pOpen_DBX*(path: cstring; dbx: ptr DAZZ_DBX; preload: bool): cint {.cdecl,
    importc: "Open_DBX", header: "DBX.h".}

proc Open_DBX*(path: string, dbx: ptr DAZZ_DBX; preload: bool) =
  let ret = pOpen_DBX(path, dbx, preload)
  if ret != 0:
    os.raiseOSError(os.osLastError(), "Could not open DAZZ_DB '" & path & "'")

proc Load_ReadX*(dbx: ptr DAZZ_DBX; i: cint; read: cstring; ascii: cint): cint {.cdecl,
    importc: "Load_ReadX", header: "DBX.h", discardable.}
proc Close_DBX*(dbx: ptr DAZZ_DBX) {.cdecl, importc: "Close_DBX", header: "DBX.h".}
