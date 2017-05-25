##  vim: set et ts=2 sts=2 sw=2 :
## 
##   Utility for displaying the overlaps in a .las file in a variety of ways including
##     a minimal listing of intervals, a cartoon, and a full out alignment.
## 
##   Author:    Gene Myers
##   Creation:  July 2013
##   Last Mod:  Jan 2015
## 
##   Based on the original LAshow.c, this code is modified by Jason Chin to support generating
##     consensus sequences from daligner output
#import parseopt as popt
from parseopt as popt import nil
from strutils import `%`
from os import nil
import
  dazz_db/DB, dazz_db/DBX, dazz_db/align, dazz_db/license_myers,
  dazz_db/common

const
  MAX_OVERLAPS* = 50000

#template MIN*(X, Y: untyped): untyped =
#  if ((X) < (Y)): (X) else: (Y)

var GROUP*: bool = false

##  Allows us to group overlaps between a pair of a/b reads as a unit, one per
##  direction (if applicable).  beg/end will point to the same overlap when
##  only one overlap found.

type
  OverlapGroup* = object
    beg*: Overlap
    `end`*: Overlap
    score*: cint
    blen*: cint

common.usePtr[OverlapGroup]()
common.usePtr[HITS_READ]()
common.usePtr[char]()

var ovlgrps*: ptr OverlapGroup

proc compare_ovlgrps*(grp1: pointer; grp2: pointer): cint {.cdecl.} =
  return (cast[ptr OverlapGroup](grp2)).score -
      (cast[ptr OverlapGroup](grp1)).score

proc belongs*(grp: ptr OverlapGroup; ovl: ptr Overlap): bool =
  var prev: ptr Overlap = addr(grp.`end`)
  return prev.flags == ovl.flags and (ovl.path.abpos > prev.path.aepos) and
      (ovl.path.bbpos > prev.path.bepos) and
      (ovl.path.abpos - prev.path.aepos) < 251

##  Add a new overlap to a new or existing overlap group. Always adds when group
##  flag is false, effectively greating groups of 1.
##  Returns 1 if added as a new overlap group, otherwise 0.
##  caller keeps track of count

proc add_overlap*(aln: ptr Alignment; ovl: ptr Overlap; count: cint): bool =
  var added: bool = false
  ##  we assume breads are in order
  if not GROUP or count < 0 or ovlgrps[count].beg.bread != ovl.bread:
    ##  Haven't seen this bread yet (or we're not grouping), move to new overlap group
    var next: ptr OverlapGroup = addr(ovlgrps[count + 1])
    next.beg = ovl[]
    next.`end` = ovl[]
    next.blen = aln.blen
    var p: ptr Path = addr(ovl.path)
    var olen: cint = p.bepos - p.bbpos
    var hlen: cint = (min(p.abpos, p.bbpos)) +
        (min(aln.alen - p.aepos, aln.blen - p.bepos))
    next.score = olen - hlen
    added = true
  else:
    var curr: ptr OverlapGroup = addr(ovlgrps[count])
    ##  Seen, should we combine it with the previous overlap group or move
    ##  on to the next?
    if belongs(curr, ovl):
      curr.`end` = ovl[]
      ##  rescore
      var beg: ptr Overlap = addr(curr.beg)
      var `end`: ptr Overlap = addr(curr.`end`)
      var olen: cint = `end`.path.bepos - beg.path.bbpos
      var hlen: cint = (min(beg.path.abpos, beg.path.bbpos)) +
          (min(aln.alen - `end`.path.aepos, aln.blen - `end`.path.bepos))
      curr.score = olen - hlen
    else:
      var next: ptr OverlapGroup = addr(ovlgrps[count + 1])
      next.beg = ovl[]
      next.`end` = ovl[]
      next.blen = aln.blen
      var p: ptr Path = addr(ovl.path)
      var olen: cint = p.bepos - p.bbpos
      var hlen: cint = (min(p.abpos, p.bbpos)) +
          (min(aln.alen - p.aepos, aln.blen - p.bepos))
      next.score = olen - hlen
      added = true
  return added

type
  qsort_cmp = proc(a,b: pointer): cint {.cdecl.}

proc qsort(p: pointer, nel, width: csize, compar: qsort_cmp) {.cdecl, importc, header:"<stdlib.h>".}

proc print_hits*(hit_count: cint; dbx2: ptr HITS_DBX; bbuffer: ptr char,
                print_scratch: var string; MAX_HIT_COUNT: cint) =
  var tmp_idx: cint
  let bsize = 1 shl 17 #capacity(print_scratch)
  qsort(ovlgrps, (hit_count + 1), sizeof((OverlapGroup)), compare_ovlgrps)
  tmp_idx = 0
  while tmp_idx < (hit_count + 1) and tmp_idx < MAX_HIT_COUNT:
    var grp: ptr OverlapGroup = addr(ovlgrps[tmp_idx])
    ## Load_ReadX assuming db2 == db1 is true
    Load_ReadX(dbx2, grp.`end`.bread, bbuffer, 0)
    if COMP(grp.`end`.flags): Complement_Seq(bbuffer, grp.blen)
    Upper_Read(bbuffer)
    var rlen: int64 = (int64)(grp.`end`.path.bepos) - (int64)(grp.beg.path.bbpos)
    if rlen < bsize:
      print_scratch.setLen(rlen)
      copyMem(addr print_scratch[0], bbuffer + grp.beg.path.bbpos, rlen)
      #strncpy(print_scratch, bbuffer + grp.beg.path.bbpos, rlen)
      printf("%08d %s\x0A", grp.`end`.bread, print_scratch)
    else:
      fprintf(stderr,
              "[WARNING]Skipping super-long read %08d, len=%lld, buf=%lld\x0A",
              grp.`end`.bread, rlen, bsize)
    inc(tmp_idx)
  printf("+ +\x0A")

var Usage* = ["[-mfsocargUFM] [-i<int(4)>] [-w<int(100)>] [-b<int(10)>] ", "    <src1:db|dam> [ <src2:db|dam> ] <align:las> [ <reads:FILE> | <reads:range> ... ]"]

const
  LAST_READ_SYMBOL* = '$'

proc ORDER*(left: pointer; rite: pointer): cint {.cdecl.} =
  var x: cint = (cast[ptr int32](left))[]
  var y: cint = (cast[ptr int32](rite))[]
  return x - y

proc SYSTEM_ERROR() =
  os.raiseOSError(os.osLastError(), "Exiting")

let Prog_Name: string = os.paramStr(0)

proc exitHelp() =
  var msg = """$# $# $#
where
-i=<Indent>
-w=<Alignment width>
-b=<Alignment border>
-H=<seed threshold (in bp)>
-n=<max numer of supporting read ouput (used for FALCON consensus. default 400, max: 2000)>

-h   help
-U   uppercase
-a   align
-r   reference
-c   cartoon
-F   flip
-M   map
-f   falcon
-s   skip
-g   group
-P   preload ##  Preload DB reads, if possible.
"""
  msg = msg % [Prog_Name, Usage[0], Usage[1]]
  stderr.write(msg)
  system.quit(system.QuitSuccess)
proc ARG_POSITIVE(flag: string, val: string): int =
  echo "flag:", flag, ", val:", repr(val)
  doAssert(len(flag) == 1)
  result = strutils.parseInt(val)
  doAssert(result > 0, "-$#=$# integer must be positive" % [flag, val])
proc ARG_NON_NEGATIVE(flag: string, val: string): int =
  doAssert(len(flag) == 1)
  result = strutils.parseInt(val)
  doAssert(result >= 0, "-$#=$# integer cannot be negative" % [flag, val])
proc main*(): int =
  var
    dbx1: HITS_DBX
  var
    dbx2: HITS_DBX
  var db1: ptr HITS_DB = addr(dbx1.db)
  var db2: ptr HITS_DB = addr(dbx2.db)
  var
    v_ovl: Overlap
    ovl: ptr Overlap = addr(v_ovl)
  var
    v_aln: Alignment
    aln: ptr Alignment = addr(v_aln)
  var input: FILE
  var novl: int64
  var
    tspace: cint
    tbytes: cint
    small: bool
  var
    reps: cint
    pts: seq[cint]
  var input_pts: bool
  var
    ALIGN: bool
    CARTOON: bool
    REFERENCE: bool
    FLIP: bool
  var
    INDENT: cint
    WIDTH: cint
    BORDER: cint
    UPPERCASE: cint
  var ISTWO: cint
  var MAP: bool
  var
    FALCON: bool
    OVERLAP: bool
    M4OVL: bool
  ##  XXX: MAX_HIT_COUNT should be renamed
  var
    SEED_MIN: cint
    MAX_HIT_COUNT: cint
    SKIP: bool
  var PRELOAD: bool
  ##   Process options
  var
    i: cint
    j: cint
    k: cint
  var flags: array[128, cint]
  var eptr: cstring
  #ARG_INIT("LA4Falcon")
  INDENT = 4
  WIDTH = 100
  BORDER = 10
  FALCON = false
  M4OVL = false
  SEED_MIN = 8000
  SKIP = false
  ALIGN = false
  REFERENCE = false
  CARTOON = false
  FLIP = false
  MAX_HIT_COUNT = 400
  var argv = @[os.paramStr(0)]
  # os.getAppFilename() to resolve

  discard popt.initOptParser()
  for kind, key, val in popt.getopt():
    case kind
    of popt.cmdArgument:
      argv.add(key)
    of popt.cmdLongOption:
      case key
      of "help": exitHelp()
      #of "version": writeVersion()
    of popt.cmdShortOption:
      case key
      of "h": exitHelp()
      of "i": INDENT = ARG_NON_NEGATIVE(key, val).cint
      of "w": WIDTH = ARG_POSITIVE(key, val).cint
      of "b": BORDER = ARG_NON_NEGATIVE(key, val).cint
      of "H": SEED_MIN = ARG_POSITIVE(key, val).cint
      of "n":
        MAX_HIT_COUNT = ARG_POSITIVE(key, val).cint
        if MAX_HIT_COUNT > 2000: MAX_HIT_COUNT = 2000
      of "U": UPPERCASE = 1
      of "a": ALIGN = true
      of "r": REFERENCE = true
      of "c": CARTOON = true
      of "F": FLIP = true
      of "M": MAP = true
      of "o": OVERLAP = true
      of "m": M4OVL = true
      of "f": FALCON = true
      of "s": SKIP = true
      of "g": GROUP = true
      of "P": PRELOAD = true
      else:
        echo "NADA!"
    of popt.cmdEnd: assert(false) # cannot happen
  echo "f=", FALCON, " P=", PRELOAD, ", o=", OVERLAP
  echo "ARGV:", repr(argv)
  let argc = len(argv)
  doAssert(argc > 2, "We need at least 2 args after prog-name and flags.")
  var dbname = argv[1]
  ##   Open trimmed DB or DB pair
  var
    pwd: string
    root: string
  ISTWO = 0.cint
  Open_DBX(dbname, addr dbx1, PRELOAD)
  if db1.part > 0:
    fprintf(stderr, "%s: Cannot be called on a block: %s\x0A", Prog_Name, dbname.cstring)
    system.quit(system.QuitFailure)
  if argc > 3:
    pwd = PathTo(argv[3])
    root = Root(argv[3], ".las")
    input = fopen(Catenate(pwd, "/", root, ".las").cstring, "r")
    if input != nil:
      ISTWO = 1
      fclose(input)
      var db2name = argv[2]
      Open_DBX(db2name, addr dbx2, PRELOAD)
      if db2.part > 0:
        fprintf(stderr, "%s: Cannot be called on a block: %s\x0A", Prog_Name,
                db2name.cstring)
        system.quit(system.QuitFailure)
      Trim_DB(db2)
    else:
      dbx2 = dbx1
      db2 = db1
  else:
    dbx2 = dbx1
    db2 = db1
  Trim_DB(db1)
  ##   Process read index arguments into a sorted list of read ranges
  input_pts = false
  if argc == ISTWO + 4:
    if argv[ISTWO + 3][0] != LAST_READ_SYMBOL or argv[ISTWO + 3][1] != '\0':
      var
        eptr: ptr char
        fptr: ptr char
      var
        b: clong
        e: clong
      var arg = argv[ISTWO + 3]
      var bptr = addr(arg[0])
      b = strtol(bptr, addr(eptr), 10)
      if eptr > bptr and b > 0:
        if eptr[] == '-':
          if eptr[1] != LAST_READ_SYMBOL or eptr[2] != '\0':
            e = strtol(eptr + 1, addr(fptr), 10)
            input_pts = (fptr <= eptr + 1 or fptr[] != '\0' or e <= 0)
        else:
          input_pts = (eptr[] != '\0')
      else:
        input_pts = true
  if input_pts:
    var
      x: cint
    var input: FILE
    #input = Fopen(argv[ISTWO + 3], "r")
    input = fopen(argv[ISTWO + 3], "r")
    doAssert(input != nil)
    reps = 0
    while true:
      let v = fscanf(input, " %d", addr(x))
      if v == -1.cint: break
      if v == 0:
        fprintf(stderr, "%s: %d\'th item of input file %s is not an integer\x0A",
                Prog_Name, reps + 1, argv[2])
        system.quit(system.QuitFailure)
      else:
        inc(reps, 1)
    reps = reps * 2
    # "Allocating read parameters"
    newSeq(pts, reps)
    rewind(input)
    var v: cint = 0
    while v < reps:
      fscanf(input, " %d", addr(x))
      pts[v] = x
      pts[v + 1] = x
      inc(v, 2)
    fclose(input)
  else:
    # "Allocating read parameters"
    newSeq(pts, 2 * argc)
    reps = 0
    if argc > 3 + ISTWO:
      var
        c: cint
        b: cint
        e: cint
      var
        eptr: ptr char
        fptr: ptr char
      c = 3 + ISTWO
      while c < argc:
        var argvc = argv[c]
        var bptr = addr argvc[0]
        if argvc[0] == LAST_READ_SYMBOL:
          b = db1.nreads
          eptr = addr(argvc[0]) + 1
        else:
          b = strtol(bptr, addr(eptr), 10).cint
        if eptr > bptr:
          if b <= 0:
            fprintf(stderr, "%s: %d is not a valid index\x0A", Prog_Name, b)
            system.quit(system.QuitFailure)
          if eptr[] == '\0':
            pts[reps] = b
            inc(reps)
            pts[reps] = b
            inc(reps)
            continue
          elif eptr[] == '-':
            if eptr[1] == LAST_READ_SYMBOL:
              e = int32.high
              fptr = eptr + 2
            else:
              e = strtol(eptr + 1, addr(fptr), 10).cint
            if fptr > eptr + 1 and fptr[] == 0.char and e > 0:
              pts[reps] = b
              inc(reps)
              pts[reps] = e
              inc(reps)
              if b > e:
                fprintf(stderr, "%s: Empty range \'%s\'\x0A", Prog_Name, argv[c])
                system.quit(system.QuitFailure)
              continue
        fprintf(stderr, "%s: argument \'%s\' is not an integer range\x0A",
                Prog_Name, argv[c])
        system.quit(system.QuitFailure)
        inc(c)
      qsort(addr(pts[0]), reps div 2, 2*sizeof(cint), ORDER)
      b = 0
      c = 0
      while c < reps:
        if b > 0 and pts[b - 1] >= pts[c] - 1:
          if pts[c + 1] > pts[b - 1]: pts[b - 1] = pts[c + 1]
        else:
          pts[b] = pts[c]
          inc(b)
          pts[b] = pts[c + 1]
          inc(b)
        inc(c, 2)
      pts[b] = int32.high
      inc(b)
      reps = b
    else:
      pts[reps] = 1
      inc(reps)
      pts[reps] = int32.high
      inc(reps)
  ##   Initiate file reading and read (novl, tspace) header
  var
    over: string
    #pwd: cstring
    #root: cstring
  pwd = PathTo(argv[2 + ISTWO])
  root = Root(argv[2 + ISTWO], ".las")
  over = Catenate(pwd, "/", root, ".las")
  #input = Fopen(over, "r")
  input = fopen(over, "r")
  doAssert(input != nil)
  if fread(addr(novl), sizeof(novl), 1, input) != 1: SYSTEM_ERROR()
  if fread(addr(tspace), sizeof(tspace), 1, input) != 1: SYSTEM_ERROR()
  if tspace == 0:
    printf("\x0ACRITICAL ERROR: tspace=0 in \'%s\'", root)
    system.quit(system.QuitFailure)
  if tspace <= TRACE_XOVR:
    small = true
    tbytes = sizeof(uint8).cint
  else:
    small = false
    tbytes = sizeof(uint16).cint
  if not (FALCON or M4OVL):
    printf("\x0A%s: ", root)
    Print_Number(novl, 0, stdout)
    printf(" records\x0A")
  #free(pwd)
  #free(root)
  ##   Read the file and display selected records
  #var j: cint
  var trace: ptr uint16
  var work: ptr Work_Data
  var tmax: cint
  var
    inside: bool
    npt: cint
    idx: cint
    ar: cint
  var tps: int64
  var p_aread: int64 = - 1
  var print_scratch = newStringOfCap(131072)
  var skip_rest: cint = 0
  var
    abuffer: ptr char
    bbuffer: ptr char
  var
    ar_wide: cint
    br_wide: cint
  var
    ai_wide: cint
    bi_wide: cint
  var
    mn_wide: cint
    mx_wide: cint
  var tp_wide: cint
  var
    blast: cint
    match: bool
    seen: bool
    lhalf: bool
    rhalf: bool
  var hit_count: cint
  aln.path = addr((ovl.path))
  if ALIGN or REFERENCE or FALCON:
    work = New_Work_Data()
    abuffer = New_Read_Buffer(db1)
    bbuffer = New_Read_Buffer(db2)
    if FALCON:
      ovlgrps = cast[ptr OverlapGroup](calloc(sizeof((OverlapGroup)), MAX_OVERLAPS + 1))
      hit_count = - 1
  else:
    abuffer = nil
    bbuffer = nil
    work = nil
  tmax = 1000
  trace = cast[ptr uint16](Malloc(sizeof(uint16) * tmax, "Allocating trace vector"))
  doAssert(trace != nil)
  inside = false
  npt = pts[0]
  idx = 1
  ar_wide = Number_Digits(cast[int64](db1.nreads))
  br_wide = Number_Digits(cast[int64](db2.nreads))
  ai_wide = Number_Digits(cast[int64](db1.maxlen))
  bi_wide = Number_Digits(cast[int64](db2.maxlen))
  if db1.maxlen < db2.maxlen:
    mn_wide = ai_wide
    mx_wide = bi_wide
    tp_wide = Number_Digits(cast[int64](db1.maxlen div tspace) + 2)
  else:
    mn_wide = bi_wide
    mx_wide = ai_wide
    tp_wide = Number_Digits(cast[int64](db2.maxlen div tspace) + 2)
  inc(ar_wide, (ar_wide - 1) div 3)
  inc(br_wide, (br_wide - 1) div 3)
  inc(ai_wide, (ai_wide - 1) div 3)
  inc(bi_wide, (bi_wide - 1) div 3)
  inc(mn_wide, (mn_wide - 1) div 3)
  inc(tp_wide, (tp_wide - 1) div 3)
  if FLIP:
    var x: cint
    x = ar_wide
    ar_wide = br_wide
    br_wide = x
    x = ai_wide
    ai_wide = bi_wide
    bi_wide = x
  blast = - 1
  match = false
  seen = false
  lhalf = false
  rhalf = false
  j = 0
  while j < novl:
    Read_Overlap(input, ovl)
    if ovl.path.tlen > tmax:
      tmax = (cint(1.2 * ovl.path.tlen.float)) + 100
      trace = cast[ptr uint16](Realloc(trace, sizeof(uint16) * tmax,
                                    "Allocating trace vector"))
      doAssert(trace != nil)
    ovl.path.trace = cast[pointer](trace)
    Read_Trace(input, ovl, tbytes)
    ##   Determine if it should be displayed
    ar = ovl.aread + 1
    if inside:
      while ar > npt:
        npt = pts[idx]
        inc(idx)
        if ar < npt:
          inside = false
          break
        npt = pts[idx]
        inc(idx)
    else:
      while ar >= npt:
        npt = pts[idx]
        inc(idx)
        if ar <= npt:
          inside = true
          break
        npt = pts[idx]
        inc(idx)
    if not inside:
      continue
    aln.alen = db1.reads[ovl.aread].rlen
    aln.blen = db2.reads[ovl.bread].rlen
    aln.flags = ovl.flags
    tps = ((ovl.path.aepos - 1) div tspace - ovl.path.abpos div tspace)
    if OVERLAP and not FALCON:
      if ovl.path.abpos != 0 and ovl.path.bbpos != 0: continue
      if ovl.path.aepos != aln.alen and ovl.path.bepos != aln.blen: continue
    if MAP:
      while ovl.bread != blast:
        if not match and seen and not (lhalf and rhalf):
          printf("Missing ")
          Print_Number(cast[int64](blast) + 1, br_wide + 1, stdout)
          printf(" %d ->%lld\x0A", db2.reads[blast].rlen, db2.reads[blast].coff)
        match = false
        seen = false
        lhalf = false
        rhalf = false
        inc(blast, 1)
      seen = true
      if ovl.path.abpos == 0: rhalf = true
      if ovl.path.aepos == aln.alen: lhalf = true
      if ovl.path.bbpos != 0 or ovl.path.bepos != aln.blen: continue
      match = true
    if not (FALCON or M4OVL):
      if ALIGN or CARTOON or REFERENCE: printf("\x0A")
      if FLIP:
        Flip_Alignment(aln, 0)
        Print_Number(cast[int64](ovl.bread) + 1, ar_wide + 1, stdout)
        printf("  ")
        Print_Number(cast[int64](ovl.aread) + 1, br_wide + 1, stdout)
      else:
        Print_Number(cast[int64](ovl.aread) + 1, ar_wide + 1, stdout)
        printf("  ")
        Print_Number(cast[int64](ovl.bread) + 1, br_wide + 1, stdout)
      if COMP(ovl.flags): printf(" c")
      else: printf(" n")
      printf("   [")
      Print_Number(cast[int64](ovl.path.abpos), ai_wide, stdout)
      printf("..")
      Print_Number(cast[int64](ovl.path.aepos), ai_wide, stdout)
      printf("] x [")
      Print_Number(cast[int64](ovl.path.bbpos), bi_wide, stdout)
      printf("..")
      Print_Number(cast[int64](ovl.path.bepos), bi_wide, stdout)
      printf("]")
    if M4OVL:
      var
        bbpos: int64
        bepos: int64
      var acc: cdouble
      if COMP(ovl.flags):
        bbpos = cast[int64](aln.blen) - cast[int64](ovl.path.bepos)
        bepos = cast[int64](aln.blen) - cast[int64](ovl.path.bbpos)
      else:
        bbpos = cast[int64](ovl.path.bbpos)
        bepos = cast[int64](ovl.path.bepos)
      acc = 100 -
          (200.0 * ovl.path.diffs.float) /
          cdouble(ovl.path.aepos - ovl.path.abpos + ovl.path.bepos - ovl.path.bbpos)
      printf("%09lld %09lld %lld %5.2f ", cast[int64](ovl.aread),
             cast[int64](ovl.bread), cast[int64](bbpos) - cast[int64](bepos), acc)
      printf("0 %lld %lld %lld ", cast[int64](ovl.path.abpos),
             cast[int64](ovl.path.aepos), cast[int64](aln.alen))
      printf("%d %lld %lld %lld ", COMP(ovl.flags), bbpos, bepos,
             cast[int64](aln.blen))
      if (cast[int64](aln.blen) < cast[int64](aln.alen)) and
          (cast[int64](ovl.path.bbpos) < 1) and
          (cast[int64](aln.blen) - cast[int64](ovl.path.bepos) < 1):
        printf("contains\x0A")
      elif (cast[int64](aln.alen) < cast[int64](aln.blen)) and
          (cast[int64](ovl.path.abpos) < 1) and
          (cast[int64](aln.alen) - cast[int64](ovl.path.aepos) < 1):
        printf("contained\x0A")
      else:
        printf("overlap\x0A")
    if FALCON:
      if p_aread == - 1:
        Load_ReadX(addr dbx1, ovl.aread, abuffer, 2)
        printf("%08d %s\x0A", ovl.aread, abuffer)
        p_aread = ovl.aread
        skip_rest = 0
      if p_aread != ovl.aread:
        print_hits(hit_count, addr dbx2, bbuffer, print_scratch,
                   MAX_HIT_COUNT)
        hit_count = - 1
        Load_ReadX(addr dbx1, ovl.aread, abuffer, 2)
        printf("%08d %s\x0A", ovl.aread, abuffer)
        p_aread = ovl.aread
        skip_rest = 0
      if skip_rest == 0:
        if add_overlap(aln, ovl, hit_count): inc(hit_count)
        if (hit_count + 1) > MAX_OVERLAPS: skip_rest = 1
        #[
        if false:
          tps = ((ovl.path.aepos - 1) div tspace - ovl.path.abpos div tspace)
          if small: Decompress_TraceTo16(ovl)
          Load_ReadX(addr dbx1, ovl.aread, abuffer, 0)
          Load_ReadX(addr dbx2, ovl.bread, bbuffer, 0)
          if COMP(aln.flags): Complement_Seq(bbuffer, aln.blen)
          Compute_Trace_PTS(aln, work, tspace) #?
          var tlen: cint = aln.path.tlen
          var trace: ptr cint = aln.path.trace
          var u: cint
          printf(" ")
          u = 0
          while u < tlen:
            printf("%d,", int16(trace[u]))
            inc(u)
        ]#
        if SKIP:
          if (int64(aln.alen) < int64(aln.blen)) and
              (int64(ovl.path.abpos) < 1) and
              (int64(aln.alen) - int64(ovl.path.aepos) < 1):
            printf("* *\x0A")
            skip_rest = 1
        ## if not SKIP, then skip_rest is always 0
    if ALIGN or CARTOON or REFERENCE:
      if ALIGN or REFERENCE:
        var
          aseq: ptr char
          bseq: ptr char
        var
          amin: cint
          amax: cint
        var
          bmin: cint
          bmax: cint
        if FLIP: Flip_Alignment(aln, 0)
        if small: Decompress_TraceTo16(ovl)
        amin = ovl.path.abpos - BORDER
        if amin < 0: amin = 0
        amax = ovl.path.aepos + BORDER
        if amax > aln.alen: amax = aln.alen
        if COMP(aln.flags):
          bmin = (aln.blen - ovl.path.bepos) - BORDER
          if bmin < 0: bmin = 0
          bmax = (aln.blen - ovl.path.bbpos) + BORDER
          if bmax > aln.blen: bmax = aln.blen
        else:
          bmin = ovl.path.bbpos - BORDER
          if bmin < 0: bmin = 0
          bmax = ovl.path.bepos + BORDER
          if bmax > aln.blen: bmax = aln.blen
        aseq = Load_Subread(db1, ovl.aread, amin, amax, abuffer, 0)
        bseq = Load_Subread(db2, ovl.bread, bmin, bmax, bbuffer, 0)
        aln.aseq = aseq - amin
        if COMP(aln.flags):
          Complement_Seq(bseq, bmax - bmin)
          aln.bseq = bseq - (aln.blen - bmax)
        else:
          aln.bseq = bseq - bmin
        Compute_Trace_PTS(aln, work, tspace, GREEDIEST)
        if FLIP:
          if COMP(aln.flags):
            Complement_Seq(aseq, amax - amin)
            Complement_Seq(bseq, bmax - bmin)
            aln.aseq = aseq - (aln.alen - amax)
            aln.bseq = bseq - bmin
          Flip_Alignment(aln, 1)
      if CARTOON:
        printf("  (")
        Print_Number(tps, tp_wide, stdout)
        printf(" trace pts)\x0A\x0A")
        Alignment_Cartoon(stdout, aln, INDENT, mx_wide)
      else:
        printf(" :   = ")
        Print_Number(cast[int64](ovl.path.diffs), mn_wide, stdout)
        printf(" diffs  (")
        Print_Number(tps, tp_wide, stdout)
        printf(" trace pts)\x0A")
      if REFERENCE:
        Print_Reference(stdout, aln, work, INDENT, WIDTH, BORDER, UPPERCASE, mx_wide)
      if ALIGN:
        Print_Alignment(stdout, aln, work, INDENT, WIDTH, BORDER, UPPERCASE, mx_wide)
    elif not (FALCON or M4OVL):
      printf(" :   < ")
      Print_Number(cast[int64](ovl.path.diffs), mn_wide, stdout)
      printf(" diffs  (")
      Print_Number(tps, tp_wide, stdout)
      printf(" trace pts)\x0A")
    inc(j)                    ##   Read it in
  if FALCON and hit_count != - 1:
    print_hits(hit_count, addr dbx2, bbuffer, print_scratch,
               MAX_HIT_COUNT)
    printf("- -\x0A")
    free(ovlgrps)
  free(trace)
  if ALIGN or FALCON:
    free(bbuffer - 1)
    free(abuffer - 1)
    Free_Work_Data(work)
  Close_DBX(addr dbx1)
  if ISTWO != 0: Close_DBX(addr dbx2)
  system.quit(system.QuitSuccess)

if isMainModule:
  system.programResult = main()
