# vim: sw=2 ts=2 sts=2 tw=80 et:
import math

proc throw*(msg: string) =
    raise newException(Exception, msg)

# https://akehrer.github.io/nim/2015/01/14/getting-started-with-nim-pt2.html
#const
#  Nan = 0.0/0.0 # floating point not a number (NaN)

proc cIsNaN(x: float): cint {.importc: "isnan", header: "<math.h>".}
  ## returns non-zero if x is not a number

proc cIsInf(x: float): cint {.importc: "isinf", header: "<math.h>".}
  ## returns non-zero if x is infinity

proc isNaN*(x: float): bool =
  ## converts the integer result from cIsNaN to a boolean
  if cIsNaN(x) != 0.cint:
    true
  else:
    false

proc isInf*(x: float): bool =
  ## converts the integer result from cIsInf to a boolean
  if cIsInf(x) != 0.cint:
    true
  else:
    false

# For ptr arithmetic: https://forum.nim-lang.org/t/1188#7366
template usePtr*[T] =
  template `+`(p: ptr T, off: SomeInteger): ptr T =
    cast[ptr type(p[])](cast[ByteAddress](p) +% int(off) * sizeof(p[]))

  template `+=`(p: ptr T, off: SomeInteger) =
    p = p + off

  template `-`(p: ptr T, off: SomeInteger): ptr T =
    cast[ptr type(p[])](cast[ByteAddress](p) -% int(off) * sizeof(p[]))

  template `-`(p: ptr T, off: ptr T): ByteAddress =
    (cast[ByteAddress](p) -% cast[ByteAddress](off))

  template `-=`(p: ptr T, off: SomeInteger) =
    p = p - int(off)

  template `[]`(p: ptr T, off: SomeInteger): T =
    (p + int(off))[]

  template `[]=`(p: ptr T, off: SomeInteger, val: T) =
    (p + off)[] = val

# https://forum.nim-lang.org/t/2943
template asarray*[T](p: pointer): auto =
  type A {.unchecked.} = array[0..0, T]
  cast[ptr A](p)

# https://forum.nim-lang.org/t/1353
# var buf = newString(8)
# buf.sprintf(format, ...)
# buf.setLen(buf.cstring.len)
#proc sprintf*(buf, format: cstring) {.header: "<stdio.h>", importc, varargs, noSideEffect.}

proc printf*(format: cstring) {.header: "<stdio.h>", importc, varargs.}
proc fprintf*(stream: File, formatstr: cstring) {.header: "<stdio.h>", varargs.}
proc fopen*(filename, mode: cstring): File {.header: "<stdio.h>", importc.}
proc fclose*(f: File): cint {.header: "<stdio.h>", importc, discardable.}
proc fread*(p: pointer, size: csize, nitems: csize, f: File): csize {.header: "<stdio.h>", importc.}
proc rewind*(c: File) {.importc, header: "<stdio.h>".}
proc free*(p: pointer) {.header: "<stdlib.h>", importc.}
proc calloc*(count, size: csize): pointer {.header: "<stdlib.h>", importc.}
proc strtol*(str: ptr char, endptr: ptr ptr char, base: int): clong {.cdecl, importc, header: "inttypes.h".}
proc fscanf*(c: File, frmt: cstring): cint {.varargs, importc, header: "<stdio.h>", discardable.}
