# vim: sw=2 ts=2 sts=2 tw=80 et:
import daligner/dazz_db
import daligner/db
import daligner/dbx

from os import execShellCmd

echo "DAZZ_DB test running ..."

var mydb: db.DAZZ_DB = db.DAZZ_DB()
block:
  var badpath = "/badpath"
  doAssert -1 == db.Open_DB(badpath, addr mydb)
let goodpath = "tmpdb.db"
block:
  doAssert 0 == os.execShellCmd("fasta2DB " & goodpath & " -ifoo < /dev/null")
block:
  doAssert 0 == db.Open_DB(goodpath, addr mydb)
  db.Close_DB(addr mydb)
block:
  var mydbx: dbx.DAZZ_DBX = dbx.DAZZ_DBX()
  dbx.Open_DBX(goodpath, addr mydbx, true)
  defer: dbx.Close_DBX(addr mydbx)
block:
  doAssert 0 == os.execShellCmd("DBrm " & goodpath)

echo "DAZZ_DB test passed!"
