# Package

version       = "0.0.0"
author        = "Christopher Dunn"
description   = "Nim wrappers/translations of DALIGNER (Gene Myers), with additions by PacBio (Jason Chin, et al.)"
license       = "BSD 3 Clear"

# Dependencies

requires "nim >= 0.17.0"

srcDir = "./src"
installDirs = @["daligner/"]
bin = @["daligner/LA4Falcon"]

task test, "Test daligner wrapper":
    withDir("tests"):
        exec("make")
