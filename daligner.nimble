# Package

version       = "0.0.0"
author        = "Christopher Dunn"
description   = "Nim wrappers/translations of DALIGNER (Gene Myers), with additions by PacBio (Jason Chin, et al.)"
license       = "BSD 3 Clear"

# Dependencies

requires "nim >= 0.17.0"

srcDir = "./src"

if not fileExists("repos/DALIGNER/align.h"):
    let msg = "git submodule update --init"
    echo msg
    exec(msg)

task test, "Test daligner wrapper":
    withDir("tests"):
        exec("make")
