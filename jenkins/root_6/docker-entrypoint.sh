#!/bin/bash

function executeCommand {
    $@
    rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
    echo "Exit code[" $@ "]: $rc"
}

executeCommand "mkdir -p build"
executeCommand "cd build"
executeCommand "export CMAKE_LIBRARY_PATH=$CMAKE_LIBRARY_PATH:/framework-dependencies/lib"
executeCommand "export CMAKE_INCLUDE_PATH=$CMAKE_INCLUDE_PATH:/framework-dependencies/include"
executeCommand "export CLHEP_DIR=/clhep-install"
executeCommand "export CLHEP_INCLUDE_DIR=/clhep-install/include"
executeCommand "export CLHEP_LIBRARY=/clhep-install/lib"
executeCommand "export LD_LIBRARY_PATH=${CLHEP_LIBRARY}:${LD_LIBRARY_PATH}"
executeCommand "export PATH=/clhep-install/bin:$PATH"
executeCommand "source /geant4-install/bin/geant4.sh"
executeCommand "source /root-system/bin/thisroot.sh"
executeCommand "export LD_LIBRARY_PATH=/CADMesh-install/lib:${LD_LIBRARY_PATH}"
executeCommand "export cadmesh_DIR=/CADMesh-install"
executeCommand "cmake .."
executeCommand "make"

