#!/bin/bash

# build.sh
# 1 - determine host, load modules on supported hosts; proceed w/o otherwise
# 2 - configure; build; install

set -eu
CMAKE_OPTS=()

log() { echo "[build.sh] $*"; }

# Check for required tools
for tool in cmake make; do
  if ! command -v $tool >/dev/null 2>&1; then
    echo "Error: $tool not found in PATH. Please install $tool." >&2
    exit 2
  fi
done

# Print tool versions
log "cmake version: $(cmake --version | head -1)"
log "make version: $(make --version | head -1)"

# Get script directory
log "Start ... $(date)"
dir_root="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Check for required scripts
if [[ ! -f $dir_root/ush/detect_machine.sh ]]; then
  echo "Error: $dir_root/ush/detect_machine.sh not found." >&2
  exit 2
fi
source $dir_root/ush/detect_machine.sh

# ==============================================================================
usage() {
  set +x
  echo
  echo "Usage: $0 [options]"
  echo
  echo "  -p <prefix>   Installation prefix (default: $dir_root/install)"
  echo "  -t <target>   Target to build for (default: $MACHINE_ID)"
  echo "  -c <opt>      Additional CMake option (can be used multiple times)"
  echo "  -v            Build with verbose output"
  echo "  -f            Force a clean build"
  echo "  --clean       Only clean the build directory and exit"
  echo "  --install     Run 'make install' after build"
  echo "  --no-modules  Skip module-setup.sh"
  echo "  -h, --help    Display this message and quit"
  echo
  echo "Examples:"
  echo "  $0 -c '-DUSE_NETCDF=OFF' -f"
  echo "  $0 --clean"
  exit 1
}

# ==============================================================================

# Defaults:
INSTALL_PREFIX="$dir_root/install"
BUILD_TARGET="${MACHINE_ID:-localhost}"
BUILD_VERBOSE="NO"
CLEAN_BUILD="NO"
DO_CLEAN_ONLY="NO"
DO_INSTALL="NO"
COMPILER="${COMPILER:-intel}"
WORKFLOW_BUILD=${WORKFLOW_BUILD:-"OFF"}
NO_MODULES="NO"

# Detect macOS or generic Linux if MACHINE_ID is not set
if [[ -z "${MACHINE_ID:-}" ]]; then
  case "$(uname -s)" in
    Darwin)
      BUILD_TARGET="macos" ;;
    Linux)
      BUILD_TARGET="linux" ;;
    *)
      BUILD_TARGET="${MACHINE_ID:-localhost}" ;;
  esac
fi

# Parse arguments (support long options and multiple -c)
while [[ $# -gt 0 ]]; do
  case $1 in
    -p)
      INSTALL_PREFIX="$2"; shift 2;;
    -t)
      BUILD_TARGET="$2"; shift 2;;
    -c)
      CMAKE_OPTS+=("$2"); shift 2;;
    -v)
      BUILD_VERBOSE="YES"; shift;;
    -f)
      CLEAN_BUILD="YES"; shift;;
    --clean)
      DO_CLEAN_ONLY="YES"; shift;;
    --install)
      DO_INSTALL="YES"; shift;;
    --no-modules)
      NO_MODULES="YES"; shift;;
    -h|--help|\?)
      usage;;
    *)
      echo "Unknown option: $1"; usage;;
  esac
done

# Print configuration summary
log "Configuration:"
log "  INSTALL_PREFIX: $INSTALL_PREFIX"
log "  BUILD_TARGET:   $BUILD_TARGET"
log "  CLEAN_BUILD:    $CLEAN_BUILD"
log "  CMAKE_OPTS:     ${CMAKE_OPTS[@]:-}"
log "  DO_INSTALL:     $DO_INSTALL"

# Host-specific module setup
if [[ $NO_MODULES == "NO" ]]; then
  case ${BUILD_TARGET} in
    hera | orion | hercules | wcoss2 | noaacloud | gaeac5 | gaeac6 | ursa )
      log "Building on $BUILD_TARGET"
      if [[ ! -f $dir_root/ush/module-setup.sh ]]; then
        echo "Error: $dir_root/ush/module-setup.sh not found." >&2
        exit 2
      fi
      source $dir_root/ush/module-setup.sh
      module use $dir_root/modulefiles
      module load ufs_$BUILD_TARGET.$COMPILER
      module list
      ;;
    $(hostname))
      log "Building on $BUILD_TARGET"
      ;;
    macos)
      log "Detected macOS host. Skipping module setup."
      ;;
    linux)
      log "Detected generic Linux host. Skipping module setup."
      ;;
    *)
      log "Building on unknown target: $BUILD_TARGET"
      ;;
  esac
else
  log "Skipping module-setup.sh due to --no-modules option."
fi

BUILD_DIR=${BUILD_DIR:-$dir_root/build}
if [[ $CLEAN_BUILD == 'YES' || $DO_CLEAN_ONLY == 'YES' ]]; then
  log "Cleaning build directory: $BUILD_DIR"
  [[ -d ${BUILD_DIR} ]] && rm -rf ${BUILD_DIR}
  if [[ $DO_CLEAN_ONLY == 'YES' ]]; then
    log "Clean complete. Exiting."
    exit 0
  fi
fi
mkdir -p ${BUILD_DIR} && cd ${BUILD_DIR}

# Configure
log "Configuring ... $(date)"
set -x
if (( ${#CMAKE_OPTS[@]} )); then
  cmake -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" "${CMAKE_OPTS[@]}" $dir_root
else
  cmake -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" $dir_root
fi
set +x

# Build
log "Building ... $(date)"
set -x
make -j ${BUILD_JOBS:-8}
set +x

# Install if requested
if [[ $DO_INSTALL == 'YES' ]]; then
  log "Installing ... $(date)"
  set -x
  make install
  set +x
fi

# Print location of built executable(s)
log "Build finished. Executables in: $dir_root, $BUILD_DIR, or $INSTALL_PREFIX/bin (if installed)"
log "Finish ... $(date)"
exit 0