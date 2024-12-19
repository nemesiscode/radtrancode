#!/usr/bin/env bash

# TODO:
# * Add command line arguments for inputs
#   - INSTALL_DEST
#   - REPO_DIR
#   - FCOMP
#   - DEBUG_EXEC
#   - LOG_FILES
#   - ISYS
# * Remove C compiler information, it is not used anywhere
#

# Use strict mode
set -o nounset -o pipefail -o errexit

echo "#######################################################"
echo "# COMPILE ON AN INTERACTIVE NODE TO AVOID STRANGENESS #"
echo "#######################################################"

## FUNCTIONS ##

get_abspath() {
	for item in $@; do
		echo "$(readlink -fn ${item})"
	done
}

get_dir() {
	abspaths=($(get_abspath $@))
	for item in "${abspaths[@]}"; do
		echo "$(dirname ${item})"
	done
}


## SCRIPT VARIABLES ##

# Constants

TRUE=0
FALSE=1
DEBUG_EXEC=$TRUE


# Variables

SCRIPT_DIR="$(get_dir $0)"

# In this instance, the script is at the top-level directory of the repository
REPO_DIR="${SCRIPT_DIR}" 
RADREPO="${REPO_DIR}"

# Files to store build logs in
declare -rA LOG_FILES=(
	[LIB]="${REPO_DIR}/make_all_lib.log"
	[BIN]="${REPO_DIR}/make_all_bin.log"
)

# Folder to install to
INSTALL_DEST="${HOME}/.local"

# build script variables
BIN=${INSTALL_DEST}/bin
#OBJ=${INSTALL_DEST}/obj/radtrancode
LIB=${INSTALL_DEST}/lib/radtrancode
RADSRC=${REPO_DIR}/radtran

# Setting for value of ${REPO_DIR}/radtran/rtm_util/isys.f
# Swap between 1 and 4. Generally use 4 for `gfortran` and 1 for `ifort`
ISYS_FILE="${REPO_DIR}/radtran/rtm_util/isys.f"
ISYS=4

# Compiler selection and options

## C compiler
# NOTE: This is not used and can be removed
CC="gcc"
CFLAGS=""

## Flags to `ar`
ARFLAGS="rvU"

## FORTRAN compiler
# FCOMP - Used in makefiles to set the fortran compiler (FC) variable
# STATIC_FLAG - Used in all makefiles to set common fortran compiler flags
# FCFLAGS1 - Used in makefiles to set the fortran compiler flags (FCFLAGS) for everything that is not fovgreg
# FCFLAGS1_FOVGREG - Used in the fovgreg makefile to set the fortran compiler flags
# DEBUG_FCFLAGS - Added to STATIC_FLAG when debugging is enabled

FCOMP="gfortran"
STATIC_FLAG=""
FCFLAGS1=""
FCFLAGS1_FOVGREG=""
DEBUG_FCFLAGS=""



case "${FCOMP}" in
	
	*ifort*)
		STATIC_FLAG="-mcmode=large -O3"
		FCFLAGS1="-i-dynamic -cm -zero -w"
		FCLAGS1_FOVGREG="-i-dynamic -cm -zero -w -w90 -w95 -FR"
		;;
	
	*gfortran*)
		STATIC_FLAG="-mcmodel=large -O3 -fallow-argument-mismatch -fallow-invalid-boz -std=legacy"
		FCFLAGS1=""
		FCFLAGS1_FOVGREG="-ffree-form"

		# NOTE: Using "-fcheck=all" means that the compiler will check the calls to functions and procedures
		# no matter if they are invoked or not. Which means compilation will fail as some files are not linked
		# properly.
		DEBUG_FCFLAGS="-g -fcheck=all"

		FC_VERSION="$(${FCOMP} --version | head -n 1 | cut -d ' ' -f 5)"
		echo "gfortran version '${FC_VERSION}'"

		# Add additional flags depending on `gfortran` version.
		case "${FC_VERSION}" in
	
			5.2.* | 7.4.* | 8.4.* )
				;;
			
			10.3.* | 13.* )
				;;
			
			*)
				echo "Unknown gfortran version '${FC_VERSION}', exiting..."
				exit 1
				;;
		
		esac

		;;
	
	*)
		echo "Unknown FORTRAN compiler '${FCOMP}', exiting..."
		exit 1
		;;
esac

if [ ${DEBUG_EXEC} -eq ${TRUE} ]; then
	echo "Debugging executable enabled."
	STATIC_FLAG="${STATIC_FLAG} ${DEBUG_FCFLAGS}"
fi



# Export variables required by makefiles etc.
export RADREPO RADSRC LIB BIN CC FCOMP CFLAGS STATIC_FLAG FCFLAGS1 FCFLAGS1_FOVGREG


echo "###########################################"
echo "Compiler information:"
echo "  CC=${CC}"
echo "    CFLAGS=${CFLAGS}"
echo "  FCOMP=${FCOMP}"
echo "    STATIC_FLAG=${STATIC_FLAG}" 
echo "    FCFLAGS1=${FCFLAGS1}"
echo "    FCFLAGS1_FOVGREG=${FCFLAGS1_FOVGREG}"
echo "  ARFLAGS=${ARFLAGS}"
echo ""
echo "Using compiler versions:"
echo "---------------    C    ---------------"
$CC --version | head -n 1
echo "---------------------------------------"
echo "--------------- FORTRAN ---------------"
$FCOMP --version | head -n 1
echo "---------------------------------------"
echo "###########################################"




## Set up environment for building

echo "Setting unlimited stack size"
ulimit -s unlimited




echo "Setting \`ar\` options to '${ARFLAGS}' for all makefiles"
ALL_MAKEFILES=(${REPO_DIR}/*/makefile ${REPO_DIR}/*/*/makefile)
for MAKEFILE in "${ALL_MAKEFILES[@]}"; do
	sed -i -r -e "s/\\bar(\\s+)((-|--)?(\\w)+)/ar ${ARFLAGS}/g" ${MAKEFILE}
done

# Fix any problems with files not included in makefiles
declare -rA MISSING_FILE_LOCATIONS=(
	[scloud12wavex.f]="${REPO_DIR}/radtran/scatter/makefile"
)

for MISSING_FILE in "${!MISSING_FILE_LOCATIONS[@]}"; do
	if ! grep -q "${MISSING_FILE}" ${MISSING_FILE_LOCATIONS[${MISSING_FILE}]}; then
		echo "Fixing file '${MISSING_FILE}' not included in makefile '${MISSING_FILE_LOCATIONS[${MISSING_FILE}]}'"
		sed -i -r -e "s/F_SOURCES =/F_SOURCES = \\\\\n\t${MISSING_FILE}/g" ${MISSING_FILE_LOCATIONS[${MISSING_FILE}]}
	fi
done



# Set ISYS value
echo "Setting ISYS=${ISYS} in '${ISYS_FILE}'"
sed -i -r -e "s/ISYS=./ISYS=${ISYS}/g" ${ISYS_FILE}


## Declare variables to hold all of the static data

declare -ra LIB_SRC_DIRS=(
	FOVgreg
	frecipes
	radtran/rtm_util
	radtran/spec_data
	radtran/path
	radtran/scatter
	radtran/radtran
	radtran/ciatable
	radtran/cirsrad
	radtran/cirsradg
	radtran/matrices
	nemesis
)

declare -ra BIN_SRC_DIRS=(
	radtran/rtm_util
	radtran/spec_data
	radtran/path
	radtran/scatter
	radtran/radtran
	radtran/ciatable
	radtran/cirsrad
	radtran/cirsradg
	nemesis
)

declare -ra LIB_FILES=(
	ciatable.a 
	cirsrad.a 
	cirsradg.a 
	lgregfov.a 
	libfrecipies.a 
	lnemesis.a 
	matricies.a 
	monteck.a 
	path.a 
	radtran.a 
	rtm_util.a 
	scatter.a 
	spec_data.a
)

declare -ra BIN_FILES=(
	Aband_ktable
	Aband_ktablec
	Aground
	Appenddb
	Ave_table
	Calc_chisq
	Calc_fnktablec
	Calc_fnktablec_dp
	Calc_lbltable
	Channel_ave_table
	CIRSdrvg_wave
	CIRSdrvg_wavePY
	CIRSdrv_wave
	CIRSdrv_wavePY
	Combi_ktable
	Combi_ktable1
	Concat_lbl_table
	Concat_table
	Concat_table_temperature
	Convertitertomre
	Convertlines
	Convert_prf
	Convert_table
	Conv_spec
	Cp_lines
	Cp_lines_seq
	Cp_lines_seq_modwid
	Cp_lines_seq_modwidJ
	Cut_table
	Direct_ktable
	Direct_ktablec
	Dust_profile
	Extract_exo_diag
	GenerateMCSspx
	Generatespx
	GeneratespxL
	Interp_prf
	Inv_scan
	Lbldrv_wave
	Li_lines
	Li_spec
	Makeband
	Makeciatable
	Makedb
	Makedbloop
	Makefptable
	Makefptable_allcia
	Makeisotable
	Makephase
	Merge
	Mergeband
	Merge_multi_ascending
	Modifyband
	Mod_profile
	Mod_table
	Mod_table1
	Mod_table1a
	Mod_table1b
	Mod_table1c
	Mod_table2
	Mod_tableco2
	Nemesis
	Nemesisdisc
	NemesisL
	Nemesisprofile
	NemesisPT
	NemesisX
	Normxsc
	Parah2_profile
	Par_ktable
	Path
	Pl_elines
	Pl_lines
	Pl_spec
	Profile
	Read_lbltable
	Readphase
	Read_table
	Scan
	Select
	Selecttemp
	Selecttemploop
	Selecttempseq
	Selecttempseqcp
	Selecttempseqcploop
	Selecttempseqtop
	Selecttempseqtopbuff
	Selecttempseqtopbuffloop
	Selecttempseqtoploop
	Sortlines
	Sort_lines_seq
	Stitchseqtable
	Summary
	Table_trans
	Testack
	Testackx
	Testbrown
	Testfracpara
	Testh2h2v2s
	Test_pmat6
	Testray
	Testsrom
	Write_path
	Write_xsec
	Zerotable
)

declare -rA DIRS=(
	[LIB]="${LIB}"
	[BIN]="${BIN}"
)


## Ensure expected directories are set up correctly

for KEY in ${!DIRS[@]}; do
	DIR="${DIRS[${KEY}]}"
	if [ -z "${DIR}" ]; then 
		echo "ERROR: Path for ${KEY} is defined as an empty string. Exiting..."
		exit 1
	elif [ -e "${DIR}" ] && [ ! -d "${DIR}" ]; then
		echo "ERROR: Path for ${KEY}, '${DIR}' exists but is not a directory. Exiting..."
		exit 1
	elif [ ! -e "${DIR}" ]; then
		echo "Path for ${KEY} does not exist, creating directory..."
		mkdir -p "${DIR}"
	else
		echo "Path for ${KEY} exists and is a directory at '${DIR}'"
		case ${KEY} in
			LIB)
				echo "Removing previous build products..."
				for FILE in "${LIB_FILES[@]}"; do
					if [ -n ${FILE} ] && [ -e "${DIR}/${FILE}" ]; then
						#echo "  ${FILE}"
						rm -f ${DIR}/${FILE}
					fi
				done
				echo "Previous build products removed."
				;;
			BIN)
				echo "Removing previous build products..."
				for FILE in "${BIN_FILES[@]}"; do
					if [ -n ${FILE} ] && [ -e "${DIR}/${FILE}" ]; then
						#echo "  ${FILE}"
						rm -f ${DIR}/${FILE}
					fi
				done
				echo "Previous build products removed."
				;;
			*)
				echo "ERROR: Unknown path name '${KEY}', exiting..."
				exit 1
				;;
		esac

	fi
done

## Rotate log files if we already have some

for KEY in ${!LOG_FILES[@]}; do
	VAL="${LOG_FILES[${KEY}]}"
	if [ -e "${VAL}" ]; then
		echo "Rotating ${KEY} log file at '${VAL}'"
		mv "${VAL}" "${VAL}.prev"
		touch "${VAL}"
	fi
done

# Turn off parts of strict mode, we want the process to complete and
# print out errors if they happen.
set +o nounset -o pipefail +o errexit


## Begin building NEMESIS

echo "Building NEMESIS..."


for DIR in "${LIB_SRC_DIRS[@]}"; do
	echo "  making libraries from '${REPO_DIR}/${DIR}'"
	echo -e "\n\n#\n# IN '${REPO_DIR}/${DIR}'\n#" >> ${LOG_FILES[LIB]}
	
	case ${DIR} in
		radtran/*)
			(
				cd "${REPO_DIR}/${DIR}"

				make clean >> ${LOG_FILES[LIB]} 2>&1
				touch *.f
				make lib >> ${LOG_FILES[LIB]} 2>&1
			)
			;;
		*)
			(	
				cd "${REPO_DIR}/${DIR}"

				make clean >> ${LOG_FILES[LIB]} 2>&1
				make lib >> ${LOG_FILES[LIB]} 2>&1;
			)
			;;
	esac
done

for DIR in "${BIN_SRC_DIRS[@]}"; do
	echo "  making binaries from '${REPO_DIR}/${DIR}'"
	echo -e "\n\n#\n# IN '${REPO_DIR}/${DIR}'\n#" >> ${LOG_FILES[BIN]}
	
	(	
		cd "${REPO_DIR}/${DIR}"

		make clean >> ${LOG_FILES[BIN]} 2>&1
		make bin >> ${LOG_FILES[BIN]} 2>&1;
	)
done

echo "### START WARNINGS ###"
grep -nHi 'warning' ${LOG_FILES[@]}
echo "###  END WARNINGS  ###"
echo ""
echo "### START ERRORS ###"
grep -nHi 'error' ${LOG_FILES[@]}
echo "###  END ERRORS  ###"


