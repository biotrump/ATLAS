#!/bin/bash
ATLAS_OUT=${ATLAS_OUT:-`pwd`/build}
ATLAS_SRC=${ATLAS_SRC:-`pwd`}
ATLAS_BRANCH=${ATLAS_BRANCH:-Rev1531}

export ATLAS_OUT
export ATLAS_SRC
export ATLAS_BRANCH
echo $ATLAS_OUT
echo $ATLAS_SRC
#exit
#if [ ! -d ${BIOTRUMP_OUT} ]; then
#  echo "${BIOTRUMP_OUT} does not exist. mkdir"
#  mkdir ${BIOTRUMP_OUT}
#fi
if [ ! -d ${ATLAS_OUT} ]; then
  echo "${ATLAS_OUT} does not exist. mkdir"
  mkdir ${ATLAS_OUT}
fi
#if [ ! -d ${ATLAS_OUT} ]; then
#  echo "${ATLAS_OUT} does not exis. mkdir"
#  mkdir ${ATLAS_OUT}
#fi
#if [[ "$#" -eq 0 || "$#" -eq 1 && "$1" == "-j"* ]]; then

case `uname` in
"Darwin")
	# Should also work on other BSDs
	CORE_COUNT=`sysctl -n hw.ncpu`
	;;
"Linux")
	CORE_COUNT=`grep processor /proc/cpuinfo | wc -l`
	;;
CYGWIN*)
	CORE_COUNT=`grep processor /proc/cpuinfo | wc -l`
	;;
*)
	echo Unsupported platform: `uname`
	exit -1
esac

#config or build
case $1 in
	"config" )

	pushd ${ATLAS_OUT}
	#x86 target : speed up to top speed,
	#so ATALAS can tune it for x86.
	#x86 corei3,i5,i7 with Intel AVX2 and FMA
	#ARM target : configure -Si archdef 0 -Fa al -mtune=native
	#http://www.vesperix.com/arm/atlas-arm/source/
	case `uname` in
	"Linux")
		if [[ "$#" -eq 1 && `uname -m` =~ "x86" ]]; then
		#default build is x86: convert KHz to MHz
			MAX_FREQ=`cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_max_freq | awk '{printf "%d\n", $1/1000}'`
			DMAX_SPEED="-D c -DPentiumCPS=${MAX_FREQ}"
			MACHINE_FLAG="-Fa alg -fPIC --nof77 -b 64 ${DMAX_SPEED}"
			echo MACHINE_FLAG=$MACHINE_FLAG
			#disabling cpu freq to top speed
			for((i=0;i<${CORE_COUNT};i=i+1))
			do
				sudo cpufreq-set -c $i -g performance
			done
		else
			if [[ "$#" -eq 2 && "$2" =~ "arm" ]]; then
			MACHINE_FLAG="-Si archdef 0 -Fa al -mtune=native -fPIC --nof77"
			fi
		fi
		;;

#	"Darwin")
#	CYGWIN*)

	*)
		echo Unsupported platform: `uname`
		exit -1
	;;
	esac

	# ATLAS needs GCC 4.7, GFortran 4.7 and above!
	#, but specify newer gcc-4.9 and gfortran-4.9 to override 4.7
	#, you may reduce performance.
	# http://math-atlas.sourceforge.net/atlas_install/node18.html

	#   mkdir my_build_dir ; cd my_build_dir
	#   /path/to/ATLAS/configure [flags]
	# PentiumCPS=3600 => 3.6Ghz CPU
	#../configure -D c -DPentiumCPS=3600 -C acg /usr/bin/gcc-4.9 -C if /usr/bin/gfortran-4.9
	#gcc 4.7 seems better performance than the other version.
	#--with-netlib-lapack=/home/thomas/build/biotrump-cv/out/lapack/lib
	#--with-netlib-lapack-tarfile=/home/thomas/build/lapack.tar
	#${ATLAS_SRC}/configure --with-netlib-lapack=/home/thomas/build/lapack/liblapack.a \
	#--nof77 ${DMAX_SPEED}
	#-fpie -fPIE for  position independent code can be only linked into executables
	#-fpic -fPIc for  position independent code
	#-Fa alg: If you use non-gnu compilers, you'll need to use -Fa
	#pass the correct flag(s) to append to force position independent code for
	#each compiler (don't forget the gcc compiler used in the index files).
	#-b 64 : which says to perform a 64 bit compile
	#echo $ATLAS_SRC

	${ATLAS_SRC}/configure ${MACHINE_FLAG}  \
	--with-netlib-lapack-git=https://github.com/biotrump/lapack.git,${ATLAS_BRANCH}:remotes/origin/${ATLAS_BRANCH}

#	${ATLAS_SRC}/configure -Fa alg -fPIC --nof77 ${DMAX_SPEED} -b 64 \
#	--with-netlib-lapack-git=https://github.com/biotrump/lapack.git,${ATLAS_BRANCH}:remotes/origin/${ATLAS_BRANCH}

	#${ATLAS_SRC}/configure --with-netlib-lapack-tarfile=/home/thomas/build/lapack.tar \
	#-Fa alg -fPIC --nof77 ${DMAX_SPEED} -b 64

	#${ATLAS_SRC}/configure ${DMAX_SPEED}

	#back to ondemand mode to save power
	#
	case `uname` in
	"Darwin")
		# Should also work on other BSDs
		;;
	"Linux")
		for((i=0;i<CORE_COUNT;i=i+1))
		do
			sudo cpufreq-set -c $i -g ondemand
		done
		;;
	CYGWIN*)
		;;
	*)
		echo Unsupported platform: `uname`
		exit -1
	esac
	popd
	;;

	"build" )
	pushd ${ATLAS_OUT}
	pushd ${ATLAS_OUT}/src/lapack/reference
	git pull
	popd
	#   make              ! tune and compile library
	#   make check        ! perform sanity tests
	#   make ptcheck      ! checks of threaded code for multiprocessor systems
	#   make time         ! provide performance summary as % of clock rate
	#   make install      ! Copy library and include files to other directories
	#gcc-4.9
	#Reference clock rate=3500Mhz, new rate=3601Mhz
	#   Refrenc : % of clock rate achieved by reference install
	#   Present : % of clock rate achieved by present ATLAS install

	#                    single precision                  double precision
	#            ********************************   *******************************
	#                  real           complex           real           complex
	#            ---------------  ---------------  ---------------  ---------------
	#Benchmark   Refrenc Present  Refrenc Present  Refrenc Present  Refrenc Present
	#=========   ======= =======  ======= =======  ======= =======  ======= =======
	#  kSelMM     1439.6  1601.2   1294.6  1434.9    757.1   839.1    700.6   775.6
	#  kGenMM      351.6   374.5    349.9   385.9    344.9   374.7    343.4   377.3
	#  kMM_NT      349.2   377.4    339.6   368.4    334.0   370.8    324.4   372.3
	#  kMM_TN      345.1   384.8    324.9   356.3    315.7   355.5    314.1   358.5
	#  BIG_MM     1366.7  1506.8   1393.7  1523.0    719.6   795.3    723.5   778.0
	#   kMV_N      235.6   210.3    484.5   470.5    112.4   116.0    216.2   207.0
	#   kMV_T      234.3   217.0    488.0   477.9    120.0   113.0    239.5   234.3
	#    kGER      166.9   155.7    320.9   300.2     78.5    75.3    160.5   151.8
	make -j ${CORE_COUNT}
	#we need lapacke to export c api
	ln -s ${ATLAS_OUT}/src/lapack/reference/liblapacke.a ${ATLAS_OUT}/lib
	make check
	make time
	popd
	;;

	* )
	atlas.sh config/build
	;;
esac