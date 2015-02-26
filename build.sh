#!/bin/bash
if [ ! -d ${BIOTRUMP_OUT} ]; then
  echo "${BIOTRUMP_OUT} does not exist. mkdir"
  mkdir ${BIOTRUMP_OUT}
fi
if [ ! -d ${ATLAS_OUT} ]; then
  echo "${ATLAS_OUT} does not exist. mkdir"
  mkdir ${ATLAS_OUT}
fi
if [ ! -d ${ATLAS_OUT} ]; then
  echo "${ATLAS_OUT} does not exis. mkdir"
  mkdir ${ATLAS_OUT}
fi
pushd ${ATLAS_OUT}

case `uname` in
"Darwin")
	# Should also work on other BSDs
	CORE_COUNT=`sysctl -n hw.ncpu`
	;;
"Linux")
	CORE_COUNT=`grep processor /proc/cpuinfo | wc -l`
	if [[ `uname -m` =~ "x86" ]]; then
	#convert KHz to MHz
		MAX_FREQ=`cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_max_freq | awk '{printf "%d\n", $1/1000}'`
		DMAX_SPEED="-D c -DPentiumCPS=${MAX_FREQ}"
	fi
	#disabling cpu freq to top speed
	for((i=0;i<CORE_COUNT;i=i+1))
	do
		sudo cpufreq-set -c $i -g performance
	done

	;;
CYGWIN*)
	CORE_COUNT=`grep processor /proc/cpuinfo | wc -l`
	;;
*)
	echo Unsupported platform: `uname`
	exit -1
esac

#x86 corei3,i5,i7 with Intel AVX2 and FMA
# GCC 4.7, GFortran 4.7 and above!

#specify new gcc-4.9 and gfortran-4.9 to override 4.7
#http://math-atlas.sourceforge.net/atlas_install/node18.html
#If you use a different gcc than 4.7.0, you may reduce performance.
#Therefore once your build is finished, you should make sure to
#compare your achieved performance against what ATLAS's architectural defaults achieved.
#   mkdir my_build_dir ; cd my_build_dir
#   /path/to/ATLAS/configure [flags]
# PentiumCPS=3600 => 3.6Ghz CPU
#../configure -D c -DPentiumCPS=3600 -C acg /usr/bin/gcc-4.9 -C if /usr/bin/gfortran-4.9
#gcc 4.7 seems better performance than the other version.
#--with-netlib-lapack=/home/thomas/build/biotrump-cv/out/lapack/lib
${ATLAS_SRC}/configure ${DMAX_SPEED}
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
make check
make time

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
