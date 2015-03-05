#!/bin/sh
pre=d
la=_al
bl=_ab
N0=0
efs=0
EF=
lrts="getrf potrf geqrf"
while [ $# -gt 0 ]
do
   case $1 in
   -p)  pre=$2
        shift 2
        ;;
   -B)  bl=$2
        shift 2
        ;;
   -L)  la=$2
        shift 2
        ;;
   -E)  EF=$2
        efs=`echo "${EF}" | ${chksum} | cut -d ' ' -f 1`
        shift 2
        ;;
   -N)  N0=$2
        NN=$3
        Ni=$4
        reps=$5
        shift 5
        ;;
   *)   echo "USAGE: $0 [-N N0 NN Ni reps] [-p <pre>] [-L <lapack>] -B [<blas>]"
        echo "   N params & reps are for timings > 6000, by default not done"
        echo "   lapack : _al or _sl"
        echo "   blas   : _ab or _sb"
        echo "   pre    : z, c, d or s"
        exit 1
        ;;
   esac
done

for rt in ${lrts}
do
   echo "\nRunning 200-2000 ${rt} with 20 reps"
   ./tvrun_rng.sh x${pre}tlatime${la}${bl} ${rt} 200 2000 200 20 -L _${rt}
   echo "\nRunning 3000-6000 ${rt} with 4 reps"
   ./tvrun_rng.sh x${pre}tlatime${la}${bl} ${rt} 3000 6000 1000 4 -L _${rt}
   if [ $N0 -gt 6000 ]
   then
      echo "\nRunning ${N0}-${NN} ${rt} with ${reps} reps"
      ./tvrun_rng.sh x${pre}tlatime${la}${bl} ${rt} $N0 $NN $Ni $reps -L _${rt}
   fi
done
rt=gemm
echo "\nRunning 200-2000 ${rt} with 20 reps"
./tvrun_rng.sh x${pre}tl3time${bl} ${rt} 200 2000 200 20 -L _${rt}
echo "\nRunning 3000-6000 ${rt} with 4 reps"
./tvrun_rng.sh x${pre}tl3time${bl} ${rt} 3000 6000 1000 4 -L _${rt}
if [ $N0 -gt 6000 ]
then
   echo "\nRunning ${N0}-${NN} ${rt} with ${reps} reps"
   ./tvrun_rng.sh x${pre}tl3time${bl} ${rt} $N0 $NN $Ni $reps -L _${rt}
fi
