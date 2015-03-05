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
BINdir=../bin
t0=ttpk0.tv
t1=ttpk1.tv
rm -f $t1 $t2
if [ $N0 -gt 6000 ]
then
   eb=res/t${pre}gemm${bl}_${N0}_${NN}_${Ni}.out 
   el="res/t${pre}getrf${la}${bl}_${N0}_${NN}_${Ni}.out \
      res/t${pre}potrf${la}${bl}_${N0}_${NN}_${Ni}.out \
      res/t${pre}geqrf${la}${bl}_${N0}_${NN}_${Ni}.out"
else
   eb=
   el=
fi
cat res/t${pre}gemm${bl}_200_2000_200_20_${efs}.out \
    res/t${pre}gemm${bl}_3000_6000_1000_4_${efs}.out ${eb} \
    | $BINdir/tvecget -H 2 N MFLOP \
    | $BINdir/tvecreduce -R 1 N_gemm -C 1 MFLOP_gemm -o $t0
cat res/t${pre}getrf${la}${bl}_200_2000_200_20_${efs}.out \
    res/t${pre}potrf${la}${bl}_200_2000_200_20_${efs}.out \
    res/t${pre}geqrf${la}${bl}_200_2000_200_20_${efs}.out \
    res/t${pre}getrf${la}${bl}_3000_6000_1000_4_${efs}.out \
    res/t${pre}potrf${la}${bl}_3000_6000_1000_4_${efs}.out \
    res/t${pre}geqrf${la}${bl}_3000_6000_1000_4_${efs}.out ${el} \
    | $BINdir/tvecget -H 2 N MFLOP \
    | $BINdir/tvecreduce -C 3 MFLOP_getrf MFLOP_potrf MFLOP_geqrf -o $t1 
cat $t0 $t1 | $BINdir/tvecselect -# 2 \
    -S 2 N_gemm MFLOP_gemm_avg -R 2 N_gemm N MFLOP_gemm_avg MFLOP_gemm \
    -S 3 MFLOP_getrf_avg MFLOP_potrf_avg MFLOP_geqrf_avg \
    -R 3 MFLOP_getrf_avg MFLOP_getrf MFLOP_potrf_avg MFLOP_potrf \
         MFLOP_geqrf_avg MFLOP_geqrf \
    | $BINdir/tvecunify -r 1 \
    | $BINdir/tvecscale -m 1.0 -b MFLOP_gemm -C 1 MFLOP_gemm \
    -R 3 MFLOP_getrf MFLOP_potrf MFLOP_geqrf \
    | $BINdir/tvecrename \
    -R 3 MFLOP_getrf GETRFpcMM MFLOP_potrf POTRFpcMM MFLOP_geqrf GEQRFpcMM \
    | $BINdir/tvecprint -h 1 \
    -C 5 N MFLOP_gemm GETRFpcMM POTRFpcMM GEQRFpcMM
