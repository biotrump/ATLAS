#!/bin/sh
atldir=\$1
shift
atlexe=\$1
shift
#
# Replace this if with your setup structure
#
if [ ! -e "/tmp/evm-input.fifo" ]
then
   pushd ~/runscript
   ./mmsetup
   popd
fi
start=\$((\$(date +%s%N)/1000000))
echo \"[\`date -R\`] \$atldir/\$atlexe.out \$*\" >> ATLrun.log
echo \"\$atldir/\$atlexe.out \$*\" >> ATLrun.log

#
# Replace this line with your accelerator runscript
#
~/runscript/mmrun \$atldir/\$atlexe.out \$*

finish=\$((\$(date +%s%N)/1000000))
echo \"[\`date -R\`] \$atldir/\$atlexe.out \$*\" >> ATLrun.log
let z=finish-start
duration=\`bc <<< \"scale=3; \$z/1000\"\`
echo \"[\`date -R\`] command finished in \$duration seconds\" >> ATLrun.log
