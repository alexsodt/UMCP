mpirun -np 2 ../../../optimized/hd.opt constant_e.inp > constant_e.out

set green = '\033[0;32m'
set red = '\033[0;31m'
set off = '\033[0m'

set EFIN = `grep "^t: 1.000000e-02 ns" constant_e.out | awk '{ print $11; }'`
set RESULT = `echo $EFIN | awk '{ del = ($1 - 1.30913883543231e+01 ); if( del < 1e-8 && del > -1e-8 ) { print "PASS"; } else { print "FAIL "; } }'`


if( $RESULT == "PASS" ) then
	echo ${green}
	set str = "Energy conservation test (in parallel) PASSED."
else
	echo ${red}
	set str = "Energy conservation test (in parallel) FAILED."
endif

echo $str
echo ${off}
echo $str >> ../test.log
