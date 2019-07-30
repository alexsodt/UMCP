rm -f will_save.save

mpirun -np 2 ../../../optimized/hd.opt save.inp > save.out

if( ! -e will_save.save ) then
	echo "Load/save test FAILED. Didn't make 'will_save.save'" >> ../test.log
	exit
endif

set e1 = `grep -i "energy at save" save.out | awk '{ print $4; }'`
mpirun -np 2 ../../../optimized/hd.opt load.inp > load.out
set e2 = `grep -i "^t: " load.out | head -1 | awk '{ print $9; }'`

set RESULT = `echo $e1 $e2 | awk '{ del = ($1 - $2); if( del < 1e-6 && del > -1e-6 ) { print "PASS"; } else { print "FAIL "; } }'`


if( $RESULT == "PASS" ) then
	set str = "Load/save test PASSED (in parallel)."
else
	set str = "Load/save test FAILED with an energy violation (in parallel)."
endif

echo $str
echo $str >> ../test.log
