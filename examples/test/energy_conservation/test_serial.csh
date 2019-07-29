../../../optimized/hd.opt constant_e.inp > constant_e.out


set EFIN = `grep "^t: 1.000000e-02 ns" constant_e.out | awk '{ print $11; }'`
set RESULT = `echo $EFIN | awk '{ del = ($1 - 1.30490035596984e+01); if( del < 1e-8 && del > -1e-8 ) { print "PASS"; } else { print "FAIL ", del; } }'`
if( $RESULT == "PASS" ) then
	set str = "Energy conservation test (in parallel) PASSED."
else
	set str = "Energy conservation test (in parallel) FAILED."
endif

echo $str
echo $str >> ../test.log
