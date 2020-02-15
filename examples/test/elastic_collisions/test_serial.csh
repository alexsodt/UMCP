../../../optimized/hd.opt collisions.inp > collisions.out

set nl = `grep NCol collisions.out | wc -l`
set elast = `grep "Energy at save" collisions.out | awk '{ print $4; }'`

if( $nl != "9" ) then
	set str = "Collision check FAILED to read nine collision lines"
	echo $str
	echo $str > ../test.log
	exit
endif

set RESULT = `echo $elast | awk '{ del = ($1 - 3.18582334129301e+02 ); if( del < 1e-8 && del > -1e-8 ) { print "PASS"; } else { print "FAIL "; } }'`

if( $RESULT == "PASS" ) then
	set str = "Collision test (in serial) PASSED."
else
	set str = "Collision test FAILED with a bad energy (in serial)."
endif

echo $str
echo $str >> ../test.log
