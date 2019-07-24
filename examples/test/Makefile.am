TESTS = energy_conservation.check

all: $(TESTS)

%.check:
	cd $* ; csh setup.csh
if PARALLEL
	cd $* ; csh test_parallel.csh 
else
	cd $* ; csh test_serial.csh
endif		
#energy_conservation:
#	cd energy_conservation
#	csh setup.csh
##if PARALLEL
#	csh test_parallel.csh 
#	if( $status == "0" ) then
#		echo "Energy conservation (parallel) PASSED" >> test.log
#	else
#		echo "Energy conservation (parallel) FAILED" >> test.log
##else
#	csh test_serial.csh
##endif	
