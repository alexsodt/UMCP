TESTS = energy_conservation.check load_save.check elastic_collisions.check

all: 
	@rm -f test.log
	@make $(TESTS)
	@echo "---FINAL TEST LOG---"
	@cat test.log
	@echo "---FINAL TEST LOG---"

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
