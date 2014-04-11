SHELL = /bin/sh
CFLAGS = -O2 -DNDEBUG -I/Users/doyen/install/boost/include -L/Users/doyen/install/boost/lib
VERSION = `grep 'define BEAND_VERSION_STRING' beand.cpp | sed 's/.*"\(.*\)"/\1/'`

bns: beand.cpp lognum.hpp
	g++ $(CFLAGS) beand.cpp -o beandprior -lboost_program_options

clean:
	rm -f beandprior beandprior.tar.gz gmon.out

recompile: clean beandprior

package: clean
	tar -cpzf BEANDiscoPrior-$(VERSION).tgz *

