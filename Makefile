include ../c_flags

###########################################################################
#
#        Options
#   Examples of depending on libaries
#        LDLIBS =  -lthing -lregression 
#        libs_used =  ../lib/libthing.a ../lib/libregression.a 
# 
###########################################################################


PROJECT_NAME = bellman

# OPT = -O3 -std=c++0x -DNDEBUG

OPT = -O3  -std=c++0x


USES = utils random

level_1 = distribution.o
level_2 = wealth.o
level_3 = utility.o
level_4 = bellman.o

############################################################################
#
#            INCLUDING RULES AND DEFINITIONS
#
###########################################################################

# TAGS 
# find . | grep ".*\.\(h\|cc\)" | xargs etags -f TAGS


# find utility for simulated mean process
simulate_risk: simulate
	./simulate

# -------------------------------------------------------------------
# bellman recursion for competitive value
#
# constrained:  univ univ 2 0.05   50   0.5 7     -0.015669 0.015669 0.015669
# uncons      uncon g50 2 0.05   50   0.05 7     37.1289 37.7115 0.291287
# -------------------------------------------------------------------

compiler: bellman
	gcc --version

constrained_test: bellman
	./bellman --gamma 2 --rounds 50 --constrain --oracleprob 0.5 --bidderprob 0.0 --write    # geometric oracle
	./bellman --gamma 2 --rounds 50 --constrain --oracleprob 0.0 --bidderprob 0.5 --write    # univ oracle
	./bellman --gamma 2 --rounds 50             --oracleprob 0.5 --bidderprob 0.0 --write    # warning message

unconstrained_test: bellman
	./bellman --gamma 2 --rounds 50                              --bidderprob 0.5 --write   



bellman_main.o: bellman_main.cc

bellman: bellman.o wealth.o utility.o distribution.o bellman_main.o
	$(GCC) $^ $(LDLIBS) -o  $@

simulate: bellman.o wealth.o utility.o distribution.o bellman_simulator.o
	$(GCC) $^ $(LDLIBS) -o  $@

# --- using latest code, angular style with sine and cosine
# g01000 univ(1) 1 2       153.435 0.05   200   0.05 10     -1.35152 -0.622174 -0.364673                                                                                   

bellman_test: bellman 
	./bellman --risk --omega 0.5 --angle 210 --rounds 200  --oracleprob 0.05 --bidderprob 0 --scale 1 --write   # unconstrained

# uncon(0.05) scaled_univ(1) UTIL: * Warning *  Bid beta too large; reduced to 0.99 UTIL: Message limit reached. 
# 210 0.5   200   0.05 20     427.449 -488.064 -211.792

#	./bellman --risk --angle 153.435   --rounds 200 --constrain --oracleprob 0.01 --bidderprob 0

#	./bellman --gamma 2   --rounds 400 --constrain --oracleprob 0 --bidderprob 0.01 --write
#	./bellman --gamma .2   --rounds 500 --constrain --oracleprob 0.01 --bidderprob 0 --write
#	./bellman --gamma 100   --rounds 500 --constrain --oracleprob 0.01 --bidderprob 0 --write=
#	./bellman --gamma 2     --rounds 500 --constrain --oracleprob 0.01 --bidderprob 0 --write

# Unconstrained 0 2.5 0.05 7 0.5 1.5 6.5 -0.0691835 0.068553 0.0550946
#   gamma style: g01000 univ1                          2   0.05   7   0.05 10     19.4462 -20.0157 -19.731
#   angle      : g01000 univ(1) 0.447214 -0.894427 153.435 0.05   7   0.05 10     8.50491 -19.82   -19.4188

risk_check: bellman
	./bellman --risk --omega 0.5 --angle 153.434949  --rounds 7  --constrain --oracleprob 0.01 --bidderprob 0 --write

risk_test: bellman
	./bellman --gamma 100  --rounds 100  --constrain --oracleprob 0 --bidderprob 0.05

reject_check: bellman
	./bellman --reject  --angle 0  --rounds 7  --constrain --oracleprob 0 --bidderprob 0.1 --write

# ---  $^ are prereq    $@ is target    $* is stem
#      change n to change path, file names, and the length of run;  gp is path
#      Once run, cat combines these lines to show all of the results.

# define these constants, then use a command like  (use uruns for unconstrained)
#    make -j lots  -k runs/summary.reject_psi0090_n100
# or
#    make -k -j lots  runs/summary.risk_psi0010_n250
#    make -k -j lots uruns/summary.risk_alpha0500_scale2_n200
# with these values chosen to match (don't know how to pick them from make input
# so you have to define the constants here and match them in the make command.
# Builds a directory in runs for these results, then files for each.
n = 200

# define geometric expert by rate
psi = 0.0500
ptxt=   0500

# define uncontrained expert by alpha level
alpha = 0.1000
atxt=     1000

# define expert by uniform n (one more than n)
# psi =   251
# ptxt=   251

# criterion should be risk or reject (and make it so in the C++ code)
goal = risk

# multiplier for unconstrained universal code
scale = 2

#--------------------------------------------------------------------------------------------
#  below here is automagic, building output in runs/   
#--------------------------------------------------------------------------------------------

# -----  unconstrained -----

# define path within uruns subdirectory for each alpha (oracle) and n combination
up = uruns/$(goal)_alpha$(atxt)_scale$(scale)_n$(n)

$(up)/.directory_built: 
	echo Building directory for unconstrained runs $(up)
	mkdir $(up)
	touch $@

# main unconstrained target with parameters that identify angle over tasks
uruns/summary.risk_alpha$(atxt)_scale$(scale)_n$(n): bellman bellman.sh $(up)/0 $(up)/15 $(up)/30 $(up)/45 $(up)/60 $(up)/75 $(up)/90 $(up)/105 $(up)/120 $(up)/135 $(up)/150 $(up)/165 $(up)/180 $(up)/195 $(up)/210 $(up)/225 $(up)/240 $(up)/255 $(up)/270 $(up)/285 $(up)/290 $(up)/295  $(up)/300 $(up)/315 $(up)/330 $(up)/345
	rm -f $@
	cat $(filter $(up)/%,$^) >> $@

# actual run command for unconstrained solution
$(up)/%: bellman bellman.sh $(up)/.directory_built
	./bellman --$(goal) --omega 0.5 --angle $* --oracleprob $(alpha) --bidderprob 0 --scale $(scale)  --rounds $(n) >  $@

# -----  constrained -----

# define path within runs subdirectory for each psi (oracle) and n combination; 0 ids universal
pp = runs/$(goal)_psi$(ptxt)_n$(n)

$(pp)/.directory_built: 
	echo Building directory for contrained runs $(pp)
	mkdir $(pp)
	touch $@

# main constrained target
runs/summary.risk_psi$(ptxt)_n$(n): bellman bellman.sh $(pp)/0 $(pp)/15 $(pp)/30 $(pp)/45 $(pp)/60 $(pp)/75 $(pp)/90 $(pp)/105 $(pp)/120 $(pp)/135 $(pp)/150 $(pp)/165 $(pp)/180 $(pp)/195 $(pp)/210 $(pp)/225 $(pp)/240 $(pp)/255 $(pp)/270 $(pp)/285 $(pp)/300 $(pp)/315 $(pp)/330 $(pp)/345
	rm -f $@
	cat $(filter $(pp)/%,$^) >> $@

# actual run command for constrained solution, with univ and geometric
$(pp)/%: bellman bellman.sh $(pp)/.directory_built
	./bellman --$(goal) --angle $* --constrain --oracleprob $(psi) --bidderprob 0      --rounds $(n) >  $@



###########################################################################
include ../rules_for_makefiles

