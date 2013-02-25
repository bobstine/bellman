include ../c_flags

###########################################################################
#
#   Note on special forms
#      $^ are prereq    $@ is target    $* is stem
#
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
level_5 = bellman_main.o bellman_calculator.o

############################################################################
#
#            INCLUDING RULES AND DEFINITIONS
#
###########################################################################

# TAGS 
# find . | grep ".*\.\(h\|cc\)" | xargs etags -f TAGS

#--------------------------------------------------------------------


compiler: bellman
	gcc --version



# -------------------------------------------------------------------------------------------------------------
# Find risk for specific mixture of signal and probabilities
# Need to prefix variables like 'n' to avoid conflicts
#
#    Needs to have subdirectory 'risk' created
#    then run <c><c> calc_risk
#
# All symbols in this section have c_ prefix

c_n         = 500
c_scale     = 2     # use small scale for omega 0.05
c_scale_txt = 2
c_omega     = 0.5
c_omega_txt = 50

# define path within risk subdirectory
c_path = risk/omega$(c_omega_txt)_scale$(c_scale_txt)_n$(c_n)

$(c_path)/.directory_built: 
	echo "Building directory for risk output:" $(c_path)
	mkdir $(c_path)
	touch $@


# main command to run
calc_risk:  $(c_path)/.directory_built $(c_path)/m0.5 $(c_path)/m1.0 $(c_path)/m1.5 $(c_path)/m2.0 $(c_path)/m2.5 $(c_path)/m3.0 $(c_path)/m3.5 $(c_path)/m4.0
	echo "Computed risks in directory: " $(c_path)

# actual run command for unconstrained solution needs m to keep unique
$(c_path)/m%: calculate
	echo $(c_path)
	./calculate --signal $* --omega $(c_omega) --scale $(c_scale) --rounds $(c_n) > $@


# executable
calculate: bellman.o wealth.o utility.o distribution.o bellman_calculator.o
	$(GCC) $^ $(LDLIBS) -o  $@


# use these to write means, probs to file for simulation (which is done in simulate_means.R)
sim_details/.directory_built: 
	echo "Building directory for holding simulation details."
	mkdir sim_details
	touch $@

sim_gen: bellman sim_details/.directory_built
	./bellman --risk --omega 0.5 --angle 334 --rounds 500  --oracleprob 0.05 --bidderprob 0 --scale 2 --write   # unconstrained



# -------------------------------------------------------------------------------------------------------------
# bellman recursion for competitive value
#

bellman_main.o: bellman_main.cc

bellman: bellman.o wealth.o utility.o distribution.o bellman_main.o
	$(GCC) $^ $(LDLIBS) -o  $@

#  Bellman options (see code for details)
#       constrain  : oracle has state (2 dim)
#	oracle prob: 0	LS, 1 RI, otherwise testimator
#	bidder prob: 0 for univ, alpha for geometric
#       omega      : 0 for fixed bidder

reject_check: bellman
	./bellman --reject           --angle 0  --rounds 7     --constrain --oracleprob 0    --bidderprob 0.1 --write

risk_check: bellman
	./bellman --risk --omega 0.05 --angle 0 --oracleprob 1 --bidderprob 0.05  --rounds 200

#	./bellman --risk --omega 0   --scale 0.5 --angle 150 --rounds 200             --oracleprob 1   --bidderprob 1 --write   # fixed wealth bidder
#	./bellman --risk --omega 0.5 --scale 0.5 --angle 150 --rounds 200             --oracleprob 1   --bidderprob 0 --write
#	./bellman --risk --omega 0.5 --scale 0.5 --angle 330 --rounds 200 --constrain --oracleprob 0.1 --bidderprob 0 --write

# define the constants n, omega, alpha, beta, and scale below.  Then use the command
#    make -j lots  -k runs/summary.reject_psi0090_n100
# or for unconstrained
#    make -k -j lots uruns/summary.risk_alpha5_beta100_omega50_scale2_n200
# with these values chosen to match (don't know how to pick them from make input
# so you have to define the constants here and match them in the make command.
# Builds a directory and summary in runs or uruns for these results, then files for each.

# define unconstrained expert by alpha level (0 is ls, 1 is risk inf)
alpha =   1
atxt=     100

# define the bidder by beta (0 for universal or rate for geometric; if omega=0, divided by n for bonferroni)
beta = 0
btxt = 0

# omega=0 indicates a fixed wealth bonferroni bidder
omega = 0.75
otxt  =   75

# multiplier for unconstrained universal code
scale = 2
stxt  = 2

# number of tests
n = 200

# criterion should be risk or reject (and make it so in the C++ code)
goal = risk


#--------------------------------------------------------------------------------------------
#  setup for paper figures as simple targets with one command
#--------------------------------------------------------------------------------------------

# Figure 1: comparisons to risk inflation oracle (make directory

f1a = figures/f1/risk_alpha100_beta10_omega05_scale2_n200
f1b = figures/f1/risk_alpha100_beta10_omega25_scale2_n200
f1c = figures/f1/risk_alpha100_beta10_omega50_scale2_n200
f1d = figures/f1/risk_alpha100_beta10_omega75_scale2_n200

figure1: figures/f1/.directory_built figures/f1/summary.$(f1a) figures/f1/summary.$(f1b) figures/f1/summary.$(f1c)
	echo "Making Figure 1..."
	ls figures/f1

figures/f1/.directory_built:
	echo "Making Figure 1..."
	mkdir figures/f1
	touch $@


# -- F1A
$(f1a)/.directory_built: 
	echo "Building directory for Figure 1 runs" $(f1a)
	mkdir $(f1a)
	touch $@

$(f1a)/%: bellman bellman.sh  $(f1a)/.directory_built
	./bellman --risk --oracleprob 1 --bidderprob 0.1 --omega 0.50 --angle $*  --scale 2  --rounds 200 >  $@

figures/f1/summary.$(f1a): bellman bellman.sh $(f1a)/0 $(f1a)/15 $(f1a)/30 $(f1a)/45 $(f1a)/60 $(f1a)/75 $(f1a)/90 $(f1a)/105 $(f1a)/120 $(f1a)/135 $(f1a)/150 $(f1a)/165 $(f1a)/180 $(f1a)/195 $(f1a)/210 $(f1a)/225 $(f1a)/240 $(f1a)/255 $(f1a)/270 $(f1a)/285 $(f1a)/290 $(f1a)/295  $(f1a)/300 $(f1a)/315 $(f1a)/330 $(f1a)/345
	rm -f $@
	cat $(filter $(f1a)/%,$^) >> $@

# -- F1B
$(f1a)/.directory_built: 
	echo "Building directory for Figure 1 runs" $(f1a)
	mkdir $(f1a)
	touch $@

$(f1a)/%: bellman bellman.sh  $(f1a)/.directory_built
	./bellman --risk --oracleprob 1 --bidderprob 0.1 --omega 0.50 --angle $*  --scale 2  --rounds 200 >  $@

figures/f1/summary.$(f1a): bellman bellman.sh $(f1a)/0 $(f1a)/15 $(f1a)/30 $(f1a)/45 $(f1a)/60 $(f1a)/75 $(f1a)/90 $(f1a)/105 $(f1a)/120 $(f1a)/135 $(f1a)/150 $(f1a)/165 $(f1a)/180 $(f1a)/195 $(f1a)/210 $(f1a)/225 $(f1a)/240 $(f1a)/255 $(f1a)/270 $(f1a)/285 $(f1a)/290 $(f1a)/295  $(f1a)/300 $(f1a)/315 $(f1a)/330 $(f1a)/345
	rm -f $@
	cat $(filter $(f1a)/%,$^) >> $@

# -- F1C
$(f1b)/.directory_built: 
	echo "Building directory for Figure 1 runs" $(f1b)
	mkdir $(f1b)
	touch $@

$(f1b)/%: bellman bellman.sh  $(f1b)/.directory_built
	./bellman --risk --oracleprob 1 --bidderprob 0.1 --omega 0.50 --angle $*  --scale 2  --rounds 200 >  $@

figures/f1/summary.$(f1b): bellman bellman.sh $(f1b)/0 $(f1b)/15 $(f1b)/30 $(f1b)/45 $(f1b)/60 $(f1b)/75 $(f1b)/90 $(f1b)/105 $(f1b)/120 $(f1b)/135 $(f1b)/150 $(f1b)/165 $(f1b)/180 $(f1b)/195 $(f1b)/210 $(f1b)/225 $(f1b)/240 $(f1b)/255 $(f1b)/270 $(f1b)/285 $(f1b)/290 $(f1b)/295  $(f1b)/300 $(f1b)/315 $(f1b)/330 $(f1b)/345
	rm -f $@
	cat $(filter $(f1b)/%,$^) >> $@

# -- F1D
$(f1c)/.directory_built: 
	echo "Building directory for Figure 1 runs" $(f1c)
	mkdir $(f1c)
	touch $@

$(f1c)/%: bellman bellman.sh  $(f1c)/.directory_built
	./bellman --risk --oracleprob 1 --bidderprob 0.1 --omega 0.50 --angle $*  --scale 2  --rounds 200 >  $@

figures/f1/summary.$(f1c): bellman bellman.sh $(f1c)/0 $(f1c)/15 $(f1c)/30 $(f1c)/45 $(f1c)/60 $(f1c)/75 $(f1c)/90 $(f1c)/105 $(f1c)/120 $(f1c)/135 $(f1c)/150 $(f1c)/165 $(f1c)/180 $(f1c)/195 $(f1c)/210 $(f1c)/225 $(f1c)/240 $(f1c)/255 $(f1c)/270 $(f1c)/285 $(f1c)/290 $(f1c)/295  $(f1c)/300 $(f1c)/315 $(f1c)/330 $(f1c)/345
	rm -f $@
	cat $(filter $(f1c)/%,$^) >> $@



#--------------------------------------------------------------------------------------------
#  below here is automagic, building output in uruns/  and runs/
#--------------------------------------------------------------------------------------------

# -----  unconstrained -----

# define path within uruns subdirectory for each alpha (oracle) and n combination
up = uruns/$(goal)_alpha$(atxt)_beta$(btxt)_omega$(otxt)_scale$(stxt)_n$(n)

$(up)/.directory_built: 
	echo "Building directory for unconstrained runs" $(up)
	mkdir $(up)
	touch $@

# actual run command for unconstrained solution
$(up)/%: bellman bellman.sh  $(up)/.directory_built
	./bellman --$(goal) --omega $(omega) --angle $* --oracleprob $(alpha) --bidderprob $(beta) --scale $(scale)  --rounds $(n) >  $@

# main unconstrained target with parameters that identify angle over tasks
# coarse spacing  [change name for each version]
uruns/summary.risk_alpha$(atxt)_beta$(btxt)_omega$(otxt)_scale$(stxt)_n$(n): bellman bellman.sh $(up)/0 $(up)/15 $(up)/30 $(up)/45 $(up)/60 $(up)/75 $(up)/90 $(up)/105 $(up)/120 $(up)/135 $(up)/150 $(up)/165 $(up)/180 $(up)/195 $(up)/210 $(up)/225 $(up)/240 $(up)/255 $(up)/270 $(up)/285 $(up)/290 $(up)/295  $(up)/300 $(up)/315 $(up)/330 $(up)/345
	rm -f $@
	cat $(filter $(up)/%,$^) >> $@

# fine spacing  [change name for each version]
fine_uruns/summary.risk_alpha$(atxt)_beta$(btxt)_omega$(otxt)_scale$(stxt)_n$(n): bellman bellman.sh $(up)/0 $(up)/15 $(up)/30 $(up)/45 $(up)/60 $(up)/65 $(up)/70 $(up)/75 $(up)/80 $(up)/85 $(up)/90 $(up)/91 $(up)/92 $(up)/93 $(up)/93.1 $(up)/93.2 $(up)/93.3 $(up)/93.4 $(up)/93.5 $(up)/93.7 $(up)/93.8 $(up)/93.9 $(up)/94 $(up)/94.1 $(up)/94.2 $(up)/94.3 $(up)/94.4 $(up)/94.5 $(up)/94.6 $(up)/94.7 $(up)/94.8 $(up)/94.9 $(up)/95 $(up)/95.1 $(up)/95.2 $(up)/95.3 $(up)/95.4 $(up)/95.5 $(up)/95.6 $(up)/95.7 $(up)/95.8 $(up)/95.9 $(up)/95.95 $(up)/96 $(up)/96.1 $(up)/96.2  $(up)/96.5 $(up)/97 $(up)/97.5 $(up)/98 $(up)/98.5 $(up)/99 $(up)/100 $(up)/100 $(up)/102 $(up)/104 $(up)/106 $(up)/108 $(up)/110 $(up)/115 $(up)/120 $(up)/135 $(up)/150 $(up)/165 $(up)/180 $(up)/195 $(up)/210 $(up)/225 $(up)/240 $(up)/255 $(up)/270 $(up)/285 $(up)/290 $(up)/295  $(up)/300 $(up)/315 $(up)/330 $(up)/345
	rm -f $@
	cat $(filter $(up)/%,$^) >> $@


# -----  constrained -----

# define path within druns subdirectory for each psi (oracle) and n combination; 0 ids universal
pp = druns/$(goal)_omega$(otxt)_scale$(stxt)_psi$(ptxt)_n$(n)

$(pp)/.directory_built: 
	echo Building directory for constrained runs $(pp)
	mkdir $(pp)
	touch $@

# main constrained target
druns/summary.risk_omega$(otxt)_scale$(stxt)_psi$(ptxt)_n$(n): bellman bellman.sh $(pp)/0 $(pp)/15 $(pp)/30 $(pp)/45 $(pp)/60 $(pp)/75 $(pp)/90 $(pp)/105 $(pp)/120 $(pp)/135 $(pp)/150 $(pp)/165 $(pp)/180 $(pp)/195 $(pp)/210 $(pp)/225 $(pp)/240 $(pp)/255 $(pp)/270 $(pp)/285 $(pp)/300 $(pp)/315 $(pp)/330 $(pp)/345
	rm -f $@
	cat $(filter $(pp)/%,$^) >> $@

# actual run command for constrained solution, with univ and geometric
$(pp)/%: bellman bellman.sh  $(pp)/.directory_built
	./bellman --$(goal) --angle $* --scale $(scale) --omega $(omega) --constrain --oracleprob $(psi) --bidderprob 0      --rounds $(n) >  $@



###########################################################################
include ../rules_for_makefiles

