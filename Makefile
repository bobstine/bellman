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

# Wealth array scaled_univ(1)[dim=78] with wealth vector beginning W[0]=3.38774 and at zero index W[7]=0.482787
# uncon(0) scaled_univ(1) 153.435 0.5   70   0.05 20     103.028 -2.2802e-07 -115.189
# Dual Wealth Array 'Universal'[dim=37] with wealth W[0]=<3.38774,2.08137> and at zero index W[29]=<0.487736,0.0292789>
# uncon(0) Universal      153.435 0.5   70   0.05 20     104.155 -6.44542e-07 -116.449

risk_check: bellman
	./bellman --risk --omega 0.5 --scale 0.5 --angle 330 --rounds 200 --constrain --oracleprob 0.1 --bidderprob 0 --write

reject_check: bellman
	./bellman --reject           --angle 0  --rounds 7     --constrain --oracleprob 0    --bidderprob 0.1 --write

# define these constants, then use a command like  (use uruns for unconstrained)
#    make -j lots  -k runs/summary.reject_psi0090_n100
# or
#    make -k -j lots  runs/summary.risk_psi0010_n250
#    make -k -j lots uruns/summary.risk_alpha5_omega50_scale2_n200
# with these values chosen to match (don't know how to pick them from make input
# so you have to define the constants here and match them in the make command.
# Builds a directory in runs for these results, then files for each.

n = 500

omega = 0.50
otxt  =   50

# define uncontrained expert by alpha level
alpha = 0.05
atxt=      5

# criterion should be risk or reject (and make it so in the C++ code)
goal = risk

# multiplier for unconstrained universal code
scale = 2
stxt  = 2

# define expert by uniform n (one more than n)
# psi =   251
# ptxt=   251

#--------------------------------------------------------------------------------------------
#  below here is automagic, building output in runs/   
#--------------------------------------------------------------------------------------------

# -----  unconstrained -----

# define path within uruns subdirectory for each alpha (oracle) and n combination
up = uruns/$(goal)_alpha$(atxt)_omega$(otxt)_scale$(stxt)_n$(n)

$(up)/.directory_built: 
	echo "Building directory for unconstrained runs" $(up)
	mkdir $(up)
	touch $@

# main unconstrained target with parameters that identify angle over tasks
uruns/summary.risk_alpha$(atxt)_omega$(otxt)_scale$(stxt)_n$(n): bellman bellman.sh $(up)/0 $(up)/15 $(up)/30 $(up)/45 $(up)/60 $(up)/75 $(up)/90 $(up)/105 $(up)/120 $(up)/135 $(up)/150 $(up)/165 $(up)/180 $(up)/195 $(up)/210 $(up)/225 $(up)/240 $(up)/255 $(up)/270 $(up)/285 $(up)/290 $(up)/295  $(up)/300 $(up)/315 $(up)/330 $(up)/333 $(up)/334 $(up)/335 $(up)/336 $(up)/337 $(up)/345
	rm -f $@
	cat $(filter $(up)/%,$^) >> $@

# actual run command for unconstrained solution
$(up)/%: bellman bellman.sh  $(up)/.directory_built
	./bellman --$(goal) --omega $(omega) --angle $* --oracleprob $(alpha) --bidderprob 0 --scale $(scale)  --rounds $(n) >  $@

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
$(pp)/%: bellman bellman.sh # $(pp)/.directory_built
	./bellman --$(goal) --angle $* --constrain --oracleprob $(psi) --bidderprob 0      --rounds $(n) >  $@



###########################################################################
include ../rules_for_makefiles

