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
	./bellman --reject --angle 0 --rounds 7     --oracle_omega 0.05 --oracle_prob 0  --bidder_omega 0.5  --bidder_prob 0.1 --write

risk_check: bellman
	./bellman --risk --angle 0 --rounds 100 --oracle_omega 0.5  --oracle_prob 0  --bidder_omega 0.5 --bidder_prob 0.10  # constrained

#	./bellman --risk --angle 0 --rounds 100 --oracle_omega 1    --oracle_prob 1  --bidder_omega 0.5 --bidder_prob 0.10  # unconstrained
#	./bellman --risk --angle 0 --rounds 100 --oracle_omega 1    --oracle_prob 1  --bidder_omega 0   --bidder_prob 0.05  # fixed alpha bidder, unconstrained


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
	./bellman --$(goal) --angle $* --oracle_omega 1 --oracle_prob $(alpha) --bidder_omega $(omega) --bidder_prob $(beta) --scale $(scale)  --rounds $(n) >  $@

# main unconstrained target with parameters that identify angle over tasks
# coarse spacing  [change name for each version]
uruns/summary.risk_alpha$(atxt)_beta$(btxt)_omega$(otxt)_scale$(stxt)_n$(n): bellman bellman.sh $(up)/0 $(up)/15 $(up)/30 $(up)/45 $(up)/60 $(up)/75 $(up)/90 $(up)/105 $(up)/120 $(up)/135 $(up)/150 $(up)/165 $(up)/180 $(up)/195 $(up)/210 $(up)/225 $(up)/240 $(up)/255 $(up)/270 $(up)/285 $(up)/290 $(up)/295  $(up)/300 $(up)/315 $(up)/330 $(up)/345
	rm -f $@
	cat $(filter $(up)/%,$^) >> $@

# fine spacing
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



#--------------------------------------------------------------------------------------------
#  setup for paper figures as simple targets with one command (files in figures directory)
#--------------------------------------------------------------------------------------------

base_angles = 0 15 30 45 60 65 70 75 80 85 90 95 100 105 110 115 120 125 135 150 165 180 195 210 225 240 255 270 285 300 315 330 345 359

extra_angles = 91 92 93 94 96 97 98 99

############################################################################################################
#
# Figure 1: comparisons of geometric spending testimator to risk inflation oracle, with varying payouts
#
############################################################################################################

f1_n = 500
f1_angles := $(base_angles) $(extra_angles) 93.2 93.4 93.6 93.8 94.2 94.4 94.6 94.8 95.2 95.4 95.6 95.8 96.2 96.4 96.6 96.8

f1a := risk_alpha100_100_beta10_05_scale2_n$(f1_n)
f1b := risk_alpha100_100_beta10_25_scale2_n$(f1_n)
f1c := risk_alpha100_100_beta10_50_scale2_n$(f1_n)
f1d := risk_alpha100_100_beta10_75_scale2_n$(f1_n)

# --- F1A

f1a_dir := figures/f1/$(f1a)

f1a_obj := $(addprefix $(f1a_dir)/, $(f1_angles))

f1a_sum := figures/f1/summary_$(f1a)

$(f1a_sum): $(f1a_obj)
	rm -f $@
	cat $(f1a_obj) > $@

$(f1a_dir)/%: bellman bellman.sh $(f1a_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 1 --oracle_omega 1 --bidder_prob 0.10 --bidder_omega 0.1 --scale 2  --rounds $(f1_n)  >  $@

$(f1a_dir)/.dir_created :
	mkdir $(f1a_dir)
	touch $@

# --- F1B

f1b_dir := figures/f1/$(f1b)

f1b_obj := $(addprefix $(f1b_dir)/, $(f1_angles))

f1b_sum := figures/f1/summary_$(f1b)

$(f1b_sum): $(f1b_obj)
	rm -f $@
	cat $(f1b_obj) > $@

$(f1b_dir)/%: bellman bellman.sh $(f1b_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 1 --oracle_omega 1 --bidder_prob 0.10 --bidder_omega 0.25 --scale 2  --rounds $(f1_n)  >  $@

$(f1b_dir)/.dir_created :
	mkdir $(f1b_dir)
	touch $@


# --- F1C

f1c_dir := figures/f1/$(f1c)

f1c_obj := $(addprefix $(f1c_dir)/, $(f1_angles))

f1c_sum := figures/f1/summary_$(f1c)

$(f1c_sum): $(f1c_obj)
	rm -f $@
	cat $(f1c_obj) > $@

$(f1c_dir)/%: bellman bellman.sh $(f1c_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 1 --oracle_omega 1 --bidder_prob 0.10 --bidder_omega 0.50  --scale 2  --rounds $(f1_n)  >  $@

$(f1c_dir)/.dir_created :
	mkdir $(f1c_dir)
	touch $@


# --- F1D

f1d_dir := figures/f1/$(f1d)

f1d_obj := $(addprefix $(f1d_dir)/, $(f1_angles))

f1d_sum := figures/f1/summary_$(f1d)

$(f1d_sum): $(f1d_obj)
	rm -f $@
	cat $(f1d_obj) > $@

$(f1d_dir)/%: bellman bellman.sh $(f1d_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 1 --oracle_omega 1 --bidder_prob 0.10 --bidder_omega 0.75 --scale 2  --rounds $(f1_n)  >  $@

$(f1d_dir)/.dir_created:
	mkdir $(f1d_dir)
	touch $@

figure1: $(f1a_sum) $(f1b_sum) $(f1c_sum) $(f1d_sum)


############################################################################################################
#
# Figure 2: comparisons of geometric spending testimator to fixed Bayes oracle (testimator) with fixed alpha
#
############################################################################################################


f2_n = 500
f2_angles := $(base_angles) $(extra_angles)

f2a := risk_alpha10_100_beta10_05_scale2_n$(f2_n)
f2b := risk_alpha10_100_beta10_25_scale2_n$(f2_n)
f2c := risk_alpha10_100_beta10_50_scale2_n$(f2_n)
f2d := risk_alpha10_100_beta10_75_scale2_n$(f2_n)

# --- F2A

f2a_dir := figures/f2/$(f2a)

f2a_obj := $(addprefix $(f2a_dir)/, $(f2_angles))

f2a_sum := figures/f2/summary_$(f2a)

$(f2a_sum): $(f2a_obj)
	rm -f $@
	cat $(f2a_obj) > $@

$(f2a_dir)/%: bellman bellman.sh $(f2a_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0.1 --oracle_omega 1 --bidder_prob 0.1 --bidder_omega 0.05 --scale 2  --rounds $(f2_n)  >  $@

$(f2a_dir)/.dir_created :
	mkdir $(f2a_dir)
	touch $@


# --- F2B

f2b_dir := figures/f2/$(f2b)

f2b_obj := $(addprefix $(f2b_dir)/, $(f2_angles))

f2b_sum := figures/f2/summary_$(f2b)

$(f2b_sum): $(f2b_obj)
	rm -f $@
	cat $(f2b_obj) > $@

$(f2b_dir)/%: bellman bellman.sh $(f2b_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0.1 --oracle_omega 1 --bidder_prob 0.1 --bidder_omega 0.25 --scale 2  --rounds $(f2_n)  >  $@

$(f2b_dir)/.dir_created :
	mkdir $(f2b_dir)
	touch $@

# --- F2C

f2c_dir := figures/f2/$(f2c)

f2c_obj := $(addprefix $(f2c_dir)/, $(f2_angles))

f2c_sum := figures/f2/summary_$(f2c)

$(f2c_sum): $(f2c_obj)
	rm -f $@
	cat $(f2c_obj) > $@

$(f2c_dir)/%: bellman bellman.sh $(f2c_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0.1 --oracle_omega 1 --bidder_prob 0.1 --bidder_omega 0.50 --scale 2  --rounds $(f2_n)  >  $@


$(f2c_dir)/.dir_created :
	mkdir $(f2c_dir)
	touch $@

# --- F2D

f2d_dir := figures/f2/$(f2d)

f2d_obj := $(addprefix $(f2d_dir)/, $(f2_angles))

f2d_sum := figures/f2/summary_$(f2d)

$(f2d_sum): $(f2d_obj)
	rm -f $@
	cat $(f2d_obj) > $@

$(f2d_dir)/%: bellman bellman.sh $(f2d_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0.1 --oracle_omega 1 --bidder_prob 0.1 --bidder_omega 0.75 --scale 2  --rounds $(f2_n)  >  $@

$(f2d_dir)/.dir_created :
	mkdir $(f2d_dir)
	touch $@

figure2: $(f2a_sum) $(f2b_sum) $(f2c_sum) $(f2d_sum)


############################################################################################################
#
# Figure 3: risks of paths within the feasible set for geometric vs testimator oracle
#
############################################################################################################

f3_alpha     = 0.1
f3_atxt      = 10
f3_beta      = 0.1
f3_btxt      = 10
f3_omega     = 0.5
f3_omega_txt = 50
f3_scale     = 2     # use small scale for omega 0.05
f3_scale_txt = 2
f3_n         = 500

# define dir within risk subdirectory
f3_dir = figures/f3/risk_alpha$(f3_atxt)_beta$(f3_btxt)_omega$(f3_omega_txt)_scale$(f3_scale_txt)_n$(f3_n)

$(f3_dir)/.directory_built: 
	echo "Building directory for risk output:" $(f3_dir)
	mkdir $(f3_dir)
	touch $@

# main command to run
f3_mu := 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0

f3_obj := $(addprefix $(f3_dir)/m, $(f3_mu))

figure3:  $(f3_dir)/.directory_built $(f3_obj)
	echo $(f3_obj)
	echo "Computed risks in directory: " $(f3_dir)

$(f3_dir)/m%: calculate
	./calculate --signal $* --alpha $(f3_alpha) --beta $(f3_beta) --omega $(f3_omega) --scale $(f3_scale) --rounds $(f3_n) > $@

# executable
calculate: bellman.o wealth.o utility.o distribution.o bellman_calculator.o
	$(GCC) $^ $(LDLIBS) -o  $@


############################################################################################################
#
# Figure 4: Comparison of geometric bidders with different omega payouts
#
############################################################################################################

f4_n = 200
f4_angles := $(base_angles)

f4a := risk_alpha10_05_beta10_50_scale2_n$(f4_n)
f4b := risk_alpha01_05_beta10_50_scale2_n$(f4_n)
f4c := risk_alpha01_05_beta01_50_scale2_n$(f4_n)
f4d := risk_alpha10_05_beta01_50_scale2_n$(f4_n)

# --- F4A

f4a_dir := figures/f4/$(f4a)

f4a_obj := $(addprefix $(f4a_dir)/, $(f4_angles))

f4a_sum := figures/f4/summary_$(f4a)

$(f4a_sum): $(f4a_obj)
	rm -f $@
	cat $(f4a_obj) > $@

$(f4a_dir)/%: bellman bellman.sh $(f4a_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0.1 --oracle_omega 0.05  --bidder_prob 0.1 --bidder_omega 0.50 --scale 2  --rounds $(f4_n)  >  $@


$(f4a_dir)/.dir_created :
	mkdir $(f4a_dir)
	touch $@

figure4: $(f4a_sum)


# --- F4B

f4b_dir := figures/f4/$(f4b)

f4b_obj := $(addprefix $(f4b_dir)/, $(f4_angles))

f4b_sum := figures/f4/summary_$(f4b)

$(f4b_sum): $(f4b_obj)
	rm -f $@
	cat $(f4b_obj) > $@

$(f4b_dir)/%: bellman bellman.sh $(f4b_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0.01 --oracle_omega 0.05  --bidder_prob 0.1 --bidder_omega 0.50 --scale 2  --rounds $(f4_n)  >  $@


$(f4b_dir)/.dir_created :
	mkdir $(f4b_dir)
	touch $@

figure4: $(f4b_sum)


# --- F4C
f4c_dir := figures/f4/$(f4c)

f4c_obj := $(addprefix $(f4c_dir)/, $(f4_angles))

f4c_sum := figures/f4/summary_$(f4c)

$(f4c_sum): $(f4c_obj)
	rm -f $@
	cat $(f4c_obj) > $@

$(f4c_dir)/%: bellman bellman.sh $(f4c_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0.01 --oracle_omega 0.05  --bidder_prob 0.01 --bidder_omega 0.50 --scale 2  --rounds $(f4_n)  >  $@


$(f4c_dir)/.dir_created :
	mkdir $(f4c_dir)
	touch $@

figure4: $(f4c_sum)


# --- F4D
f4d_dir := figures/f4/$(f4d)

f4d_obj := $(addprefix $(f4d_dir)/, $(f4_angles))

f4d_sum := figures/f4/summary_$(f4d)

$(f4d_sum): $(f4d_obj)
	rm -f $@
	cat $(f4d_obj) > $@

$(f4d_dir)/%: bellman bellman.sh $(f4d_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0.10 --oracle_omega 0.05  --bidder_prob 0.01 --bidder_omega 0.50 --scale 2  --rounds $(f4_n)  >  $@


$(f4d_dir)/.dir_created :
	mkdir $(f4d_dir)
	touch $@

figure4: $(f4a_sum) $(f4b_sum) $(f4c_sum) $(f4d_sum)

###########################################################################
include ../rules_for_makefiles

