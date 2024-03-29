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

OPT = -O3 -std=c++0x -DNDEBUG

# OPT = -O3

USES = utils random

level_1 = distribution.o spending_rule.o
level_2 = wealth.o
level_3 = utility.o
level_4 = bellman.o
level_5 = bellman_main.o bellman_calculator.o bellman_optimize.o

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
	./bellman --risk --angle 296.565 --rounds 200  --oracle_omega 0.25 --oracle_prob 0 --bidder_omega 0.25 --bidder_prob 0.001 --write



# -------------------------------------------------------------------------------------------------------------
# bellman recursion for competitive value
#

bellman_main.o: bellman_main.cc

bellman: bellman.o wealth.o utility.o bellman_main.o spending_rule.o
	$(GCC) $^ $(LDLIBS) -o  $@

bellman_optimize.o: bellman_optimize.cc

optimize: bellman.o wealth.o utility.o spending_rule.o bellman_optimize.o
	$(GCC) $^ $(LDLIBS) -o  $@


#  Bellman options (see code for details)
#       constrain  : oracle has state (2 dim)
#	oracle prob: 0 LS, 1 RI, otherwise testimator
#	bidder prob: 0 for univ, alpha for geometric
#       omega      : 0 for fixed bidder

reject_check: bellman
	./bellman --reject --angle 0       --rounds   7  --oracle_omega 0.05 --oracle_prob 0   --bidder_omega 0.5  --bidder_prob 0.10 --write

risk_check: bellman
	./bellman --risk --angle 296.565 --rounds 100 --oracle_omega 0.25 --oracle_prob 0   --bidder_omega 0.25 --bidder_prob 0.001 --write

risk_inflation: optimize
	./optimize

#  risk check results 14 Aug 2014  (before finer wealth grid)
#	165        250.25   200   12.53072    229.18925   903.76172                                                   
#       296.565    250.25   200    3.0054817   50.305794   21.792614
#                     18 Aug 2014  (finer grid with zero point at 93 and after patched line search, upper limit at 10)
#       165        250.25   200    7.7545943   63.933151  268.5632
#       165        250.25   400   10.038836    71.703323  306.38754
#       296.565    250.25   200    3.0122957   49.897739   21.580969
#       296.565    250.25   400    3.1305635   65.771912   29.385805

#	./bellman --risk --angle 0 --rounds 100 --oracle_omega 0.5  --oracle_prob 0  --bidder_omega 0.5 --bidder_prob 0.10  # constrained
#	./bellman --risk --angle 0 --rounds 100 --oracle_omega 1    --oracle_prob 1  --bidder_omega 0.5 --bidder_prob 0.10  # unconstrained
#	./bellman --risk --angle 0 --rounds 100 --oracle_omega 1    --oracle_prob 1  --bidder_omega 0   --bidder_prob 0.05  # fixed alpha bid,uncon


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



#-----------------------------------------------------------------------------------------------
#  Figures     Figures     Figures     Figures     Figures     Figures     Figures     Figures     
#
#     alpha: x-axis in figures (often an oracle)
#      beta: y-axis
#
#     alpha(i)(p)(o)beta(i)(p)(o)   the values of i,p,o indicate W0, prob, and omega
#
#-----------------------------------------------------------------------------------------------

base_angles = 0 15 30 45 60 75 90 105 120 135 140 145 150 155 160 165 170 175 180 195 210 225 240 255 270 275 280 285 290 295 300 305 310 315 330 345 359 

extra_angles = 171 172 173 173.5 174 174.5   176 177 178 179   306 307 308 309   311 312 313 314 314.5

more_angles = 


############################################################################################################
#
# Figure 0: Universal vs risk inflation oracle
#
############################################################################################################

f0_n = 1000
f0_more_angles :=  314.6 314.7 314.8 314.81 314.82 314.83 314.84 314.85 314.86 314.87 314.873  314.874
f0_even_more = 175.1 175.11 175.12 175.13 175.14 175.15 175.151 175.152 175.153 175.154 175.156 175.157 175.158 175.159
f0_yet_more = 175.2 175.3 175.4 175.41 175.42 175.43 175.44 175.45 175.46 175.47 175.48 175.49 175.5 
f0_angles := $(base_angles) $(extra_angles) $(f0_more_angles)                $(f0_yet_more) 

f0a := risk_alpha100_100_100_beta50_00_50_scale2_n$(f0_n)


# --- F0A

f0a_dir := figures/f0/$(f0a)

f0a_obj := $(addprefix $(f0a_dir)/, $(f0_angles))

f0a_sum := figures/f0/summary_$(f0a)

$(f0a_sum): $(f0a_obj)
	rm -f $@
	cat $(f0a_obj) > $@

$(f0a_dir)/%: bellman bellman.sh $(f0a_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 1 --oracle_omega 1  --bidder_prob 0 --bidder_omega 0.50 --scale 2  --rounds $(f0_n)  >  $@

$(f0a_dir)/.dir_created :
	mkdir $(f0a_dir)
	touch $@

figure0: $(f0a_sum)


test0: bellman Makefile
	rm -f test0
	./bellman --risk --angle 310 --oracle_prob 1 --oracle_omega 1  --bidder_prob 0 --bidder_omega 0.50 --scale 2  --rounds 1000 > $@
	cat test0



############################################################################################################
#
# Figure 0A: risks of paths within the feasible set for universal vs risk inflation oracle
#
############################################################################################################

f0A_alpha     = 1
f0A_atxt      = 100
f0A_beta      = 0
f0A_btxt      = 00
f0A_omega     = 0.5
f0A_omega_txt = 50
f0A_scale     = 2  
f0A_scale_txt = 2
f0A_n         = 1000

# define dir within risk subdirectory
f0A_dir = figures/f0/risk_alpha$(f0A_atxt)_beta$(f0A_btxt)_omega$(f0A_omega_txt)_scale$(f0A_scale_txt)_n$(f0A_n)

$(f0A_dir)/.directory_built: 
	echo "Building directory for risk output:" $(f0A_dir)
	mkdir $(f0A_dir)
	touch $@

# main command to run
f0A_mu := 0.5 0.75 1.0 1.5 2.0 2.5 3.0 3.5 4.0 6.0

f0A_obj := $(addprefix $(f0A_dir)/m, $(f0A_mu))

figure0A:  $(f0A_dir)/.directory_built $(f0A_obj)
	echo $(f0A_obj)
	echo "Computed risks in directory: " $(f0A_dir)

$(f0A_dir)/m%: calculate
	./calculate --signal $* --alpha $(f0A_alpha) --beta $(f0A_beta) --omega $(f0A_omega) --scale $(f0A_scale) --rounds $(f0A_n) > $@

# executable
calculate: bellman.o wealth.o utility.o spending_rule.o bellman_calculator.o
	$(GCC) $^ $(LDLIBS) -o  $@


##############################################################################################################################
#
# Figure 11: comparisons of universal spending testimator to risk inflation oracle, with varying payouts
#
#############################################################################################################################

f11_n = 1000

# for RI oracle
f11_alpha     = 1
f11_alpha_txt = 100

# for Testimator oracle
# f11_alpha     = 0.10
# f11_alpha_txt =   10

f11_beta      = 0
f11_beta_txt  = 00
f11_angles_05 =  177.1 177.2 177.25 177.275 177.276 177.277 177.278 177.279 177.27967 177.27968  177.280 
f11_angles_25 =  174.6 174.65 174.67 174.68 174.69 174.7    314.6 314.7 314.8 314.85 314.9 
f11_angles := $(base_angles) $(extra_angles) $(f11_angles_05) $(f11_angles_25) 175.1 175.15 175.155 175.158 175.2 175.5 175.7 176.2 176.5

f11a := risk_alpha100_$(f11_alpha_txt)_100_beta05_$(f11_beta_txt)_05_scale2_n$(f11_n)
f11b := risk_alpha100_$(f11_alpha_txt)_100_beta25_$(f11_beta_txt)_25_scale2_n$(f11_n)
f11c := risk_alpha100_$(f11_alpha_txt)_100_beta50_$(f11_beta_txt)_50_scale2_n$(f11_n)

# --- F11A

f11a_dir := figures/f11/$(f11a)

f11a_obj := $(addprefix $(f11a_dir)/, $(f11_angles))

f11a_sum := figures/f11/summary_$(f11a)

$(f11a_sum): $(f11a_obj)
	rm -f $@
	cat $(f11a_obj) > $@

$(f11a_dir)/%: bellman bellman.sh $(f11a_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f11_alpha) --oracle_omega 1 --bidder_prob $(f11_beta) --bidder_omega 0.05 --scale 2  --rounds $(f11_n)  >$@

$(f11a_dir)/.dir_created :
	mkdir $(f11a_dir)
	touch $@

# --- F11B

f11b_dir := figures/f11/$(f11b)

f11b_obj := $(addprefix $(f11b_dir)/, $(f11_angles))

f11b_sum := figures/f11/summary_$(f11b)

$(f11b_sum): $(f11b_obj)
	rm -f $@
	cat $(f11b_obj) > $@

$(f11b_dir)/%: bellman bellman.sh $(f11b_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f11_alpha) --oracle_omega 1 --bidder_prob $(f11_beta) --bidder_omega 0.25 --scale 2  --rounds $(f11_n)  >$@

$(f11b_dir)/.dir_created :
	mkdir $(f11b_dir)
	touch $@

# --- F11C

f11c_dir := figures/f11/$(f11c)

f11c_obj := $(addprefix $(f11c_dir)/, $(f11_angles))

f11c_sum := figures/f11/summary_$(f11c)

$(f11c_sum): $(f11c_obj)
	rm -f $@
	cat $(f11c_obj) > $@

$(f11c_dir)/%: bellman bellman.sh $(f11c_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f11_alpha) --oracle_omega 1 --bidder_prob $(f11_beta) --bidder_omega 0.50 --scale 2  --rounds $(f11_n)  >$@

$(f11c_dir)/.dir_created :
	mkdir $(f11c_dir)
	touch $@

figure11: $(f11a_sum) $(f11b_sum) $(f11c_sum)

##############################################################################################################################
#
# Figure 1: comparisons of geometric spending testimator to risk inflation oracle, with varying rates and payouts
#
#############################################################################################################################

f1_n = 500

# for RI oracle
f1_alpha     = 1
f1_alpha_txt = 100

# for Testimator oracle
# f1_alpha     = 0.10
# f1_alpha_txt =   10

f1_beta      = 0.10
f1_beta_txt  =   10

f1_angles := $(base_angles)  $(extra_angles)  0.05 1 2 4 6 8 10 12

f1a := risk_alpha100_$(f1_alpha_txt)_100_beta25_001_25_scale2_n$(f1_n)
f1aa:= risk_alpha100_$(f1_alpha_txt)_100_beta25_005_25_scale2_n$(f1_n)
f1b := risk_alpha100_$(f1_alpha_txt)_100_beta25_01_25_scale2_n$(f1_n)
f1c := risk_alpha100_$(f1_alpha_txt)_100_beta25_05_25_scale2_n$(f1_n)

# f1a := risk_alpha100_$(f1_alpha_txt)_100_beta05_$(f1_beta_txt)_05_scale2_n$(f1_n)
# f1aa:= risk_alpha100_$(f1_alpha_txt)_100_beta10_$(f1_beta_txt)_10_scale2_n$(f1_n)
# f1b := risk_alpha100_$(f1_alpha_txt)_100_beta25_$(f1_beta_txt)_25_scale2_n$(f1_n)
# f1c := risk_alpha100_$(f1_alpha_txt)_100_beta50_$(f1_beta_txt)_50_scale2_n$(f1_n)
f1d := risk_alpha100_$(f1_alpha_txt)_100_beta75_$(f1_beta_txt)_75_scale2_n$(f1_n)
# fix W0 = 0.05
f1ea:= risk_alpha100_$(f1_alpha_txt)_100_beta05_$(f1_beta_txt)_10_scale2_n$(f1_n)
f1e := risk_alpha100_$(f1_alpha_txt)_100_beta05_$(f1_beta_txt)_25_scale2_n$(f1_n)
f1f := risk_alpha100_$(f1_alpha_txt)_100_beta05_$(f1_beta_txt)_50_scale2_n$(f1_n)
f1g := risk_alpha100_$(f1_alpha_txt)_100_beta05_$(f1_beta_txt)_75_scale2_n$(f1_n)
# fix W0 = 0.50
f1h := risk_alpha100_$(f1_alpha_txt)_100_beta50_$(f1_beta_txt)_05_scale2_n$(f1_n)
f1ha:= risk_alpha100_$(f1_alpha_txt)_100_beta50_$(f1_beta_txt)_10_scale2_n$(f1_n)
f1i := risk_alpha100_$(f1_alpha_txt)_100_beta50_$(f1_beta_txt)_25_scale2_n$(f1_n)
f1j := risk_alpha100_$(f1_alpha_txt)_100_beta50_$(f1_beta_txt)_75_scale2_n$(f1_n)


# --- F1A

f1a_dir := figures/f1/$(f1a)

f1a_obj := $(addprefix $(f1a_dir)/, $(f1_angles))

f1a_sum := figures/f1/summary_$(f1a)

$(f1a_sum): $(f1a_obj)
	rm -f $@
	cat $(f1a_obj) > $@

$(f1a_dir)/%: bellman bellman.sh $(f1a_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f1_alpha) --oracle_omega 1 --bidder_prob .001 --bidder_omega 0.25 --scale 2  --rounds $(f1_n)  >$@

$(f1a_dir)/.dir_created :
	mkdir $(f1a_dir)
	touch $@


# --- F1AA

f1aa_dir := figures/f1/$(f1aa)

f1aa_obj := $(addprefix $(f1aa_dir)/, $(f1_angles))

f1aa_sum := figures/f1/summary_$(f1aa)

$(f1aa_sum): $(f1aa_obj)
	rm -f $@
	cat $(f1aa_obj) > $@

$(f1aa_dir)/%: bellman bellman.sh $(f1aa_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f1_alpha) --oracle_omega 1 --bidder_prob .005 --bidder_omega 0.25 --scale 2  --rounds $(f1_n)  >$@

$(f1aa_dir)/.dir_created :
	mkdir $(f1aa_dir)
	touch $@

# --- F1B

f1b_dir := figures/f1/$(f1b)

f1b_obj := $(addprefix $(f1b_dir)/, $(f1_angles))

f1b_sum := figures/f1/summary_$(f1b)

$(f1b_sum): $(f1b_obj)
	rm -f $@
	cat $(f1b_obj) > $@

$(f1b_dir)/%: bellman bellman.sh $(f1b_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f1_alpha) --oracle_omega 1 --bidder_prob 0.01 --bidder_omega 0.25 --scale 2  --rounds $(f1_n)  >$@

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
	./bellman --risk --angle $* --oracle_prob $(f1_alpha) --oracle_omega 1 --bidder_prob 0.05 --bidder_omega 0.25  --scale 2  --rounds $(f1_n)  >$@

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
	./bellman --risk --angle $* --oracle_prob $(f1_alpha) --oracle_omega 1 --bidder_prob $(f1_beta) --bidder_omega 0.75 --scale 2  --rounds $(f1_n)  >$@

$(f1d_dir)/.dir_created:
	mkdir $(f1d_dir)
	touch $@

# --- F1EA

f1ea_dir := figures/f1/$(f1ea)

f1ea_obj := $(addprefix $(f1ea_dir)/, $(f1_angles))

f1ea_sum := figures/f1/summary_$(f1ea)

$(f1ea_sum): $(f1ea_obj)
	rm -f $@
	cat $(f1ea_obj) > $@

$(f1ea_dir)/%: bellman bellman.sh $(f1ea_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f1_alpha) --oracle_omega 1 --bidder_w0 0.05 --bidder_prob $(f1_beta) --bidder_omega 0.10 --scale 2  --rounds $(f1_n)  >$@

$(f1ea_dir)/.dir_created:
	mkdir $(f1ea_dir)
	touch $@

# --- F1E

f1e_dir := figures/f1/$(f1e)

f1e_obj := $(addprefix $(f1e_dir)/, $(f1_angles))

f1e_sum := figures/f1/summary_$(f1e)

$(f1e_sum): $(f1e_obj)
	rm -f $@
	cat $(f1e_obj) > $@

$(f1e_dir)/%: bellman bellman.sh $(f1e_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f1_alpha) --oracle_omega 1 --bidder_w0 0.05 --bidder_prob $(f1_beta) --bidder_omega 0.25 --scale 2  --rounds $(f1_n)  >$@

$(f1e_dir)/.dir_created:
	mkdir $(f1e_dir)
	touch $@

# --- F1F

f1f_dir := figures/f1/$(f1f)

f1f_obj := $(addprefix $(f1f_dir)/, $(f1_angles))

f1f_sum := figures/f1/summary_$(f1f)

$(f1f_sum): $(f1f_obj)
	rm -f $@
	cat $(f1f_obj) > $@

$(f1f_dir)/%: bellman bellman.sh $(f1f_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f1_alpha) --oracle_omega 1 --bidder_w0 0.05 --bidder_prob $(f1_beta) --bidder_omega 0.50 --scale 2  --rounds $(f1_n)  >$@

$(f1f_dir)/.dir_created:
	mkdir $(f1f_dir)
	touch $@

# --- F1G

f1g_dir := figures/f1/$(f1g)

f1g_obj := $(addprefix $(f1g_dir)/, $(f1_angles))

f1g_sum := figures/f1/summary_$(f1g)

$(f1g_sum): $(f1g_obj)
	rm -f $@
	cat $(f1g_obj) > $@

$(f1g_dir)/%: bellman bellman.sh $(f1g_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f1_alpha) --oracle_omega 1 --bidder_w0 0.05 --bidder_prob $(f1_beta) --bidder_omega 0.75 --scale 2  --rounds $(f1_n)  >$@

$(f1g_dir)/.dir_created:
	mkdir $(f1g_dir)
	touch $@

# --- F1H

f1h_dir := figures/f1/$(f1h)

f1h_obj := $(addprefix $(f1h_dir)/, $(f1_angles))

f1h_sum := figures/f1/summary_$(f1h)

$(f1h_sum): $(f1h_obj)
	rm -f $@
	cat $(f1h_obj) > $@

$(f1h_dir)/%: bellman bellman.sh $(f1h_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f1_alpha) --oracle_omega 1 --bidder_w0 0.50 --bidder_prob $(f1_beta) --bidder_omega 0.05 --scale 2  --rounds $(f1_n)  >$@

$(f1h_dir)/.dir_created:
	mkdir $(f1h_dir)
	touch $@

# --- F1HA

f1ha_dir := figures/f1/$(f1ha)

f1ha_obj := $(addprefix $(f1ha_dir)/, $(f1_angles))

f1ha_sum := figures/f1/summary_$(f1ha)

$(f1ha_sum): $(f1ha_obj)
	rm -f $@
	cat $(f1ha_obj) > $@

$(f1ha_dir)/%: bellman bellman.sh $(f1ha_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f1_alpha) --oracle_omega 1 --bidder_w0 0.50 --bidder_prob $(f1_beta) --bidder_omega 0.10 --scale 2  --rounds $(f1_n)  >$@

$(f1ha_dir)/.dir_created:
	mkdir $(f1ha_dir)
	touch $@

# --- F1I

f1i_dir := figures/f1/$(f1i)

f1i_obj := $(addprefix $(f1i_dir)/, $(f1_angles))

f1i_sum := figures/f1/summary_$(f1i)

$(f1i_sum): $(f1i_obj)
	rm -f $@
	cat $(f1i_obj) > $@

$(f1i_dir)/%: bellman bellman.sh $(f1i_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f1_alpha) --oracle_omega 1 --bidder_w0 0.50 --bidder_prob $(f1_beta) --bidder_omega 0.25 --scale 2  --rounds $(f1_n)  >$@

$(f1i_dir)/.dir_created:
	mkdir $(f1i_dir)
	touch $@

# --- F1J

f1j_dir := figures/f1/$(f1j)

f1j_obj := $(addprefix $(f1j_dir)/, $(f1_angles))

f1j_sum := figures/f1/summary_$(f1j)

$(f1j_sum): $(f1j_obj)
	rm -f $@
	cat $(f1j_obj) > $@

$(f1j_dir)/%: bellman bellman.sh $(f1j_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f1_alpha) --oracle_omega 1 --bidder_w0 0.50 --bidder_prob $(f1_beta) --bidder_omega 0.75 --scale 2  --rounds $(f1_n)  >$@

$(f1j_dir)/.dir_created:
	mkdir $(f1j_dir)
	touch $@


#                               w0 = omega                                           w0 = 0.05                                         w0=0.50
# figure1: $(f1a_sum) $(f1aa_sum) $(f1b_sum) $(f1c_sum) $(f1d_sum)    $(f1ea_sum) $(f1e_sum) $(f1f_sum) $(f1g_sum)    $(f1h_sum) $(f1ha_sum) $(f1i_sum) $(f1j_sum)
figure1: $(f1a_sum) $(f1aa_sum) $(f1b_sum) $(f1c_sum)



##############################################################################################################################
#
# Figure 2: comparisons of geometric spending testimator to testimator oracle, with varying payouts
#            (same as figure 1, just different oracle)
#
#############################################################################################################################

f2_n = 1000

# for RI oracle
# f2_alpha     = 1
# f2_alpha_txt = 100

# for Testimator oracle
f2_alpha     = 0.10
f2_alpha_txt =   10

f2_beta      = 0.10
f2_beta_txt  =   10

f2_angles := $(base_angles)  $(extra_angles) 

f2a := risk_alpha100_$(f2_alpha_txt)_100_beta05_$(f2_beta_txt)_05_scale2_n$(f2_n)
f2aa:= risk_alpha100_$(f2_alpha_txt)_100_beta10_$(f2_beta_txt)_10_scale2_n$(f2_n)
f2b := risk_alpha100_$(f2_alpha_txt)_100_beta25_$(f2_beta_txt)_25_scale2_n$(f2_n)
f2c := risk_alpha100_$(f2_alpha_txt)_100_beta50_$(f2_beta_txt)_50_scale2_n$(f2_n)
f2d := risk_alpha100_$(f2_alpha_txt)_100_beta75_$(f2_beta_txt)_75_scale2_n$(f2_n)
# fix W0 = 0.05
f2ea:= risk_alpha100_$(f2_alpha_txt)_100_beta05_$(f2_beta_txt)_10_scale2_n$(f2_n)
f2e := risk_alpha100_$(f2_alpha_txt)_100_beta05_$(f2_beta_txt)_25_scale2_n$(f2_n)
f2f := risk_alpha100_$(f2_alpha_txt)_100_beta05_$(f2_beta_txt)_50_scale2_n$(f2_n)
f2g := risk_alpha100_$(f2_alpha_txt)_100_beta05_$(f2_beta_txt)_75_scale2_n$(f2_n)
# fix W0 = 0.50
f2h := risk_alpha100_$(f2_alpha_txt)_100_beta50_$(f2_beta_txt)_05_scale2_n$(f2_n)
f2ha:= risk_alpha100_$(f2_alpha_txt)_100_beta50_$(f2_beta_txt)_10_scale2_n$(f2_n)
f2i := risk_alpha100_$(f2_alpha_txt)_100_beta50_$(f2_beta_txt)_25_scale2_n$(f2_n)
f2j := risk_alpha100_$(f2_alpha_txt)_100_beta50_$(f2_beta_txt)_75_scale2_n$(f2_n)


# --- F2A

f2a_dir := figures/f2/$(f2a)

f2a_obj := $(addprefix $(f2a_dir)/, $(f2_angles))

f2a_sum := figures/f2/summary_$(f2a)

$(f2a_sum): $(f2a_obj)
	rm -f $@
	cat $(f2a_obj) > $@

$(f2a_dir)/%: bellman bellman.sh $(f2a_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f2_alpha) --oracle_omega 1 --bidder_prob $(f2_beta) --bidder_omega 0.05 --scale 2  --rounds $(f2_n)  >$@

$(f2a_dir)/.dir_created :
	mkdir $(f2a_dir)
	touch $@


# --- F2AA

f2aa_dir := figures/f2/$(f2aa)

f2aa_obj := $(addprefix $(f2aa_dir)/, $(f2_angles))

f2aa_sum := figures/f2/summary_$(f2aa)

$(f2aa_sum): $(f2aa_obj)
	rm -f $@
	cat $(f2aa_obj) > $@

$(f2aa_dir)/%: bellman bellman.sh $(f2aa_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f2_alpha) --oracle_omega 1 --bidder_prob $(f2_beta) --bidder_omega 0.10 --scale 2  --rounds $(f2_n)  >$@

$(f2aa_dir)/.dir_created :
	mkdir $(f2aa_dir)
	touch $@

# --- F2B

f2b_dir := figures/f2/$(f2b)

f2b_obj := $(addprefix $(f2b_dir)/, $(f2_angles))

f2b_sum := figures/f2/summary_$(f2b)

$(f2b_sum): $(f2b_obj)
	rm -f $@
	cat $(f2b_obj) > $@

$(f2b_dir)/%: bellman bellman.sh $(f2b_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f2_alpha) --oracle_omega 1 --bidder_prob $(f2_beta) --bidder_omega 0.25 --scale 2  --rounds $(f2_n)  >$@

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
	./bellman --risk --angle $* --oracle_prob $(f2_alpha) --oracle_omega 1 --bidder_prob $(f2_beta) --bidder_omega 0.50  --scale 2  --rounds $(f2_n)  >$@

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
	./bellman --risk --angle $* --oracle_prob $(f2_alpha) --oracle_omega 1 --bidder_prob $(f2_beta) --bidder_omega 0.75 --scale 2  --rounds $(f2_n)  >$@

$(f2d_dir)/.dir_created:
	mkdir $(f2d_dir)
	touch $@

# --- F2EA

f2ea_dir := figures/f2/$(f2ea)

f2ea_obj := $(addprefix $(f2ea_dir)/, $(f2_angles))

f2ea_sum := figures/f2/summary_$(f2ea)

$(f2ea_sum): $(f2ea_obj)
	rm -f $@
	cat $(f2ea_obj) > $@

$(f2ea_dir)/%: bellman bellman.sh $(f2ea_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f2_alpha) --oracle_omega 1 --bidder_w0 0.05 --bidder_prob $(f2_beta) --bidder_omega 0.10 --scale 2  --rounds $(f2_n)  >$@

$(f2ea_dir)/.dir_created:
	mkdir $(f2ea_dir)
	touch $@

# --- F2E

f2e_dir := figures/f2/$(f2e)

f2e_obj := $(addprefix $(f2e_dir)/, $(f2_angles))

f2e_sum := figures/f2/summary_$(f2e)

$(f2e_sum): $(f2e_obj)
	rm -f $@
	cat $(f2e_obj) > $@

$(f2e_dir)/%: bellman bellman.sh $(f2e_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f2_alpha) --oracle_omega 1 --bidder_w0 0.05 --bidder_prob $(f2_beta) --bidder_omega 0.25 --scale 2  --rounds $(f2_n)  >$@

$(f2e_dir)/.dir_created:
	mkdir $(f2e_dir)
	touch $@

# --- F2F

f2f_dir := figures/f2/$(f2f)

f2f_obj := $(addprefix $(f2f_dir)/, $(f2_angles))

f2f_sum := figures/f2/summary_$(f2f)

$(f2f_sum): $(f2f_obj)
	rm -f $@
	cat $(f2f_obj) > $@

$(f2f_dir)/%: bellman bellman.sh $(f2f_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f2_alpha) --oracle_omega 1 --bidder_w0 0.05 --bidder_prob $(f2_beta) --bidder_omega 0.50 --scale 2  --rounds $(f2_n)  >$@

$(f2f_dir)/.dir_created:
	mkdir $(f2f_dir)
	touch $@

# --- F2G

f2g_dir := figures/f2/$(f2g)

f2g_obj := $(addprefix $(f2g_dir)/, $(f2_angles))

f2g_sum := figures/f2/summary_$(f2g)

$(f2g_sum): $(f2g_obj)
	rm -f $@
	cat $(f2g_obj) > $@

$(f2g_dir)/%: bellman bellman.sh $(f2g_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f2_alpha) --oracle_omega 1 --bidder_w0 0.05 --bidder_prob $(f2_beta) --bidder_omega 0.75 --scale 2  --rounds $(f2_n)  >$@

$(f2g_dir)/.dir_created:
	mkdir $(f2g_dir)
	touch $@

# --- F2H

f2h_dir := figures/f2/$(f2h)

f2h_obj := $(addprefix $(f2h_dir)/, $(f2_angles))

f2h_sum := figures/f2/summary_$(f2h)

$(f2h_sum): $(f2h_obj)
	rm -f $@
	cat $(f2h_obj) > $@

$(f2h_dir)/%: bellman bellman.sh $(f2h_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f2_alpha) --oracle_omega 1 --bidder_w0 0.50 --bidder_prob $(f2_beta) --bidder_omega 0.05 --scale 2  --rounds $(f2_n)  >$@

$(f2h_dir)/.dir_created:
	mkdir $(f2h_dir)
	touch $@

# --- F2HA

f2ha_dir := figures/f2/$(f2ha)

f2ha_obj := $(addprefix $(f2ha_dir)/, $(f2_angles))

f2ha_sum := figures/f2/summary_$(f2ha)

$(f2ha_sum): $(f2ha_obj)
	rm -f $@
	cat $(f2ha_obj) > $@

$(f2ha_dir)/%: bellman bellman.sh $(f2ha_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f2_alpha) --oracle_omega 1 --bidder_w0 0.50 --bidder_prob $(f2_beta) --bidder_omega 0.10 --scale 2  --rounds $(f2_n)  >$@

$(f2ha_dir)/.dir_created:
	mkdir $(f2ha_dir)
	touch $@

# --- F2I

f2i_dir := figures/f2/$(f2i)

f2i_obj := $(addprefix $(f2i_dir)/, $(f2_angles))

f2i_sum := figures/f2/summary_$(f2i)

$(f2i_sum): $(f2i_obj)
	rm -f $@
	cat $(f2i_obj) > $@

$(f2i_dir)/%: bellman bellman.sh $(f2i_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f2_alpha) --oracle_omega 1 --bidder_w0 0.50 --bidder_prob $(f2_beta) --bidder_omega 0.25 --scale 2  --rounds $(f2_n)  >$@

$(f2i_dir)/.dir_created:
	mkdir $(f2i_dir)
	touch $@

# --- F2J

f2j_dir := figures/f2/$(f2j)

f2j_obj := $(addprefix $(f2j_dir)/, $(f2_angles))

f2j_sum := figures/f2/summary_$(f2j)

$(f2j_sum): $(f2j_obj)
	rm -f $@
	cat $(f2j_obj) > $@

$(f2j_dir)/%: bellman bellman.sh $(f2j_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob $(f2_alpha) --oracle_omega 1 --bidder_w0 0.50 --bidder_prob $(f2_beta) --bidder_omega 0.75 --scale 2  --rounds $(f2_n)  >$@

$(f2j_dir)/.dir_created:
	mkdir $(f2j_dir)
	touch $@


#                               w0 = omega                                           w0 = 0.05                                         w0=0.50
figure2: $(f2a_sum) $(f2aa_sum) $(f2b_sum) $(f2c_sum) $(f2d_sum)    $(f2ea_sum) $(f2e_sum) $(f2f_sum) $(f2g_sum)    $(f2h_sum) $(f2ha_sum) $(f2i_sum) $(f2j_sum)




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

$(f3_dir)/m%: calculateGeo
	./calculate --signal $* --alpha $(f3_alpha) --beta $(f3_beta) --omega $(f3_omega) --scale $(f3_scale) --rounds $(f3_n) > $@

# executable
calculateGeo: bellman.o wealth.o utility.o distribution.o bellman_calculator.o
	$(GCC) $^ $(LDLIBS) -o  $@


############################################################################################################
#
# Figure 4: Comparison of geometric bidders with different omega payouts
#
############################################################################################################

f4_n = 20
f4_angles := $(base_angles)

f4a := risk_alpha05_10_05_beta50_10_50_scale2_n$(f4_n)
f4b := risk_alpha05_01_05_beta50_10_50_scale2_n$(f4_n)
f4c := risk_alpha05_10_05_beta50_01_50_scale2_n$(f4_n)
f4d := risk_alpha05_01_05_beta50_01_50_scale2_n$(f4_n)

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


############################################################################################################
#
# Figure 5: Universal vs geometric bidder  (two separate batches, first batch for paper)
#
############################################################################################################

f5_n = 1000

f5_xtra_angles = 138   141 142 143 144  146 147 148 149  151 152 153 154  157   281 282 283 284 286 287 296 297

f5_angles := $(base_angles) $(extra_angles) $(f5_xtra_angles)

f5a := risk_alpha25_001_25_beta25_00_25_scale2_n$(f5_n)
f5b := risk_alpha25_002_25_beta25_00_25_scale2_n$(f5_n)
f5c := risk_alpha25_005_25_beta25_00_25_scale2_n$(f5_n)

# --- F5A

f5a_dir := figures/f5/$(f5a)

f5a_obj := $(addprefix $(f5a_dir)/, $(f5_angles))

f5a_sum := figures/f5/summary_$(f5a)

$(f5a_sum): $(f5a_obj)
	rm -f $@
	cat $(f5a_obj) > $@

$(f5a_dir)/%: bellman bellman.sh $(f5a_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0.001 --oracle_omega 0.25  --bidder_prob 0 --bidder_omega 0.25 --scale 2  --rounds $(f5_n)  >  $@

$(f5a_dir)/.dir_created:
	mkdir $(f5a_dir)
	touch $@

# --- F5B

f5b_dir := figures/f5/$(f5b)

f5b_obj := $(addprefix $(f5b_dir)/, $(f5_angles))

f5b_sum := figures/f5/summary_$(f5b)

$(f5b_sum): $(f5b_obj)
	rm -f $@
	cat $(f5b_obj) > $@

$(f5b_dir)/%: bellman bellman.sh $(f5b_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0.002 --oracle_omega 0.25  --bidder_prob 0 --bidder_omega 0.25 --scale 2  --rounds $(f5_n)  >  $@

$(f5b_dir)/.dir_created:
	mkdir $(f5b_dir)
	touch $@

# --- F5C
f5c_dir := figures/f5/$(f5c)

f5c_obj := $(addprefix $(f5c_dir)/, $(f5_angles))

f5c_sum := figures/f5/summary_$(f5c)

$(f5c_sum): $(f5c_obj)
	rm -f $@
	cat $(f5c_obj) > $@

$(f5c_dir)/%: bellman bellman.sh $(f5c_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0.005 --oracle_omega 0.25  --bidder_prob 0 --bidder_omega 0.25 --scale 2  --rounds $(f5_n)  >  $@

$(f5c_dir)/.dir_created:
	mkdir $(f5c_dir)
	touch $@

figure5: $(f5a_sum) $(f5b_sum) $(f5c_sum) 


# --- F51 series flips axes to place universal on x with geo on y, and explores wider range of geometrics

f51_n = 200

f51_xtra_angles = 291.801 296.565 # 138   141 142 143 144  146 147 148 149  151 152 153 154  157   281 282 283 284 286 287 296 297

f51_angles := $(base_angles) $(extra_angles) $(f51_xtra_angles)

f51a := risk_alpha25_00_25_beta25_001_25_scale2_n$(f51_n)
f51b := risk_alpha25_00_25_beta25_01_25_scale2_n$(f51_n)
f51c := risk_alpha25_00_25_beta25_10_25_scale2_n$(f51_n)
f51d := risk_alpha25_00_25_beta25_25_25_scale2_n$(f51_n)
f51e := risk_alpha25_00_25_beta25_35_25_scale2_n$(f51_n)
f51f := risk_alpha25_00_25_beta25_0001_25_scale2_n$(f51_n)

# --- F51A

f51a_dir := figures/f51/$(f51a)

f51a_obj := $(addprefix $(f51a_dir)/, $(f51_angles))

f51a_sum := figures/f51/summary_$(f51a)

$(f51a_sum): $(f51a_obj)
	rm -f $@
	cat $(f51a_obj) > $@

$(f51a_dir)/%: bellman bellman.sh $(f51a_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25  --bidder_prob 0.001 --bidder_omega 0.25 --scale 2  --rounds $(f51_n)  >  $@

$(f51a_dir)/.dir_created :
	mkdir $(f51a_dir)
	touch $@

# --- F51B

f51b_dir := figures/f51/$(f51b)

f51b_obj := $(addprefix $(f51b_dir)/, $(f51_angles))

f51b_sum := figures/f51/summary_$(f51b)

$(f51b_sum): $(f51b_obj)
	rm -f $@
	cat $(f51b_obj) > $@

$(f51b_dir)/%: bellman bellman.sh $(f51b_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25  --bidder_prob 0.01 --bidder_omega 0.25 --scale 2  --rounds $(f51_n)  >  $@

$(f51b_dir)/.dir_created :
	mkdir $(f51b_dir)
	touch $@

# --- F51C

f51c_dir := figures/f51/$(f51c)

f51c_obj := $(addprefix $(f51c_dir)/, $(f51_angles))

f51c_sum := figures/f51/summary_$(f51c)

$(f51c_sum): $(f51c_obj)
	rm -f $@
	cat $(f51c_obj) > $@

$(f51c_dir)/%: bellman bellman.sh $(f51c_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25  --bidder_prob 0.10 --bidder_omega 0.25 --scale 2  --rounds $(f51_n)  >  $@

$(f51c_dir)/.dir_created:
	mkdir $(f51c_dir)
	touch $@

# --- F51D

f51d_dir := figures/f51/$(f51d)

f51d_obj := $(addprefix $(f51d_dir)/, $(f51_angles))

f51d_sum := figures/f51/summary_$(f51d)

$(f51d_sum): $(f51d_obj)
	rm -f $@
	cat $(f51d_obj) > $@

$(f51d_dir)/%: bellman bellman.sh $(f51d_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25  --bidder_prob 0.25 --bidder_omega 0.25 --scale 2  --rounds $(f51_n)  >  $@

$(f51d_dir)/.dir_created :
	mkdir $(f51d_dir)
	touch $@

# --- F51E

f51e_dir := figures/f51/$(f51e)

f51e_obj := $(addprefix $(f51e_dir)/, $(f51_angles))

f51e_sum := figures/f51/summary_$(f51e)

$(f51e_sum): $(f51e_obj)
	rm -f $@
	cat $(f51e_obj) > $@

$(f51e_dir)/%: bellman bellman.sh $(f51e_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25  --bidder_prob 0.35 --bidder_omega 0.25 --scale 2  --rounds $(f51_n)  >  $@

$(f51e_dir)/.dir_created :
	mkdir $(f51e_dir)
	touch $@

# --- F51F

f51f_dir := figures/f51/$(f51f)

f51f_obj := $(addprefix $(f51f_dir)/, $(f51_angles))

f51f_sum := figures/f51/summary_$(f51f)

$(f51f_sum): $(f51f_obj)
	rm -f $@
	cat $(f51f_obj) > $@

$(f51f_dir)/%: bellman bellman.sh $(f51f_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25  --bidder_prob 0.0001 --bidder_omega 0.25 --scale 2  --rounds $(f51_n)  >  $@

$(f51f_dir)/.dir_created :
	mkdir $(f51f_dir)
	touch $@


figure51: $(f51a_sum) $(f51b_sum) $(f51c_sum) $(f51d_sum) $(f51e_sum) $(f51f_sum)



##############################################################################################################################
#
# Figure 6:  Comparison of two testimators (both unconstrained)
#
#############################################################################################################################

f6_n = 500

f6_alpha     = 0.05
f6_alpha_txt =   05

f6_angles := $(base_angles)  138 139 139.5 140 140.5 141 142 143 145 146 147 148 149 150 155

f6a := risk_alpha$(f6_alpha_txt)_100_100_beta10_100_100_scale2_n$(f6_n)
f6b := risk_alpha$(f6_alpha_txt)_100_100_beta15_100_100_scale2_n$(f6_n)
f6c := risk_alpha$(f6_alpha_txt)_100_100_beta20_100_100_scale2_n$(f6_n)


# --- F6A

f6a_dir := figures/f6/$(f6a)

f6a_obj := $(addprefix $(f6a_dir)/, $(f6_angles))

f6a_sum := figures/f6/summary_$(f6a)

$(f6a_sum): $(f6a_obj)
	rm -f $@
	cat $(f6a_obj) > $@

$(f6a_dir)/%: bellman bellman.sh $(f6a_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_w0 $(f6_alpha) --oracle_omega 1 --bidder_w0 0.10  --bidder_omega 1 --scale 2  --rounds $(f6_n)  >$@

$(f6a_dir)/.dir_created :
	mkdir $(f6a_dir)
	touch $@


# --- F6B

f6b_dir := figures/f6/$(f6b)

f6b_obj := $(addprefix $(f6b_dir)/, $(f6_angles))

f6b_sum := figures/f6/summary_$(f6b)

$(f6b_sum): $(f6b_obj)
	rm -f $@
	cat $(f6b_obj) > $@

$(f6b_dir)/%: bellman bellman.sh $(f6b_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_w0 $(f6_alpha) --oracle_omega 1 --bidder_w0 0.15  --bidder_omega 1 --scale 2  --rounds $(f6_n)  >$@

$(f6b_dir)/.dir_created :
	mkdir $(f6b_dir)
	touch $@


# --- F6C

f6c_dir := figures/f6/$(f6c)

f6c_obj := $(addprefix $(f6c_dir)/, $(f6_angles))

f6c_sum := figures/f6/summary_$(f6c)

$(f6c_sum): $(f6c_obj)
	rm -f $@
	cat $(f6c_obj) > $@

$(f6c_dir)/%: bellman bellman.sh $(f6c_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_w0 $(f6_alpha) --oracle_omega 1 --bidder_w0 0.20  --bidder_omega 1 --scale 2  --rounds $(f6_n)  >$@

$(f6c_dir)/.dir_created :
	mkdir $(f6c_dir)
	touch $@


figure6: $(f6a_sum) $(f6b_sum) $(f6c_sum)




############################################################################################################
#
# Figure X,Y: Universal vs Universal, with varying initial wealths
#
#     Manual... need to make directory figures/fX
#
############################################################################################################

fX_n = 200

fX_angles =  2.5 5 10    97 112 142 147 152 157 162 167 171 172 173 175.5 176 176.33 176.67 177    292 302 307 312 317 322

fX_angles := $(base_angles) $(fX_angles)

fXa := risk_alpha50_00_50_beta25_00_25_scale2_n$(fX_n)
fXb := risk_alpha50_00_50_beta10_00_10_scale2_n$(fX_n)
fXc := risk_alpha50_00_50_beta35_00_35_scale2_n$(fX_n)
fXd := risk_alpha50_00_50_beta45_00_45_scale2_n$(fX_n)
fXe := risk_alpha50_00_50_beta05_00_05_scale2_n$(fX_n)


# --- FXA

fXa_dir := figures/fX/$(fXa)

fXa_obj := $(addprefix $(fXa_dir)/, $(fX_angles))

fXa_sum := figures/fX/summary_$(fXa)

$(fXa_sum): $(fXa_obj)
	rm -f $@
	cat $(fXa_obj) > $@

$(fXa_dir)/%: bellman bellman.sh $(fXa_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.50 --bidder_prob 0 --bidder_omega 0.25 --rounds $(fX_n)  >  $@

$(fXa_dir)/.dir_created :
	mkdir $(fXa_dir)
	touch $@


# --- FXB

fXb_dir := figures/fX/$(fXb)

fXb_obj := $(addprefix $(fXb_dir)/, $(fX_angles))

fXb_sum := figures/fX/summary_$(fXb)

$(fXb_sum): $(fXb_obj)
	rm -f $@
	cat $(fXb_obj) > $@

$(fXb_dir)/%: bellman bellman.sh $(fXb_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.50 --bidder_prob 0 --bidder_omega 0.10 --rounds $(fX_n)  >  $@

$(fXb_dir)/.dir_created :
	mkdir $(fXb_dir)
	touch $@


# --- FXC

fXc_dir := figures/fX/$(fXc)

fXc_obj := $(addprefix $(fXc_dir)/, $(fX_angles))

fXc_sum := figures/fX/summary_$(fXc)

$(fXc_sum): $(fXc_obj)
	rm -f $@
	cat $(fXc_obj) > $@

$(fXc_dir)/%: bellman bellman.sh $(fXc_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.50 --bidder_prob 0 --bidder_omega 0.35 --rounds $(fX_n)  >  $@

$(fXc_dir)/.dir_created :
	mkdir $(fXc_dir)
	touch $@

# --- FXD

fXd_dir := figures/fX/$(fXd)

fXd_obj := $(addprefix $(fXd_dir)/, $(fX_angles))

fXd_sum := figures/fX/summary_$(fXd)

$(fXd_sum): $(fXd_obj)
	rm -f $@
	cat $(fXd_obj) > $@

$(fXd_dir)/%: bellman bellman.sh $(fXd_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.50 --bidder_prob 0 --bidder_omega 0.45 --rounds $(fX_n)  >  $@

$(fXd_dir)/.dir_created :
	mkdir $(fXd_dir)
	touch $@

# --- FXE

fXe_dir := figures/fX/$(fXe)

fXe_obj := $(addprefix $(fXe_dir)/, $(fX_angles))

fXe_sum := figures/fX/summary_$(fXe)

$(fXe_sum): $(fXe_obj)
	rm -f $@
	cat $(fXe_obj) > $@

$(fXe_dir)/%: bellman bellman.sh $(fXe_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.50 --bidder_prob 0 --bidder_omega 0.05 --rounds $(fX_n)  >  $@

$(fXe_dir)/.dir_created :
	mkdir $(fXe_dir)
	touch $@

# ---

figureX: $(fXa_sum) $(fXb_sum)  $(fXc_sum) $(fXd_sum)  $(fXe_sum) 

# virtual machine...
# cp figures/fX/summary_risk_alpha50_00_50_beta25_00_25_scale2_n50 /mnt/hgfs/bob/C/bellman/figures/fX


# ---------------------------------------  FY  -----------------------------------------------------

fY_n = 200

fY_angles =  97 112 142 147 152 157 162 167 171 172 173 175.5 176 176.33 176.67 177 178 179.71 179.72 179.73 179.74 179.75 179.8  \
            292 302 307 312 317 322     0.25 0.26 0.27 0.4 2.5 5 10  

fY_angles := $(base_angles)  $(fY_angles)

fY0 := risk_alpha25_00_25_beta02_00_02_scale2_n$(fY_n)
fYa := risk_alpha25_00_25_beta05_00_05_scale2_n$(fY_n)
fYb := risk_alpha25_00_25_beta10_00_10_scale2_n$(fY_n)
fYc := risk_alpha25_00_25_beta15_00_15_scale2_n$(fY_n)
fYd := risk_alpha25_00_25_beta20_00_20_scale2_n$(fY_n)
fYe := risk_alpha25_00_25_beta30_00_30_scale2_n$(fY_n)
fYf := risk_alpha25_00_25_beta35_00_35_scale2_n$(fY_n)
fYg := risk_alpha25_00_25_beta40_00_40_scale2_n$(fY_n)
fYh := risk_alpha25_00_25_beta45_00_45_scale2_n$(fY_n)
fYi := risk_alpha25_00_25_beta50_00_50_scale2_n$(fY_n)
fYj := risk_alpha25_00_25_beta60_00_60_scale2_n$(fY_n)
fYk := risk_alpha25_00_25_beta70_00_70_scale2_n$(fY_n)
fYm := risk_alpha25_00_25_beta90_00_90_scale2_n$(fY_n)

# --- FYA

fYa_dir := figures/fY/$(fYa)

fYa_obj := $(addprefix $(fYa_dir)/, $(fY_angles))

fYa_sum := figures/fY/summary_$(fYa)

$(fYa_sum): $(fYa_obj)
	rm -f $@
	cat $(fYa_obj) > $@

$(fYa_dir)/%: bellman bellman.sh $(fYa_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25 --bidder_prob 0 --bidder_omega 0.05 --rounds $(fY_n)  >  $@

$(fYa_dir)/.dir_created :
	mkdir $(fYa_dir)
	touch $@


# --- FYB

fYb_dir := figures/fY/$(fYb)

fYb_obj := $(addprefix $(fYb_dir)/, $(fY_angles))

fYb_sum := figures/fY/summary_$(fYb)

$(fYb_sum): $(fYb_obj)
	rm -f $@
	cat $(fYb_obj) > $@

$(fYb_dir)/%: bellman bellman.sh $(fYb_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25 --bidder_prob 0 --bidder_omega 0.10 --rounds $(fY_n)  >  $@

$(fYb_dir)/.dir_created :
	mkdir $(fYb_dir)
	touch $@


# --- FYC

fYc_dir := figures/fY/$(fYc)

fYc_obj := $(addprefix $(fYc_dir)/, $(fY_angles))

fYc_sum := figures/fY/summary_$(fYc)

$(fYc_sum): $(fYc_obj)
	rm -f $@
	cat $(fYc_obj) > $@

$(fYc_dir)/%: bellman bellman.sh $(fYc_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25 --bidder_prob 0 --bidder_omega 0.15 --rounds $(fY_n)  >  $@

$(fYc_dir)/.dir_created :
	mkdir $(fYc_dir)
	touch $@


# --- FYD

fYd_dir := figures/fY/$(fYd)

fYd_obj := $(addprefix $(fYd_dir)/, $(fY_angles))

fYd_sum := figures/fY/summary_$(fYd)

$(fYd_sum): $(fYd_obj)
	rm -f $@
	cat $(fYd_obj) > $@

$(fYd_dir)/%: bellman bellman.sh $(fYd_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25 --bidder_prob 0 --bidder_omega 0.20 --rounds $(fY_n)  >  $@

$(fYd_dir)/.dir_created :
	mkdir $(fYd_dir)
	touch $@

# --- FYE

fYe_dir := figures/fY/$(fYe)

fYe_obj := $(addprefix $(fYe_dir)/, $(fY_angles))

fYe_sum := figures/fY/summary_$(fYe)

$(fYe_sum): $(fYe_obj)
	rm -f $@
	cat $(fYe_obj) > $@

$(fYe_dir)/%: bellman bellman.sh $(fYe_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25 --bidder_prob 0 --bidder_omega 0.30 --rounds $(fY_n)  >  $@

$(fYe_dir)/.dir_created :
	mkdir $(fYe_dir)
	touch $@

# --- FYF

fYf_dir := figures/fY/$(fYf)

fYf_obj := $(addprefix $(fYf_dir)/, $(fY_angles))

fYf_sum := figures/fY/summary_$(fYf)

$(fYf_sum): $(fYf_obj)
	rm -f $@
	cat $(fYf_obj) > $@

$(fYf_dir)/%: bellman bellman.sh $(fYf_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25 --bidder_prob 0 --bidder_omega 0.35 --rounds $(fY_n)  >  $@

$(fYf_dir)/.dir_created :
	mkdir $(fYf_dir)
	touch $@


# --- FYG

fYg_dir := figures/fY/$(fYg)

fYg_obj := $(addprefix $(fYg_dir)/, $(fY_angles))

fYg_sum := figures/fY/summary_$(fYg)

$(fYg_sum): $(fYg_obj)
	rm -f $@
	cat $(fYg_obj) > $@

$(fYg_dir)/%: bellman bellman.sh $(fYg_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25 --bidder_prob 0 --bidder_omega 0.40 --rounds $(fY_n)  >  $@

$(fYg_dir)/.dir_created :
	mkdir $(fYg_dir)
	touch $@

# --- FYH

fYh_dir := figures/fY/$(fYh)

fYh_obj := $(addprefix $(fYh_dir)/, $(fY_angles))

fYh_sum := figures/fY/summary_$(fYh)

$(fYh_sum): $(fYh_obj)
	rm -f $@
	cat $(fYh_obj) > $@

$(fYh_dir)/%: bellman bellman.sh $(fYh_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25 --bidder_prob 0 --bidder_omega 0.45 --rounds $(fY_n)  >  $@

$(fYh_dir)/.dir_created :
	mkdir $(fYh_dir)
	touch $@

# --- FYI

fYi_dir := figures/fY/$(fYi)

fYi_obj := $(addprefix $(fYi_dir)/, $(fY_angles))

fYi_sum := figures/fY/summary_$(fYi)

$(fYi_sum): $(fYi_obj)
	rm -f $@
	cat $(fYi_obj) > $@

$(fYi_dir)/%: bellman bellman.sh $(fYi_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25 --bidder_prob 0 --bidder_omega 0.50 --rounds $(fY_n)  >  $@

$(fYi_dir)/.dir_created :
	mkdir $(fYi_dir)
	touch $@

# --- FYJ

fYj_dir := figures/fY/$(fYj)

fYj_obj := $(addprefix $(fYj_dir)/, $(fY_angles))

fYj_sum := figures/fY/summary_$(fYj)

$(fYj_sum): $(fYj_obj)
	rm -f $@
	cat $(fYj_obj) > $@

$(fYj_dir)/%: bellman bellman.sh $(fYj_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25 --bidder_prob 0 --bidder_omega 0.60 --rounds $(fY_n)  >  $@

$(fYj_dir)/.dir_created :
	mkdir $(fYj_dir)
	touch $@


# --- FYK

fYk_dir := figures/fY/$(fYk)

fYk_obj := $(addprefix $(fYk_dir)/, $(fY_angles))

fYk_sum := figures/fY/summary_$(fYk)

$(fYk_sum): $(fYk_obj)
	rm -f $@
	cat $(fYk_obj) > $@

$(fYk_dir)/%: bellman bellman.sh $(fYk_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25 --bidder_prob 0 --bidder_omega 0.70 --rounds $(fY_n)  >  $@

$(fYk_dir)/.dir_created :
	mkdir $(fYk_dir)
	touch $@

# --- FY0

fY0_dir := figures/fY/$(fY0)

fY0_obj := $(addprefix $(fY0_dir)/, $(fY_angles))

fY0_sum := figures/fY/summary_$(fY0)

$(fY0_sum): $(fY0_obj)
	rm -f $@
	cat $(fY0_obj) > $@

$(fY0_dir)/%: bellman bellman.sh $(fY0_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25 --bidder_prob 0 --bidder_omega 0.02 --rounds $(fY_n)  >  $@

$(fY0_dir)/.dir_created :
	mkdir $(fY0_dir)
	touch $@

# --- FYM

fYm_dir := figures/fY/$(fYm)

fYm_obj := $(addprefix $(fYm_dir)/, $(fY_angles))

fYm_sum := figures/fY/summary_$(fYm)

$(fYm_sum): $(fYm_obj)
	rm -f $@
	cat $(fYm_obj) > $@

$(fYm_dir)/%: bellman bellman.sh $(fYm_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25 --bidder_prob 0 --bidder_omega 0.90 --rounds $(fY_n)  >  $@

$(fYm_dir)/.dir_created :
	mkdir $(fYm_dir)
	touch $@

# ---

figureY: $(fYa_sum) $(fYb_sum) $(fYc_sum) $(fYd_sum) $(fYe_sum) $(fYf_sum) $(fYg_sum) $(fYh_sum) $(fYi_sum) $(fYj_sum) $(fYk_sum) $(fY0_sum) $(fYm_sum)

# virtual machine...
# cp figures/fY/summary_risk_alpha50_00_50_beta25_00_25_scale2_n50 /mnt/hgfs/bob/C/bellman/figures/fY


# ---------------------------------------  FZ  -----------------------------------------------------

fZ_n = 200

fZ_angles =  97 112 142 147 152 157 162 167 171 172 173 175.5 176 176.33 176.67 177 178 179.71 179.72 179.73 179.74 179.75 179.8  \
            292 302 307 312 317 322     0.25 0.26 0.27 0.4 2.5 5 10  

fZ_angles := $(base_angles) # $(fZ_angles)

fZa := risk_alpha25_00_25_beta05_00_25_scale2_n$(fZ_n)
fZb := risk_alpha25_00_25_beta10_00_25_scale2_n$(fZ_n)
fZc := risk_alpha25_00_25_beta20_00_25_scale2_n$(fZ_n)
fZd := risk_alpha25_00_25_beta30_00_25_scale2_n$(fZ_n)
fZe := risk_alpha25_00_25_beta40_00_25_scale2_n$(fZ_n)
fZf := risk_alpha25_00_25_beta50_00_25_scale2_n$(fZ_n)

# --- FZa

fZa_dir := figures/fZ/$(fZa)

fZa_obj := $(addprefix $(fZa_dir)/, $(fZ_angles))

fZa_sum := figures/fZ/summary_$(fZa)

$(fZa_sum): $(fZa_obj)
	rm -f $@
	cat $(fZa_obj) > $@

$(fZa_dir)/%: bellman bellman.sh $(fZa_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25 --bidder_w0 0.05 --bidder_prob 0 --bidder_omega 0.25 --rounds $(fZ_n)  >  $@

$(fZa_dir)/.dir_created :
	mkdir $(fZa_dir)
	touch $@

# --- FZb

fZb_dir := figures/fZ/$(fZb)

fZb_obj := $(addprefix $(fZb_dir)/, $(fZ_angles))

fZb_sum := figures/fZ/summary_$(fZb)

$(fZb_sum): $(fZb_obj)
	rm -f $@
	cat $(fZb_obj) > $@

$(fZb_dir)/%: bellman bellman.sh $(fZb_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25 --bidder_w0 0.10 --bidder_prob 0 --bidder_omega 0.25 --rounds $(fZ_n)  >  $@

$(fZb_dir)/.dir_created :
	mkdir $(fZb_dir)
	touch $@


# --- FZc

fZc_dir := figures/fZ/$(fZc)

fZc_obj := $(addprefix $(fZc_dir)/, $(fZ_angles))

fZc_sum := figures/fZ/summary_$(fZc)

$(fZc_sum): $(fZc_obj)
	rm -f $@
	cat $(fZc_obj) > $@

$(fZc_dir)/%: bellman bellman.sh $(fZc_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25 --bidder_w0 0.20 --bidder_prob 0 --bidder_omega 0.25 --rounds $(fZ_n)  >  $@

$(fZc_dir)/.dir_created :
	mkdir $(fZc_dir)
	touch $@


# --- FZd

fZd_dir := figures/fZ/$(fZd)

fZd_obj := $(addprefix $(fZd_dir)/, $(fZ_angles))

fZd_sum := figures/fZ/summary_$(fZd)

$(fZd_sum): $(fZd_obj)
	rm -f $@
	cat $(fZd_obj) > $@

$(fZd_dir)/%: bellman bellman.sh $(fZd_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25 --bidder_w0 0.30 --bidder_prob 0 --bidder_omega 0.25 --rounds $(fZ_n)  >  $@

$(fZd_dir)/.dir_created :
	mkdir $(fZd_dir)
	touch $@


# --- FZe

fZe_dir := figures/fZ/$(fZe)

fZe_obj := $(addprefix $(fZe_dir)/, $(fZ_angles))

fZe_sum := figures/fZ/summary_$(fZe)

$(fZe_sum): $(fZe_obj)
	rm -f $@
	cat $(fZe_obj) > $@

$(fZe_dir)/%: bellman bellman.sh $(fZe_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25 --bidder_w0 0.40 --bidder_prob 0 --bidder_omega 0.25 --rounds $(fZ_n)  >  $@

$(fZe_dir)/.dir_created :
	mkdir $(fZe_dir)
	touch $@

# --- FZf

fZf_dir := figures/fZ/$(fZf)

fZf_obj := $(addprefix $(fZf_dir)/, $(fZ_angles))

fZf_sum := figures/fZ/summary_$(fZf)

$(fZf_sum): $(fZf_obj)
	rm -f $@
	cat $(fZf_obj) > $@

$(fZf_dir)/%: bellman bellman.sh $(fZf_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0 --oracle_omega 0.25 --bidder_w0 0.50 --bidder_prob 0 --bidder_omega 0.25 --rounds $(fZ_n)  >  $@

$(fZf_dir)/.dir_created :
	mkdir $(fZf_dir)
	touch $@
# ---

figureZ: $(fZa_sum) $(fZb_sum) $(fZc_sum) $(fZd_sum) $(fZe_sum) $(fZf_sum) # $(fZg_sum) $(fZh_sum) $(fZi_sum) $(fZj_sum) $(fZk_sum) 




###########################################################################
include ../rules_for_makefiles

