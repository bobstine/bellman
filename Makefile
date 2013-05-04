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

level_1 = distribution.o spending_rule.o
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

bellman: bellman.o wealth.o utility.o bellman_main.o spending_rule.o
	$(GCC) $^ $(LDLIBS) -o  $@

#  Bellman options (see code for details)
#       constrain  : oracle has state (2 dim)
#	oracle prob: 0	LS, 1 RI, otherwise testimator
#	bidder prob: 0 for univ, alpha for geometric
#       omega      : 0 for fixed bidder

reject_check: bellman
	./bellman --reject --angle 0 --rounds 7  --oracle_omega 0.05 --oracle_prob 0  --bidder_omega 0.5  --bidder_prob 0.1 --write

risk_check: bellman
	./bellman --risk --angle 125 --rounds 10 --oracle_w0 0.10 --oracle_omega 1 --bidder_w0 0.10 --bidder_omega 1  --write  # both unconstrained, constant mean

#	./bellman --risk --angle 0 --rounds 100 --oracle_omega 0.5  --oracle_prob 0  --bidder_omega 0.5 --bidder_prob 0.10  # constrained
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
f0_angles := $(base_angles) $(extra_angles) $(f0_more_angles) 175.1 175.11 175.12 175.13 175.14 175.15 175.151 175.152 175.153 175.154 175.156 175.157 175.158 175.159

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
f0A_mu := 0.5 0.75 1.0 1.5 2.0 2.5 3.0 3.5 4.0

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
# Figure 1: comparisons of geometric spending testimator to risk inflation oracle, with varying payouts
#
#############################################################################################################################

f1_n = 1000

# for RI oracle
f1_alpha     = 1
f1_alpha_txt = 100

# for Testimator oracle
# f1_alpha     = 0.10
# f1_alpha_txt =   10

f1_beta      = 0.10
f1_beta_txt  =   10

f1_angles := $(base_angles)  $(extra_angles) 

f1a := risk_alpha100_$(f1_alpha_txt)_100_beta05_$(f1_beta_txt)_05_scale2_n$(f1_n)
f1aa:= risk_alpha100_$(f1_alpha_txt)_100_beta10_$(f1_beta_txt)_10_scale2_n$(f1_n)
f1b := risk_alpha100_$(f1_alpha_txt)_100_beta25_$(f1_beta_txt)_25_scale2_n$(f1_n)
f1c := risk_alpha100_$(f1_alpha_txt)_100_beta50_$(f1_beta_txt)_50_scale2_n$(f1_n)
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
	./bellman --risk --angle $* --oracle_prob $(f1_alpha) --oracle_omega 1 --bidder_prob $(f1_beta) --bidder_omega 0.05 --scale 2  --rounds $(f1_n)  >$@

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
	./bellman --risk --angle $* --oracle_prob $(f1_alpha) --oracle_omega 1 --bidder_prob $(f1_beta) --bidder_omega 0.10 --scale 2  --rounds $(f1_n)  >$@

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
	./bellman --risk --angle $* --oracle_prob $(f1_alpha) --oracle_omega 1 --bidder_prob $(f1_beta) --bidder_omega 0.25 --scale 2  --rounds $(f1_n)  >$@

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
	./bellman --risk --angle $* --oracle_prob $(f1_alpha) --oracle_omega 1 --bidder_prob $(f1_beta) --bidder_omega 0.50  --scale 2  --rounds $(f1_n)  >$@

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
figure1: $(f1a_sum) $(f1aa_sum) $(f1b_sum) $(f1c_sum) $(f1d_sum)    $(f1ea_sum) $(f1e_sum) $(f1f_sum) $(f1g_sum)    $(f1h_sum) $(f1ha_sum) $(f1i_sum) $(f1j_sum)




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
# Figure 5: Universal vs geometric bidder
#
############################################################################################################

f5_n = 100
f5_angles := $(base_angles) $(extra_angles)

f5a := risk_alpha25_01_25_beta25_00_25_scale2_n$(f5_n)
f5b := risk_alpha25_05_25_beta25_00_25_scale2_n$(f5_n)
f5c := risk_alpha25_10_25_beta25_00_25_scale2_n$(f5_n)

# --- F5A

f5a_dir := figures/f5/$(f5a)

f5a_obj := $(addprefix $(f5a_dir)/, $(f5_angles))

f5a_sum := figures/f5/summary_$(f5a)

$(f5a_sum): $(f5a_obj)
	rm -f $@
	cat $(f5a_obj) > $@

$(f5a_dir)/%: bellman bellman.sh $(f5a_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0.01 --oracle_omega 0.25  --bidder_prob 0 --bidder_omega 0.25 --scale 2  --rounds $(f5_n)  >  $@

$(f5a_dir)/.dir_created :
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
	./bellman --risk --angle $* --oracle_prob 0.05 --oracle_omega 0.25  --bidder_prob 0 --bidder_omega 0.25 --scale 2  --rounds $(f5_n)  >  $@

$(f5b_dir)/.dir_created :
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
	./bellman --risk --angle $* --oracle_prob 0.10 --oracle_omega 0.25  --bidder_prob 0 --bidder_omega 0.25 --scale 2  --rounds $(f5_n)  >  $@


$(f5c_dir)/.dir_created :
	mkdir $(f5c_dir)
	touch $@

# --- F5D
f5d_dir := figures/f5/$(f5d)

f5d_obj := $(addprefix $(f5d_dir)/, $(f5_angles))

f5d_sum := figures/f5/summary_$(f5d)

$(f5d_sum): $(f5d_obj)
	rm -f $@
	cat $(f5d_obj) > $@

$(f5d_dir)/%: bellman bellman.sh $(f5d_dir)/.dir_created 
	./bellman --risk --angle $* --oracle_prob 0.10 --oracle_omega 0.05  --bidder_prob 0.01 --bidder_omega 0.50 --scale 2  --rounds $(f5_n)  >  $@


$(f5d_dir)/.dir_created :
	mkdir $(f5d_dir)
	touch $@


figure5: $(f5a_sum) $(f5b_sum) $(f5c_sum)

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


###########################################################################
include ../rules_for_makefiles

