# This Makefile requires BioMake
# https://github.com/evoldoers/biomake

prolog
set_len(T,Y,Trials,InitLen,MaxGap) :-
	atom_concat('0',T,Tz),
	atom_number(Tz,Tn),
	atom_concat('0',Y,Yz),
	atom_number(Yz,Yn),
	((Tn >= 1; Yn > 0.6)
	 -> Trials = 100,
	    InitLen = 100000,
	    MaxGap = 1000;
	 InitLen = 1000,
	 MaxGap = 100,
	 (Tn >= 0.0625
		  -> Trials = 100000;
		  (Tn >= 0.03125
			   -> Trials = 1000000;
			   Trials = 10000000))).
endprolog

TL = ../trajeclike -v2
RELENT = ../relent.js

all: relents stats times

times:
	../trajeclike --maxlen 40 --benchmark --trials 1000 >$@

relents: scan-gammas
	cat *results | perl -e 'print"lambda mu x y t maxgap seqlen n method D\n";while(<>){@m=/^D l(\S+)-m(\S+)-x(\S+)-y(\S+)-t(\S+)-g(\S+)-s(\S+)-n(\S+)-sim \S+-(\S+) (\S+)/;if(@m){print"@m\n"}}' >$@

stats: scan-gammas
	cat *results | perl -e 'print"stat lambda mu x y t maxgap-seqlen-n method value\n";while(<>){@m=/^([^D]\S*) l(\S+)-m(\S+)-x(\S+)-y(\S+)-t(\S+)-g(\S+)-(\S+) (\S+)/;if(@m){print"@m\n"}}' >$@

scan-gammas: x.5-y.5-scan-t y.5-t.5-scan-x x.5-t.5-scan-y t.5-scan-xy x.5-y.5-t.5-scan-lambda x.5-y.5-t.5-scan-mu

TIME_RANGE = $(shell perl -e '@m=grep(s/^0\././ || 1,map(2**($$_/8),-56..8));print"@m"')
X_RANGE = $(shell perl -e '@m=grep(s/^0\././ || 1,map($$_/100,0..70));print"@m"')
Y_RANGE = $(shell perl -e '@m=grep(s/^0\././ || 1,map($$_/100,0..65));print"@m"')
XY_RANGE = $(shell perl -e '@m=grep(s/^0\././ || 1,map($$_/100,0..95));print"@m"')
LAMBDA_RANGE = $(shell perl -e '@m=grep(s/^0\././ || 1,map($$_/50,1..50));print"@m"')
MU_RANGE = $(shell perl -e '@m=grep(s/^0\././ || 1,map($$_/50,1..50));print"@m"')

x$(X)-y$(Y)-scan-t: $(addprefix l1-m1-x$(X)-y$(Y)-t,$(addsuffix -set-seqlen,$(TIME_RANGE)))
	echo $^

y$(Y)-t$(T)-scan-x: $(addprefix l1-m1-x,$(addsuffix -y$(Y)-t$(T)-set-seqlen,$(X_RANGE)))
	echo $^

x$(X)-t$(T)-scan-y: $(addprefix l1-m1-x$(X)-y,$(addsuffix -t$(T)-set-seqlen,$(Y_RANGE)))
	echo $^

t$(T)-scan-xy: $(addprefix x,$(addsuffix -t$(T)-set-y,$(XY_RANGE)))
	echo $^

x$(X)-y$(Y)-t$(T)-scan-lambda: $(addprefix l,$(addsuffix -m1-x$(X)-y$(Y)-t$(T)-set-seqlen,$(LAMBDA_RANGE)))
	echo $^

x$(X)-y$(Y)-t$(T)-scan-mu: $(addprefix l1-m,$(addsuffix -x$(X)-y$(Y)-t$(T)-set-seqlen,$(MU_RANGE)))
	echo $^

x$(X)-t$(T)-set-y: l1-m1-x$(X)-y$(X)-t$(T)-set-seqlen
	echo $^

l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-set-seqlen { set_len(T,Y,Trials,InitLen,MaxGap) }: l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-g$(MaxGap)-s$(InitLen)-n$(Trials)-results

$(PREFIX)-s$(S)-n$(N)-results: $(PREFIX)-s$(S)-n$(N)-sim $(PREFIX)-moments $(PREFIX)-traj $(PREFIX)-cim $(PREFIX)-tkf91 $(PREFIX)-tkf92 $(PREFIX)-rs07 $(PREFIX)-prank
	$(RELENT) --normalize $^ >$@

l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-g$(G)-s$(S)-n$(N)-sim:
	(time $(TL) --lambda $(L) --mu $(M) --x $(X) --y $(Y) --time $(T) --maxlen $(G) --seedtime --simulate --trials $(N) --initlen $(S) >$@) >&$@.log
	grep -v "Beginning simulation trial" $@.log >$@.log.short
	rm $@.log

# Override gaplen for trajectory method, which is sloooow
MAX_TRAJ_GAPLEN = 30
l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-g$(G)-traj:
	(time $(TL) --lambda $(L) --mu $(M) --x $(X) --y $(Y) --time $(T) --maxlen $(MAX_TRAJ_GAPLEN) >$@) >&$@.log
	grep -v "Likelihood for" $@.log >$@.log.short
	rm $@.log

l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-g$(G)-cim: { atom_concat('0',T,Tz), atom_number(Tz,Tn), (Tn > 0.1 -> DT = 0.0001; DT = 0.00001) }
	(time $(TL) --lambda $(L) --mu $(M) --x $(X) --y $(Y) --time $(T) --maxlen $(G) --cim --dt $(DT) >$@) >&$@.log

l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-g$(G)-moments: { atom_concat('0',T,Tz), atom_number(Tz,Tn), (Tn > 0.1 -> DT = 0.0001; DT = 0.00001) }
	(time $(TL) --lambda $(L) --mu $(M) --x $(X) --y $(Y) --time $(T) --maxlen $(G) --moments --dt $(DT) >$@) >&$@.log

l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-g$(G)-tkf91:
	(time $(TL) --lambda $(L) --mu $(M) --x $(X) --y $(Y) --time $(T) --maxlen $(G) --tkf91 >$@) >&$@.log

l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-g$(G)-tkf92:
	(time $(TL) --lambda $(L) --mu $(M) --x $(X) --y $(Y) --time $(T) --maxlen $(G) --tkf92 >$@) >&$@.log

l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-g$(G)-rs07:
	(time $(TL) --lambda $(L) --mu $(M) --x $(X) --y $(Y) --time $(T) --maxlen $(G) --rs07 >$@) >&$@.log

# Override maxlen for PRANK, which has a high covariance
l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-g$(G)-prank:
	(time $(TL) --lambda $(L) --mu $(M) --x $(X) --y $(Y) --time $(T) --maxlen 1000 --prank >$@) >&$@.log

shorten-logs: $(addsuffix .short,$(wildcard *-traj.log) $(wildcard *-sim.log))

%-traj.log.short: %-traj.log
	grep -v "Likelihood for" $< >$@

%-sim.log.short: %-sim.log
	grep -v "Beginning simulation trial" $< >$@


.SECONDARY:
