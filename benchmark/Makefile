
TL = ../trajeclike -v2
RELENT = ../relent.js

all: relents stats times

times:
	../trajeclike --maxlen 40 --benchmark --trials 1000 >$@

# Transform (gamma,mu,rins/gamma,rdel) of trajeclike into (lambda,mu,x,y) of general indel model
relents: scan-gammas
	cat *results | perl -e 'print"lambda mu x y t maxgap seqlen n method D reversible\n";while(<>){@m=/^D l(\S+)-m(\S+)-x(\S+)-y(\S+)-t(\S+)-g(\S+)-s(\S+)-n(\S+)-sim \S+-(\S+) (\S+)/;if(@m){print"@m\n"}}' >$@

stats: scan-gammas
	cat *results | perl -e 'print"stat lambda mu x y t maxgap seqlen n method value reversible\n";while(<>){@m=/^([^D]\S*) l(\S+)-m(\S+)-x(\S+)-y(\S+)-t(\S+)-g(\S+)-s(\S+)-n(\S+)-(\S+) (\S+)/;if(@m){print"@m\n"}}' >$@

scan-gammas: x.5-y.5-scan-t y.5-t.5-scan-x x.5-t.5-scan-y t.5-scan-gaplen

TIME_RANGE = $(shell perl -e '@m=grep(s/^0\././ || 1,map(2**($$_/8),-48..24));print"@m"')
X_RANGE = $(shell perl -e '@m=grep(s/^0\././ || 1,map($$_/100,40..60));print"@m"')
Y_RANGE = $(shell perl -e '@m=grep(s/^0\././ || 1,map($$_/100,40..60));print"@m"')
GAPLEN_RANGE = $(shell perl -e '@m=grep(s/^0\././ || 1,map($$_/10,10..40));print"@m"')

x$(X)-y$(Y)-scan-t: $(addprefix l1-m1-x$(X)-y$(Y)-t,$(addsuffix -set-seqlen,$(TIME_RANGE)))
	echo $^

y$(Y)-t$(T)-scan-x: $(addprefix l1-m1-x,$(addsuffix -y$(Y)-t$(T)-set-seqlen,$(X_RANGE)))
	echo $^

x$(X)-t$(T)-scan-y: $(addprefix l1-m1-x$(X)-y,$(addsuffix -t$(T)-set-seqlen,$(Y_RANGE)))
	echo $^

t$(T)-scan-gaplen: $(addprefix gaplen,$(addsuffix -t$(T)-set-xy,$(GAPLEN_RANGE)))
	echo $^

gaplen$(L)-t$(T)-set-xy { atom_number(L,Ln), Y is 1 - (1 / Ln) }: l1-m1-x$(Y)-y$(Y)-t$(T)-set-seqlen
	echo $^

l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-set-seqlen { atom_concat('0',T,Tz), atom_number(Tz,Tn), (Tn => 1 -> Trials = 100, InitLen = 100000, MaxGap = 1000; InitLen = 1000, MaxGap = 100, (Tn => 0.0625 -> Trials = 100000; (Tn => 0.03125 -> Trials = 1000000; (Trials = 10000000)))) }: l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-g$(MaxGap)-s$(InitLen)-n$(Trials)-results

l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-g$(G)-s$(S)-n$(N)-results: l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-g$(G)-s$(S)-n$(N)-sim l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-g$(G)-methods
	$(RELENT) --normalize $^ >$@

$(PREFIX)-methods: $(PREFIX)-moments $(PREFIX)-traj $(PREFIX)-cim $(PREFIX)-tkf91 $(PREFIX)-tkf92 $(PREFIX)-rs07 $(PREFIX)-prank

l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-g$(G)-s$(S)-n$(N)-sim:
	(time $(TL) --ggi --lambda $(L) --mu $(M) --x $(X) --y $(Y) --time $(T) --maxlen $(G) --seedtime --simulate --trials $(N) --initlen $(S) >$@) >&$@.log
	grep -v "Beginning simulation trial" $@.log >$@.log.short
	rm $@.log

# Override gaplen for trajectory method, which is sloooow
MAX_TRAJ_GAPLEN = 30
l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-g$(G)-traj:
	(time $(TL) --ggi --lambda $(L) --mu $(M) --x $(X) --y $(Y) --time $(T) --maxlen $(MAX_TRAJ_GAPLEN) >$@) >&$@.log
	grep -v "Likelihood for" $@.log >$@.log.short
	rm $@.log

l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-g$(G)-cim: { atom_concat('0',T,Tz), atom_number(Tz,Tn), (Tn > 0.1 -> DT = 0.0001; DT = 0.00001) }
	(time $(TL) --ggi --lambda $(L) --mu $(M) --x $(X) --y $(Y) --time $(T) --maxlen $(G) --cim --dt $(DT) >$@) >&$@.log

l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-g$(G)-moments: { atom_concat('0',T,Tz), atom_number(Tz,Tn), (Tn > 0.1 -> DT = 0.0001; DT = 0.00001) }
	(time $(TL) --ggi --lambda $(L) --mu $(M) --x $(X) --y $(Y) --time $(T) --maxlen $(G) --moments --dt $(DT) >$@) >&$@.log

l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-g$(G)-tkf91:
	(time $(TL) --ggi --lambda $(L) --mu $(M) --x $(X) --y $(Y) --time $(T) --maxlen $(G) --tkf91 >$@) >&$@.log

l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-g$(G)-tkf92:
	(time $(TL) --ggi --lambda $(L) --mu $(M) --x $(X) --y $(Y) --time $(T) --maxlen $(G) --tkf92 >$@) >&$@.log

l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-g$(G)-rs07:
	(time $(TL) --ggi --lambda $(L) --mu $(M) --x $(X) --y $(Y) --time $(T) --maxlen $(G) --rs07 >$@) >&$@.log

l$(L)-m$(M)-x$(X)-y$(Y)-t$(T)-g$(G)-prank:
	(time $(TL) --ggi --lambda $(L) --mu $(M) --x $(X) --y $(Y) --time $(T) --maxlen $(G) --prank >$@) >&$@.log

shorten-logs: $(addsuffix .short,$(wildcard *-traj.log) $(wildcard *-sim.log))

%-traj.log.short: %-traj.log
	grep -v "Likelihood for" $< >$@

%-sim.log.short: %-sim.log
	grep -v "Beginning simulation trial" $< >$@


.SECONDARY:
