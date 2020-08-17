
TL = ../trajeclike -v2
RELENT = ../relent.js

all: relents stats times

times:
	../trajeclike --maxlen 40 --benchmark --trials 1000 >$@

# Transform (gamma,mu,rins/gamma,rdel) of trajeclike into (lambda,mu,x,y) of general indel model
relents: scan-gammas
	cat *results | perl -e 'print"gamma lambda mu rdel rins t method D reversible\n";while(<>){@m=/^D g(\S+)-m(\S+)-r(\S+)-i(\S+)-t(\S+)-sim \S+-(\S+) (\S+)/;if(@m){push@m,($$m[3]==$$m[2]?"yes":"no");$$m[3]*=$$m[0];$$m[1]*=1-$$m[2];unshift@m,$$m[0]*$$m[1]*(1-$$m[2])/(1-$$m[3]);@m[0,1]=@m[1,0];print"@m\n"}}' >$@

stats: scan-gammas
	cat *results | perl -e 'print"stat gamma lambda mu rdel rins t method value reversible\n";while(<>){@m=/^([^D]\S*) g(\S+)-m(\S+)-r(\S+)-i(\S+)-t(\S+)-(\S+) (\S+)/;if(@m){push@m,($$m[4]==$$m[3]?"yes":"no");$$m[4]*=$$m[1];$$m[2]*=1-$$m[3];unshift@m,$$m[1]*$$m[2]*(1-$$m[3])/(1-$$m[4]);@m[0..2]=@m[1,2,0];print"@m\n"}}' >$@

#scan-gammas: g1-scan-mu g.99-scan-mu g.9-scan-mu g1-m.5-r.5-i.5-scan-time-fine g.8-m.5-r.5-i.5-scan-time-fine g1-m1-r.5-t.25-scan-rins-fine
#scan-gammas: g1-m.5-r.5-i.5-scan-time-fine g.8-m.5-r.5-i.5-scan-time-fine g1-m1-r.5-t.25-scan-rins-fine g.99-m1-r.5-t.25-scan-rins-fine
scan-gammas: g1-m.5-r.5-i.5-scan-time-fine g.8-m.5-r.5-i.5-scan-time-fine g1-m1-r.5-t.25-scan-rins-fine

$(PREFIX)-scan-mu: $(PREFIX)-m.05-scan-rdel $(PREFIX)-m.5-scan-rdel

$(PREFIX)-scan-rdel: $(PREFIX)-r.5-scan-rins $(PREFIX)-r.75-scan-rins $(PREFIX)-r.9-scan-rins

$(PREFIX)-r$(R)-scan-rins: $(PREFIX)-r$(R)-i$(R)-scan-time $(PREFIX)-r$(R)-i.45-scan-time $(PREFIX)-r$(R)-i.8-scan-time-short $(PREFIX)-r$(R)-i.95-scan-time-tiny

$(PREFIX)-scan-time: $(PREFIX)-scan-time-short $(PREFIX)-t1-results $(PREFIX)-t2-results $(PREFIX)-t4-results $(PREFIX)-t8-results

$(PREFIX)-scan-time-short: $(PREFIX)-scan-time-tiny $(PREFIX)-t.125-results $(PREFIX)-t.25-results $(PREFIX)-t.5-results

$(PREFIX)-scan-time-tiny: $(PREFIX)-t.015625-results $(PREFIX)-t.03125-results $(PREFIX)-t.0625-results

TIME_RANGE = $(shell perl -e '@m=grep(s/^0\././ || 1,map(2**($$_/8),-48..24));print"@m"')

$(PREFIX)-scan-time-fine: $(addprefix $(PREFIX)-t,$(addsuffix -results,$(TIME_RANGE)))

RINS_RANGE = $(shell perl -e '@m=grep(s/^0\././ || 1,map($$_/100,20..80));print"@m"')

g$(G)-m$(M)-r$(R)-t$(T)-scan-rins-fine: $(addprefix g$(G)-m$(M)-r$(R)-i,$(addsuffix -t$(T)-results,$(RINS_RANGE)))

$(PREFIX)-results: $(PREFIX)-sim $(PREFIX)-moments $(PREFIX)-traj $(PREFIX)-cim $(PREFIX)-tkf91 $(PREFIX)-tkf92 $(PREFIX)-rs07 $(PREFIX)-prank
	$(RELENT) --normalize $^ >$@

g$(G)-m$(M)-r$(R)-i$(I)-t$(T)-sim: { atom_concat('0',T,Tz), atom_number(Tz,Tn), Tn =< 0.0625 }
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --rins $(I) --time $(T) --seedtime --simulate --trials 100000000 --initlen 1000 --maxlen 100 >$@) >&$@.log
	grep -v "Beginning simulation trial" $@.log >$@.log.short
	rm $@.log

g$(G)-m$(M)-r$(R)-i$(I)-t$(T)-sim: { atom_concat('0',T,Tz), atom_number(Tz,Tn), Tn > 0.0625, Tn =< 0.1 }
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --rins $(I) --time $(T) --seedtime --simulate --trials 10000000 --initlen 1000 --maxlen 100 >$@) >&$@.log
	grep -v "Beginning simulation trial" $@.log >$@.log.short
	rm $@.log

g$(G)-m$(M)-r$(R)-i$(I)-t$(T)-sim: { atom_concat('0',T,Tz), atom_number(Tz,Tn), Tn > 0.1, Tn < 1 }
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --rins $(I) --time $(T) --seedtime --simulate --trials 1000000 --initlen 1000 --maxlen 100 >$@) >&$@.log
	grep -v "Beginning simulation trial" $@.log >$@.log.short
	rm $@.log

g$(G)-m$(M)-r$(R)-i$(I)-t$(T)-sim: { atom_concat('0',T,Tz), atom_number(Tz,Tn), Tn >= 1 }
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --rins $(I) --time $(T) --seedtime --simulate --trials 100 --initlen 100000 --maxlen 1000 >$@) >&$@.log
	grep -v "Beginning simulation trial" $@.log >$@.log.short
	rm $@.log

g$(G)-m$(M)-r$(R)-i$(I)-t$(T)-traj:
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --rins $(I) --time $(T) --maxlen 30 >$@) >&$@.log
	grep -v "Likelihood for" $@.log >$@.log.short
	rm $@.log


g$(G)-m$(M)-r$(R)-i$(I)-t$(T)-cim: { atom_concat('0',T,Tz), atom_number(Tz,Tn), Tn =< 0.1 }
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --rins $(I) --time $(T) --cim --maxlen 100 --dt .00001 >$@) >&$@.log

g$(G)-m$(M)-r$(R)-i$(I)-t$(T)-cim: { atom_concat('0',T,Tz), atom_number(Tz,Tn), Tn > 0.1, (Tn < 1 -> MaxLen = 100; MaxLen = 1000) }
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --rins $(I) --time $(T) --cim --maxlen $(MaxLen) --dt .0001 >$@) >&$@.log


g$(G)-m$(M)-r$(R)-i$(I)-t$(T)-moments: { atom_concat('0',T,Tz), atom_number(Tz,Tn), Tn =< 0.1 }
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --rins $(I) --time $(T) --moments --maxlen 100 --dt .00001 >$@) >&$@.log

g$(G)-m$(M)-r$(R)-i$(I)-t$(T)-moments: { atom_concat('0',T,Tz), atom_number(Tz,Tn), Tn > 0.1, (Tn < 1 -> MaxLen = 100; MaxLen = 1000) }
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --rins $(I) --time $(T) --moments --maxlen $(MaxLen) --dt .0001 >$@) >&$@.log


g$(G)-m$(M)-r$(R)-i$(I)-t$(T)-tkf91: { atom_concat('0',T,Tz), atom_number(Tz,Tn), (Tn < 1 -> MaxLen = 100; MaxLen = 1000) }
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --rins $(I) --time $(T) --tkf91 --maxlen $(MaxLen) >$@) >&$@.log

g$(G)-m$(M)-r$(R)-i$(I)-t$(T)-tkf92: { atom_concat('0',T,Tz), atom_number(Tz,Tn), (Tn < 1 -> MaxLen = 100; MaxLen = 1000) }
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --rins $(I) --time $(T) --tkf92 --maxlen $(MaxLen) >$@) >&$@.log

g$(G)-m$(M)-r$(R)-i$(I)-t$(T)-rs07: { atom_concat('0',T,Tz), atom_number(Tz,Tn), (Tn < 1 -> MaxLen = 100; MaxLen = 1000) }
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --rins $(I) --time $(T) --rs07 --maxlen $(MaxLen) >$@) >&$@.log

g$(G)-m$(M)-r$(R)-i$(I)-t$(T)-prank: { atom_concat('0',T,Tz), atom_number(Tz,Tn), (Tn < 1 -> MaxLen = 100; MaxLen = 1000) }
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --rins $(I) --time $(T) --prank --maxlen $(MaxLen) >$@) >&$@.log


shorten-logs: $(addsuffix .short,$(wildcard *-traj.log) $(wildcard *-sim.log))

%-traj.log.short: %-traj.log
	grep -v "Likelihood for" $< >$@

%-sim.log.short: %-sim.log
	grep -v "Beginning simulation trial" $< >$@


.SECONDARY:
