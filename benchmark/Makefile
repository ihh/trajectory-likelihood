
TL = ../trajeclike -v2
RELENT = ../relent.js

all: relents stats times

times:
	../trajeclike --maxlen 40 --benchmark --trials 1000 >$@

relents: scan-gammas
	cat *results | perl -e 'print"gamma mu rdel rins t method D reversible\n";while(<>){@m=/^D g(\S+)-m(\S+)-r(\S+)-i(\S+)-t(\S+)-sim \S+-(\S+) (\S+)/;if(@m){push@m,($$m[3]==$$m[2]?"yes":"no");$$m[3]*=$$m[0];print"@m\n"}}' >$@

stats: scan-gammas
	cat *results | perl -e 'print"stat gamma mu rdel rins t method value reversible\n";while(<>){@m=/^([^D]\S*) g(\S+)-m(\S+)-r(\S+)-i(\S+)-t(\S+)-(\S+) (\S+)/;if(@m){push@m,($$m[4]==$$m[3]?"yes":"no");$$m[4]*=$$m[1];print"@m\n"}}' >$@

scan-gammas: g1-scan-mu g.99-scan-mu g.9-scan-mu

$(PREFIX)-scan-mu: $(PREFIX)-m.05-scan-rdel $(PREFIX)-m.5-scan-rdel

$(PREFIX)-scan-rdel: $(PREFIX)-r.5-scan-rins $(PREFIX)-r.75-scan-rins $(PREFIX)-r.9-scan-rins

$(PREFIX)-r$(R)-scan-rins: $(PREFIX)-r$(R)-i$(R)-scan-time $(PREFIX)-r$(R)-i.45-scan-time $(PREFIX)-r$(R)-i.8-scan-time-short $(PREFIX)-r$(R)-i.95-scan-time-tiny

$(PREFIX)-scan-time: $(PREFIX)-scan-time-short $(PREFIX)-t1-results $(PREFIX)-t2-results $(PREFIX)-t4-results $(PREFIX)-t8-results

#$(PREFIX)-scan-time-short: $(PREFIX)-scan-time-tiny $(PREFIX)-t.125-results $(PREFIX)-t.176776695296637-results $(PREFIX)-t.25-results $(PREFIX)-t.353553390593274-results $(PREFIX)-t.5-results
$(PREFIX)-scan-time-short: $(PREFIX)-scan-time-tiny $(PREFIX)-t.125-results $(PREFIX)-t.25-results $(PREFIX)-t.5-results

#$(PREFIX)-scan-time-tiny: $(PREFIX)-t.015625-results $(PREFIX)-t.0220970869120796-results $(PREFIX)-t.03125-results $(PREFIX)-t.0441941738241592-results $(PREFIX)-t.0625-results $(PREFIX)-t.0883883476483184-results
$(PREFIX)-scan-time-tiny: $(PREFIX)-t.015625-results $(PREFIX)-t.03125-results $(PREFIX)-t.0625-results

g1-$(PREFIX)-results: g1-$(PREFIX)-sim g1-$(PREFIX)-moments g1-$(PREFIX)-traj g1-$(PREFIX)-cim
	$(RELENT) --normalize $^ >$@

g.9$(PREFIX)-results: g.9$(PREFIX)-sim g.9$(PREFIX)-moments g.9$(PREFIX)-traj
	$(RELENT) --normalize $^ >$@

g$(G)-m$(M)-r$(R)-i$(I)-t$(T)-sim: { atom_concat('0',T,Tz), atom_number(Tz,Tn), Tn =< 1 }
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --rins $(I) --time $(T) --seedtime --simulate --trials 1000000 --initlen 1000 --maxlen 100 >$@) >&$@.log
	grep -v "Beginning simulation trial" $@.log >$@.log.short
	rm $@.log

g$(G)-m$(M)-r$(R)-i$(I)-t$(T)-sim: { atom_concat('0',T,Tz), atom_number(Tz,Tn), Tn > 1, Tn =< 8 }
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --rins $(I) --time $(T) --seedtime --simulate --trials 100000 --initlen 1000 --maxlen 100 >$@) >&$@.log
	grep -v "Beginning simulation trial" $@.log >$@.log.short
	rm $@.log

g$(G)-m$(M)-r$(R)-i$(I)-t$(T)-sim: { atom_concat('0',T,Tz), atom_number(Tz,Tn), Tn > 8 }
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --rins $(I) --time $(T) --seedtime --simulate --trials 10000 --initlen 10000 --maxlen 100 >$@) >&$@.log
	grep -v "Beginning simulation trial" $@.log >$@.log.short
	rm $@.log

g$(G)-m$(M)-r$(R)-i$(I)-t$(T)-traj:
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --rins $(I) --time $(T) --maxlen 30 >$@) >&$@.log
	grep -v "Likelihood for" $@.log >$@.log.short
	rm $@.log

g$(G)-m$(M)-r$(R)-i$(I)-t$(T)-cim:
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --rins $(I) --time $(T) --cim --maxlen 100 --dt .0001 >$@) >&$@.log

g$(G)-m$(M)-r$(R)-i$(I)-t$(T)-moments:
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --rins $(I) --time $(T) --moments --maxlen 100 --dt .0001 >$@) >&$@.log

shorten-logs: $(addsuffix .short,$(wildcard *-traj.log) $(wildcard *-sim.log))

%-traj.log.short: %-traj.log
	grep -v "Likelihood for" $< >$@

%-sim.log.short: %-sim.log
	grep -v "Beginning simulation trial" $< >$@


.SECONDARY:
