
TL = ../trajeclike -v2
RELENT = ../relent.js

all: relents

relents: # scan-gammas
	cat *results | perl -e 'print"gamma mu r t method D\n";while(<>){@m=/^D g(\S+)-m(\S+)-r(\S+)-t(\S+)-sim \S+-(\S+) (\S+)/;if(@m){print"@m\n"}}' >$@

stats: # scan-gammas
	cat *results | perl -e 'print"stat gamma mu r t method value\n";while(<>){@m=/^([^D]\S*) g(\S+)-m(\S+)-r(\S+)-t(\S+)-(\S+) (\S+)/;if(@m){print"@m\n"}}' >$@

scan-gammas: g1-scan-mu g.99-scan-mu g.9-scan-mu

$(PREFIX)-scan-mu: $(PREFIX)-m.05-scan-r $(PREFIX)-m.5-scan-r

$(PREFIX)-scan-r: $(PREFIX)-r.5-scan-time $(PREFIX)-r.75-scan-time $(PREFIX)-r.9-scan-time

$(PREFIX)-scan-time: $(PREFIX)-t.015625-results $(PREFIX)-t.03125-results $(PREFIX)-t.03125-results $(PREFIX)-t.0625-results $(PREFIX)-t.125-results $(PREFIX)-t.25-results $(PREFIX)-t.5-results $(PREFIX)-t1-results $(PREFIX)-t2-results $(PREFIX)-t4-results $(PREFIX)-t8-results

g1-$(PREFIX)-results: g1-$(PREFIX)-sim g1-$(PREFIX)-moments g1-$(PREFIX)-traj g1-$(PREFIX)-cim
	$(RELENT) --normalize $^ >$@

g.9$(PREFIX)-results: g.9$(PREFIX)-sim g.9$(PREFIX)-moments g.9$(PREFIX)-traj
	$(RELENT) --normalize $^ >$@

g$(G)-m$(M)-r$(R)-t$(T)-sim:
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --time $(T) --seedtime --simulate --maxlen 100 >$@) >&$@.log

g$(G)-m$(M)-r$(R)-t$(T)-traj:
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --time $(T) --maxlen 30 >$@) >&$@.log

g$(G)-m$(M)-r$(R)-t$(T)-cim:
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --time $(T) --cim --maxlen 100 >$@) >&$@.log

g$(G)-m$(M)-r$(R)-t$(T)-moments:
	(time $(TL) --gamma $(G) --mu $(M) --r $(R) --time $(T) --moments --maxlen 100 >$@) >&$@.log

