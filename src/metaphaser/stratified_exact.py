"""Exact whole-read test from stratified analysis v1.0.0. See docs/PROVENANCE.json."""
import bisect
import collections
import itertools
import math
from fractions import Fraction

class ExactLimit(RuntimeError):pass

def check(ok,msg):
 if not ok:raise ValueError(msg)

def choose_distribution(values,k,limit):
 counts=tuple(sorted(collections.Counter(values).items()))
 if len(counts)==1:return {counts[0][0]*k:math.comb(len(values),k)}
 dp={(0,0):1};remaining=len(values)
 for value,freq in counts:
  remaining-=freq;nxt=collections.defaultdict(int)
  for (used,score),weight in dp.items():
   lo=max(0,k-used-remaining);hi=min(freq,k-used)
   for take in range(lo,hi+1):nxt[used+take,score+take*value]+=weight*math.comb(freq,take)
  if len(nxt)>limit:raise ExactLimit(f'Exact subset-sum distribution needs >{limit} states')
  dp=nxt
 result={score:weight for (used,score),weight in dp.items() if used==k}
 check(sum(result.values())==math.comb(len(values),k),'Exact distribution does not normalize')
 return result

def stratified(strata,limit=2000000):
 # strata contain integer high-call counts on one fixed per-stratum position set.
 weights=[];observed=[];denominators=[]
 for s in strata:
  a,b=s['baseline_counts'],s['state_counts'];na,nb=len(a),len(b);D=s['shared_positions'];N=na+nb;T=sum(a)+sum(b)
  check(na>0 and nb>0 and D>0,'Undefined stratum')
  denominator=na*nb*D
  observed.append(Fraction(na*sum(b)-nb*sum(a),denominator))
  denominators.append(denominator)
 K=len(strata);effect=sum(observed)/K
 total_space=math.prod(math.comb(len(s['baseline_counts'])+len(s['state_counts']),len(s['baseline_counts'])) for s in strata)
 if effect==0:return dict(p_value=1.0,p_numerator=total_space,p_denominator=total_space,delta_pp=0.0)
 scale=math.lcm(*denominators)
 for s,den in zip(strata,denominators):
  a,b=s['baseline_counts'],s['state_counts'];na,nb=len(a),len(b);N=na+nb;T=sum(a)+sum(b)
  take_baseline=na<=nb;k=min(na,nb)
  raw=choose_distribution(a+b,k,limit)
  scaled={((na*T-N*v) if take_baseline else (N*v-nb*T))*(scale//den):w for v,w in raw.items()}
  weights.append(scaled)
 cutoff=abs(sum(observed)*scale);check(cutoff.denominator==1,'Noninteger cutoff');cutoff=int(cutoff)
 if K==1:
  extreme=sum(w for x,w in weights[0].items() if abs(x)>=cutoff)
 elif K==2:
  v=sorted(weights[1]);prefix=[0]
  for x in v:prefix.append(prefix[-1]+weights[1][x])
  interior=0
  for x,w in weights[0].items():
   left=bisect.bisect_right(v,-cutoff-x);right=bisect.bisect_left(v,cutoff-x)
   interior+=w*(prefix[right]-prefix[left])
  extreme=total_space-interior
 else:raise ValueError('Only forward and reverse orientations are expected')
 check(0<=extreme<=total_space,'Exact p outside [0,1]')
 return dict(p_value=extreme/total_space,p_numerator=extreme,p_denominator=total_space,delta_pp=float(effect*100))

def brute(strata):
 all_diffs=[];observed=[]
 for s in strata:
  a,b=s['baseline_counts'],s['state_counts'];na,nb=len(a),len(b);D=s['shared_positions'];v=a+b
  observed.append(Fraction(sum(b),nb*D)-Fraction(sum(a),na*D))
  ds=[]
  for ids in itertools.combinations(range(len(v)),na):
   selected=sum(v[i] for i in ids)
   ds.append(Fraction(sum(v)-selected,nb*D)-Fraction(selected,na*D))
  all_diffs.append(ds)
 obs=sum(observed)
 extreme=total=0
 for ds in itertools.product(*all_diffs):extreme+=abs(sum(ds))>=abs(obs);total+=1
 return extreme/total

def controls():
 fixtures=[
  [dict(baseline_counts=[0,1,2],state_counts=[0,0,1],shared_positions=3)],
  [dict(baseline_counts=[18,15],state_counts=[10,10],shared_positions=92),dict(baseline_counts=[38],state_counts=[9,9,10,11],shared_positions=144)],
  [dict(baseline_counts=[0,1],state_counts=[1],shared_positions=4),dict(baseline_counts=[3],state_counts=[0,2],shared_positions=5)],
  [dict(baseline_counts=[0,0],state_counts=[0],shared_positions=1),dict(baseline_counts=[1],state_counts=[0,0],shared_positions=1)]
 ]
 for ss in fixtures:
  a=stratified(ss);check(abs(a['p_value']-brute(ss))<1e-14,'Brute-force control mismatch')
  swapped=[dict(s,baseline_counts=s['state_counts'],state_counts=s['baseline_counts']) for s in ss]
  b=stratified(swapped);check(abs(a['p_value']-b['p_value'])<1e-14 and abs(a['delta_pp']+b['delta_pp'])<1e-12,'State-swap invariance failed')
  reversed_profiles=[dict(s,baseline_counts=list(reversed(s['baseline_counts'])),state_counts=list(reversed(s['state_counts']))) for s in reversed(ss)]
  c=stratified(reversed_profiles);check(a['p_value']==c['p_value'],'Input-order invariance failed')
 check(stratified(fixtures[1])['p_value']==1/30,'NOSH exact regression failed')
 check(abs(stratified([dict(baseline_counts=[0]*20+[1]*2,state_counts=[0]*21,shared_positions=3)])['p_value']-0.4883720930232558)<1e-14,'5mC SETUN regression')
 check(abs(stratified([dict(baseline_counts=[0]*15+[1]*7,state_counts=[0]*21,shared_positions=3)])['p_value']-0.008900911907151272)<1e-14,'5hmC SETUN regression')
 return dict(status='passed',brute_force_cases=len(fixtures),state_swap_invariance=True,input_order_invariance=True,NOSH_1_of_30=True,SETUN_selected_p_checks=2)
