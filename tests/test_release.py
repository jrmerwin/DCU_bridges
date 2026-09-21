"""Independent structural/arithmetic checks; no external data or network."""
from pathlib import Path
import itertools, json, math, re, sys, unittest
from collections import Counter
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from reproduce import load_json,saved,definitions

class ReleaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.info=saved(38);cls.inverse=saved(39)
        env={'__name__':'independent_test'}
        definitions(ROOT/'code/registry/DCU_Rosetta_Calibration_F_Registry_Anatomy.py',env)
        cls.registry=env['rosetta_dag7_reference']()

    def test_registry_incidence_and_counts(self):
        r=self.registry;v=self.info['registry'];self.assertEqual(len(r['parents']),173)
        self.assertEqual(Counter(x['sector'] for x in r['rows']),Counter(S=81,I=40,G=16))
        hist=Counter();edges=0
        for row in v['pairs']:
            i,j=row['i'],row['j'];a=r['ancestors'][r['rows'][i]['object_id']];b=r['ancestors'][r['rows'][j]['object_id']]
            c=len(a&b)-3;hist[c]+=1;edges+=c>0
            self.assertEqual(c,row['conditional_shared_trits'])
            self.assertEqual(c>0,j in r['adjacency'][i])
        self.assertEqual(hist,{0:3969,1:3774,2:1323,3:250});self.assertEqual(edges,5347)
        self.assertEqual(Counter(map(len,r['adjacency'])),{73:126,136:11})

    def test_register_bundle_aliases(self):
        r=self.registry;sets=[frozenset(r['ancestors'][x['object_id']]-{0,1,x['object_id']}) for x in r['rows']]
        self.assertEqual(Counter(Counter(sets).values()),{5:25,1:12})
        self.assertTrue(all(len(x)==4 for x in sets))
        v=self.info['registry'];branches=v['branches']
        unions=[set().union(*(sets[i] for i in b)) for b in branches]
        self.assertEqual(list(map(len,unions)),[17,17]);self.assertEqual(unions[0]&unions[1],{r['pairs'][(0,1)]})

    def test_conditional_information_from_finite_field_rank(self):
        count=0
        # Fixed entry reflections are affine/invertible. The rank of the selected
        # coordinate projection equals the number of distinct registers.
        for env in self.info['panels'].values():
            for panel in env['lengths'].values():
                supports=[{z for z,e in path[1]} for path in panel['paths']]
                for row in panel['selected_triples']:
                    a,b,c=(supports[i] for i in row['indices'])
                    rank_entropy=len(a|c)+len(b|c)-len(c)-len(a|b|c)
                    self.assertAlmostEqual(row['CMI_bits'],rank_entropy*math.log2(3),places=11);count+=1
        self.assertEqual(count,46)

    def test_one_step_sampler_exact_enumeration(self):
        for state in self.inverse['states']:
            n,F,R=state['n'],state['F'],state['R'];tot=math.comb(R,3);weights=[0,0,0]
            for k in range(4):
                if not 0<=3-k<=F:continue
                for subset in itertools.combinations(range(n),k):
                    hit=int(3 in subset)+int(4 in subset);weights[hit]+=math.comb(F,3-k)
            self.assertEqual(sum(weights),tot)
            p=3/R;p2=6/(R*(R-1));truth=[1-2*p+p2,2*(p-p2),p2]
            np.testing.assert_allclose(np.array(weights)/tot,truth,rtol=1e-13,atol=1e-15)

    def test_outcome_counts_and_estimation_statistics(self):
        for depth in self.inverse['depths'].values():
            for row in depth['rows']:
                counts=np.asarray(row['counts']);self.assertEqual(counts.shape,(64,9,3))
                self.assertTrue(np.all(counts.sum(axis=2)==32768))
                estimates=np.asarray(row['estimates_p']);p=row['true_p']
                self.assertAlmostEqual(float(np.median(estimates)),row['median_estimated_p'],places=13)
                self.assertAlmostEqual(float(np.sqrt(np.mean((estimates-p)**2))),row['RMSE_p'],places=13)
                self.assertTrue(np.all((estimates>=0)&(estimates<=.6)))

    def test_native_events_independent_charge_reconstruction(self):
        # Reconstruct exact ancestry and first/repeat charges, without the Cell39
        # simulator implementation. Every enabled absent pair must be present.
        productive=0;zero=0
        for H in (0,6):
            parents=[None,None,(0,1),(0,2),(1,2)];anc=[{0},{1},{0,1,2},{0,1,2,3},{0,1,2,4}]
            ch=[1,1,2,3,3];recorded={0,1,2};pairset={(0,1),(0,2),(1,2)}
            for ev in (e for e in self.inverse['native_events'] if e['H']==H):
                self.assertEqual(len(parents),ev['n'])
                needed=[tuple(p) for p in itertools.combinations(sorted(ev['served']),2) if tuple(p) not in pairset]
                self.assertEqual(needed,[tuple(x) for x in ev['pairs']])
                h=Counter(z for a,b in needed for z in anc[a]|anc[b]);new=set(h)-recorded
                work=2*sum((2*ch[z]-2)*m for z,m in h.items())+9*sum(2*ch[z]-2 for z in new)
                self.assertEqual(work,ev['cost']);self.assertEqual(sorted(new-{0,1}),ev['newly_recorded_nonprimitive'])
                if needed:productive+=1;zero+=bool(work and not(new-{0,1}))
                recorded.update(h)
                for a,b in needed:
                    z=len(parents);parents.append((a,b));anc.append(anc[a]|anc[b]|{z});ch.append(ch[a]+ch[b]);pairset.add((a,b))
        self.assertEqual(productive,21);self.assertEqual(zero,3)

    def test_figures_and_citations_resolve(self):
        text=(ROOT/'paper/main.tex').read_text();bib=(ROOT/'paper/references.bib').read_text()
        keys=set(re.findall(r'@\w+\{([^,]+),',bib));cited=set()
        for group in re.findall(r'\\cite(?:\[[^]]*\])?\{([^}]+)\}',text):cited.update(group.split(','))
        self.assertTrue(cited<=keys,cited-keys)
        figs=re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}',text)
        self.assertEqual(len(figs),7)
        for f in figs:self.assertTrue((ROOT/'paper'/f).is_file(),f)
        self.assertEqual(text.count(r'\begin{table}'),5)
        labels=re.findall(r'\\label\{([^}]+)\}',text);self.assertEqual(len(labels),len(set(labels)))

if __name__=='__main__':unittest.main()
