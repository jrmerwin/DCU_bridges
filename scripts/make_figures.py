#!/usr/bin/env python3
"""Regenerate every paper figure from frozen result files. No scientific fits.

Each figure is a separate matplotlib plot. PDF/SVG are vector outputs; PNG is a
convenience preview. Color defaults are left to matplotlib. SOURCE_MAP.json
records the exact inputs and transformations for each figure.
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from collections import Counter
from fractions import Fraction
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from reproduce import ROOT, saved, load_json

plt.rcParams.update({'font.size':10, 'axes.labelsize':10, 'axes.titlesize':11,
                     'legend.fontsize':9, 'pdf.fonttype':42, 'ps.fonttype':42,
                     'svg.fonttype':'none', 'svg.hashsalt':'DCU_BRIDGES_0_1', 'savefig.dpi':180})

def number(x):
    if isinstance(x,dict) and 'numerator' in x:return float(Fraction(x['numerator'],x['denominator']))
    if isinstance(x,str):return float(Fraction(x))
    return float(x)

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,default=ROOT/'paper/figures');args=ap.parse_args()
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=True);mapping={}
    def finish(fig,name,sources,transformation):
        fig.tight_layout(pad=1.0)
        for ext in ('pdf','svg','png'):
            kw={'metadata':{'CreationDate':None,'ModDate':None}} if ext=='pdf' else ({'metadata':{'Date':None}} if ext=='svg' else {})
            fig.savefig(out/f'{name}.{ext}',bbox_inches='tight',**kw)
        plt.close(fig)
        mapping[name]={'sources':sources,'transformation':transformation}

    # Figure 1: keep dependent transfers and different approximation levels visible.
    a,b,c=saved(17),saved(18),saved(19)
    labels=[r'$\tau\to e$ partial rate',r'$\tau\to\mu$ partial rate','Pion-decay recoil','H–muonium optical shift','Pion ratio: leading order','Pion ratio: + pointlike QED']
    fig,ax=plt.subplots(figsize=(6.7,3.9))
    for recipe,mark,offset,title in [('structural_control','o',-.11,'Structural-only control'),('screened_primary','s',.11,'Screened prescription')]:
        vals=[100*next(r['fractional_error'] for r in a['comparisons'] if r['recipe']==recipe and r['channel']==ch) for ch in ('electron','muon')]
        vals += [100*number(next(r['fractional_error'] for r in b['comparisons'] if r['recipe']==recipe and r['quantity']==q)) for q in ('muon_recoil_MeV_c','H_minus_Mu_1S2S_Hz','pion_e_mu_ratio_LO')]
        vals += [100*c['pion_comparisons'][recipe]['corrected_error']]
        ax.plot(vals,np.arange(6)+offset,mark,linestyle='none',label=title)
    ax.axvline(0,linestyle=':',linewidth=1)
    ax.set_yticks(np.arange(6),labels);ax.invert_yaxis();ax.set_xlabel('Signed residual relative to the saved benchmark (%)')
    ax.set_xlim(-2.15,6.35);ax.grid(axis='x',alpha=.22);ax.legend(loc='center right',bbox_to_anchor=(1,.56),frameon=False)
    finish(fig,'horizontal_residuals',[17,18,19],'100 * saved fractional residuals; no combination into a score; no uncertainty bars invented.')

    # Figure 2: released profile statistic at one fixed source amplitude.
    a=saved(22);fig,ax=plt.subplots(figsize=(6.7,3.65))
    curve=a['curve'];ck=next(k for k in curve if k!='R');print('Daya Bay curve field:',ck)
    ax.plot([number(v) for v in curve['R']],[number(v) for v in curve[ck]],label='Released-surface interpolation')
    for r,mark in [(32,'o'),(33,'*'),(34,'s')]:
        y=number(a['likelihood_tests'][str(r)]['delta_chi2'])
        ax.plot([r],[y],marker=mark,markersize=10 if r==33 else 6,linestyle='none',label=f'R = {r}: {y:.3f}')
    ax.axhline(2.30,linestyle='--',linewidth=1,label='Two-parameter reference: 2.30')
    ax.axhline(6.18,linestyle=':',linewidth=1,label='Two-parameter reference: 6.18')
    ax.set(xlim=(30,36),ylim=(0,12),xlabel=r'$R=\Delta m_{32}^{2}/\Delta m_{21}^{2}$ (solar scale fixed)',ylabel=r'Released $\Delta\chi^{2}$')
    ax.legend(loc='lower center',bbox_to_anchor=(.5,1.01),frameon=False,ncol=2,fontsize=8);ax.grid(alpha=.2)
    finish(fig,'neutrino_profile',[22],'Fixed-amplitude slice of released NO surface; primary bilinear values, no new minimum subtraction or fit.')

    # Figure 3: independent per-shot model values, but dependent states/worlds as captioned.
    a=load_json(ROOT/'experiments/cell35/replication/replication_results.json')
    fig,ax=plt.subplots(figsize=(6.7,3.6));ages=[128,4194304]
    for h,mode,marker,ls in [(0,'global8192','o','-'),(6,'global8192','s','-'),(0,'local32','^','--'),(6,'local32','D','--')]:
        mean=[];se=[]
        for age in ages:
            vs=[r['correct_word']['27'] for task in a['tasks'] if task['H']==h and task['epoch']==age for r in task['readouts'] if r['mode']==mode and r['complete']]
            assert len(vs)==128
            mean.append(np.mean(vs));se.append(np.std(vs,ddof=1)/np.sqrt(len(vs)))
        label=f'H={h}: '+('8,192 global ticks' if mode=='global8192' else '32 support services')
        ax.errorbar(ages,mean,yerr=se,marker=marker,linestyle=ls,capsize=3,label=label)
    ax.set_xscale('log');ax.set_xticks(ages,['128','4,194,304'])
    ax.set(xlabel='Maintenance age at preparation',ylabel='Probability of correct decoding',ylim=(0.2,1.05),xlim=(80,6700000))
    ax.legend(loc='upper left',frameon=False,fontsize=8.5,ncol=1);ax.grid(alpha=.2)
    finish(fig,'epoch_activity',['experiments/cell35/replication/replication_results.json'],'Per-arm mean and sample SE of 128 expected readouts, conditional on four worlds. Points only at two saved ages; no interpolation claim.')

    # Figure 4: entropy bias is retained, not set to the theoretical value.
    a=saved(38);rows=a['panels']['no_relief']['lengths']['3']['finite_sampling']
    fig,ax=plt.subplots(figsize=(6.7,3.4))
    x=np.arange(4);ax.plot(x,x*np.log2(3),linestyle='--',label='Exact shared-register law')
    for count,marker in [(4096,'o'),(65536,'s')]:
        rr=sorted((r for r in rows if r['samples']==count),key=lambda r:r['shared'])
        ax.plot([r['shared'] for r in rr],[r['plugin_MI_bits'] for r in rr],marker,markersize=6,linestyle='none',label=f'{count:,} first-write realizations')
    ax.set(xticks=x,xlabel='Shared registers in two three-register paths',ylabel='Mutual information (bits)',ylim=(-.12,5.2))
    ax.legend(frameon=False,loc='upper left');ax.grid(alpha=.2)
    finish(fig,'record_information',[38],'Exact identity and saved plug-in MI estimates for all four H0 length-3 representatives; no new Monte Carlo.')

    # Figure 5: exact conditional MI, reordered but not clustered/fitted.
    a=saved(38)['registry'];order=a['branches'][0]+a['universal_indices']+a['branches'][1]
    assert sorted(order)==list(range(137))
    matrix=np.zeros((137,137))
    for r in a['pairs']:matrix[r['i'],r['j']]=matrix[r['j'],r['i']]=r['conditional_shared_trits']
    matrix=matrix[np.ix_(order,order)];np.fill_diagonal(matrix,np.nan)
    fig,ax=plt.subplots(figsize=(6.3,5.1))
    im=ax.imshow(np.ma.masked_invalid(matrix),vmin=0,vmax=3,interpolation='none',rasterized=False)
    ax.set_xticks([31,68,105],[r'$B_a$ (63)',r'$V_0$ (11)',r'$B_b$ (63)'])
    ax.set_yticks([31,68,105],[r'$B_a$ (63)',r'$V_0$ (11)',r'$B_b$ (63)'])
    for v in (62.5,73.5):ax.axvline(v,linestyle='--',linewidth=.65);ax.axhline(v,linestyle='--',linewidth=.65)
    ax.set_xlabel('Registry root, structurally ordered');ax.set_ylabel('Registry root, structurally ordered')
    bar=fig.colorbar(im,ax=ax,ticks=[0,1,2,3],fraction=.045,pad=.04);bar.set_label(r'$I(Y_i;Y_j\mid g_c)/\log_2 3$ (shared trits)')
    finish(fig,'registry_information',[38],'Exact 9,316 pair table ordered by inherited branches/junction. Diagonal masked (formula scoped to distinct roots). S/I/G labels not inferred from matrix.')

    # Figure 6: quantify synthetic estimation spread, not a fabricated CI.
    a=saved(39);fig,ax=plt.subplots(figsize=(6.3,4.05))
    lim=[.0035,.58];ax.plot(lim,lim,linestyle=':',linewidth=1,label='Exact recovery')
    for q,marker in [('27','o'),('81','s')]:
        rows=[r for r in a['depths'][q]['rows'] if r['label']!='preparation_control']
        x=np.array([r['true_p'] for r in rows]);y=np.array([r['median_estimated_p'] for r in rows]);interval=np.array([r['p_empirical_95percent_range'] for r in rows])
        err=np.stack([y-interval[:,0],interval[:,1]-y])
        ax.errorbar(x,y,yerr=err,fmt=marker,capsize=3,label=f'Phase increment 1/{q}',markersize=5)
    ax.set_xscale('log');ax.set_yscale('log');ax.set(xlim=lim,ylim=lim,xlabel='Hidden native service probability p',ylabel='Estimated p from outcome counts')
    ax.legend(loc='upper left',frameon=False);ax.grid(alpha=.2,which='both')
    finish(fig,'inverse_activity',[39],'Six new prepared states; median and empirical central 95% range across 64 synthetic datasets, not confidence intervals for a single dataset.')

    # Figure 7: every productive event, duplicate coordinates explicitly sized.
    a=saved(39);fig,ax=plt.subplots(figsize=(6.3,3.5))
    for h,marker in [(0,'o'),(6,'s')]:
        counts=Counter((e['new_record_information_trits'],e['cost']) for e in a['native_events'] if e['H']==h and e['born'])
        ax.scatter([x for x,y in counts],[y for x,y in counts],s=[36*n for n in counts.values()],marker=marker,label=f'H={h}',alpha=.8)
        for (x,y),n in counts.items():
            if n>1:ax.annotate(f'×{n}',(x,y),xytext=(7,4),textcoords='offset points',fontsize=8)
    ax.set(xlabel='New first-written trits in the complete burst',ylabel='Native recording obligations',xticks=[0,1,2],xlim=(-.25,2.4),ylim=(0,440));ax.grid(axis='y',alpha=.2);ax.legend(frameon=False)
    finish(fig,'work_information',[39],'All 21 productive events. Marker area counts repeated coordinates within regime. Charges are work obligations, not heat or energy.')
    (out/'SOURCE_MAP.json').write_text(json.dumps(mapping,indent=2)+'\n',encoding='utf-8')
    print('Generated',len(mapping),'figures, each PDF/SVG/PNG, from archived results.')
if __name__=='__main__':main()
