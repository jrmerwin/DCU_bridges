# Optional conversion helper: overwrites paper/main.tex and references.bib.
from pathlib import Path
import re,subprocess,json,shutil,zipfile
R=Path(__file__).resolve().parents[1];P=R/'paper'
s=(R/'docs/source/DCU_BRIDGES_MANUSCRIPT_v0_1.md').read_text()
# Preserve the original source verbatim under docs/source. Editorial edition only.
abstract=s.split('## Abstract\n\n',1)[1].split('\n\n**Keywords:',1)[0]
body=s.split('## 1. Introduction',1)[1].split('## References and archived evidence',1)[0]
body='## Introduction'+body
body=re.sub(r'^(#{2,4}) \d+(?:\.\d+)*\.? ',r'\1 ',body,flags=re.M)
# Archive citation labels become proper, linked LaTeX citations. No scientific changes.
def citations(m):
 text=m.group(1)
 if not re.match(r'^[DREN]\d',text):return m.group(0)
 notes='';parts=text.split(';');keys=[]
 for part in parts:
  part=part.strip()
  section=re.search(r',\s*(Section[^;]+)',part)
  if section:notes=section.group(1);part=part[:section.start()]
  for tok in re.split(r',\s*',part):
   ma=re.fullmatch(r'([DREN])(\d+)[–\-]([DREN])?(\d+)',tok)
   if ma:keys.extend(f'{ma[1]}{i}' for i in range(int(ma[2]),int(ma[4])+1))
   elif re.fullmatch(r'[DREN]\d+',tok):keys.append(tok)
   else:raise ValueError(('citation',tok))
 return '\\cite'+('['+notes+']' if notes else '')+'{'+','.join(keys)+'}'
body=re.sub(r'\[([^\]\n]+)\]',citations,body)
# Add cross-references and exact data-driven figures without changing the claims.
figures=[
 ('### Recoil and reduced-mass spectroscopy','horizontal_residuals','horizontal',
  'Signed residuals for the frozen structural-only and screened prescriptions, extracted from Cells 17--19. Rows share mass inputs and sometimes observations; they are not independent confirmations. The two pion-ratio rows differ in their borrowed radiative approximation. No common theoretical uncertainty or combined score is assigned.',.97),
 ('### Thresholds, mixed sectors and the present boundary','neutrino_profile','neutrino',
  'The fixed-amplitude slice through the released Daya Bay 3158-day normal-ordering surface, using the unchanged source solar scale. Points mark the already declared R=32, 33 and 34 cases; R=33 remains primary. Dashed reference levels concern the original two-parameter surface, not a newly fitted one-parameter confidence interval or a probability of the theory. The collaboration solar and nuisance treatment is retained.',.96),
 ('## New information and inverse-measurement bridges','epoch_activity','epoch',
  'Activity-matched persistence of the same early-selected record pair in early and late native states (Cell 35). The phase increment is 1/27. Error bars are sample standard errors for 128 expected readouts in each regime/age/window, conditional on four source worlds. Lines connect the two measured preparation ages only. Global-tick longevity improves while retention after matched support activity does not show the same gain.',.95),
 ('### The registry as a conditional-information graph','record_information','information',
  'Mutual information of the four frozen H0 three-register representatives (Cell 38). Each sample is an independent first-write valuation at a fixed graph, not a new trit drawn at every maintenance service. All four shared-register classes are retained. The positive plug-in estimate at zero shared registers is finite-sample bias; the exact model value is zero.',.95),
 ('### Inferring native activity without exposing hidden counts','registry_information','registry',
  'The exact conditional-information matrix for all distinct pairs of registry roots, in units of shared trits. Roots are ordered by the inherited branch--junction--branch decomposition, not by fitted clustering. The diagonal is excluded. Positive entries recover the 5,347-edge overlap graph; the two 63-root branches share no remaining information after conditioning on the genesis register. This is a dependence structure of the declared record variables, not a spatial map or an independent derivation of sector identities.',.83),
 ('### Recording workload and information are not interchangeable','inverse_activity','inverse',
  'Recovery of the hidden native service probability from finite synthetic outcome counts (Cell 39). Markers are medians; bars show empirical central 95-percent ranges across 64 datasets at each of six new prepared states. They are not confidence intervals from a single measurement. Each dataset contains 32,768 ideal shots at each of nine settings. The smaller phase increment worsens sensitivity at the same shot budget.',.91),
 ('## Implementation and identification boundaries','work_information','work',
  'Native recording obligations versus newly written trits for all 21 productive events in the two new continuation windows. Marker area reflects coincident event multiplicity within each regime, with repeated coordinates explicitly annotated. Equal trit increments can incur different work, and positive work can occur without a new trit. No temperature, heat or energy unit is assigned.',.91)
]
for marker,name,label,caption,width in figures:
 block=f'''\n\n\\begin{{figure}}[tbp]
\\centering
\\includegraphics[width={width}\\linewidth]{{figures/{name}.pdf}}
\\caption{{{caption}}}\\label{{fig:{label}}}
\\end{{figure}}\n\n'''
 if marker not in body:raise ValueError(marker)
 body=body.replace(marker,block+marker,1)
# Small inline-math typography repairs to variable tokens in prose, not formulas.
for a,b in [('root r_i,','root \\(r_i\\),'),('Each B_i contains','Each \\(B_i\\) contains'),('to r_i.','to \\(r_i\\).'),('Let Y_i be','Let \\(Y_i\\) be'),('composite c={a,b}.','composite \\(c=\\{a,b\\}\\).'),('If W_new is','If \\(W_{\\mathrm{new}}\\) is'),('Only R=F+n is identified.','Only \\(R=F+n\\) is identified.')]:
 body=body.replace(a,b)
# Brief linking sentences make placement deliberate and references searchable.
links={
 '### Recoil and reduced-mass spectroscopy':'Figure \\ref{fig:horizontal} compares these rate results with the subsequent recoil, optical, and pion-ratio transfers.\n\n',
 '### Thresholds, mixed sectors and the present boundary':'Figure \\ref{fig:neutrino} displays the frozen spectral-surface comparison.\n\n',
 '## New information and inverse-measurement bridges':'Figure \\ref{fig:epoch} contrasts the global-age and local-activity observation windows.\n\n',
 '### The registry as a conditional-information graph':'Figure \\ref{fig:information} retains the exact sharing law and its finite-sample estimates.\n\n',
 '### Inferring native activity without exposing hidden counts':'Figure \\ref{fig:registry} displays the complete pairwise conditional-information table.\n\n',
 '### Recording workload and information are not interchangeable':'Figure \\ref{fig:inverse} shows the finite-outcome recovery and its sampling spread.\n\n',
 '## Implementation and identification boundaries':'Figure \\ref{fig:work} keeps all productive events in this work--information comparison.\n\n'
}
for marker,text in links.items():body=body.replace(marker,text+marker,1)
# Formatting only for bare symbols used as variables in Markdown prose.
body=body.replace('Let Y_A denote','Let \\(Y_A\\) denote').replace(' set A and',' set \\(A\\) and')
body=body.replace(' K=a+b is sufficient',' \\(K=a+b\\) is sufficient').replace(' K of two tagged',' \\(K\\) of two tagged')
body=body.replace('of first-written nonprimitive registers','of first-written nonprimitive registers')
# State actual packaging scope separately from the historical manuscript-assembly scope.
body += r'''

## Distribution and outstanding work

This LaTeX edition preserves the scientific results of working draft 0.1. Figures are generated directly from the frozen numerical files; the original Markdown remains in the repository. The compact reproduction command regenerates the horizontal Cells 14 and 17--24, re-evaluates Cells 32--34 from archived histories, executes the instrument calculations in Cells 36--37, reruns the new Cells 38--39, and checks the saved-record analysis of Cell 35. Exact rational fields and discrete records are checked separately from floating-point tolerances. The earlier native-history suite and the million-tick epoch generator are separate opt-in routes. The late neutron-search corpus is not reproduced by this distribution.

The two bounded horizontal additions---pointlike electromagnetic production and beta-endpoint effective mass---remain proposed rather than executed. The historical comparisons keep their original releases and conventions. Final primary-reference reconciliation and author review remain prerequisites for submission; this typeset edition does not convert those pending tasks into completed results.
'''
(P/'manuscript.md').write_text(body,encoding='utf-8')
# Pandoc is a conversion-time convenience; users can edit/compile main.tex directly.
def tex(text):
 p=subprocess.run(['pandoc','--from','markdown+tex_math_single_backslash+raw_tex','--to','latex','--top-level-division=section','--shift-heading-level-by=-1','--wrap=none'],input=text,text=True,capture_output=True,check=True)
 return p.stdout
bodytex=tex(body)
# Give the five source tables a caption and label, keeping every row.
captions=[('Frozen mass prescriptions used by the executed horizontal calculations. Entries are in electron-mass units.','dictionary'),('Exact registry-pair distribution after conditioning on the common genesis register.','registry_pairs'),('Finite-outcome activity estimates for all six new prepared states. Relative RMSE is computed across 64 synthetic datasets per state at the primary phase increment.','inverse'),('All productive events in the new continuation windows, grouped by new trit count.','work'),('Categorical reach of the executed bridges and the limits of each correspondence.','reach')]
idx=0
def tablecap(m):
 global idx
 cap,lab=captions[idx];idx+=1
 return m.group(1)+f'\n\\caption{{{cap}}}\\label{{tab:{lab}}}\\\\\n'+m.group(2)
bodytex=re.sub(r'(\\begin\{longtable\}.*?)(\\toprule)',tablecap,bodytex,flags=re.S)
assert idx==5,idx

# These five short tables fit on one page; use real floats rather than a
# multipage longtable head that can collide with queued figures.
def short_table(m):
    text=m.group(1)
    cap=re.search(r'\\caption\{.*?\}\\label\{tab:[^}]+\}\\\\',text,re.S).group(0)
    spec=text[:text.index('\\caption')].strip()
    spec=re.sub(r'^\\begin\{longtable\}(?:\[[^]]*\])?', '',spec).strip()
    rest=text[text.index(cap)+len(cap):]
    head,tail=rest.split('\\endhead',1)
    tail=tail.split('\\endlastfoot',1)[1]
    if 'tab:inverse' in cap:
        spec=r'{@{}lrrrr@{}}'
        head=r'\toprule State & True $R$ & True $p$ & Median $\widehat p$ & Relative RMSE \\\midrule'
    out=(r'\begin{table}[htbp]\centering'+'\n'+cap[:-2]+'\n'+
         r'\small\setlength{\tabcolsep}{5pt}\renewcommand{\arraystretch}{1.15}'+'\n'+
         r'\begin{tabular}'+spec+'\n'+head+tail+r'\bottomrule'+'\n'+
         r'\end{tabular}'+'\n'+r'\end{table}')
    return out
bodytex=re.sub(r'(\\begin\{longtable\}.*?)(?:\\end\{longtable\})',short_table,bodytex,flags=re.S)
# Make display equations individually referenceable; content untouched.
e=0
def equations(m):
 global e
 e+=1
 return '\\begin{equation}\\label{eq:'+str(e)+'}\n'+m.group(1)+'\n\\end{equation}'
bodytex=re.sub(r'\\\[\n?(.*?)\n?\\\]',equations,bodytex,flags=re.S)
preamble=r'''% LaTeX/figure edition of the supplied working manuscript 0.1.
% Scientific source: docs/source/DCU_BRIDGES_MANUSCRIPT_v0_1.md.
% Compile from paper/: latexmk -lualatex -interaction=nonstopmode -halt-on-error main.tex
\documentclass[11pt,a4paper]{article}
\usepackage[margin=25mm,headheight=15pt]{geometry}
\usepackage{fontspec}
\setmainfont{Latin Modern Roman}
\setsansfont{Latin Modern Sans}
\setmonofont[Scale=MatchLowercase]{Latin Modern Mono}
\usepackage{amsmath,amssymb,amsthm,mathtools}
\usepackage{microtype}
\usepackage{graphicx,booktabs,longtable,array,calc}
\usepackage{caption}
\captionsetup{font=small,labelfont=bf,skip=7pt}
\usepackage{fancyhdr}
\usepackage{enumitem}
\usepackage{etoolbox}
\usepackage{url}
\usepackage{xurl}
\usepackage[hidelinks,unicode]{hyperref}
\usepackage[backend=biber,style=numeric,sorting=none,maxnames=3,giveninits=true,url=true,doi=true]{biblatex}
\addbibresource{references.bib}
\setlength{\bibitemsep}{4pt}
\renewcommand*{\bibfont}{\small}
\setlength{\emergencystretch}{2em}
\setlength{\parskip}{3pt}
\setlength{\parindent}{1em}
\setlength{\textfloatsep}{17pt plus 3pt minus 3pt}
\renewcommand{\topfraction}{0.88}
\renewcommand{\bottomfraction}{0.75}
\renewcommand{\textfraction}{0.10}
\renewcommand{\floatpagefraction}{0.72}
\setcounter{topnumber}{2}
\providecommand{\tightlist}{\setlength{\itemsep}{0pt}\setlength{\parskip}{0pt}}
\AtBeginEnvironment{longtable}{\small\setlength{\tabcolsep}{5pt}\renewcommand{\arraystretch}{1.14}}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small Registry structure and operational bridges}
\fancyhead[R]{\small J. Merwin}
\fancyfoot[C]{\thepage}
\renewcommand{\headrulewidth}{0.3pt}
\title{\vspace{-1.1em}\textbf{From Registry Structure to Physical Observables}\\[0.45em]\large Horizontal Transfers and Operational Bridges in the\\ Distinction Combinatorial Universe}
\author{Jason Merwin\\\small Independent Researcher}
\date{\small September 2026\quad\textbullet\quad Working manuscript 0.1, LaTeX edition}
\hypersetup{pdftitle={From Registry Structure to Physical Observables},pdfauthor={Jason Merwin},pdfsubject={Horizontal transfers and operational bridges in the DCU}}
\begin{document}
\maketitle
\thispagestyle{plain}
\begin{abstract}
'''
status=r'''
\end{abstract}
\noindent\textbf{Keywords:} distinction registry; operational observables; parameter transfer; shared records; quantum channels; information; finite service.

\medskip
\noindent\begin{minipage}{\linewidth}\small\textit{Working-draft status.} This is a typeset and illustrated edition of the integrated first draft, not a submission-ready paper. No new bridge experiment is added here. The two proposed horizontal additions and final primary-reference reconciliation remain pending. Reproduction coverage is stated at the end of the paper and in the accompanying repository.\end{minipage}
\medskip
'''
(P/'main.tex').write_text(preamble+tex(abstract)+status+bodytex+'\n\\printbibliography[title={References and archived evidence}]\n\\end{document}\n',encoding='utf-8')
# Bibliography preserves supplied archive identities, without inventing publication DOI/URLs.
entries=[]
def bib(kind,key,**fields):
 entries.append('@'+kind+'{'+key+',\n'+',\n'.join('  '+k+' = {'+str(v)+'}' for k,v in fields.items())+'\n}\n')
D={
'D1':('Structure, Symmetry, and Motif Embedding of the 137-Object Distinction Registry','Working manuscript 0.1, September 2026; structural predecessor and its original reproducibility archive'),
'D2':('Self-Limited Growth, Intrinsic Time, and Record Geometry in the Distinction Combinatorial Universe','Revised manuscript 0.5, September 2026; supplied native-law source'),
'D3':('Relational Mathematical Realism: Registry Architecture Predicts Lepton, Baryon, and Strange Baryon Mass Spectra','March 2026; interpreted with the exact executed-branch reconciliation in this paper'),
'D4':('Thermodynamic Screening Corrections in Relational Mathematical Realism: Mass Predictions from Graph Topology and a Universal Screening Unit','March 2026; printed rounded bases and precision claims are not substituted for executed results'),
'D5':('Structural Taxonomy of Fundamental Constants from a 137-Bit Relational Registry','Unpublished manuscript, February 2026; selected neutrino branch. Integer-grammar significance claims are not adopted')}
for key,(title,note) in D.items():bib('unpublished',key,author=('Merwin, Jason R.' if key=='D4' else 'Merwin, Jason'),title=title,year='2026',note=note)
bib('article','E1',author='Shannon, Claude E.',title='A Mathematical Theory of Communication',journaltitle='Bell System Technical Journal',volume='27',pages='379--423 and 623--656',year='1948',doi='10.1002/j.1538-7305.1948.tb01338.x',note='Second installment: DOI 10.1002/j.1538-7305.1948.tb00917.x')
titles={17:'Frozen-mass weak partial-rate transfer',18:'Recoil, optical interval, and leading-order pion branching',19:'Pion radiative correction, mixing angles, and hyperon limitations',20:'Coherent vacuum electron-neutrino survival',21:'KamLAND 2005 approximate detector-folded transfer',22:'Daya Bay released normal-ordering spectral surface',23:'Production thresholds, operator ambiguity, and finite character algebra',24:'Finite record phases and energy-map selection',25:'Native preparation and service-driven phases',26:'Record-supported clocks and exact overlap response',27:'Native paired-response audit',28:'Native neutron-search scope',29:'Native recording back-action',30:'Local records and communication',31:'Moving-carrier identity and clock limitations',32:'Joint-history quantum channel and filtered-state transfer',33:'Fixed echo and temporal-information controls',34:'Protected-subspace classification',35:'Activity-matched epochs and collective four-path protection',36:'Measurement state updates and coherent probe',37:'Nondemolition obstruction and code-restricted renewal',38:'Record information and registry conditional-dependence graph',39:'Inverse activity measurement and work--information audit'}
for n,title in titles.items():
 bib('misc',f'R{n}',author='Merwin, Jason',title=title+f' (Cell {n})',year='2026',howpublished=r'Executable research record: \path{code/cells/DCU_Mass_Cell_'+str(n)+'.py}',note=r'Frozen outputs in \path{data/reference/}; methods and provenance in \path{docs/notes/}. '+('Collaboration source: arXiv:2211.14988; Physical Review Letters 130, 161802 (2023). The ancillary-column convention and partial transcription are retained.' if n==22 else ''))
for key,title,name in [('N1','Fixed K3 trial-reference screen','SCREEN_SUMMARY.md'),('N2','Contained-K3 successor investigation','SUCCESSOR_REVIEW(1).md')]:
 bib('misc',key,author='Merwin, Jason',title=title,year='2026',note=r'Archived saved-data report: \path{docs/notes/'+name+'}. Structural-work proxies, not particle identifications.')
(P/'references.bib').write_text('\n'.join(entries),encoding='utf-8')
print('LaTeX source:',len(bodytex),'characters,',e,'numbered displays;',idx,'tables; 7 figures')
