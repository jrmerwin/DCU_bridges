# CELL 28 — portable state-sensitive neutron-search bootstrap.
# Unzip DCU_NEUTRON_SEARCH_v1 beside this notebook, or set PACKAGE_ROOT below.
# Standard library only; NO previous notebook cell is required for execution.
# The bundled original data/model are used, and no prior dictionary is overwritten.
from pathlib import Path
import sys

PACKAGE_ROOT = Path('DCU_NEUTRON_SEARCH_v1').resolve()
if not (PACKAGE_ROOT / 'neutron_search.py').is_file():
    if (Path.cwd() / 'neutron_search.py').is_file():
        PACKAGE_ROOT = Path.cwd()
    else:
        raise FileNotFoundError('Unzip the complete DCU_NEUTRON_SEARCH_v1 folder beside the notebook, '
                                'or edit PACKAGE_ROOT to its location.')
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))
import dcu_search
if Path(dcu_search.__file__).resolve().parent.parent != PACKAGE_ROOT:
    raise RuntimeError('A different search-package version is already loaded; restart the kernel before mixing versions.')
from dcu_search.search import bootstrap
from dcu_search.serialization import dump

dcu_mass_28 = bootstrap()
c = dcu_mass_28['counts']
print('CELL 28 — LOCAL RECORD STATES AND CONSTRAINED STRUCTURAL PARTNERS')
print('No proton/neutron identification, charge calibration, or new mass formula.')
print(f"{c['prior_matched_structures']} inherited supports; {c['state_assignments']:,} complete internal state assignments;")
print(f"{c['local_state_reflections']:,} local-reflection checks; total internal entry-odd current is always zero.")
print('The current is a declared diagnostic, not electric charge. Passive record values do not alter native work.')
print(f"{c['K3_archive']:,} archived K3s; {c['raw']} local shared-parent alternatives;")
print(f"  {c['matched']} match support size/grade; {c['selected']} occur in common archived hosts.")
print('Sign patterns over OLD 32 chains:',dcu_mass_28['patterns']['chain32'])
print('Sign patterns over 137 registry-root validation contexts:',dcu_mass_28['patterns']['registry137'])
print('Eight one-sided (including zeros) chain contrasts reverse sign under the broader validation contexts.')
print(f"All {c['old_Q_checks']:,} old native Q values replay exactly; {c['robust_work_leads']} robust-work alerts.")
print('No numerical mass target participated. A work alert alone would still NOT identify a neutron.')
path = PACKAGE_ROOT/'runs/notebook/Cell28_RESULTS.json.gz'
dump(path,dcu_mass_28)
print('Saved:',path)
print('The separate terminal runner can continue on fresh exact hosts and native histories.')
