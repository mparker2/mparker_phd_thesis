import sys
import os
from glob import glob
from operator import itemgetter

af_path = os.path.join(
    os.path.split(os.path.realpath(__file__))[0],
    'acronym-finder'
)
print(af_path)
sys.path.append(af_path)
import acronym_finding as acro

if __name__ == '__main__':
    all_acros = {}
    for fn in glob('*/text.md'):
        with open(fn) as f:
            text = f.read()
            all_acros.update(acro.find_all_acronyms(text))
    with open('acronyms_unlabelled.csv', 'w') as f:
        for ac, val in sorted(all_acros.items(), key=itemgetter(0)):
            f.write('{},{}\n'.format(ac, val[0]))