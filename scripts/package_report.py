#!/usr/bin/env python3
"""Bundle the actual report sources and generated figures; verify an isolated build."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import subprocess
import time
import zipfile
from common import save_json
from inventory import GROUPS

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    root=Path(__file__).resolve().parent.parent;os.chdir(str(root))
    def expanded(path):
        return re.sub(r'\\input\{([^}]+)\}',lambda match: expanded(
            Path('report',match.group(1) if match.group(1).endswith('.tex') else match.group(1)+'.tex')),
            Path(path).read_text())
    Path('report/reproduction-complete.tex').write_text(expanded('report/reproduction.tex'))
    paths=sorted(Path('report').glob('*.tex'))+[Path('figures',name+'.pdf') for name in GROUPS]
    assert Path('report/reproduction.pdf').read_bytes()[:4]==b'%PDF'
    paths.append(Path('report/reproduction.pdf'))
    archive=Path('report/latex-source.zip')
    manifest={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    with zipfile.ZipFile(str(archive),'w',compression=zipfile.ZIP_DEFLATED) as z:
        for p in paths:z.write(str(p),str(p))
        z.writestr('README.txt','Independent fault-valving report.\n\nCompile from report/:\n'
            '  latexmk -pdf -interaction=nonstopmode -halt-on-error reproduction-complete.tex\n\n'
            'reproduction-complete.tex is a single expanded LaTeX source.\n'
            'The modular reproduction.tex source and fragments are also included.\n'
            'Numerical panels use data from the independent C++ implementation.\n'
            'Figure 1a includes a newly drawn geometry schematic and assumed friction profiles;\n'
            'the schematic is not a reproduced simulation.\n'
            'Code, panel data, and provenance: https://github.com/bradlipovsky/zhu2020-fault-valving\n')
        z.writestr('source-manifest.json',json.dumps(manifest,indent=2)+'\n')
    verification=dict(archive=str(archive),bytes=archive.stat().st_size,
        sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),isolated_compile=False)
    if args.verify:
        test=root/'.tmp'/('report-package-'+str(int(time.time()*1e9)));test.mkdir(parents=True)
        with zipfile.ZipFile(str(archive)) as z:z.extractall(str(test))
        # Remove the supplied PDF to ensure the verification actually compiles it.
        (test/'report/reproduction.pdf').unlink()
        env=dict(os.environ,TMPDIR=str(root/'.tmp'))
        command=['latexmk','-pdf','-interaction=nonstopmode','-halt-on-error','reproduction-complete.tex']
        with open('data/report_bundle_build.log','w') as log:
            subprocess.run(command,cwd=str(test/'report'),env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        assert (test/'report/reproduction-complete.pdf').read_bytes()[:4]==b'%PDF'
        verification.update(isolated_compile=True,command=command,source_files=manifest)
    save_json('data/report_bundle_verification.json',verification)
    print(str(archive),verification)

if __name__=='__main__':main()
