"""Download the authors' OSF archive; verify every file against OSF SHA-256."""
import concurrent.futures
import hashlib
import json
from pathlib import Path
import requests

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'reference' / 'data'
DEST.mkdir(parents=True, exist_ok=True)

def download(item):
    a = item['attributes']
    path = DEST / a['name']
    expected = a['extra']['hashes']['sha256']
    def digest(p):
        h = hashlib.sha256()
        with p.open('rb') as f:
            for block in iter(lambda: f.read(2**20), b''):
                h.update(block)
        return h.hexdigest()
    if path.exists() and digest(path) == expected:
        return str(path)
    print('Downloading', a['name'], a['size'], flush=True)
    tmp = path.with_suffix(path.suffix + '.part')
    with requests.get(item['links']['download'], stream=True, timeout=(30, 180)) as r:
        r.raise_for_status()
        with tmp.open('wb') as f:
            for block in r.iter_content(2**20):
                f.write(block)
    if digest(tmp) != expected:
        raise ValueError('SHA-256 mismatch: ' + str(path))
    tmp.rename(path)
    print('Verified', a['name'], flush=True)
    return str(path)

if __name__ == '__main__':
    url = 'https://api.osf.io/v2/nodes/9ygrp/files/osfstorage/'
    items = []
    while url:
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        page = r.json()
        items.extend(page['data'])
        url = page['links']['next']
    (ROOT / 'reference' / 'osf_manifest.json').write_text(json.dumps(items, indent=2) + '\n')
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(download, items))
