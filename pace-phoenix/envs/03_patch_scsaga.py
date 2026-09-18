#!/usr/bin/env python
"""Clone the AluruLab scSAGA checkout and apply the SAVE-H patch.

WHY THIS IS NECESSARY
---------------------
Upstream scSAGA (github.com/AluruLab/scSAGA, and its fork) computes the 30-dim
joint embedding inside `Saga.run_multi()` and then DISCARDS it -- the released
code never writes `joint_embedding_H.npy`. The reverse-imputeKNN step in this
pipeline is built entirely on that file, so the published runs used a locally
modified copy of `scmint/scsaga.py` that persists it.

The patch is exactly the one used for the published results: after the existing
`save_joint_embedding_plot(...)` call, write

    aligned_<name>.npy / aligned_<name>_barcodes.txt   for every dataset
    joint_embedding_H.npy                              (anchor block first)
    joint_embedding_barcodes.txt

The patched module is installed under BOTH names, so either import works:
    from scmint.scsaga import main            # patched in place
    from scmint.scsaga_saveH import main      # copy kept for compatibility

Usage:
  python envs/03_patch_scsaga.py                    # clone + patch into $SCAGA_REPO
  SCAGA_REPO=~/tools/scSAGA python envs/03_patch_scsaga.py
  python envs/03_patch_scsaga.py --check            # verify, do not modify

Exit code 0 only if the final checkout contains a working SAVE-H patch.
"""
import argparse
import hashlib
import os
import shutil
import subprocess
import sys

UPSTREAM = 'https://github.com/AluruLab/scSAGA.git'
ANCHOR_CALL = '    save_joint_embedding_plot(aligned_by_name=aligned, outdir=outdir)\n'

PATCH_BODY = '''    # --- [SAVE-H PATCH] Persist the full joint embedding H directly.
    # --- run_multi() computed the 30-dim aligned embeddings and the original code
    # --- discarded them here, which is why upstream needs a separate
    # --- extract_embedding.py step.  This pipeline's reverse-imputeKNN reads
    # --- joint_embedding_H.npy, so we write it out.
    anchor_bc = scdata.barcodes_by_name.get(scdata.anchor)
    for name, emb in aligned.items():
        np.save(os.path.join(outdir, f"aligned_{name}.npy"), emb)
        bc = scdata.barcodes_by_name.get(name)
        if bc is not None:
            with open(os.path.join(outdir, f"aligned_{name}_barcodes.txt"), "w") as f:
                f.write("\\n".join(bc) + "\\n")
    H = np.vstack([aligned[scdata.anchor]] +
                  [aligned[n] for n in aligned if n != scdata.anchor]).astype(np.float32)
    np.save(os.path.join(outdir, "joint_embedding_H.npy"), H)
    if anchor_bc is not None:
        parts = [anchor_bc] + [scdata.barcodes_by_name[n]
                               for n in aligned if n != scdata.anchor]
        with open(os.path.join(outdir, "joint_embedding_barcodes.txt"), "w") as f:
            f.write("\\n".join(np.concatenate(parts)) + "\\n")
    print(f"[SAVE-H PATCH] wrote aligned_*.npy + joint_embedding_H.npy ({H.shape})")
'''

MARKER = 'SAVE-H PATCH'


def run(cmd, cwd=None):
    print('  $', ' '.join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def git_head(repo):
    try:
        out = subprocess.run(['git', '-C', repo, 'rev-parse', 'HEAD'],
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return 'unknown'


def already_patched(repo):
    p = os.path.join(repo, 'scmint', 'scsaga.py')
    if not os.path.exists(p):
        return False
    return MARKER in open(p).read()


def check(repo):
    """Return (ok, message)."""
    p = os.path.join(repo, 'scmint', 'scsaga.py')
    if not os.path.exists(p):
        return False, f'{p} not found -- scSAGA not cloned here'
    src = open(p).read()
    if MARKER not in src:
        return False, ('SAVE-H patch NOT applied: scmint/scsaga.py never writes '
                       'joint_embedding_H.npy')
    for needle in ['joint_embedding_H.npy', 'aligned_{name}.npy',
                   'joint_embedding_barcodes.txt']:
        if needle not in src:
            return False, f'patch incomplete: {needle} missing'
    # must still compile
    try:
        import py_compile
        import tempfile
        with tempfile.NamedTemporaryFile(suffix='.pyc', delete=True) as tf:
            py_compile.compile(p, cfile=tf.name, doraise=True)
    except Exception as e:
        return False, f'scsaga.py does not compile after patching: {e}'
    if not os.path.exists(os.path.join(repo, 'scmint', 'scsaga_saveH.py')):
        return False, 'scmint/scsaga_saveH.py (compatibility module) missing'
    return True, f'SAVE-H patch present and compiling (HEAD {git_head(repo)})'


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.environ.get('PACE_ROOT', os.path.dirname(here))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--repo', default=os.environ.get('SCAGA_REPO', f'{root}/tools/scSAGA'))
    ap.add_argument('--check', action='store_true', help='verify only, do not modify')
    ap.add_argument('--force', action='store_true',
                    help='re-clone from scratch even if the checkout exists')
    args = ap.parse_args()
    repo = os.path.abspath(os.path.expanduser(args.repo))

    if args.check:
        ok, msg = check(repo)
        print(('OK   ' if ok else 'FAIL ') + msg)
        return 0 if ok else 1

    if args.force and os.path.isdir(repo):
        print(f'removing {repo}')
        shutil.rmtree(repo)

    os.makedirs(os.path.dirname(repo), exist_ok=True)
    if not os.path.isdir(os.path.join(repo, '.git')):
        print(f'=== cloning {UPSTREAM} -> {repo}')
        run(['git', 'clone', UPSTREAM, repo])
    else:
        print(f'=== reusing existing checkout {repo} (HEAD {git_head(repo)})')

    if already_patched(repo):
        print('=== already patched -- nothing to do')
    else:
        target = os.path.join(repo, 'scmint', 'scsaga.py')
        src = open(target).read()
        n = src.count(ANCHOR_CALL)
        if n != 1:
            print(f'ERROR: anchor line found {n} times in {target}; upstream may '
                  f'have changed. Patch by hand using PATCH_BODY in this script.',
                  file=sys.stderr)
            return 2
        open(target, 'w').write(src.replace(ANCHOR_CALL, ANCHOR_CALL + PATCH_BODY))
        print('=== applied SAVE-H patch to scmint/scsaga.py')

        # compatibility copy: pipelines may import scmint.scsaga_saveH
        dest = os.path.join(repo, 'scmint', 'scsaga_saveH.py')
        shutil.copyfile(target, dest)
        print(f'=== wrote compatibility module {dest}')

    ok, msg = check(repo)
    print(('OK   ' if ok else 'FAIL ') + msg)
    if ok:
        print(f'\nscSAGA ready at {repo}. Export SCAGA_REPO={repo} before running '
              f'the integration step.')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
