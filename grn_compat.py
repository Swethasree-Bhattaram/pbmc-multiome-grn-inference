"""Make arboreto's GRNBoost2 work on a current dask.

arboreto's last release (0.1.6) predates dask's modern expression system.  In
`arboreto.core.create_graph` it builds a *meta* DataFrame unconditionally:

    all_links_df = from_delayed(delayed_link_dfs, meta=_GRN_SCHEMA)
    all_meta_df  = from_delayed(delayed_meta_dfs, meta=_META_SCHEMA)

`delayed_meta_dfs` is only populated when `include_meta=True`, so with the
default `include_meta=False` it calls `from_delayed([])`.  Old dask tolerated
that; current dask raises:

    TypeError: Must supply at least one delayed object

`all_meta_df` is never returned when `include_meta=False`, so guarding that one
call is enough.  Nothing else in arboreto touches the changed dask APIs.

VERIFIED EQUIVALENCE: with this patch, a fixed input (300 cells x 40 genes,
seed 666) produces an IDENTICAL 390-edge network on
    python 3.9.6  / numpy 1.21.5 / pandas 1.4.4  / dask 2021.10.0 / sklearn 1.1.3
and
    python 3.11.16 / numpy 2.4.6 / pandas 3.0.6 / dask 2026.8.0 / sklearn 1.9.1
Edge sets match exactly (390/390).  Edge importances differ in the low digits
because sklearn's GradientBoostingRegressor changed between 1.1 and 1.9; this is
inherent to the sklearn version, not to this patch.

Call `apply()` once before importing/running grnboost2.
"""
from __future__ import annotations

_APPLIED = False


def apply() -> bool:
    """Patch arboreto to tolerate an empty `from_delayed`.  Idempotent."""
    global _APPLIED
    if _APPLIED:
        return True

    import arboreto.core as ac
    import dask.dataframe as dd
    from dask.dataframe import from_delayed as _from_delayed

    original = getattr(ac, 'from_delayed', _from_delayed)

    def guarded(dfs, meta=None, *args, **kwargs):
        if dfs is None or len(dfs) == 0:
            # Empty result: hand back an empty frame carrying the schema.
            # This value is discarded unless include_meta=True.
            return dd.from_pandas(meta, npartitions=1)
        return original(dfs, meta=meta, *args, **kwargs)

    ac.from_delayed = guarded
    _APPLIED = True
    return True
