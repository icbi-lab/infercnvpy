"""`use`/`out` arguments like scanpy>=1.13: :mod:`anndata.acc` accessors, or strings resolved by `A.resolve`."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import scanpy as sc
from anndata.acc import A, AdRef, GraphAcc, LayerAcc, MetaAcc, MultiAcc
from scverse_misc import Deprecation

if TYPE_CHECKING:
    from anndata import AnnData

DEPR_USE_REP = Deprecation("0.7.0", "Use `use` instead, e.g. `use=A.obsm['X_cnv']`.")
DEPR_NEIGHBORS_KEY = Deprecation("0.7.0", "Use `use` instead, e.g. `use=A.obsp['cnv_neighbors']`.")
DEPR_LAYER = Deprecation("0.7.0", "Use `use` instead, e.g. `use=A.layers['counts']`.")
DEPR_KEY_ADDED = Deprecation("0.7.0", "Use `out` instead.")
DEPR_INPLACE = Deprecation("0.7.0", "Use `out=True` to return the result instead of writing it.")
DEPR_INPLACE_COPY = Deprecation("0.7.0", "Copy `adata` before calling this function instead.")


def _v2_preset() -> bool:
    return sc.settings.preset == "scanpy-v2-preview"


def _obsm_key(name: str) -> str:
    """Default `.obsm` key: `X_{name}`, or `name` under scanpy's 2.0 preview preset."""
    return name if _v2_preset() else f"X_{name}"


def _dest(spec: Any) -> tuple[str, str | None]:
    """Slot and key an accessor or spec string points to, e.g. `("obsm", "X_cnv")` or `("X", None)`."""
    acc = A.resolve(spec) if isinstance(spec, str) else spec
    match acc:
        case MultiAcc():
            return f"{acc.dim}m", acc.k
        case GraphAcc():
            return f"{acc.dim}p", acc.k
        case LayerAcc():
            return ("X", None) if acc.k is None else ("layers", acc.k)
        case AdRef(acc=MetaAcc() as meta) if isinstance(acc.idx, str):
            return meta.dim, acc.idx
    raise TypeError(f"Expected an anndata accessor such as `A.obsm['X_cnv']` or a spec string, got {spec!r}.")


def _key(spec: Any, *slots: str) -> str:
    """Key `spec` points to, which has to be in one of `slots` (that all have keys, i.e. not `X`)."""
    slot, key = _dest(spec)
    if slot not in slots or key is None:
        raise TypeError(f"Expected an accessor into `.{'`/`.'.join(slots)}`, got {spec!r}.")
    return key


def _legacy(new: Any, old: Any, translated: Any, *, name: str, default: Any = None) -> Any:
    """Translate a deprecated argument `old` (`None` if not passed) to its replacement `new`."""
    if old is None:
        return new
    if new is not default:
        raise TypeError(f"Pass either the new argument or `{name}`, not both.")
    return translated


def _uns_entry(adata: AnnData, obsm_key: str) -> dict:
    """The `.uns` entry written along `.obsm[obsm_key]`: named after it, without an `X_` prefix."""
    try:
        return adata.uns[obsm_key.removeprefix("X_")]
    except KeyError:
        raise KeyError(f"`adata.uns[{obsm_key.removeprefix('X_')!r}]` not found. Did you run `tl.infercnv`?") from None


def _write_obsm(adata: AnnData, val: Any, out: Any, default: str) -> Any:
    """Return `val` for `out=True`, else write it to `.obsm` (`default` for `out=False`)."""
    if out is True:
        return val
    adata.obsm[_obsm_key(default) if out is False else _key(out, "obsm")] = val
    return None


def _spec(use: Any) -> str:
    """Spec string of an accessor for `.uns[...]["params"]`, e.g. `"obsm.X_cnv"`."""
    slot, key = _dest(use)
    return slot if key is None else f"{slot}.{key}"


def _get(adata: AnnData, use: Any) -> Any:
    """The array `use` points to; strings are resolved by `A.resolve`."""
    return adata[A.resolve(use, vec=False) if isinstance(use, str) else use]
