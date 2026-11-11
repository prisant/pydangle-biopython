"""Alternate-conformation (altloc) selection.

Two modes, chosen with ``--altloc``:

``first`` (default)
    Measure the first conformer.  "First" is local: among the alternates
    of the atoms concerned (a residue, a phi/psi window, one atom), the
    label with the lowest atom serial number, i.e. the one written first
    in the file (normally A).  Ranking labels over the whole file instead
    picks B everywhere in an entry whose first alternate residue has only
    B and C (3zeu A1), where the Richardson Lab measured A (3zeu B168).
    BioPython lists alternates alphabetically, so the serial numbers, not
    that list, give the file order.  Before measuring, :func:`prepare_model`

    * resolves microheterogeneity -- two residue types at one position --
      by keeping one residue type: the first conformer whose type is a
      standard residue (the 20 amino acids or MSE, :data:`STANDARD`), or
      the first conformer if none is.  So where conformer A is a modified
      residue and B the standard one (5k0a A157 CSO/CYS, 1xpm D111
      SCY/CYS), the standard residue is kept and measured, as in the
      Richardson Lab's Top2018 database, whose phi/psi there are those of
      the standard conformer.  BioPython represents microheterogeneity
      either as a ``DisorderedResidue`` (both types written as ATOM, e.g.
      2zs0 D31 Ile/Val) or as two separate residues with the same number
      and insertion code and different hetero flags (one type written as
      HETATM, e.g. 1din A123 CSD/Cys); both forms are handled.  (Where the
      first type has an atom without altloc label, BioPython's parser
      keeps only that type and discards the others, e.g. 1ejg A22 Pro/Ser
      keeps the Pro; this happens in both modes);
    * points every disordered atom at its first conformer.

    Each residue is then measured in one conformer: the first conformer
    shared by every atom of its phi/psi window -- C(i-1), N, CA, C, N(i+1)
    -- (:func:`set_residue`).  Every measurement of the residue whose
    disordered atoms all have that conformer uses it (:func:`vectors`);
    otherwise the measurement uses the first conformer its own disordered
    atoms share, and failing that each atom's own first conformer.  So a
    dihedral is never built from two conformers when a common one exists,
    a residue present only as conformer B is measured with its neighbours'
    B atoms, and so is the residue after it (6b6w B301 and B302, 5j1s B432:
    the values in the Richardson Lab's Top2018 database).

``occupancy``
    BioPython's default and pydangle's behaviour before v0.6: for each atom
    separately, the alternate with the highest occupancy (ties: the first
    listed).  Nothing is changed before measuring.  Microheterogeneity
    positions give one record per residue type.

The peptide-bond test (``measure._PeptideBondBuilder``) does not depend on
the mode: it accepts any compatible altloc pair.

Decisions 2026-11-09 (first conformer) and 2026-11-10 (per-residue order,
standard residue type at microheterogeneity).
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterator
from typing import Any

MODES = ("first", "occupancy")

# Residue types preferred at microheterogeneity positions: the 20 standard
# amino acids and selenomethionine, the residue types with Ramachandran
# values in the Richardson Lab's Top2018 database.
STANDARD = frozenset(
    {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "MSE",
    }
)

_state: dict[str, Any] = {"mode": "first", "residue": None}


def set_mode(mode: str) -> None:
    """Set the altloc mode used by :func:`vectors` (``first`` or ``occupancy``)."""
    if mode not in MODES:
        raise ValueError(f"unknown altloc mode {mode!r}; expected one of {MODES}")
    _state["mode"] = mode


def get_mode() -> str:
    """Return the current altloc mode."""
    return str(_state["mode"])


def _residue_children(residue: Any) -> list[Any]:
    """Return the residue, or every residue type of a ``DisorderedResidue``."""
    if residue.is_disordered() == 2:
        return list(residue.disordered_get_list())
    return [residue]


def _atom_children(atom: Any) -> list[Any]:
    return list(atom.disordered_get_list()) if atom.is_disordered() else [atom]


def _alt_atoms(residue: Any) -> Iterator[Any]:
    """Every alternate of every atom of every residue type of *residue*."""
    for child in _residue_children(residue):
        for atom in child:
            yield from _atom_children(atom)


def _labels(residue: Any) -> set[str]:
    """Non-blank altloc labels over every atom of every residue type."""
    return {a.get_altloc() for a in _alt_atoms(residue) if a.get_altloc() != " "}


def _all_labelled(residue: Any) -> bool:
    """True if every atom of every residue type carries a non-blank altloc."""
    return all(a.get_altloc() != " " for a in _alt_atoms(residue))


def _key(label: str, atoms: list[Any]) -> tuple[int, str]:
    """Sort key of *label*: its lowest atom serial among *atoms*' alternates."""
    serials = [
        alt.get_serial_number()
        for atom in atoms
        for alt in _atom_children(atom)
        if alt.get_altloc() == label and isinstance(alt.get_serial_number(), int)
    ]
    return (min(serials) if serials else 10**12, label)


def _first_label(labels: set[str], atoms: list[Any]) -> str:
    """The label of *labels* written first among *atoms*' alternates."""
    return min(labels, key=lambda label: _key(label, atoms))


def _preferred(candidates: list[tuple[str, tuple[int, str]]]) -> int:
    """Index of the residue type to keep: (resname, first-label key) pairs.

    The first conformer among the standard residue types, or the first
    conformer if none is standard.
    """
    standard = [i for i, (name, _) in enumerate(candidates) if name in STANDARD]
    pool = standard or list(range(len(candidates)))
    return min(pool, key=lambda i: candidates[i][1])


def _resolve_split_microheterogeneity(chain: Any) -> None:
    """Keep one residue where two share a number but differ in hetero flag.

    The residue kept is chosen by :func:`_preferred`.  Only when each of
    them consists entirely of altloc-labelled atoms and their labels do not
    overlap (alternates of one position, e.g. 1din A123: HETATM CSD altloc
    A, ATOM CYS altloc B).  A ligand that merely shares a residue number is
    left alone.
    """
    groups: dict[tuple[int, str], list[Any]] = defaultdict(list)
    for residue in chain:
        groups[(residue.id[1], residue.id[2])].append(residue)
    for members in groups.values():
        if len(members) < 2 or not all(_all_labelled(r) for r in members):
            continue
        label_sets = [_labels(r) for r in members]
        if any(a & b for i, a in enumerate(label_sets) for b in label_sets[i + 1:]):
            continue
        candidates = []
        for residue, labels in zip(members, label_sets, strict=True):
            atoms = [a for child in _residue_children(residue) for a in child]
            candidates.append(
                (residue.get_resname(), _key(_first_label(labels, atoms), atoms))
            )
        keep = _preferred(candidates)
        for i, residue in enumerate(members):
            if i != keep:
                chain.detach_child(residue.id)


def _select_residue_type(residue: Any) -> None:
    """In a ``DisorderedResidue``, select the type chosen by :func:`_preferred`."""
    names = []
    candidates = []
    for name in residue.disordered_get_id_list():
        child = residue.disordered_get(name)
        labels = _labels(child)
        if labels:
            atoms = list(child)
            names.append(name)
            candidates.append((name, _key(_first_label(labels, atoms), atoms)))
    if candidates:
        residue.disordered_select(names[_preferred(candidates)])


def prepare_model(model: Any) -> None:
    """In ``first`` mode, select the conformers to measure.

    Changes the model in place: at split microheterogeneity positions
    detaches the residues not kept, selects the kept residue type in each
    ``DisorderedResidue`` and the first conformer of each disordered atom.
    Does nothing in ``occupancy`` mode.
    """
    if _state["mode"] != "first":
        return
    for chain in model:
        _resolve_split_microheterogeneity(chain)
        for residue in chain:
            if residue.is_disordered() == 2:
                _select_residue_type(residue)
            for child in _residue_children(residue):
                for atom in child:
                    if atom.is_disordered():
                        atom.disordered_select(
                            _first_label(set(atom.disordered_get_id_list()), [atom])
                        )


def _common_labels(atoms: list[Any]) -> set[str] | None:
    """Altloc labels shared by all disordered *atoms* (None if none is disordered)."""
    common: set[str] | None = None
    for atom in atoms:
        if atom.is_disordered():
            ids = set(atom.disordered_get_id_list())
            common = ids if common is None else common & ids
    return common


def set_residue(residue_list: list[Any] | None, index: int = 0) -> None:
    """Fix the conformer for the residue at *index* (``first`` mode).

    The first conformer shared by every disordered atom of the phi/psi
    window C(i-1), N, CA, C, N(i+1) (atoms that are missing are skipped).
    ``set_residue(None)`` clears it.
    """
    _state["residue"] = None
    if residue_list is None or _state["mode"] != "first":
        return
    window = []
    if index > 0 and residue_list[index - 1].has_id("C"):
        window.append(residue_list[index - 1]["C"])
    residue = residue_list[index]
    window += [residue[name] for name in ("N", "CA", "C") if residue.has_id(name)]
    if index + 1 < len(residue_list) and residue_list[index + 1].has_id("N"):
        window.append(residue_list[index + 1]["N"])
    common = _common_labels(window)
    if common:
        _state["residue"] = _first_label(common, window)


def vectors(atoms: list[Any]) -> list[Any]:
    """Coordinates of *atoms* for one measurement.

    ``occupancy`` mode: each atom's current selection.  ``first`` mode: the
    residue's conformer (:func:`set_residue`) if every disordered atom has
    it; otherwise the first conformer the disordered atoms share; if they
    share none, each atom's own first conformer (its selection after
    :func:`prepare_model`).
    """
    if _state["mode"] != "first":
        return [a.get_vector() for a in atoms]
    common = _common_labels(atoms)
    if not common:
        return [a.get_vector() for a in atoms]
    preferred = _state["residue"]
    label = preferred if preferred in common else _first_label(common, atoms)
    return [
        a.disordered_get(label).get_vector() if a.is_disordered() else a.get_vector()
        for a in atoms
    ]
