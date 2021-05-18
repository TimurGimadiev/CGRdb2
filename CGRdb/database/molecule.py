from pickle import dumps, loads
from CGRtools.containers import MoleculeContainer
from pony.orm import PrimaryKey, Required, Set, IntArray
from functools import cached_property
from typing import Union, Optional, List
from StructureFingerprint import LinearFingerprint
from . import config
from .config import Entity


class Molecule(Entity):
    id = PrimaryKey(int, auto=True)
    structures = Set('MoleculeStructure')
    substances = Set('SubstanceStructure')

    def __init__(self, mol: Optional[MoleculeContainer] = None, /):
        super().__init__()
        if mol is not None:
            MoleculeStructure(mol, self)

    @cached_property
    def canonic_structure(self) -> 'MoleculeStructure':
        return MoleculeStructure.get(molecule=self, is_canonic=True)


class MoleculeStructure(Entity):
    id = PrimaryKey(int, auto=True)
    molecule = Required('Molecule')
    signature = Required(bytes, unique=True)
    smiles = Required(str)
    fingerprint = Required(IntArray)
    is_canonic = Required(bool)
    _structure = Required(bytes)
    _fast_mapping = Required(IntArray)

    def __init__(self, mol: Union[bytes, MoleculeContainer], molecule: 'Molecule', /, is_canonic: bool = True, *,
                 fingerprint: Optional[List[int]] = None, signature: Optional[bytes] = None,
                 smiles: Optional[str] = None, fast_mapping: Optional[List[int]] = None):
        if isinstance(mol, bytes):  # low-level. errors not controlled
            if signature is None or fingerprint is None or smiles is None or fast_mapping is None:
                raise ValueError('for dumps signature, fingerprint, smiles and fast_mapping is required')
        elif not isinstance(mol, MoleculeContainer):
            raise TypeError('CGRtools.MoleculeContainer container or pickle dumped CGRtools.MoleculeContainer expected')
        else:
            if molecule.id is None:
                # This was done in optimistic expectation that new Molecule obj has no id before commit
                if not is_canonic:
                    # potential problem:
                    # more than one MoleculeStructure objects creation for single Molecule in same session
                    raise ValueError('first structure should be canonic')
            elif is_canonic:
                raise ValueError('only one structure should be canonic')
            elif {(n, a.atomic_number, a.isotope) for n, a in molecule.canonic_structure.structure.atoms()} != \
                    {(n, a.atomic_number, a.isotope) for n, a in mol.atoms()}:
                raise ValueError('atoms order in CGRtools.MoleculeContainer of new structure '
                                 'not corresponds to canonic structure')
            if fingerprint is None:
                fingerprint = LinearFingerprint(**config.fingerprint).transform_hashes([mol])[0]
            if signature is None:
                signature = bytes(mol)
            if smiles is None or fast_mapping is None:
                fast_mapping = mol.smiles_atoms_order
                smiles = str(mol)
            mol = dumps(mol)
        super().__init__(_structure=mol, fingerprint=fingerprint, molecule=molecule, signature=signature,
                         is_canonic=is_canonic, smiles=smiles, _fast_mapping=fast_mapping)

    @cached_property
    def structure(self) -> MoleculeContainer:
        return loads(self._structure)


class NonOrganic(Entity):
    id = PrimaryKey(int, auto=True)
    signature = Required(bytes, unique=True)
    smiles = Required(str)
    _structure = Required(bytes)
    _fast_mapping = Required(IntArray)
    substances = Set('SubstanceStructure')

    def __init__(self, mol: Union[bytes, MoleculeContainer], /, *, signature: Optional[bytes] = None,
                 smiles: Optional[str] = None, fast_mapping: Optional[List[int]] = None):
        if isinstance(mol, bytes):
            if signature is None or smiles is None or fast_mapping is None:
                raise ValueError('for dumps signature is required')
        elif not isinstance(mol, MoleculeContainer):
            raise TypeError('CGRtools.MoleculeContainer container or pickle dumped CGRtools.MoleculeContainer expected')
        else:
            if signature is None:
                signature = bytes(mol)
            if smiles is None or fast_mapping is None:
                fast_mapping = mol.smiles_atoms_order
                smiles = str(mol)
            mol = dumps(mol)
        super().__init__(_structure=mol, signature=signature, smiles=smiles, _fast_mapping=fast_mapping)

    @cached_property
    def structure(self) -> MoleculeContainer:
        return loads(self._structure)


__all__ = ['Molecule', 'MoleculeStructure', 'NonOrganic']
