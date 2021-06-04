from functools import cached_property
from pickle import dumps, loads
from typing import Union, Optional, List, Callable, Any

from CGRtools.containers import MoleculeContainer
from StructureFingerprint import LinearFingerprint
from pony.orm import PrimaryKey, Required, Set, IntArray, db_session
from datasketch import MinHash
from .indexes import CursorHolder
from functools import partial

import json
from . import config
from .config import Entity, db


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

    @classmethod
    def _similarity_pre_request(cls, fingerprint: list) -> str:
        h = MinHash(num_perm=config.lsh_num_permute, hashfunc=hash)
        hashranges = json.loads(config.hashranges)
        h.update_batch(fingerprint)
        request_end = "\nOR\n".join(f"""band={n} AND key='{hash(int.from_bytes(h.hashvalues[i:j].byteswap().data,
                                                                               'big'))}'""" for n, (i, j) in
                                    enumerate(hashranges, 1))
        return request_end

    @classmethod
    def __postprocess_unordered(cls, result, *, fingerprint=None):
        struct = MoleculeStructure[result[1]]
        tanimoto = len(set(fingerprint).intersection(struct.fingerprint)) / \
                   len(set(fingerprint).union(struct.fingerprint))
        return Molecule[result[0]], struct, tanimoto

    @classmethod
    def __postprocess_ordered(cls, result):
        return Molecule[result[0]], MoleculeStructure[result[1]], result[2]

    @classmethod
    def _search_similar_unorderd_request(cls, fingerprint: list):
        request_end = cls._similarity_pre_request(fingerprint)
        request = f"""
            SELECT x.molecule m, x.id
            FROM MoleculeStructure x
            WHERE x.id IN (
            SELECT distinct unnest(records) m FROM moleculesimilarityindex WHERE {request_end})"""
        return request, partial(cls.__postprocess_unordered, fingerprint=fingerprint)

    @classmethod
    def _search_similar_orderd_request(cls, fingerprint: list) -> tuple[str, Callable]:
        request_end = cls._similarity_pre_request(fingerprint)
        request = f"""
           SELECT x.molecule m, x.id, icount(x.fingerprint & '{set(fingerprint)}') / icount(x.fingerprint |
                '{set(fingerprint)}')::float t
           FROM MoleculeStructure x
           WHERE x.id IN (
           SELECT distinct unnest(records) m FROM moleculesimilarityindex WHERE {request_end})
           ORDER BY t DESC"""
        return request, cls.__postprocess_ordered

    def similars(self, ordered=True):
        return self.similar_to_structure(self, ordered)

    @classmethod
    def similar_to_structure(cls, mol: [MoleculeContainer, "Molecule"], ordered=True):
        if isinstance(mol, Molecule):
            fingerprint = mol.canonic_structure.fingerprint
        elif isinstance(mol, MoleculeContainer):
            fingerprint = LinearFingerprint(**config.fingerprint).transform_bitset([mol])[0]
        else:
            raise ValueError(" Only CGRtools.MoleculeContainer or CGRdb.Molecule")
        if ordered:
            request, postprocess = cls._search_similar_orderd_request(fingerprint)
        else:
            request, postprocess = cls._search_similar_unorderd_request(fingerprint)
        return CursorHolder(request, cls._database_, postprocess)

    @classmethod
    def _search_substructure_ordered_request(cls, fingerprint: list):
        len_fp = len(fingerprint)
        request = f'''
        WITH table_1 AS (SELECT a.molecule, a.fingerprint, a.id, a.fingerprint_len 
        FROM MoleculeStructure a
        WHERE a.fingerprint_len BETWEEN {len_fp} AND 
        {int(len_fp // 0.1)+5} AND a.fingerprint @> '{set(fingerprint)}'),
        table_2 AS (SELECT DISTINCT ON (table_1.molecule)  table_1.molecule m, table_1.id s,
        {len_fp} / table_1.fingerprint_len::float t
        FROM table_1)
        SELECT h.m, h.s, h.t
        FROM table_2 h JOIN MoleculeStructure s ON h.s = s.id
        ORDER BY h.t DESC
        '''
        return request, cls.__postprocess_ordered

    @classmethod
    def _search_substructure_unordered_request(cls, fingerprint: list):
        request = f'''
        WITH table_1 AS (SELECT a.molecule m, a.id s FROM MoleculeStructure a
        WHERE a.fingerprint @> '{set(fingerprint)}')
        SELECT h.m, h.s
        FROM table_1 h JOIN MoleculeStructure s ON h.s = s.id
        '''
        return request, partial(cls.__postprocess__unordered, fingerprint=fingerprint)

    def clear_substructure_cursor(self):
        if hasattr(self, "substructure_cursor"):
            del self.substructure_cursor

    def substructures(self, ordered=True):
        return self.substructres_to_structure(self, ordered=ordered)

    @classmethod
    def substructres_to_structure(cls, mol: [MoleculeContainer, "Molecule"], ordered=True):
        if isinstance(mol, Molecule):
            fingerprint = mol.canonic_structure.fingerprint
        elif isinstance(mol, MoleculeContainer):
            fingerprint = LinearFingerprint(**config.fingerprint).transform_bitset([mol])[0]
        else:
            raise ValueError(" Only CGRtools.MoleculeContainer or CGRdb.Molecule")
        if ordered:
            request, postprocess = cls._search_substructure_ordered_request(fingerprint)
        else:
            request, postprocess = cls._search_substructure_unordered_request(fingerprint)
        return CursorHolder(request, cls._database_, postprocess)


class MoleculeStructure(Entity):
    id = PrimaryKey(int, auto=True)
    molecule = Required('Molecule')
    signature = Required(bytes, unique=True)
    smiles = Required(str)
    fingerprint = Required(IntArray)
    fingerprint_len = Required(int)
    is_canonic = Required(bool)
    _structure = Required(bytes)
    _fast_mapping = Required(IntArray)

    def __init__(self, mol: Union[bytes, MoleculeContainer], molecule: Molecule, /, is_canonic: bool = True, *,
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
            if smiles is None or fast_mapping is None:
                fast_mapping = mol.smiles_atoms_order
                smiles = str(mol)
            if signature is None:
                signature = bytes(mol)
            if fingerprint is None:
                fingerprint = LinearFingerprint(**config.fingerprint).transform_bitset([mol])[0]
            mol = dumps(mol)
        super().__init__(_structure=mol, fingerprint=fingerprint, fingerprint_len=len(fingerprint), molecule=molecule,
                         signature=signature, is_canonic=is_canonic, smiles=smiles, _fast_mapping=fast_mapping)

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
            if smiles is None or fast_mapping is None:
                fast_mapping = mol.smiles_atoms_order
                smiles = str(mol)
            if signature is None:
                signature = bytes(mol)
            mol = dumps(mol)
        super().__init__(_structure=mol, signature=signature, smiles=smiles, _fast_mapping=fast_mapping)

    @cached_property
    def structure(self) -> MoleculeContainer:
        return loads(self._structure)


__all__ = ['Molecule', 'MoleculeStructure', 'NonOrganic']
