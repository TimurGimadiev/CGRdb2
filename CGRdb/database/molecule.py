from functools import cached_property
from pickle import dumps, loads
from typing import Union, Optional, List

from CGRtools.containers import MoleculeContainer
from StructureFingerprint import LinearFingerprint
from pony.orm import PrimaryKey, Required, Set, IntArray
from datasketch import MinHash
from .indexes import CursorHolder, RequestPack
from . import reaction
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
    def _similarity_filtering(cls, fingerprint: list) -> str:
        h = MinHash(num_perm=config.lsh_num_permute, hashfunc=hash)
        hashranges = json.loads(config.hashranges)
        h.update_batch(fingerprint)
        request_end = "\nOR\n".join(f"""band={n} AND key='{hash(int.from_bytes(h.hashvalues[i:j].byteswap().data,
                                                                               'big'))}'""" for n, (i, j) in
                                    enumerate(hashranges, 1))
        return request_end

    @classmethod
    def __postprocess_unordered_molecules(cls, result, *, fingerprint=None):
        struct = MoleculeStructure[result[1]]
        tanimoto = len(set(fingerprint).intersection(struct.fingerprint)) / \
                   len(set(fingerprint).union(struct.fingerprint))
        return Molecule[result[0]], struct, tanimoto

    @classmethod
    def __postprocess_ordered_molecules(cls, result):
        return Molecule[result[0]], MoleculeStructure[result[1]], result[2]

    @classmethod
    def _similarity_unorderd_request(cls, fingerprint: list) -> RequestPack:
        request_end = cls._similarity_filtering(fingerprint)
        request = f"""
            SELECT x.molecule m, x.id
            FROM MoleculeStructure x
            WHERE x.id IN (
            SELECT distinct unnest(records) m FROM moleculesimilarityindex WHERE {request_end})"""
        return RequestPack(request, partial(cls.__postprocess_unordered_molecules, fingerprint=fingerprint))

    @classmethod
    def _similarity_ordered(cls, fingerprint: list) -> RequestPack:
        request_end = cls._similarity_filtering(fingerprint)
        request = f'''
           SELECT x.molecule m, x.id, icount(x.fingerprint & '{set(fingerprint)}') / icount(x.fingerprint |
                '{set(fingerprint)}')::float t
           FROM MoleculeStructure x
           WHERE x.id IN (
           SELECT distinct unnest(records) m FROM moleculesimilarityindex WHERE {request_end})
           ORDER BY t DESC
           '''
        return RequestPack(request, cls.__postprocess_ordered_molecules)

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
            request_pack = cls._similarity_ordered(fingerprint)
        else:
            request_pack = cls._similarity_unorderd_request(fingerprint)
        return CursorHolder(request_pack, cls._database_)

    @classmethod
    def _substructure_ordered(cls, fingerprint: list) -> RequestPack:
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
        return RequestPack(request, cls.__postprocess_ordered_molecules)

    @classmethod
    def _search_substructure_unordered_request(cls, fingerprint: list):
        request = f'''
        WITH table_1 AS (SELECT a.molecule m, a.id s FROM MoleculeStructure a
        WHERE a.fingerprint @> '{set(fingerprint)}')
        SELECT h.m, h.s
        FROM table_1 h JOIN MoleculeStructure s ON h.s = s.id
        '''
        return RequestPack(request, partial(cls.__postprocess__unordered, fingerprint=fingerprint))

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
            request_pack = cls._substructure_ordered(fingerprint)
        else:
            request_pack = cls._search_substructure_unordered_request(fingerprint)
        return CursorHolder(request_pack, cls._database_)

    @classmethod
    def __postprocess_unordered_reaction(cls, result, *, fingerprint=None):
        if not fingerprint:
            raise ValueError("Initial molecule fingerprint is not provided")
        mol = Molecule[result[1]]
        tanimoto = len(set(fingerprint).intersection(mol.canonic_structure.fingerprint)) / \
                   len(set(fingerprint).union(mol.canonic_structure.fingerprint))
        return reaction.Reaction[result[0]], mol, tanimoto

    @classmethod
    def __postprocess_ordered_reaction(cls, result):
        return reaction.Reaction[result[0]], Molecule[result[1]], result[2]

    @classmethod
    def _search_similar_unordered_in_reactions_request(cls, fingerprint):
        request_core, _ = cls._similarity_unorderd_request(fingerprint)
        request = f'''
        WITH
        molecules as ({request_core}),
        substances as (SELECT s.substance, molecules.m
            FROM substancestructure s JOIN molecules ON s.molecule = molecules.m)
        SELECT r.reaction as reaction_id, r_s.m as molecule_id FROM reactionsubstance r JOIN substances r_s 
        ON r.substance = r_s.substance'''
        return RequestPack(request, partial(cls.__postprocess_unordered_reaction, fingerprint=fingerprint))

    @classmethod
    def __role_filter(cls, request_pack: RequestPack, role: bool):  # True for products, False - reactants as in DB
        request, postprocess = request_pack
        if role:
            request += ' WHERE r.is_product = TRUE'
        else:
            request += ' WHERE r.is_product = FALSE'
        return RequestPack(request, postprocess)

    @classmethod
    def _search_similar_ordered_in_reactions_request(cls, fingerprint):
        request_core, _ = cls._similarity_ordered(fingerprint)
        request = f'''
            WITH
            molecules as ({request_core}),
            substances as (SELECT s.substance, molecules.m, molecules.t
            FROM substancestructure s JOIN molecules ON s.molecule = molecules.m)
            SELECT r.reaction as reaction_id, r_s.m as molecule_id, r_s.t as tanimoto
            FROM reactionsubstance r JOIN substances r_s ON r.substance = r_s.substance'''
        return RequestPack(request, cls.__postprocess_ordered_reaction)

    @classmethod
    def similar_to_structure_in_reactions(cls, mol: [MoleculeContainer, "Molecule"], ordered=True,
                                          role: Optional[str] = None):
        if isinstance(mol, Molecule):
            fingerprint = mol.canonic_structure.fingerprint
        elif isinstance(mol, MoleculeContainer):
            fingerprint = LinearFingerprint(**config.fingerprint).transform_bitset([mol])[0]
        else:
            raise ValueError("CGRtools.MoleculeContainer or CGRdb.Molecule should be provided")
        if ordered:
            request_pack = cls._search_similar_ordered_in_reactions_request(fingerprint=fingerprint)
        else:
            request_pack = cls._search_similar_unordered_in_reactions_request(fingerprint=fingerprint)
        if role is not None:
            if isinstance(role, str) and role == "product":
                cls.__role_filter(request_pack, True)
            elif isinstance(role, str) and role == "reactant":
                cls.__role_filter(request_pack, False)
            else:
                raise ValueError("role can be String (product, reactant) or None")
        return CursorHolder(request_pack, cls._database_)

    def similar_in_reaction(self, ordered=True, role: Optional[str] = None):
        return self.similar_to_structure_in_reactions(self, ordered=ordered, role=role)

    @classmethod
    def _exact_match(cls, mol: MoleculeContainer):
        request = f'''SELECT x.molecule
        FROM MoleculeStructure x
        WHERE x.signature = '\\x{bytes(mol).hex()}'::bytea'''
        return RequestPack(request, cls.__postprocess_exact_mol)

    @classmethod
    def __postprocess_exact_mol(cls, result):
        return Molecule[result]

    @classmethod
    def __postprocess_exact_reactions(cls, result):
        return reaction.Reaction[result]

    @classmethod
    def get(cls, mol: MoleculeContainer):
        if not isinstance(mol, MoleculeContainer):
            raise ValueError("CGRtools.MoleculeContainer should be provided")
        return CursorHolder(cls._exact_match(mol), cls._database_)

    @classmethod
    def _exact_match_in_reactions(cls, mol):
        request_core, postprocess = cls._exact_match(mol)
        request = f'''  WITH molecules as ({request_core}),
                        substances as (SELECT s.molecule m, s.substance 
                        FROM substancestructure s JOIN molecules ON s.molecule = molecules.molecule)
                        SELECT r.reaction as reaction_id
                        FROM reactionsubstance r JOIN substances r_s ON r.substance = r_s.substance'''
        return RequestPack(request, cls.__postprocess_exact_reactions)

    @classmethod
    def get_reactions_with(cls, mol: MoleculeContainer, role: Optional[str] = None):
        if not isinstance(mol, MoleculeContainer):
            raise ValueError("CGRtools.MoleculeContainer should be provided")
        request_pack = cls._exact_match_in_reactions(mol)
        if role is not None:
            if isinstance(role, str) and role == "product":
                request_pack = cls.__role_filter(request_pack, True)
            elif isinstance(role, str) and role == "reactant":
                request_pack = cls.__role_filter(request_pack, False)
            else:
                raise ValueError("role can be String (product, reactant) or None")
        return CursorHolder(request_pack, cls._database_)

    def in_reactions(self, role: Optional[str] = None):
        request_pack = RequestPack(request=f'''WITH substances as (SELECT s.substance
                            FROM substancestructure s WHERE s.molecule={self.id})
                            SELECT r.reaction as reaction_id
                            FROM reactionsubstance r JOIN substances r_s ON r.substance = r_s.substance''',
                                   postprocess=self.__postprocess_exact_reactions)
        if role is not None:
            if isinstance(role, str) and role == "product":
                self.__role_filter(request_pack, True)
            elif isinstance(role, str) and role == "reactant":
                self.__role_filter(request_pack, False)
            else:
                raise ValueError("role can be String (product, reactant) or None")
        return CursorHolder(request_pack, self._database_)

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

    @classmethod
    def __postprocess_exact_non_org(cls, result):
        return NonOrganic[result]

    @classmethod
    def __postprocess_exact_reactions(cls, result):
        return reaction.Reaction[result]

    @classmethod
    def _exact_match(cls, mol: MoleculeContainer):
        request = f'''SELECT x.id
            FROM NonOrganic x
            WHERE x.signature = '\\x{bytes(mol).hex()}'::bytea'''
        return RequestPack(request, cls.__postprocess_exact_non_org)

    @classmethod
    def get(cls, mol: MoleculeContainer):
        if not isinstance(mol, MoleculeContainer):
            raise ValueError("CGRtools.MoleculeContainer should be provided")
        return CursorHolder(cls._exact_match(mol), cls._database_)

    @classmethod
    def _exact_match_in_reactions(cls, mol):
        request_core, postprocess = cls._exact_match(mol)
        request = f'''  WITH non_organic as ({request_core}),
                            substances as (SELECT s.non_organic n, s.substance 
                            FROM substancestructure s JOIN non_organic ON s.non_organic = non_organic.id)
                            SELECT r.reaction as reaction_id
                            FROM reactionsubstance r JOIN substances r_s ON r.substance = r_s.substance'''
        return RequestPack(request, cls.__postprocess_exact_reactions)

    @classmethod
    def __role_filter(cls, request_pack: RequestPack, role: bool):  # True for products, False - reactants as in DB
        request, postprocess = request_pack
        if role:
            request += ' WHERE r.is_product = TRUE'
        else:
            request += ' WHERE r.is_product = FALSE'
        return RequestPack(request, postprocess)

    @classmethod
    def get_reactions_with(cls, mol: MoleculeContainer, role: Optional[str] = None):
        if not isinstance(mol, MoleculeContainer):
            raise ValueError("CGRtools.MoleculeContainer should be provided")
        request_pack = cls._exact_match_in_reactions(mol)
        if role is not None:
            if isinstance(role, str) and role == "product":
                request_pack = cls.__role_filter(request_pack, True)
            elif isinstance(role, str) and role == "reactant":
                request_pack = cls.__role_filter(request_pack, False)
            else:
                raise ValueError("role can be String (product, reactant) or None")
        return CursorHolder(request_pack, cls._database_)




__all__ = ['Molecule', 'MoleculeStructure', 'NonOrganic']
