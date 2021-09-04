# -*- coding: utf-8 -*-
#
#  Copyright 2021 Timur Gimadiev <timur.gimadiev@gmail.com>
#  Copyright 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
#  This file is part of CGRdb.
#
#  CGRdb is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
import time
from datetime import datetime
from .config import Entity
from . import config
from . import reaction
from .indexes import CursorHolder, RequestPack

from functools import cached_property
from pickle import dumps, loads
from typing import Union, Optional, List
from CGRtools.containers import MoleculeContainer
from StructureFingerprint import LinearFingerprint
from pony.orm import PrimaryKey, Required, Set, IntArray, Json, select
from datasketch import MinHash
from functools import partial
import json


class Molecule(Entity):
    id = PrimaryKey(int, auto=True)
    structures = Set('MoleculeStructure')
    substances = Set('SubstanceStructure')
    metadata = Set("MoleculeProperties")
    classes = Set("MoleculeClass")

    def __init__(self, mol: Optional[MoleculeContainer] = None, /):
        super().__init__()
        if mol is not None:
            MoleculeStructure(mol, self)

    @cached_property
    def canonic_structure(self) -> 'MoleculeStructure':
        return MoleculeStructure.get(molecule=self, is_canonic=True)

# request modifiers and response processors

    @classmethod
    def __similarity_filtering(cls, fingerprint: list) -> str:
        h = MinHash(num_perm=config.lsh_num_permute, hashfunc=hash)
        hashranges = json.loads(config.hashranges)
        h.update_batch(fingerprint)
        request_end = "\nOR\n".join(f"""band={n} AND key='{hash(int.from_bytes(h.hashvalues[i:j].byteswap().data,
                                                                               'big'))}'""" for n, (i, j) in
                                    enumerate(hashranges, 1))
        return request_end

    @classmethod
    def __postprocess_unordered_molecules(cls, result, *, fingerprint=None, no_graph=False):
        struct = MoleculeStructure[result[1]]
        tanimoto = len(set(fingerprint).intersection(struct.fingerprint)) / \
                   len(set(fingerprint).union(struct.fingerprint))
        mol = Molecule[result[0]]
        if no_graph:
            for i in mol.structures:
                i._vals_.pop(MoleculeStructure.fingerprint)
                i._vals_.pop(MoleculeStructure._structure)
                i._vals_.pop(MoleculeStructure.signature)
                i._vals_.pop(MoleculeStructure._fast_mapping)
        return mol, struct, tanimoto

    @staticmethod
    def __postprocess_ordered_molecules(result, substr=None, tanimoto_limit=None):
        mol = Molecule[result[0]]
        struct = MoleculeStructure[result[1]]
        tanimoto = result[2]
        if substr is None:
            if tanimoto > tanimoto_limit:
                return mol, struct, tanimoto
            else:
                return "stop"
        else:
            if substr <= struct.structure:
                return mol, struct, tanimoto
            # clean cache for memory keeping
            del struct._vals_[MoleculeStructure._structure]
            del struct.__dict__['structure']


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
    def __postprocess_exact_mol(cls, result):
        return Molecule[result]

    @classmethod
    def __postprocess_exact_reactions(cls, result):
        return reaction.Reaction[result]

    @classmethod
    def __role_filter(cls, request_pack: RequestPack, is_product: bool):  # True for products, False - reactants as in DB
        request, postprocess, prefetch_map = request_pack
        if is_product:
            request += ' WHERE r.is_product = TRUE'
        else:
            request += ' WHERE r.is_product = FALSE'
        return RequestPack(request, postprocess, prefetch_map)

    @classmethod
    def __reactions_from_unordered(cls, request_core: str, *, fingerprint=None):
        request = f'''
                WITH
                molecules as ({request_core}),
                substances as (SELECT s.substance, molecules.m
                    FROM substancestructure s JOIN molecules ON s.molecule = molecules.m)
                SELECT DISTINCT r.reaction as reaction_id, r_s.m as molecule_id FROM reactionsubstance r JOIN substances r_s 
                ON r.substance = r_s.substance'''
        return RequestPack(request, partial(cls.__postprocess_unordered_reaction, fingerprint=fingerprint),
                           prefetch_map=(reaction.Reaction, 0,
                                         [reaction.Reaction.structure, reaction.Reaction.fingerprint]))

    @classmethod
    def __reactions_from_ordered(cls, request_core: str):
        request = f'''
                    WITH
                    molecules as ({request_core}),
                    substances as (SELECT s.substance, molecules.m, molecules.t
                    FROM substancestructure s JOIN molecules ON s.molecule = molecules.m)
                    SELECT DISTINCT r.reaction as reaction_id, r_s.m as molecule_id, r_s.t as tanimoto
                    FROM reactionsubstance r JOIN substances r_s ON r.substance = r_s.substance'''
        return RequestPack(request, cls.__postprocess_ordered_reaction,
                           prefetch_map=(reaction.Reaction, 0, None))

#   exact search block

    @classmethod
    def _exact_match(cls, mol: MoleculeContainer):
        request = f'''SELECT x.molecule
        FROM MoleculeStructure x
        WHERE x.signature = '\\x{bytes(mol).hex()}'::bytea'''
        return RequestPack(request, cls.__postprocess_exact_mol)

    @classmethod
    def get(cls, mol: Optional[MoleculeContainer] = None, **kwargs):  # fix this
        if mol is not None:
            if not isinstance(mol, MoleculeContainer):
                raise ValueError("CGRtools.MoleculeContainer should be provided")
            return CursorHolder(cls._exact_match(mol), cls._database_)
        else:
            return super().get(**kwargs)


    @classmethod
    def _exact_match_in_reactions(cls, mol):
        request_core, postprocess = cls._exact_match(mol)
        request = f'''  WITH molecules as ({request_core}),
                        substances as (SELECT s.molecule m, s.substance 
                        FROM substancestructure s JOIN molecules ON s.molecule = molecules.molecule)
                        SELECT r.reaction as reaction_id
                        FROM reactionsubstance r JOIN substances r_s ON r.substance = r_s.substance'''
        return RequestPack(request, cls.__postprocess_exact_reactions, prefetch_map=(None, None, None))

    @classmethod
    def reactions(cls, mol: MoleculeContainer, is_product: Optional[bool] = None, request_only=False):  # Fix this
        if not isinstance(mol, MoleculeContainer):
            raise ValueError("CGRtools.MoleculeContainer should be provided")
        request_pack = cls._exact_match_in_reactions(mol)
        if is_product is not None:
            if isinstance(is_product, bool):
                request_pack = cls.__role_filter(request_pack, is_product)
            else:
                raise ValueError("is_product can be boolean (True for product, False for reactant) or None")
        if request_only:
            return request_pack.request
        return CursorHolder(request_pack, cls._database_)

    def in_reactions(self, is_product: Optional[bool] = None):  # fix this combine with reactions classmethod
        request_pack = RequestPack(request=f'''WITH substances as (SELECT s.substance
                            FROM substancestructure s WHERE s.molecule={self.id})
                            SELECT r.reaction as reaction_id
                            FROM reactionsubstance r JOIN substances r_s ON r.substance = r_s.substance''',
                                   postprocess=self.__postprocess_exact_reactions, prefetch_map=(None, None, None))
        if is_product is not None:
            if isinstance(is_product, bool):
                request_pack = self.__role_filter(request_pack, is_product)
            else:
                raise ValueError("is_product can be boolean (True for product, False for reactant) or None")
        return CursorHolder(request_pack, self._database_)

# similarity search block

    @classmethod
    def _similarity_unorderd(cls, fingerprint: list, tanimoto_limit: Optional[float]) -> RequestPack:
        request_end = cls.__similarity_filtering(fingerprint)
        request = f"""
            SELECT x.molecule m, x.id
            FROM MoleculeStructure x
            WHERE x.id IN (
            SELECT distinct unnest(records) m FROM moleculesimilarityindex WHERE {request_end})"""
        return RequestPack(request, partial(cls.__postprocess_unordered_molecules, fingerprint=fingerprint,
                                            tanimoto_limit=tanimoto_limit),
                           prefetch_map=(MoleculeStructure, 1, [MoleculeStructure.fingerprint]))

    @classmethod
    def _similarity_ordered(cls, fingerprint: list, tanimoto_limit: Optional[float]) -> RequestPack:
        request_end = cls.__similarity_filtering(fingerprint)
        request = f'''
           SELECT x.molecule m, x.id, icount(x.fingerprint & '{set(fingerprint)}') / icount(x.fingerprint |
                '{set(fingerprint)}')::float t
           FROM MoleculeStructure x
           WHERE x.id IN (
           SELECT distinct unnest(records) m FROM moleculesimilarityindex WHERE {request_end})
           ORDER BY t DESC
           '''
        return RequestPack(request, partial(cls.__postprocess_ordered_molecules, tanimoto_limit=tanimoto_limit),
                           prefetch_map=(MoleculeStructure, 1, None))

    def similars(self: Union[MoleculeContainer, "Molecule"], ordered=True, request_only=False,
                 tanimoto_limit: Optional[float] = None):
        if tanimoto_limit is None:
            tanimoto_limit = config.lsh_threshold
        elif not isinstance(tanimoto_limit, float) or 1. < tanimoto_limit < 0.:
            raise f"tanimoto_limit can be float from 0 to 1, but {tanimoto_limit, type(tanimoto_limit)} was given"
        if isinstance(self, Molecule):
            fingerprint = self.canonic_structure.fingerprint
        elif isinstance(self, MoleculeContainer):
            fingerprint = LinearFingerprint(**config.fingerprint).transform_bitset([self])[0]
            self = Molecule
        else:
            raise ValueError(" Only CGRtools.MoleculeContainer or CGRdb.Molecule")
        if ordered:
            request_pack = self._similarity_ordered(fingerprint, tanimoto_limit=tanimoto_limit)
        else:
            request_pack = self._similarity_unorderd(fingerprint, tanimoto_limit=tanimoto_limit)
        if request_only:
            return request_pack.request
        return CursorHolder(request_pack, self._database_)

    def similars_in_reactions(self: Union[MoleculeContainer, "Molecule"], ordered=True,
                              is_product: Optional[bool] = None, request_only=False):
        if isinstance(self, MoleculeContainer):
            request = Molecule.similars(self, ordered=ordered, request_only=True)
            self = Molecule
        elif isinstance(self, Molecule):
            request = self.similars(ordered=ordered, request_only=True)
        else:
            raise TypeError("Only CGRtools.MoleculeContainer or CGRdb.Molecule can be used")
        if ordered:
            request_pack = self.__reactions_from_ordered(request_core=request)
        else:
            request_pack = self.__reactions_from_unordered(request_core=request)
        if is_product is not None:
            if isinstance(is_product, bool):
                request_pack = self.__role_filter(request_pack, is_product)
            else:
                raise ValueError("is_product can be boolean (True for product, False for reactant) or None")
        if request_only:
            return request_pack.request
        return CursorHolder(request_pack, self._database_)


# substructure search block

    @classmethod
    def _substructure_ordered(cls, fingerprint: list, substr=None, prefilter: bool = True,
                              tanimoto_limit: Optional[float] = None) -> RequestPack:
        len_fp = len(fingerprint)
        if prefilter:
            prefilter = f"""WHERE a.fingerprint_len BETWEEN {len_fp} AND 
        {int(len_fp // 0.1)+5} AND"""
        else:
            prefilter = "WHERE"
        request = f'''
        WITH table_1 AS (SELECT a.molecule, a.fingerprint, a.id, a.fingerprint_len 
        FROM MoleculeStructure a
        {prefilter} a.fingerprint @> '{set(fingerprint)}'),
        table_2 AS (SELECT DISTINCT ON (table_1.molecule)  table_1.molecule m, table_1.id s,
        {len_fp} / table_1.fingerprint_len::float t
        FROM table_1)
        SELECT h.m, h.s, h.t
        FROM table_2 h JOIN MoleculeStructure s ON h.s = s.id
        ORDER BY h.t DESC
        '''
        return RequestPack(request, partial(cls.__postprocess_ordered_molecules, substr=substr,
                                            tanimoto_limit=tanimoto_limit),
                           prefetch_map=(MoleculeStructure, 1, [MoleculeStructure._structure]))

    @classmethod
    def _substructure_unordered(cls, fingerprint: list, substr=None, tanimoto_limit: Optional[float] = None):
        request = f'''
        WITH table_1 AS (SELECT a.molecule m, a.id s FROM MoleculeStructure a
        WHERE a.fingerprint @> '{set(fingerprint)}')
        SELECT h.m, h.s
        FROM table_1 h JOIN MoleculeStructure s ON h.s = s.id
        '''
        return RequestPack(request, partial(cls.__postprocess_unordered_molecules,
                                            fingerprint=fingerprint, substr=substr, tanimoto_limit=tanimoto_limit),
                           prefetch_map=(MoleculeStructure, 1, [MoleculeStructure.fingerprint, MoleculeStructure._structure]))

    def substructures(self: Union[MoleculeContainer, "Molecule"], ordered=True, request_only=False,
                     tanimoto_limit: Optional[float] = None):
        """
        :param ordered: Return substructures ordered by similarity (also similarity limit will be used)
        :param request_only: Option is for internal usage or debugging to see actual sql request
        :param tanimoto_limit: takes tanimoto limit from config table of DB, can be override, to go for
    lower similarities, that we do not recommend. This could lead to erroneous interpretation that request
    will give ALL other molecules from DB, but in fact it will give only rest of the molecules from the buckets
    of LSH algorithm. So in order to go for lower similarities reconfigure DB and re-index it
    with new similarity limit
        :return: generator like CursorHolder object
        """
        if isinstance(self, Molecule):
            fingerprint = self.canonic_structure.fingerprint
            mol = self.canonic_structure.structure
        elif isinstance(self, MoleculeContainer):
            fingerprint = LinearFingerprint(**config.fingerprint).transform_bitset([self])[0]
            mol = self
            self = Molecule
        else:
            raise ValueError(" Only CGRtools.MoleculeContainer or CGRdb.Molecule")
        if tanimoto_limit is None:
            tanimoto_limit = config.lsh_threshold
        elif not isinstance(tanimoto_limit, float) or 1. < tanimoto_limit < 0.:
            raise f"tanimoto_limit can be float from 0 to 1, but {tanimoto_limit, type(tanimoto_limit)} was given"
        if ordered:
            request_pack = self._substructure_ordered(fingerprint, substr=mol, tanimoto_limit=tanimoto_limit)
        else:
            request_pack = self._substructure_unordered(fingerprint, substr=mol, tanimoto_limit=tanimoto_limit)
        if request_only:
            return request_pack.request
        return CursorHolder(request_pack, self._database_)

    def substructures_in_reactions(self: Union[MoleculeContainer, "Molecule"], ordered=True,
                                   is_product: Optional[bool] = None, request_only=False):
        if isinstance(self, MoleculeContainer):
            request = Molecule.substructres(self, ordered=ordered, request_only=True)
            self = Molecule
        elif isinstance(self, Molecule):
            request = self.substructres(ordered=ordered, request_only=True)
        else:
            raise TypeError("Only CGRtools.MoleculeContainer or CGRdb.Molecule can be used")

        if ordered:
            request_pack = self.__reactions_from_ordered(request_core=request)
        else:
            request_pack = self.__reactions_from_unordered(request_core=request)
        if is_product is not None:
            if isinstance(is_product, bool):
                request_pack = self.__role_filter(request_pack, is_product)
            else:
                raise ValueError("is_product can be boolean (True for product, False for reactant) or None")
        if request_only:
            return request_pack.request
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
            mol = mol.pack()
        super().__init__(_structure=mol, fingerprint=fingerprint, fingerprint_len=len(fingerprint), molecule=molecule,
                         signature=signature, is_canonic=is_canonic, smiles=smiles, _fast_mapping=fast_mapping)

    @cached_property
    def structure(self) -> MoleculeContainer:
        return MoleculeContainer.unpack(self._structure)

    def unload(self):
        self._vals_.pop(MoleculeStructure._structure, None)
        self._vals_.pop(MoleculeStructure.fingerprint, None)
        self._vals_.pop(MoleculeStructure.smiles, None)
        self.__dict__.pop('structure', None)


class NonOrganic(Entity):
    id = PrimaryKey(int, auto=True)
    signature = Required(bytes, unique=True, lazy=True)
    smiles = Required(str, lazy=True)
    _structure = Required(bytes, lazy=True)
    _fast_mapping = Required(IntArray, lazy=True)
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
    def get(cls, mol: MoleculeContainer = None, **kwargs):  # cannot use get not to overwrite pony
        if mol is None:
            return type(cls).get(cls, **kwargs)
        elif not isinstance(mol, MoleculeContainer):
            raise ValueError("CGRtools.MoleculeContainer should be provided")
        return cls.get(signature=bytes(mol))

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
    def get_reactions_with(cls, mol: MoleculeContainer, is_product: Optional[bool] = None, request_only=False):
        if not isinstance(mol, MoleculeContainer):
            raise ValueError("CGRtools.MoleculeContainer should be provided")
        request_pack = cls._exact_match_in_reactions(mol)
        if is_product is not None:
            if isinstance(is_product, bool) :
                request_pack = cls.__role_filter(request_pack, is_product)
            else:
                raise ValueError("role can be String (product, reactant) or None")
        if request_only:
            return request_pack.request
        return CursorHolder(request_pack)


class MoleculeProperties(Entity):
    id = PrimaryKey(int, auto=True)
    date = Required(datetime, default=datetime.utcnow)
    structure = Required('Molecule')
    data = Required(Json)


class MoleculeClass(Entity):
    id = PrimaryKey(int, auto=True)
    name = Required(str)
    _type = Required(int, default=0, column='type')
    structures = Set('Molecule')

class MoleculeSearchCache:

    def __init__(self, cursor: CursorHolder):
        self.date = time.asctime()
        self._tanimotos = []
        self._molecules = []
        self._structures = []
        for m, s, t in cursor:
            self._molecules.append(m)
            self._structures.append(s)
            self._tanimotos.append(t)
            s.unload()
        self.size = len(self._molecules)

    def molecules(self, page=1, pagesize=100):
        if page < 1:
            raise ValueError('page should be greater or equal than 1')
        elif pagesize < 1:
            raise ValueError('pagesize should be greater or equal than 1')

        start = (page - 1) * pagesize
        end = start + pagesize
        mis = select(x for x in Molecule if x in self._molecules[start:end]).fetch()

        if not mis:
            return []

        return mis

    def tanimotos(self, page=1, pagesize=100):
        if page < 1:
            raise ValueError('page should be greater or equal than 1')
        elif pagesize < 1:
            raise ValueError('pagesize should be greater or equal than 1')

        start = (page - 1) * pagesize
        end = start + pagesize
        return self._tanimotos[start:end]

    def __len__(self):
        return self.size

__all__ = ['Molecule', 'MoleculeStructure', 'NonOrganic']
