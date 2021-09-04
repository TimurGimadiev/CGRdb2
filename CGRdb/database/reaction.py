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
from functools import cached_property
from typing import Optional, Union, List, Dict, Tuple
import time

from CGRtools.containers import ReactionContainer
from CGRtools.exceptions import MappingError
from pony.orm import PrimaryKey, Required, Optional as PonyOptional, Set, Json, select, IntArray

from .config import Entity, config
from . import substance
from ..utils import validate_molecule
from . import molecule
from .indexes import CursorHolder, RequestPack
from functools import partial
from StructureFingerprint import LinearFingerprint
from .cgr import CGR


class Reaction(Entity):
    id = PrimaryKey(int, auto=True)
    substances = Set('ReactionSubstance')
    CGR = PonyOptional("CGR")
    reaction_center = PonyOptional(str, lazy=True)
    reaction_index = PonyOptional("ReactionIndex")
    metadata = Set("ReactionConditions")
    classes = Set("ReactionClass")

    def __init__(self, reaction: Optional[ReactionContainer] = None, /, r_index=True):
        if reaction is not None:
            try:
                reaction_cgr = ~reaction
                rc = reaction_cgr.substructure(reaction_cgr.center_atoms)
                super().__init__(reaction_center=str(rc))
            except MappingError:
                super().__init__()
            if any(x.connected_components_count > 1 for x in reaction.molecules()):
                raise ValueError('Each CGRtools.MoleculeContainer should contain only one sole graph(molecule)')
            for x in reaction.reactants:
                ReactionSubstance(substance.Substance(((x, None),)), self, False)
            for x in reaction.products:
                ReactionSubstance(substance.Substance(((x, None),)), self, True)
            if cgr:
                CGR(reaction, reaction=self)

    @cached_property
    def structure(self):
        reactants = []
        products = []
        for i in self.substances:
            if i.is_product:
                products.append(i.structure)
            else:
                reactants.append(i.structure)
        reaction = ReactionContainer(reactants=reactants, products=products)
        return reaction

    @cached_property
    def fingerprint(self):
        return self.calculate_fingerprint(self.structure)

    @classmethod
    def calculate_rc_fingerprint(cls, reaction: ReactionContainer):
        if isinstance(ReactionContainer):
            cgr = ~reaction
            rc = cgr.substructure(cgr.center_atoms)
            lfp = LinearFingerprint(min_radius=1, max_radius=4)
            return lfp.transform_bitset([rc])[0]
        else:
            raise ValueError("CGRtools.ReactionContainer is needed")

    @classmethod
    def __postprocess_exact_reaction(cls, result):
        return Reaction[result]

    @classmethod
    def __postprocess_exact_reaction_mapped(cls, result, *, initial=None):
        if not initial:
            raise ValueError("Initial reaction is not provided")
        found = Reaction[result].structure
        if (~found) == (~initial):
            return Reaction[result]

    @classmethod
    def __postprocess_list_reactions(cls, result):
        reaction, mols, score = result
        return Reaction[reaction], [molecule.Molecule[x] for x in mols], score

    @classmethod
    def __postprocess_list_reactions_mapped(cls, result, *, initial_cgr=None):
        if not initial_cgr:
            raise ValueError("Initial CGR of reaction core is not provided")
        reaction, mols, score = result
        reaction_center = Reaction[reaction].reaction_center
        if reaction_center == initial_cgr:
            return Reaction[reaction], [molecule.Molecule[x] for x in mols], score

    @classmethod
    def get(cls, reaction: ReactionContainer = None, mapping: bool = False, fix_roles: bool = True,
                         request_only: bool = False, **kwargs):
        if reaction is None:
            return type(cls).get(cls, **kwargs)
        if not isinstance(reaction, ReactionContainer):
            raise TypeError("CGRtools.ReactionContainer expected as input")
        reactants = []
        products = []
        for mol in reaction.reactants:
            if validate_molecule(mol):
                reactants.append(molecule.Molecule.reactions(mol, request_only=True,
                                                             is_product=False if fix_roles else None))
            else:
                reactants.append(molecule.NonOrganic.get_reactions_with(mol, request_only=True,
                                                                        is_product=False if fix_roles else None))
        for mol in reaction.products:
            if validate_molecule(mol):
                products.append(molecule.Molecule.reactions(mol, request_only=True,
                                                            is_product=True if fix_roles else None))
            else:
                reactants.append(molecule.NonOrganic.get_reactions_with(mol, request_only=True,
                                                                        is_product=True if fix_roles else None))
        request = "("+") INTERSECT (".join(reactants + products)+")"
        if request_only:
            return request
        if not mapping:
            return CursorHolder(
                RequestPack(request, cls.__postprocess_exact_reaction,
                            prefetch_map=(Reaction, 0, [Reaction.structure])))
        else:
            return CursorHolder(
                RequestPack(request, partial(cls.__postprocess_exact_reaction_mapped,
                                             initial=reaction), prefetch_map=()))

    #def structurally_same(self, mapping=False):
    #    return self.get_by_structure(self.structure, mapping)

    def similars(self: Union[ReactionContainer, "Reaction"], ordered: bool = True, fix_roles: bool = True,
                 mapping: bool = False, request_only=False):
        if isinstance(self, ReactionContainer):
            reaction = self
            self = Reaction
        elif isinstance(self, Reaction):
            reaction = self.structure
        else:
            raise TypeError("CGRtools.ReactionContainer or CGRdb.Reaction expected as input")
        reactants = []
        for mol in reaction.reactants:
            if validate_molecule(mol):
                reactants.append(
                    molecule.Molecule.similars_in_reactions(mol, ordered=ordered,
                                                            is_product=False if fix_roles else None,
                                                            request_only=True))
        products = []
        for mol in reaction.products:
            if validate_molecule(mol):
                products.append(
                    molecule.Molecule.similars_in_reactions(mol, ordered=ordered,
                                                            is_product=True if fix_roles else None,
                                                            request_only=True))
        if not reactants and not products:
            raise ValueError("This reaction consist only from molecules with Non_organic atoms or ions,"
                             " similarity search is not available for them")
        if ordered:
            request_core = "(" + ") UNION ALL (".join(reactants + products) + ") "
            request = f"""SELECT s.reaction_id, s.x, t tanimoto
                                      FROM (SELECT a.reaction_id, array_agg(a.molecule_id) x, sum(a.tanimoto) t
                                           FROM ({request_core}) as a
                                           GROUP BY a.reaction_id) s
                                      ORDER BY tanimoto DESC;
                                   """
        else:
            request_core = "("+") UNION (".join(reactants + products)+") "
            request = f"""SELECT s.reaction_id, s.x, array_length(s.x,1) l 
                          FROM (SELECT a.reaction_id, array_agg(a.molecule_id) x 
                               FROM ({request_core}) as a
                               GROUP BY a.reaction_id) s
                          ORDER BY l DESC;
                       """
        if request_only:
            return request
        if not mapping:
            return CursorHolder(RequestPack(request, self.__postprocess_list_reactions,
                                            prefetch_map=(None, None, None)))
        else:
            reaction_cgr = ~reaction
            core = str(reaction_cgr.substructure(reaction_cgr.center_atoms))
            return CursorHolder(
                RequestPack(request, partial(self.__postprocess_list_reactions_mapped,
                                             initial_cgr=core), prefetch_map=(None, None, None)))

    def substructures(self: Union[ReactionContainer, "Reaction"], ordered: bool = True, fix_roles: bool = True,
                      mapping: bool = False, request_only=False):
        if isinstance(self, ReactionContainer):
            reaction = self
            self = Reaction
        elif isinstance(self, Reaction):
            reaction = self.structure
        else:
            raise TypeError("CGRtools.ReactionContainer or CGRdb.Reaction expected as input")
        reactants = []
        for mol in reaction.reactants:
            if validate_molecule(mol):
                reactants.append(
                    molecule.Molecule.substructures_in_reactions(mol,
                                                                 ordered=ordered,
                                                                 is_product=False if fix_roles else None,
                                                                 request_only=True))
        products = []
        for mol in reaction.products:
            if validate_molecule(mol):
                products.append(
                    molecule.Molecule.substructures_in_reactions(mol,
                                                                 ordered=ordered,
                                                                 is_product=True if fix_roles else None,
                                                                 request_only=True))
        if not reactants and not products:
            raise ValueError("This reaction consist only from molecules with Non_organic atoms or ions,"
                             " similarity search is not available for them")
        request_core = "(" + ") UNION ALL (".join(reactants + products) + ") "
        if ordered:
            request = f"""SELECT s.reaction_id, s.x, t tanimoto
                                      FROM (SELECT a.reaction_id, array_agg(a.molecule_id) x, sum(a.tanimoto) t
                                           FROM ({request_core}) as a
                                           GROUP BY a.reaction_id) s
                                      ORDER BY tanimoto DESC;
                                   """
        else:
            request = f"""SELECT s.reaction_id, s.x, array_length(s.x,1) l 
                          FROM (SELECT a.reaction_id, array_agg(a.molecule_id) x 
                               FROM ({request_core}) as a
                               GROUP BY a.reaction_id) s
                          ORDER BY l DESC;
                       """
        if request_only:
            return request
        if not mapping:
            return CursorHolder(RequestPack(request, self.__postprocess_list_reactions,
                                            prefetch_map=(None, None, None)))
        else:
            reaction_cgr = ~reaction
            core = str(reaction_cgr.substructure(reaction_cgr.center_atoms))
            return CursorHolder(RequestPack(request,
                                            partial(self.__postprocess_list_reactions_mapped,
                                                    initial_cgr=core), prefetch_map=(None, None, None)))

    @classmethod
    def find_similar(cls, structure) -> "ReactionSearchCache":
        cursor = cls.similars(structure)
        return ReactionSearchCache(cursor)


class ReactionSubstance(Entity):
    id = PrimaryKey(int, auto=True)
    is_product = Required(bool)
    _mapping = PonyOptional(Json, column='mapping', nullable=True)
    reaction = Required('Reaction')
    substance = Required('Substance')

    def __init__(self, substance: 'substance.Substance', reaction: Reaction, is_product: bool, /, *,
                 mapping: Union[Dict[int, int], List[Tuple[int, int]], None] = None):
        if isinstance(mapping, dict):
            if not all(isinstance(x, int) for x in mapping.items() for x in x):
                raise TypeError('Mapping is dict that contains following structure: '
                                '{int (atom of substance) : int (atom of Substance)}')
            mapping = list(mapping.items())
        elif isinstance(mapping, (tuple, list)):
            if not all(isinstance(x, (tuple, list)) and len(x) == 2 and all(isinstance(x, int) for x in x)
                       for x in mapping):
                raise TypeError('For the tuple pairs of integers following constructions should be applied: '
                                '[(int (atom of substance), int (atom of Substance)),...]')
        elif mapping is not None:
            raise TypeError('Mapping of Reaction to Substance expected')
        super().__init__(reaction=reaction, substance=substance, is_product=is_product, _mapping=mapping)

    @cached_property
    def structure(self):
        if self.mapping:
            return self.substance.structure.remap(self.mapping, copy=True)
        return self.substance.structure

    @cached_property
    def mapping(self):
        return dict(self._mapping) if self._mapping else {}


class ReactionSearchCache:

    def __init__(self, cursor: CursorHolder):
        self.date = time.asctime()
        self._tanimotos = []
        self._reactions = []
        self._structures = []
        for r, s, t in cursor:
            self._reactions.append(r)
            self._structures.append(s)
            self._tanimotos.append(t)
        self.size = len(self._reactions)

    def reactions(self, page=1, pagesize=100):
        if page < 1:
            raise ValueError('page should be greater or equal than 1')
        elif pagesize < 1:
            raise ValueError('pagesize should be greater or equal than 1')

        start = (page - 1) * pagesize
        end = start + pagesize
        ris = select(x for x in Reaction if x in self._reactions[start:end]).fetch()

        if not ris:
            return []
        return ris

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

    class ReactionClass(Entity):
        id = PrimaryKey(int, auto=True)
        name = Required(str)
        type = Required(int, default=0)
        structures = Set('Reaction')

    class ReactionConditions(Entity):
        id = PrimaryKey(int, auto=True)
        data = Required(Json, lazy=True)
        structure = Required("Reaction")

    class ReactionSearchCache:

        def __init__(self, cursor: CursorHolder):
            self.date = time.asctime()
            self._tanimotos = []
            self._reactions = []
            self._structures = []
            for r, s, t in cursor:
                self._reactions.append(r)
                self._structures.append(s)
                self._tanimotos.append(t)
                s.unload()
            self.size = len(self._reactions)

        def reactions(self, page=1, pagesize=100):
            if page < 1:
                raise ValueError('page should be greater or equal than 1')
            elif pagesize < 1:
                raise ValueError('pagesize should be greater or equal than 1')

            start = (page - 1) * pagesize
            end = start + pagesize
            ris = select(x for x in Reaction if x.id in self._reactions[start:end]).fetch()

            if not ris:
                return []
            return ris

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


__all__ = ['Reaction', 'ReactionSubstance']
