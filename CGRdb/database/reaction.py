from functools import cached_property
from typing import Optional, Union, List, Dict, Tuple

from CGRtools.containers import ReactionContainer
from pony.orm import PrimaryKey, Required, Optional as PonyOptional, Set, Json

from .config import Entity
from . import substance
from StructureFingerprint import LinearFingerprint
from ..utils import validate_molecule
from . import molecule
from . import config
from .indexes import CursorHolder, RequestPack
from functools import partial


class Reaction(Entity):
    id = PrimaryKey(int, auto=True)
    substances = Set('ReactionSubstance')

    def __init__(self, reaction: Optional[ReactionContainer] = None, /):
        super().__init__()
        if reaction is not None:
            if any(x.connected_components_count > 1 for x in reaction.molecules()):
                raise ValueError('Each CGRtools.MoleculeContainer should contain only one sole graph(molecule)')
            for x in reaction.reactants:
                ReactionSubstance(substance.Substance(((x, None),)), self, False)
            for x in reaction.products:
                ReactionSubstance(substance.Substance(((x, None),)), self, True)

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
    def calculate_fingerprint(cls, reaction: ReactionContainer):
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
    def get(cls, reaction: ReactionContainer, mapping=False):
        if not isinstance(reaction, ReactionContainer):
            raise TypeError("CGRtools.ReactionContainer expected as input")
        reactants = []
        products = []
        for mol in reaction.reactants:
            if validate_molecule(mol):
                reactants.append(molecule.Molecule._exact_match_in_reactions(mol).request)
            else:
                reactants.append(molecule.NonOrganic._exact_match_in_reactions(mol).request)
        for mol in reaction.products:
            if validate_molecule(mol):
                products.append(molecule.Molecule._exact_match_in_reactions(mol).request)
            else:
                products.append(molecule.NonOrganic._exact_match_in_reactions(mol).request)
        request = "("+") INTERSECT (".join(reactants + products)+")"
        if not mapping:
            return CursorHolder(RequestPack(request, cls.__postprocess_exact_reaction), cls._database_)
        else:
            return CursorHolder(RequestPack(request, partial(cls.__postprocess_exact_reaction_mapped,
                                                             initial=reaction)), cls._database_)

    def structurally_same(self, mapping=False):
        return self.get(self.structure, mapping)

    @classmethod
    def similar_to(cls, reaction, ordered=True, mapping=False):
        raise NotImplemented
        if not isinstance(reaction, ReactionContainer):
            raise TypeError("CGRtools.ReactionContainer expected as input")
        reactants = []
        for mol in reaction.reactants:
            if validate_molecule(mol):
                fingerprint = LinearFingerprint(**config.fingerprint).transform_bitset([mol])[0]
                reactants.append(molecule.Molecule._search_similar_unordered_in_reactions_request(fingerprint).request)
        products = []
        for mol in reaction.products:
            if validate_molecule(mol):
                fingerprint = LinearFingerprint(**config.fingerprint).transform_bitset([mol])[0]
                products.append(molecule.Molecule._search_similar_unordered_in_reactions_request(fingerprint).request)
        if not reactants and not products:
            raise ValueError("This reaction consist only from molecules with Non_organic atoms or ions,"
                             " similarity search is not available for them")
        request = "("+") INTERSECT (".join(reactants + products)+")"
        if not mapping:
            return CursorHolder(RequestPack(request, cls.__postprocess_exact_reaction), cls._database_)
        else:
            return CursorHolder(RequestPack(request, partial(cls.__postprocess_exact_reaction_mapped,
                                                             initial=reaction)), cls._database_)


    @classmethod
    def broad_similarity(cls, reaction: [ReactionContainer, "Reaction"]):
        raise NotImplemented

    def substructure(self):
        raise NotImplemented

    def substructure_mappless(self):
        raise NotImplemented


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


__all__ = ['Reaction', 'ReactionSubstance']
