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
from StructureFingerprint import LinearFingerprint
from .cgr import CGR


class Reaction(Entity):
    id = PrimaryKey(int, auto=True)
    substances = Set('ReactionSubstance')
    CGR = PonyOptional("CGR")

    def __init__(self, reaction: Optional[ReactionContainer] = None, /, keep_cgr=False):
        super().__init__()
        if reaction is not None:
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
    def __postprocess_list_reactions(cls, result):
        reaction, mols, score = result
        return Reaction[reaction], [molecule.Molecule[x] for x in mols], score

    @classmethod
    def __postprocess_list_reactions_mapped(cls, result, *, initial_cgr=None):
        if not initial_cgr:
            raise ValueError("Initial CGR of reaction core is not provided")
        reaction, mols, score = result
        reaction = Reaction[reaction]
        cgr = ~reaction.structure
        if (cgr) > initial_cgr:
            return reaction, [molecule.Molecule[x] for x in mols], score

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

    @classmethod
    def similars(cls, reaction: ReactionContainer, ordered: bool = True, fix_roles: bool = True,
                 mapping: bool = False, request_only=False):
        if not isinstance(reaction, ReactionContainer):
            raise TypeError("CGRtools.ReactionContainer expected as input")
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
            return CursorHolder(RequestPack(request, cls.__postprocess_list_reactions,
                                            prefetch_map=(None, None, None)))
        else:
            cgr = ~reaction
            core = cgr.substructure(cgr.center_atoms)
            return CursorHolder(
                RequestPack(request, partial(cls.__postprocess_list_reactions_mapped,
                                             initial_cgr=core), prefetch_map=(None, None, None)))

    @classmethod
    def substructures(cls, reaction: ReactionContainer, ordered: bool = True, fix_roles: bool = True,
                      mapping: bool = False, request_only=False):
        if not isinstance(reaction, ReactionContainer):
            raise TypeError("CGRtools.ReactionContainer expected as input")
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
            return CursorHolder(RequestPack(request, cls.__postprocess_list_reactions,
                                            prefetch_map=(None, None, None)))
        else:
            cgr = ~reaction
            core = cgr.substructure(cgr.center_atoms)
            return CursorHolder(RequestPack(request,
                                            partial(cls.__postprocess_list_reactions_mapped,
                                                    initial_cgr=core), prefetch_map=(None, None, None)))

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
