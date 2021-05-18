
from CachedMethods import cached_property
from CGRtools.containers import ReactionContainer, MoleculeContainer, QueryContainer
from collections import defaultdict
from datetime import datetime
from itertools import product
from pickle import dumps
from pony.orm import PrimaryKey, Required, Optional, Set, Json, select, IntArray, FloatArray, composite_key, raw_sql
from typing import Optional as tOptional
from .config import db


class Reaction(db.Entity):
    id = PrimaryKey(int, auto=True)
    Substances = Set('ReactionSubstance')
    Fingerprint = Optional(IntArray)


class ReactionSubstance(db.Entity):
    id = PrimaryKey(int, auto=True)
    mapping = Required(Json)
    is_product = Required(bool)
    reaction = Required('Reaction')
    substance = Required('Substance')


#class ReactionIndex(db.Entity):
#    id = PrimaryKey(int, auto=True)
#    reaction = Required('Reaction')
#    signature = Required(bytes, unique=True,  column='signature')
#    fingerprint = Required(IntArray, optimistic=False, index=False, volatile=True, column='fingerprint')
#    structures = Required(IntArray, optimistic=False, index=False,  volatile=True, column='structures')



__all__ = ['Reaction', 'ReactionSubstance']
