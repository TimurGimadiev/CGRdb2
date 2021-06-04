
from pony.orm import Database, PrimaryKey, Required, Set, IntArray, Optional, Json
from .config import Entity


class Config(Entity):
    _table_ = 'db_config'
    id = PrimaryKey(int, auto=True)
    name = Required(str, unique=True)
    config = Required(Json, index=False, optimistic=False)
    version = Required(str)

class ModificationStamp(Entity):
    id = PrimaryKey(int, auto=True)
    fileName = Optional(str)
    time = Required(str, unique=True)


class DataSource(Entity):
    id = PrimaryKey(int, auto=True)
    name = Required(str, unique=True)
    molecules = Set('MoleculeSourceID')
    reactions = Set('ReactionSourceID')


class MoleculeSourceID(Entity):
    source = Required('DataSource')
    molecule = Required('Molecule')
    id = PrimaryKey(source, molecule)


class ReactionSourceID(Entity):
    source = Required('DataSource')
    reaction = Required('Reaction')
    id = PrimaryKey(source, reaction)



__all__ = ['ModificationStamp', 'DataSource', 'MoleculeSourceID', 'ReactionSourceID']
