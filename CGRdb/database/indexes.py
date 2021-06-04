from .config import Entity, Config
from pony.orm import Required, PrimaryKey, IntArray, db_session
from ..utils import MinHashLSH
from datasketch import MinHash
from tqdm import tqdm
from . import config
from datasketch.lsh import _optimal_param
import json
from dataclasses import dataclass
from typing import Any
from psycopg2.errors import DuplicateCursor


class MoleculeSimilarityIndex(Entity):
    band = Required(int)
    key = Required(int, size=64)
    records = Required(IntArray)
    PrimaryKey(band, key)


def molecule_similatiryidx_create(db):
    num_permute = config.lsh_num_permute or 64
    threshold = config.lsh_threshold or 0.7
    db.execute(f"""DELETE FROM MoleculeSimilarityIndex WHERE 1=1""")
    db.commit()
    lsh = MinHashLSH(threshold=threshold, num_perm=num_permute, hashfunc=hash)
    b, r = _optimal_param(threshold, num_permute, 0.5, 0.5)
    hashranges = [(i * r, (i + 1) * r) for i in range(b)]
    Config(key="hashranges", value=json.dumps(hashranges))
    with db_session:
        db.execute(f"""DELETE FROM MoleculeSimilarityIndex WHERE 1=1""")
        db.commit()
        for idx, fingerprint in tqdm(db.execute(f'SELECT id, fingerprint FROM MoleculeStructure')):
            h = MinHash(num_perm=num_permute, hashfunc=hash)
            h.update_batch(fingerprint)
            lsh.insert(idx, h, check_duplication=False)
        for n, ht in enumerate(lsh.hashtables, 1):
            for key, value in ht._dict.items():
                db.insert(MoleculeSimilarityIndex, band=n, key=key, records=list(value))
            db.commit()


class CursorHolder:
    used_names = set()

    def __init__(self, request, db, function):
        self.function = function
        self.connection = db.get_connection()
        self.request = request
        self.name = f"{hash(self.request)}"
        self.cursor = self.connection.cursor(self.name, withhold=True)
        if self.name in self.used_names:
            self.flush()
            self.cursor = self.connection.cursor(self.name, withhold=True)
        self.used_names.add(self.name)
        self.cursor.execute(request)

    def __iter__(self):
        return self

    def __next__(self):
        res = self.cursor.fetchone()
        if res is not None:
            return self.function(res)
        else:
            self.flush()
            raise StopIteration

    def flush(self):
        self.cursor.close()
        self.used_names.remove(self.name)

# class ReactionSimilairityIndex(Entity):
#     band = Required(int)
#     key = Required(int, size=64)
#     records = Required(IntArray)
#     PrimaryKey(band, key)
#
#
# def reaction_similatiryidx_create(db):
#     num_permute = config.lsh_num_permute or 64
#     threshold = config.lsh_threshold or 0.7
#     db.execute(f"""DELETE FROM ReactionSimilairityIndex WHERE 1=1""")
#     db.commit()
#     lsh = MinHashLSH(threshold=threshold, num_perm=num_permute, hashfunc=hash)
#     b, r = _optimal_param(threshold, num_permute, 0.5, 0.5)
#     hashranges = [(i * r, (i + 1) * r) for i in range(b)]
#     Config(key="hashranges", value=json.dumps(hashranges))
#     with db_session:
#         db.execute(f"""DELETE FROM MoleculeSimilarityIndex WHERE 1=1""")
#         db.commit()
#         for idx, fingerprint in tqdm(db.execute(f'SELECT id, fingerprint FROM reaction')):
#             h = MinHash(num_perm=num_permute, hashfunc=hash)
#             h.update_batch(fingerprint)
#             lsh.insert(idx, h, check_duplication=False)
#         for n, ht in enumerate(lsh.hashtables, 1):
#             for key, value in ht._dict.items():
#                 db.insert(MoleculeSimilarityIndex, band=n, key=key, records=list(value))
#             db.commit()


__all__ = ['MoleculeSimilarityIndex', 'molecule_similatiryidx_create']