import json
from types import MethodType

from datasketch import MinHash
from datasketch.lsh import _optimal_param
from pony.orm import Required, PrimaryKey, IntArray, db_session
from tqdm import tqdm

from . import config, db
from .config import Entity, Config
from ..utils import MinHashLSH
from uuid import uuid4
from collections import namedtuple

RequestPack = namedtuple('RequestPack', ['request', 'postprocess'])


class MoleculeSimilarityIndex(Entity):
    band = Required(int)
    key = Required(int, size=64)
    records = Required(IntArray)
    PrimaryKey(band, key)


def fingerprintidx_create(self):
    print("creation of index for fingerprints")
    self.execute("CREATE INDEX idx_moleculestructure ON MoleculeStructure USING GIST(fingerprint gist__intbig_ops);")
    print("creation of index for fingerprint length")
    self.execute("CREATE INDEX idx_fingerprint_len ON moleculestructure(fingerprint_len);")

def molecule_similatiryidx_create(self):
    num_permute = config.lsh_num_permute or 64
    threshold = config.lsh_threshold or 0.7
    lsh = MinHashLSH(threshold=threshold, num_perm=num_permute, hashfunc=hash)
    b, r = _optimal_param(threshold, num_permute, 0.5, 0.5)
    hashranges = [(i * r, (i + 1) * r) for i in range(b)]
    with db_session:
        self.execute(f"""DELETE FROM MoleculeSimilarityIndex WHERE 1=1""")
        self.execute(f"""DELETE FROM Config c WHERE c.key='hashranges' """)
        self.commit()
        Config(key="hashranges", value=json.dumps(hashranges))
        print("Creation of LSH")
        for idx, fingerprint in tqdm(self.execute(f'SELECT id, fingerprint FROM MoleculeStructure')):
            h = MinHash(num_perm=num_permute, hashfunc=hash)
            h.update_batch(fingerprint)
            lsh.insert(idx, h, check_duplication=False)
        print("Uploading LSH tables into DB")
        for n, ht in tqdm(enumerate(lsh.hashtables, 1)):
            for key, value in ht._dict.items():
                self.insert(MoleculeSimilarityIndex, band=n, key=key, records=list(value))
            self.commit()


def reaction_similatiryidx_create(self):
    raise NotImplementedError

    num_permute = config.lsh_num_permute or 64
    threshold = config.lsh_threshold or 0.7
    self.execute(f"""DELETE FROM MoleculeSimilarityIndex WHERE 1=1""")
    self.commit()
    lsh = MinHashLSH(threshold=threshold, num_perm=num_permute, hashfunc=hash)
    b, r = _optimal_param(threshold, num_permute, 0.5, 0.5)
    hashranges = [(i * r, (i + 1) * r) for i in range(b)]
    pass


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
        while True:
            res = next(self.buffer)
            if res:
                result = self.function(res)
                if result:
                    return result
            else:
                raise StopIteration("No more suitable objects")

    def __del__(self):
        try:
            self.cursor.close()
        except Exception:
            pass

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


db.create_sim_index = MethodType(molecule_similatiryidx_create, db)
db.create_fing_index = MethodType(fingerprintidx_create, db)


__all__ = ['MoleculeSimilarityIndex', 'molecule_similatiryidx_create']