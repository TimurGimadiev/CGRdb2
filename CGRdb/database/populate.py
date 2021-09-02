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
from .substance import Substance
from .reaction import Reaction
from .config import db

from pony.orm import db_session
from CGRtools import smiles
from CGRtools.exceptions import InvalidAromaticRing, IncorrectSmiles, ValenceError, MappingError, EmptyMolecule
from multiprocess import Queue, Process
from tqdm import tqdm
import zipfile


def init_cgrdb(provider='postgres', user='postgres', host='localhost', password="example", database='test',
               port=5432, lsh_num_permute=64, lsh_threshold=0.7, cgr_lsh_num_permute=64, cgr_lsh_threshold=0.7,
               min_radius=1, max_radius=6, length=2048, number_active_bits=2, number_bit_pairs=4):
    from CGRdb.database.config import Config
    db.bind(provider=provider, user=user, host=host, password=password, database=database,
            port=port)
    db.generate_mapping(create_tables=True)
    db.execute("Create extension if not exists intarray;")
    Config(key="fingerprint", value={"min_radius": min_radius, "max_radius": max_radius, "length": length,
                                     "number_active_bits": number_active_bits, "number_bit_pairs": number_bit_pairs})
    Config(key="lsh_num_permute", value=lsh_num_permute)
    Config(key="lsh_threshold", value=lsh_threshold)
    Config(key="cgr_lsh_num_permute", value=cgr_lsh_num_permute)
    Config(key="cgr_lsh_threshold", value=cgr_lsh_threshold)
    db.commit()
    db.disconnect()
    db.unbind()


def upload_smi(filename, provider='postgres', user='postgres', host='localhost', password="example", database='test',
               port=5432, workers=10, reaction=False, verbose=True):

    def worker(q, reaction=reaction):
        db.bind(provider=provider, user=user, host=host, password=password, database=database,
                port=port)
        db.generate_mapping()
        while True:
            data = q.get()
            if data is None:
                break
            for i in data:
                mol, label = i
                mol = smiles(mol)
                try:
                    mol.canonicalize()
                except (InvalidAromaticRing, IncorrectSmiles, ValenceError, MappingError, EmptyMolecule):
                    if verbose:
                        print(mol, "standardization_failed")
                    continue
                if reaction:
                    for _ in range(10):
                        try:
                            with db_session():
                                Reaction(mol)
                                break
                        except Exception:
                            continue
                    else:
                        if verbose:
                            print(mol, "upload_failed")
                else:
                    subs = [(x, None) for x in mol.split()]
                    for _ in range(10):
                        try:
                            with db_session():
                                Substance(subs)
                                break
                        except Exception:
                            continue
                    else:
                        if verbose:
                            print(mol, "upload_failed")
        print("finished")

    q = Queue(maxsize=workers*2)
    num_workers = workers
    pr = [Process(target=worker,
                  args=[q]) for _ in range(num_workers)]
    [p.start() for p in pr]
    if filename.endswith(".zip"):
        with zipfile.ZipFile(filename, 'r') as zip_ref:
            zip_ref.extractall(filename[:-4])
        filename = filename[:-4]
    elif not filename.endswith(".smi"):
        raise ("""file with smiles should have .smi extension, it also can be zipped (.smi.zip). 
                        It should have following structure in each row: SMILES(required) id(optional) 
                        """)
    with open(filename) as f:
        tmp = []
        for n, row in tqdm(enumerate(f)):
            columns = row.rstrip('\n').lstrip().split(" ")
            smi = columns[0]
            idx = None
            if len(columns) > 1:
                if columns[1].startswith("|"):
                    smi += columns[1]
                    if len(columns) > 2:
                        idx = columns[2]
                else:
                    idx = columns[1]
            if not idx:
                idx = n
            tmp.append((smi, idx))
            if len(tmp) == 1000:
                q.put(tmp)
                tmp = []
        else:
            q.put(tmp)
        for i in range(num_workers):
            q.put(None)

