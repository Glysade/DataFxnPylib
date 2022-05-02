"""
===============
rdkit_utils.py
===============

Copyright (C) 2017-2022 Glysade, LLC

Utility functions for handling chemical structures using RDKit (www.rdkit.org)
"""

import base64
import gzip
import math
from contextlib import contextmanager
from enum import Enum
from io import BytesIO, StringIO
from typing import Iterable, Union, List, Optional

from rdkit import Chem
from rdkit.Chem import RWMol, SanitizeFlags
from rdkit.Chem.rdchem import Mol, KekulizeException, AtomValenceException
from rdkit.Geometry.rdGeometry import Point3D


class RDKitFormat(Enum):
    """
    Enumerated class for structure formats handled by RDkit interface

    Members:
        - sdf
        - pdb
        - smi
    """
    sdf = 1
    pdb = 2
    smi = 3
    sma = 4


def type_to_format(type: str) -> RDKitFormat:
    """
    Given a chemical mime type returns the equivalent :class:`RDKitFormat` Enumerated type

    :param type: Chemical mime type
    :return: RDKit format
    """

    if type in {'chemical/x-mdl-molfile', 'chemical/x-sdf', 'chemical/x-mdl-molfile-v3000'}:
        return RDKitFormat.sdf
    elif type == 'chemical/x-pdb':
        return RDKitFormat.pdb
    elif type in ['chemical/x-smarts', 'chemical/x-daylight-smarts']:
        return RDKitFormat.sma
    elif type in ['chemical/x-smiles', 'chemical/x-daylight-smiles']:
        return RDKitFormat.smi
    else:
        raise ValueError("Unknown chemical type {}".format(type))


def string_to_mol(type: RDKitFormat, mol_string: str) -> Optional[Mol]:
    """
    Converts a string to an RDKit molecule

    :param type: The structure format for the string as :class:`RDKitFormat`
    :param mol_string: Molecular string
    :return: RDkit molecule, class :class:`rdkit.Chem.rdchem.Mol`
    """

    if type == RDKitFormat.sdf:
        mol = sdf_to_mol(mol_string)
    elif type == RDKitFormat.pdb:
        mol = Chem.MolFromPDBBlock(mol_string)
    elif type == RDKitFormat.smi:
        mol = smiles_to_mol(mol_string)
    elif type == RDKitFormat.sma:
        mol = Chem.MolFromSmarts(mol_string)
    else:
        raise ValueError("Unable to convert type {} from block".format(type))

    return mol


def string_to_mols(type: RDKitFormat, mol_string: str) -> List[Mol]:
    """
    Converts a string to an RDKit molecule

    :param type: The structure format for the string as :class:`RDKitFormat`
    :param mol_string: Molecular string
    :return: RDkit molecule, class :class:`rdkit.Chem.rdchem.Mol`
    """

    mols = []
    if type == RDKitFormat.sdf:
        # if we want SD tags we need to use a supplier as Chem.MolFromMolBlock does not process SD tags
        sdf_in = BytesIO(mol_string.encode('UTF-8'))
        supplier = Chem.ForwardSDMolSupplier(sdf_in)
        mols = [m for m in supplier if m is not None]
        sdf_in.close()
        # mol = Chem.MolFromMolBlock(mol_string, True, False, False)
    elif type == RDKitFormat.pdb:
        pdb_string = ""
        for line in mol_string.splitlines():
            pdb_string += line
            if line.startswith('END'):
                pdb_string += line
                mol = Chem.MolFromPDBBlock(pdb_string)
                if mol:
                    mols.append(mol)
                pdb_string = ""
    elif type == RDKitFormat.smi:
        for line in mol_string.splitlines():
            if line:
                mol = smiles_to_mol(mol_string)
                if mol:
                    mols.append(mol)
    elif type == RDKitFormat.sma:
        for line in mol_string.splitlines():
            if line:
                mol = Chem.MolFromSmarts(mol_string)
                if mol:
                    mols.append(mol)
    else:
        raise ValueError("Unable to convert type {} from block".format(type))

    return mols


def sdf_to_mol(mol_string: str) -> Optional[Mol]:
    # if we want SD tags we need to use a supplier as Chem.MolFromMolBlock does not process SD tags
    sdf_in = BytesIO(mol_string.encode('UTF-8'))
    supplier = Chem.ForwardSDMolSupplier(sdf_in, sanitize=False)
    mol = next(supplier)
    try:
        sanitize_mol(mol)
    except Exception as ex:
        name = type(ex).__name__
        print('Got exception {} processing sdf input {}'.format(name, mol_string))
        return None
    finally:
        sdf_in.close()
    return mol
    # mol = Chem.MolFromMolBlock(mol_string, True, False, False)


def sanitize_mol(mol: Mol) -> None:
    try:
        Chem.SanitizeMol(mol)
    except AtomValenceException:
        print('Atom valence error for sdf input {}'.format(Chem.MolToSmiles(mol)))
        flags = SanitizeFlags.SANITIZE_ALL ^ SanitizeFlags.SANITIZE_PROPERTIES
        Chem.SanitizeMol(mol, flags)
    except KekulizeException:
        print('Kekulization error for molecule {}'.format(Chem.MolToSmiles(mol)))
        flags = SanitizeFlags.SANITIZE_ALL ^ SanitizeFlags.SANITIZE_KEKULIZE
        Chem.SanitizeMol(mol, flags)


def smiles_to_mol(smiles: str) -> Optional[Mol]:
    mol = Chem.MolFromSmiles(smiles)
    try:
        if not mol:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol:
                print('Failed to Sanitize smiles {}'.format(smiles))
                flags = SanitizeFlags.SANITIZE_ALL ^ SanitizeFlags.SANITIZE_PROPERTIES
                mol.UpdatePropertyCache(False)
                try:
                    Chem.SanitizeMol(mol, flags)
                except KekulizeException:
                    print('Kekulization error for smiles {}'.format(smiles))
                    flags = flags ^ SanitizeFlags.SANITIZE_KEKULIZE
                    Chem.SanitizeMol(mol, flags)
            if not mol:
                print("Failed to convert smiles {} to mol!".format(smiles))
        return mol
    except Exception as ex:
        name = type(ex).__name__
        print('Got exception {} processing smiles {}'.format(name, smiles))
        return None


def mol_to_string(type: RDKitFormat, mol: Mol, isomeric_smiles=True) -> str:
    """
    Converts an RDKit molecule to a molecular string

    :param type: The structure format for the string as :class:`RDKitFormat`
    :param mol: RDkit molecule, class :class:`rdkit.Chem.rdchem.Mol`
    :param isomeric_smiles: If True set smiles to be isomeric (default)
    :return: The molecular string
    """

    if not mol:
        raise ValueError("Molecule is None!")
    if type == RDKitFormat.sdf:
        # mol_string = Chem.MolToMolBlock(mol)
        sdf_in = StringIO()
        writer = Chem.SDWriter(sdf_in)
        try:
            writer.write(mol)
        except KekulizeException:
            writer.SetKekulize(False)
            writer.write(mol)
        writer.flush()
        mol_string = sdf_in.getvalue()
        sdf_in.close()
    elif type == RDKitFormat.pdb:
        mol_string = Chem.MolToPDBBlock()
    elif type == RDKitFormat.smi:
        mol_string = Chem.MolToSmiles(mol, isomericSmiles=isomeric_smiles)
    elif type == RDKitFormat.smi:
        mol_string = Chem.MolToSmarts(mol, isomericSmiles=isomeric_smiles)
    else:
        raise ValueError("Unknown chemical type {}".format(type))
    return mol_string


def encode_mol(type: RDKitFormat, mol: Mol) -> str:
    """
    Converts an RDKit molecule to a molecular string which is then gzipped and Base64 encoded.

    :param type: The structure format for the string as :class:`RDKitFormat`
    :param mol: RDkit molecule, class :class:`rdkit.Chem.rdchem.Mol`
    :return: The molecular string, gzipped then Base64 encoded
    """

    mol_str = mol_to_string(type, mol)
    mol_bytes = bytes(mol_str, 'utf-8')
    mol_gzip = gzip.compress(mol_bytes)
    mol_b64 = base64.b64encode(mol_gzip).decode('utf-8')
    return mol_b64


def file_to_format(file: str) -> RDKitFormat:
    """
    Elucidates the RDKit structure format from a filename

    :param file: the name of a file
    :return: The structure format for the compounds in the file as :class:`RDKitFormat`
    """

    if file.endswith(".gz"):
        file = file[:-3]
    if file[-4] != '.':
        raise ValueError("unable to find molecular file type for file {}".format(file))
    ext = file.lower()[-3:]
    if ext in ['sdf', 'mol']:
        return RDKitFormat.sdf
    elif ext in ['pdb', 'ent']:
        return RDKitFormat.pdb
    elif ext in ['smi']:
        return RDKitFormat.smi
    elif ext in ['sma']:
        return RDKitFormat.sma


@contextmanager
def mol_supplier(file: str, type: RDKitFormat = None) -> Iterable[Mol]:
    """
    Wrapper round the RDKit molecular suppliers.  Opens the file using the correct supplier and and yields RdKit molecules.

    :param file: The name of the file
    :param type: The structure format for the compounds in the file as :class:`RDKitFormat`.  If None (the default) the type will be inferred from the file name
    :return: yields RDkit molecules of class :class:`rdkit.Chem.rdchem.Mol`
    """

    fh = None
    if type is None:
        type = file_to_format(file)
    if type == RDKitFormat.sdf:
        if file.lower().endswith('.gz'):
            fh = gzip.open(file, 'r')
            supplier = Chem.ForwardSDMolSupplier(fh)
        else:
            supplier = Chem.SDMolSupplier(file)
    elif type == RDKitFormat.smi:
        assert not file.lower().endswith('.gz')
        supplier = Chem.SmilesMolSupplier(file, titleLine=False)
    else:
        raise ValueError("Unable to create supplier for type {}".format(type))
    yield supplier
    if fh:
        fh.close()


@contextmanager
def mol_writer(file: str, type: RDKitFormat = None) -> Union[Chem.SDWriter, Chem.SmilesWriter]:
    """
    Wrapper round the RDKit writers.  Returns the correct writer to a file for the molecular format

    :param file: The name of the file
    :param type: The structure format for the compounds in the file as :class:`RDKitFormat`.  If None (the default) the type will be inferred from the file name
    :return: yields the appropriate writer for the molecular type and file
    """

    fh = None
    if type is None:
        type = file_to_format(file)
    if type == RDKitFormat.sdf:
        if file.lower().endswith('.gz'):
            fh = gzip.open(file, 'wt')
            writer = Chem.SDWriter(fh)
        else:
            fh = open(file, 'w')
            writer = Chem.SDWriter(fh)
    elif type == RDKitFormat.smi:
        assert not file.lower().endswith('.gz')
        writer = Chem.SmilesWriter(file)
    yield writer
    writer.flush()
    writer.close()
    if fh:
        fh.close()
    # The SDWriter in a regular file does not always seem to be closed at this point


def mols_to_file(file: str, mols: Iterable[Mol], type: RDKitFormat = None) -> None:
    """
    Writes a list of RDKit molecules to a file
    :param file: The name of the file
    :param mols: an Iterator of :class:`rdkit.Chem.rdchem.Mol` RDKit molecules
    :param type:  The structure format for the compounds in the file as :class:`RDKitFormat`.  If None (the default) the type will be inferred from the file name
    """
    with mol_writer(file, type) as writer:
        for mol in mols:
            if mol:
                writer.write(mol)


def file_to_mols(file: str, type: RDKitFormat = None) -> List[Mol]:
    """
    Reads a list of RDKit molecules from a file

    :param file: The name of the file
    :param type: The structure format for the compounds in the file as :class:`RDKitFormat`. If None (the default) the type will be inferred from the file name
    :return: a List of :class:`rdkit.Chem.rdchem.Mol` RDKit molecules
    """
    with mol_supplier(file, type) as supplier:
        return [m for m in supplier]


def is_three_dimensional(mol: Mol) -> bool:
    """
    Return true is this mol was created from a three dimensional SD file

    :param mol:
    :return: True if the SD file header contains 3D
    """

    if not mol.HasProp('_MolFileInfo'):
        return False
    return mol.GetProp('_MolFileInfo')[20:22] == '3D'


def remove_atom_mappings(mol: Mol, remove_h: bool = True) -> Mol:
    """
    Removes all atom mappings from a molecule

    :param mol: the molecule
    """
    mol = RWMol(mol)
    atoms_to_remove = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() != 0:
            atom.SetIsotope(0)
        if atom.HasProp('molAtomMapNumber'):
            atom.ClearProp('molAtomMapNumber')
            if atom.GetAtomicNum() == 1:
                atoms_to_remove.append(atom)
    if remove_h:
        for atom in atoms_to_remove:
            mol.RemoveAtom(atom.GetIdx())
    return Mol(mol)


def remove_explicit_hydrogens(mol: Mol) -> Mol:
    """
    Removes all explicit hydrogens and hydrogen atoms from a molecule

    :param mol:
    :return:
    """

    # thought that Allchem.RemoveHs would remove all hydrogens, but it doesn't e.g
    # AllChem.RemoveHs(Chem.MolFromSmiles('Nc1nc(NC(CO)Cc2ccc(O)cc2)nc2c1ncn2C1OC(c2nn[nH]n2)C(O)C1O'))
    # see :func:`test.test_chem.test_rdkit.TestRDKit.test_remove_h`
    mol = RWMol(mol)
    atoms_to_remove = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            atoms_to_remove.append(atom)
        else:
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(0)
    for atom in atoms_to_remove:
        mol.RemoveAtom(atom.GetIdx())
    return Mol(mol)


def average_bond_length(mol: Mol) -> float:
    conf = mol.GetConformer()

    def bond_length(bond):
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        p1 = conf.GetAtomPosition(a1)
        p2 = conf.GetAtomPosition(a2)
        x = p1.x - p2.x
        y = p1.y - p2.y
        return math.sqrt(x * x + y * y)

    bond_lengths = [bond_length(b) for b in mol.GetBonds()]
    return sum(bond_lengths) / float(len(bond_lengths))


def rescale_bond_lengths(mol: Mol, bond_length: float = 1.5) -> None:
    avg_bond_len = average_bond_length(mol)
    scale = bond_length / avg_bond_len
    conf = mol.GetConformer()

    for a in range(0, mol.GetNumAtoms()):
        p = conf.GetAtomPosition(a)
        new_p = Point3D(p.x * scale, p.y * scale, 0.0)
        conf.SetAtomPosition(a, new_p)


def print_mol_information(mol: Mol) -> None:
    """
    Prints information about a molecule.  For debugging purposes

    :param mol:
    :return:
    """

    for atom in mol.GetAtoms():
        print('{} {} {} {} {}'.format(atom.GetIdx(), atom.GetSymbol(), atom.GetAtomicNum(), atom.GetIsotope(),
                                      atom.GetAtomMapNum()))
    print()
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        bt = str(bond.GetBondType())
        print('{} {} {}'.format(a1, a2, bt))
