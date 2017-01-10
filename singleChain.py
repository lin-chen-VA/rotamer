#!/usr/bin/python

from Bio.PDB import *
import numpy as np
import math
import os
import sys

def getDihedral(p):
    ''' Reference: http://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python'''
    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

def getChains(fileName):
    ''' Get all chains in a protein

    Args:
        name (string), mmcif file name

    Return:
        chains (generator), all chains in the protein
    '''
    # Create a mmcif parser
    parser = MMCIFParser();
    # Create a structure object
    name = os.path.splitext(os.path.basename(fileName))[0];
    structure = parser.get_structure(name, fileName);
    # Get all chains
    for chain in structure.get_chains():
        yield chain;

def getResidues(chain):
    return list(chain.get_residues());

def getSCDist(atomDict, resName):
    ''' Get the distance between CA and the mass center of the side-chain

        Args:
            residue (residue)
    '''
    massCenter = getSCMassCenter(atomDict, resName);
    if len(massCenter) == 0: # no atoms
        return -1;
    return getDist(atomDict['CA'], massCenter);

def getBlockDist(atomDict, resName):
    ''' Get the distance between CA and the mass center of the side-chain

        Args:
            residue (residue)
    '''
    massCenter = getBlockMassCenter(atomDict, resName);
    if len(massCenter) == 0: # no atoms
        return -1;
    return getDist(atomDict['CA'], massCenter);

def getAtomMass(atomName):
    t = atomName[0];
    mN = 14.0067;
    mC = 12.0107;
    mO = 15.9994;
    mS = 32.065;
    if t == 'N':
        return mN;
    if t == 'C':
        return mC;
    if t == 'O':
        return mO;
    if t == 'S':
        return mS;

def getBlockMassCenter(atomDict, resName):
    ''' Get the block center

        Ref: http://formulas.tutorvista.com/physics/center-of-mass-formula.html
    '''
    l = [];
    if resName == 'ALA':
        atomList = ['CB'];
    if resName == 'ARG':
        atomList = ['NE', 'CZ', 'NH1', 'NH2'];
    if resName == 'ASN':
        atomList = ['CG', 'OD1', 'ND2'];
    if resName == 'ASP':
        atomList = ['CG', 'OD1', 'OD2'];
    if resName == 'CYS':
        atomList = ['CB', 'SG'];
    if resName == 'GLU':
        atomList = ['CD', 'OE1', 'OE2'];
    if resName == 'GLN':
        atomList = ['CD', 'OE1', 'NE2'];
    if resName == 'GLY':
        atomList = [];
    if resName == 'HIS':
        atomList = ['CG', 'ND1', 'CD2', 'CE1', 'NE2'];
    if resName == 'ILE':
        atomList = ['CD1'];
    if resName == 'LEU':
        atomList = ['CG', 'CD1', 'CD2'];
    if resName == 'LYS':
        atomList = ['CE', 'NZ'];
    if resName == 'MET':
        atomList = ['SD', 'CE'];
    if resName == 'PHE':
        atomList = ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'];
    if resName == 'PRO':
        atomList = ['CB', 'CG', 'CD'];
    if resName == 'SER':
        atomList = ['OG'];
    if resName == 'THR':
        atomList = ['CB', 'OG1', 'CG2'];
    if resName == 'TRP':
        atomList = ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'];
    if resName == 'TYR':
        atomList = ['OH'];
    if resName == 'VAL':
        atomList = ['CB', 'CG1', 'CG2'];
    if len(atomList) == 0:
        return l;
    massTemp = 0;
    xTemp = 0;
    yTemp = 0;
    zTemp = 0;
    for a in atomList:
        coord = atomDict[a];
        m = getAtomMass(a);
        massTemp += m;
        xTemp += coord[0]*m;
        yTemp += coord[1]*m;
        zTemp += coord[2]*m;
    return [xTemp/massTemp, yTemp/massTemp, zTemp/massTemp];

def getSCMassCenter(atomDict, resName):
    ''' Get the residue mass center

        Ref: http://formulas.tutorvista.com/physics/center-of-mass-formula.html
    '''
    l = [];
    if resName == 'ALA':
        atomList = ['CB'];
    if resName == 'ARG':
        atomList = ['CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'];
    if resName == 'ASN':
        atomList = ['CB', 'CG', 'OD1', 'ND2'];
    if resName == 'ASP':
        atomList = ['CB', 'CG', 'OD1', 'OD2'];
    if resName == 'CYS':
        atomList = ['CB', 'SG'];
    if resName == 'GLU':
        atomList = ['CB', 'CG', 'CD', 'OE1', 'OE2'];
    if resName == 'GLN':
        atomList = ['CB', 'CG', 'CD', 'OE1', 'NE2'];
    if resName == 'GLY':
        atomList = [];
    if resName == 'HIS':
        atomList = ['CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'];
    if resName == 'ILE':
        atomList = ['CB', 'CG1', 'CG2', 'CD1'];
    if resName == 'LEU':
        atomList = ['CB', 'CG', 'CD1', 'CD2'];
    if resName == 'LYS':
        atomList = ['CB', 'CG', 'CD', 'CE', 'NZ'];
    if resName == 'MET':
        atomList = ['CB', 'CG', 'SD', 'CE'];
    if resName == 'PHE':
        atomList = ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'];
    if resName == 'PRO':
        atomList = ['CB', 'CG', 'CD'];
    if resName == 'SER':
        atomList = ['CB', 'OG'];
    if resName == 'THR':
        atomList = ['CB', 'OG1', 'CG2'];
    if resName == 'TRP':
        atomList = ['CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'];
    if resName == 'TYR':
        atomList = ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'];
    if resName == 'VAL':
        atomList = ['CB', 'CG1', 'CG2'];
    if len(atomList) == 0:
        return l;
    massTemp = 0;
    xTemp = 0;
    yTemp = 0;
    zTemp = 0;
    for a in atomList:
        coord = atomDict[a];
        m = getAtomMass(a);
        massTemp += m;
        xTemp += coord[0]*m;
        yTemp += coord[1]*m;
        zTemp += coord[2]*m;
    return [xTemp/massTemp, yTemp/massTemp, zTemp/massTemp];

def getDist(coord_1, coord_2):
    dist = 0;
    for i in xrange(3):
        dist += (coord_1[i] - coord_2[i])*(coord_1[i] - coord_2[i]);
    return math.sqrt(dist);

def getRotamer(residues, name):
    ''' Collect the data for each specific residue

        Args:

            residues (list), residue list
            name (string), residue name
    '''
    with open(name+'.csv', 'a') as outputFile:
        residueSize = len(residues);
        for index, residue in enumerate(residues):
            # Ignore the first residue and the last residue
            if index == 0 or index == residueSize-1:
                continue;
            # Calculate phi, psi and chi for the specific residues
            if residue.get_resname() == name:
                angles = getPhiPsiChi(residues, index);
                if len(angles) == 0:
                    continue;
                # output protein id, chain id, residue index, residue name
                line = residue.get_parent().get_parent().get_parent().get_id()+','+str(residue.get_parent().get_id())+','+str(residue.get_id()[1])+','+residue.get_resname();
                # calculate dist between CA and side-chain center, dist between CA and block center
                atomDict = getResidueDict(residue);
                resName = residue.get_resname();
                line = line+','+str(getSCDist(atomDict, resName));
                line = line+','+str(getBlockDist(atomDict, resName));
                # output phi, psi, chi_1
                for angle in angles:
                    line = line+','+str(angle);
                line = line+'\n';
                outputFile.write(line);

def getResidueDict(residue):
    ''' Convert residue to a dict, key is atom name and value is coordinates '''
    rDict = {};
    atoms = list(residue.get_atom());
    for atom in atoms:
        rDict[atom.get_id()] = np.array(atom.get_coord());
    return rDict;

def getPhiPsiChi(residues, index):
    angles = [];
    if not checkResidue(residues[index-1], residues[index-1].get_resname()):
        return angles;
    prevAtoms = getResidueDict(residues[index-1]); # get last residue
    if not checkResidue(residues[index], residues[index].get_resname()):
        return angles;
    currAtoms = getResidueDict(residues[index]); # get the current residue
    if not checkResidue(residues[index+1], residues[index+1].get_resname()):
        return angles;
    nextAtoms = getResidueDict(residues[index+1]); # get the next residue
    atom_Cp = prevAtoms['C']; # C_i-1
    atom_N = currAtoms['N']; # N_i
    atom_CA = currAtoms['CA']; # CA_i
    atom_Cc = currAtoms['C']; # C_i
    atom_Nn = nextAtoms['N']; # N_i+1
    phi = getDihedral([atom_Cp, atom_N, atom_CA, atom_Cc]); # phi
    psi = getDihedral([atom_N, atom_CA, atom_Cc, atom_Nn]); # psi
    angles.append(getAngle(phi));
    angles.append(getAngle(psi));
    resName = residues[index].get_resname();
    atomList = [];
    if resName == 'ARG': # X1, X2, X3, X4
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
        atomList.append(currAtoms['CD']);
        atomList.append(currAtoms['NE']);
        atomList.append(currAtoms['CZ']);
    if resName == 'ASN': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
    if resName == 'ASP': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
    if resName == 'CYS': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['SG']);
    if resName == 'GLU': # X1, X2
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
        atomList.append(currAtoms['CD']);
    if resName == 'GLN': # X1, X2
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
        atomList.append(currAtoms['CD']);
    if resName == 'HIS': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
    if resName == 'ILE': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG1']);
    if resName == 'LEU': # X1, X2
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
        atomList.append(currAtoms['CD1']);
    if resName == 'LYS': # X1, X2, X3, X4
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
        atomList.append(currAtoms['CD']);
        atomList.append(currAtoms['CE']);
        atomList.append(currAtoms['NZ']);
    if resName == 'MET': # X1, X2, X3
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
        atomList.append(currAtoms['SD']);
        atomList.append(currAtoms['CE']);
    if resName == 'PHE': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
    if resName == 'PRO': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
    if resName == 'SER': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['OG']);
    if resName == 'THR': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['OG1']);
    if resName == 'TRP': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
    if resName == 'TYR': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
    if resName == 'VAL': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG1']);
    angles = angles + getChi(atomList);
    return angles;

def getCoord(atom):
    ''' Convert atom coordinates to numpy array '''
    return np.array(atom.get_coord());

def getAngle(angle):
    ''' Convert -pi ~ pi to 0 to 2*pi '''
    if angle < 0:
        angle = 360 + angle;
    return angle;

def getChi(atomList):
    ''' Get all chi angles of the side-chain of a residue '''
    chi = [];
    if len(atomList) < 4:
        return chi;
    for i in xrange(len(atomList)-3):
        chi.append(getAngle(getDihedral([atomList[i], atomList[i+1], atomList[i+2], atomList[i+3]])));
    return chi;

def getAtom(atoms, name):
    for atom in atoms:
        if atom.get_id() == name:
            return atom;

def checkResidue(residue, name):
    ''' Validate residue
    Args:
        residue (Bio.PDB.Residue), residue
        name (string), residue name

    Return:
        True/False
    '''
    if not residue.has_id('N'):
        return False;
    if not residue.has_id('CA'):
        return False;
    if not residue.has_id('C'):
        return False;
    if not residue.has_id('O'):
        return False;
    if name == 'ALA':
        if not residue.has_id('CB'):
            return False;
    if name == 'ARG':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('CD'):
            return False;
        if not residue.has_id('NE'):
            return False;
        if not residue.has_id('CZ'):
            return False;
        if not residue.has_id('NH1'):
            return False;
        if not residue.has_id('NH2'):
            return False;
    if name == 'ASN':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('OD1'):
            return False;
        if not residue.has_id('ND2'):
            return False;
    if name == 'ASP':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('OD1'):
            return False;
        if not residue.has_id('OD2'):
            return False;
    if name == 'CYS':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('SG'):
            return False;
    if name == 'GLU':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('CD'):
            return False;
        if not residue.has_id('OE1'):
            return False;
        if not residue.has_id('OE2'):
            return False;
    if name == 'GLN':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('CD'):
            return False;
        if not residue.has_id('OE1'):
            return False;
        if not residue.has_id('NE2'):
            return False;
    if name == 'HIS':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('ND1'):
            return False;
        if not residue.has_id('CD2'):
            return False;
        if not residue.has_id('CE1'):
            return False;
        if not residue.has_id('NE2'):
            return False;
    if name == 'ILE':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG1'):
            return False;
        if not residue.has_id('CG2'):
            return False;
        if not residue.has_id('CD1'):
            return False;
    if name == 'LEU':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('CD1'):
            return False;
        if not residue.has_id('CD2'):
            return False;
    if name == 'LYS':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('CD'):
            return False;
        if not residue.has_id('CE'):
            return False;
        if not residue.has_id('NZ'):
            return False;
    if name == 'MET':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('SD'):
            return False;
        if not residue.has_id('CE'):
            return False;
    if name == 'PHE':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('CD1'):
            return False;
        if not residue.has_id('CD2'):
            return False;
        if not residue.has_id('CE1'):
            return False;
        if not residue.has_id('CE2'):
            return False;
        if not residue.has_id('CZ'):
            return False;
    if name == 'PRO':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('CD'):
            return False;
    if name == 'SER':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('OG'):
            return False;
    if name == 'THR':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('OG1'):
            return False;
        if not residue.has_id('CG2'):
            return False;
    if name == 'TRP':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('CD1'):
            return False;
        if not residue.has_id('CD2'):
            return False;
        if not residue.has_id('NE1'):
            return False;
        if not residue.has_id('CE2'):
            return False;
        if not residue.has_id('CE3'):
            return False;
        if not residue.has_id('CZ2'):
            return False;
        if not residue.has_id('CZ3'):
            return False;
        if not residue.has_id('CH2'):
            return False;
    if name == 'TYR':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('CD1'):
            return False;
        if not residue.has_id('CD2'):
            return False;
        if not residue.has_id('CE1'):
            return False;
        if not residue.has_id('CE2'):
            return False;
        if not residue.has_id('CZ'):
            return False;
        if not residue.has_id('OH'):
            return False;
    if name == 'VAL':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG1'):
            return False;
        if not residue.has_id('CG2'):
            return False;
    return True;

def getCifFiles(directory):
    ''' Get all cif files in a specific directory
    
        Args:
            directory (string), a folder name
            
        Return:
            cif file list (generator)
        '''
    for path, subdirs, files in os.walk(directory):
            for name in files:
                fileTemp, extension = os.path.splitext(name);
                if extension == '.cif':
                    yield os.path.join(path, name)

def getResidueList():
    ''' Set a list of target residues

        Return:
            string list, residue list
    '''
    #l = ['THR'];
    l = ['ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'];
    return l;

def clearOutputs(l):
    for f in l:
        os.remove(f+'.csv');

def processFile(fileName, residueName):
    ''' Process a specific residue in a cif file

        Args:
            fileName (string), cif file name
            residueName (string), specific residue name
    '''
    # Get the chains in a protein
    chains = getChains(fileName);
    # Process each chain in the protein
    for chain in chains:
        residues = getResidues(chain);
        getRotamer(residues, residueName);
        break; # process the first chain in a protein model

def main(residueName):
    ''' Generate dihedral angles for a specific residue 
    
        Args:
            residueName (string), target residue
    '''
    # Get all cif files in a directory
    cifFiles = getCifFiles(sys.argv[1]);
    count = 0;

    for cifFile in cifFiles:
        count += 1;
        print str(count)+'th: '+cifFile;
        try:
            processFile(cifFile, residueName);
        except Exception, e:
            print str(e);

if __name__ == '__main__':
    ''' Generate the distance between CA and side-chain centre, dihedral angles for each residue
    	python singleChain.py dirName
    '''
    l = getResidueList();# get the target residue list
    for r in l:
        main(r);

