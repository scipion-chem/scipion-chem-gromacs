# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************


"""
This module imports an already-prepared GROMACS system (topology, coordinates
[and trajectory]) into Scipion.
"""
import os, re, shutil

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from gromacs import Plugin as gromacsPlugin
from gromacs.objects import GromacsSystem

# ---------------------------------------------------------------------------
# Molecule classification tables
# ---------------------------------------------------------------------------
# These sets are intentionally generous so that classification keeps working
# for the many naming conventions produced by CHARMM-GUI, AMBER and GROMACS.

# Standard + CHARMM/AMBER protonation and terminal variants of amino acids.
STANDARD_AA = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    # CHARMM histidine / cysteine / lysine variants
    'HSD', 'HSE', 'HSP', 'HID', 'HIE', 'HIP', 'CYX', 'CYM', 'LYN', 'ASH',
    'GLH', 'ACE', 'NME', 'NMA',
}

# Water residue/molecule names.
WATER_NAMES = {
    'SOL', 'WAT', 'HOH', 'H2O', 'TIP', 'TIP3', 'TIP3P', 'TIP4', 'TIP4P',
    'TIP5', 'TIP5P', 'TP3', 'TP4', 'SPC', 'SPCE', 'SW', 'OPC', 'OPC3',
    'T3P', 'T4P',
}

# Monoatomic ion residue/molecule names (GROMACS + CHARMM + AMBER flavours).
ION_NAMES = {
    'NA', 'NA+', 'SOD', 'CL', 'CL-', 'CLA', 'K', 'K+', 'POT', 'CA', 'CA2+',
    'CAL', 'MG', 'MG2+', 'ZN', 'ZN2+', 'FE', 'FE2', 'CU', 'CU1', 'CU2',
    'CS', 'CS+', 'CES', 'RB', 'RB+', 'LI', 'LI+', 'BR', 'BR-', 'IOD', 'I',
    'I-', 'F', 'F-', 'MN', 'CO', 'NI', 'BA',
}

# Lipid residue/molecule names commonly emitted by CHARMM-GUI Membrane Builder.
LIPID_NAMES = {
    'POPC', 'POPE', 'POPS', 'POPG', 'POPA', 'POPI', 'POPI13', 'POPI14',
    'POPI15', 'POPI24', 'POPI25', 'DOPC', 'DOPE', 'DOPS', 'DOPG', 'DPPC',
    'DPPE', 'DPPG', 'DMPC', 'DMPE', 'DMPG', 'DSPC', 'DLPC', 'DLPE', 'PLPC',
    'PLPE', 'SAPC', 'SDPC', 'SOPC', 'PSM', 'SSM', 'CHL1', 'CHOL', 'ERG',
    'LPPC', 'DXPC', 'DBPC', 'DEPC', 'DDPC', 'DYPC', 'PYPC', 'OYPE', 'YOPC',
}


class GromacsImportSystem(EMProtocol):
    """
    AI Generated:

    Protocol to import an already prepared GROMACS system into Scipion.

    It takes a system built outside Scipion (e.g. with CHARMM-GUI: protein, membrane,
    ligands, solvent and ions) and creates a GromacsSystem with the structural metadata and
    the index (.ndx) file that the simulation, MMPBSA, trajectory and viewer protocols require.

    Overview
    --------
    1. Import of the input files
       - Coordinates, topology and every file it #includes, optional index and trajectory

    2. Classification of the molecules of the topology
       - Protein, ligand, lipid, ion and water

    3. Derivation of the structural metadata
       - Protein chains and lengths, primary ligand, membrane, force field and water model

    4. Index generation
       - Standard GROMACS groups plus chain, ligand, membrane and complex groups

    Input
    -----
    inputCoords:
        Coordinates of the system (.gro, or .pdb converted to .gro).

    inputTopology:
        GROMACS topology (.top). The files it #includes are resolved relative to it,
        recursively, and copied along with it.

    inputItpDir:
        Optional folder (e.g. the CHARMM-GUI "toppar" folder) where #included files not
        found next to the topology are looked for.

    inputPdb:
        Optional original .pdb used to recover the chain identifiers lost in the .gro file.

    ligandNames:
        Optional comma-separated residue/molecule names to force-classify as ligands. The
        first one becomes the primary ligand (its ID and topology are stored).

    inputIndex:
        Optional index file whose groups are merged into the generated one (e.g. the
        CHARMM-GUI index.ndx, needed by its .mdp files).

    inputTrajectory:
        Optional trajectory (.xtc, .trr or multi-model .pdb) matching the coordinates.

    Workflow
    --------
    1. Copies the input files, rewriting the topology #includes to their copied location

    2. Classifies every molecule of the "[ molecules ]" section from its name and residue
       composition (CHARMM-GUI naming conventions included); user-named ligands override
       the automatic classification

    3. Builds the protein chains (splitting molecules that merge several chains) and picks
       the primary ligand, its residue name and topology

    4. Generates the index file: make_ndx default groups, merged with the input index, plus
       per-chain groups, Membrane, one group per ligand, Protein_<ligand>, Protein_membrane
       and Protein_membrane_lig

    5. Reads the time and frame range of the trajectory, if given

    Output
    ------
    outputSystem:
        GromacsSystem containing:
        - coordinates, original structure and topology (with all its includes)
        - chain names and lengths, ligand ID and ligand topology
        - force field and water model (best-effort detection)
        - index file and, if given, the trajectory

    Key Features
    ------------
    - Direct use of CHARMM-GUI GROMACS outputs (membrane systems included)
    - Automatic molecule classification with manual ligand override
    - Complete index file for the downstream protocols

    Notes
    -----
    - The topology must contain a "[ molecules ]" section, and the forced ligand names must
      appear in it
    - The form warns about #included files that cannot be found
    - The trajectory atom count and order must match the coordinates and topology
    """
    _label = 'import system'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used. """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputCoords', params.FileParam,
                      label="Input Gromacs Coordinates (gro, pdb): ", allowsNull=False,
                      help='GROMACS coordinates file (.gro preferred, .pdb also accepted).')
        form.addParam('inputTopology', params.FileParam,
                      label="Input Gromacs Topology (top): ", allowsNull=False,
                      help='GROMACS topology file (.top). Any files it "#include"s (e.g. the '
                           '"toppar" directory generated by CHARMM-GUI) are resolved relative to '
                           'this file and copied along with it so the topology stays valid.')
        form.addParam('inputItpDir', params.FolderParam,
                      label="Topology .itp directory: ", allowsNull=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Optional directory containing the topology .itp files (e.g. the '
                           'CHARMM-GUI "toppar" folder). Only needed as a fallback: files the '
                           'topology "#include"s that are not found next to the .top are looked '
                           'up here.')
        form.addParam('inputPdb', params.FileParam,
                      label="Original PDB (chains): ", allowsNull=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Optional original .pdb structure. Chain information often lost'
                           ' in the .gro file. This file will returned lost information such as'
                           ' chains identifiers.')
        form.addParam('ligandNames', params.StringParam, allowsNull=True,
                      label="Ligand residue name(s): ",
                      help='Optional, comma-separated residue/molecule names to force-classify '
                           'as ligands (e.g. "NAD" or "NAD,HEM"). Ligands are detected '
                           'automatically, but sometimes fails. Use this to select a specific ligand to study or to pick'
                           ' which ligand is the primary one. The names must '
                           'match those in the topology "[ molecules ]" section; the first one '
                           'listed becomes the primary ligand (its ID and topology are stored).')
        form.addParam('inputIndex', params.FileParam,
                      label="Input Gromacs Index (ndx): ", allowsNull=True,
                      help='Optional index file (.ndx). Its groups are merged into the index '
                           'generated by this protocol.')
        form.addParam('inputTrajectory', params.FileParam,
                      label="Input Gromacs Trajectory (trr, xtc, pdb): ", allowsNull=True,
                      help='Optional trajectory file: GROMACS .xtc / .trr, or a multi-model '
                           '(multi-MODEL) .pdb file, e.g. an ensemble/morph of frames. Its atom '
                           'count and order must match the input coordinates/topology.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep)

    def createOutputStep(self):
        # 1) Bring every input file
        localFiles = self._importInputFiles()

        # 2) Parse the topology and classify every molecule it contains.
        forcedLigands = self._getForcedLigands()
        molecules = parseTopologyMolecules(localFiles['topology'])
        molInfo = self._collectMoleculeInfo(localFiles['topology'], localFiles['includes'])
        classified = classifyMolecules(molecules, molInfo, forcedLigands=forcedLigands)

        # 3) Derive the structural metadata (chains, ligands, membrane, FF).
        chains = self._buildChains(classified['protein'], molInfo)
        chains = self._applyPdbChainNames(chains)
        # Put any explicitly-named ligand first so it becomes the primary one.
        ligands = self._orderLigands(classified['ligand'], forcedLigands)
        lipidResnames = self._collectResnames(classified['lipid'], molInfo)
        ff, wff = detectForceFields(localFiles['topology'], classified['water'])

        self._logClassification(classified, chains, ligands, lipidResnames)

        # 4) Create the output object with all the metadata.
        outSystem = GromacsSystem(filename=os.path.relpath(localFiles['coords']))
        outSystem.setOriStructFile(os.path.relpath(localFiles['oriStruct']))
        outSystem.setTopologyFile(os.path.relpath(localFiles['topology']))
        if ff:
            outSystem.setForceField(ff)
        if wff:
            outSystem.setWaterForceField(wff)

        outSystem.setChainNames([c['id'] for c in chains])
        outSystem.setChainLengths([c['length'] for c in chains])

        ligResnames = []
        if ligands:
            primaryLig = ligands[0]
            ligResnames = self._collectResnames(ligands, molInfo)
            outSystem.setLigandID(self._pickLigandID(primaryLig, molInfo))
            ligItp = self._getMoleculeItp(primaryLig, localFiles['includes'])
            if ligItp:
                outSystem.setLigTopologyFile(os.path.relpath(ligItp))

        if localFiles.get('trajectory'):
            outSystem.setTrajectoryFile(localFiles['trajectory'])
            outSystem.readTrjInfo(protocol=self, outDir=self._getExtraPath())

        # 5) Build the index file (standard groups + chains/ligand/membrane).
        indexFile = self._generateIndexFile(outSystem, chains, lipidResnames,
                                            ligResnames, localFiles.get('index'))
        outSystem.setIndexFile(indexFile)

        self._defineOutputs(outputSystem=outSystem)

    # --------------------------- File import helpers --------------------------
    def _importInputFiles(self):
        """Copy the structure, topology (+ includes), index and trajectory into
        the protocol folder, converting the structure to .gro when needed.
        Returns a dict with the local (copied) paths."""
        files = {}
        coordsSrc = os.path.abspath(self.inputCoords.get())
        ext = os.path.splitext(coordsSrc)[1].lower()
        baseName = os.path.splitext(os.path.basename(coordsSrc))[0]

        localCoords = self._getExtraPath(os.path.basename(coordsSrc))
        shutil.copy(coordsSrc, localCoords)

        if ext == '.gro':
            files['coords'] = localCoords
            files['oriStruct'] = localCoords
        else:
            # A .gro is needed for editconf/analysis downstream; generate it.
            groFile = self._getExtraPath(f'{baseName}.gro')
            gromacsPlugin.runGromacs(self, 'gmx',
                                     f'editconf -f {os.path.abspath(localCoords)} '
                                     f'-o {os.path.abspath(groFile)}',
                                     cwd=self._getExtraPath())
            files['coords'] = groFile
            files['oriStruct'] = localCoords

        files['topology'], files['includes'] = self._importTopology()

        if self.inputIndex.get():
            localIndex = self._getExtraPath('input.ndx')
            shutil.copy(os.path.abspath(self.inputIndex.get()), localIndex)
            files['index'] = localIndex

        if self.inputTrajectory.get():
            trjSrc = os.path.abspath(self.inputTrajectory.get())
            localTrj = self._getExtraPath(os.path.basename(trjSrc))
            shutil.copy(trjSrc, localTrj)
            files['trajectory'] = localTrj

        return files

    def _importTopology(self):
        """Copy the .top file together with every file it #includes, following
        includes recursively (an included .itp can itself #include further
        files, e.g. a standard force field's forcefield.itp pulling in its own
        ffnonbonded.itp/ffbonded.itp). Preserves the relative layout (e.g. the
        CHARMM-GUI 'toppar/' directory, or a force field's own subfolder) so
        the include paths stay valid. Returns (localTopPath,
        {includeRelPath: localAbsPath}) for the topology's own (first-level)
        includes."""
        topSrc = os.path.abspath(self.inputTopology.get())
        topDir = os.path.dirname(topSrc)
        localTop = self._getExtraPath(os.path.basename(topSrc))
        shutil.copy(topSrc, localTop)

        topLevelIncludes = parseTopologyIncludes(topSrc)
        includes, seen = {}, set()
        pending = [(relInc, topDir) for relInc in topLevelIncludes]

        while pending:
            relInc, srcDir = pending.pop(0)
            if relInc in seen:
                continue
            seen.add(relInc)

            srcInc = self._resolveInclude(relInc, srcDir)
            if srcInc is None:
                self.info(f'WARNING: included topology file not found, skipping: {relInc}')
                continue
            # Keep the include path as written (relative) whenever it stays
            # inside the topology directory; otherwise flatten into the root.
            if not os.path.isabs(relInc) and not relInc.startswith('..'):
                dstInc = self._getExtraPath(relInc)
            else:
                dstInc = self._getExtraPath(os.path.basename(relInc))
            os.makedirs(os.path.dirname(dstInc), exist_ok=True)
            shutil.copy(srcInc, dstInc)

            if relInc in topLevelIncludes:
                includes[relInc] = os.path.abspath(dstInc)

            # Follow this file's own #includes, relative to its source directory.
            for nestedInc in parseTopologyIncludes(srcInc):
                pending.append((nestedInc, os.path.dirname(srcInc)))

        # Rewrite the #include paths to the absolute location of the copied files.
        rewriteTopologyIncludes(localTop, includes)

        return localTop, includes

    def _resolveInclude(self, relInc, topDir):
        """Locate an #include'd file: first next to the file that includes it,
        then (as a fallback) inside the optional user-provided .itp directory,
        then inside GROMACS's own standard force-field directory (covers
        standard force fields such as amber14sb.ff/charmm36.ff/oplsaa.ff, which
        are shipped with GROMACS itself and not next to the user's .top).
        Returns the absolute source path or None if it cannot be found."""
        candidates = [os.path.join(topDir, relInc)]
        itpDir = self.inputItpDir.get()
        if itpDir:
            candidates.append(os.path.join(itpDir, relInc))
            candidates.append(os.path.join(itpDir, os.path.basename(relInc)))
        candidates.append(os.path.join(gromacsPlugin.getTopDir(), relInc))
        for cand in candidates:
            cand = os.path.normpath(cand)
            if os.path.exists(cand):
                return cand
        return None

    # --------------------------- Metadata helpers ----------------------------
    def _collectMoleculeInfo(self, topFile, includes):
        """Parse every #include'd .itp (from the copied-includes map) and build
        {moleculeName: info}. Also picks up moleculetypes defined inline in the
        .top itself."""
        molInfo = {}
        for itpPath in includes.values():
            if os.path.exists(itpPath):
                molInfo.update(parseItpMolecules(itpPath))
        for name, info in parseItpMolecules(topFile).items():
            molInfo.setdefault(name, info)
        return molInfo

    def _buildChains(self, proteinMols, molInfo):
        """Expand the protein molecules (name, count) into individual chains,
        naming them A, B, C ..."""
        chains, letters = [], self._chainLetters()
        for name, count in proteinMols:
            info = molInfo.get(name, {})
            segments = info.get('chains') or [info.get('nres', 0)]
            for _ in range(count):
                for segLength in segments:
                    cid = next(letters)
                    chains.append({'id': cid, 'length': segLength, 'molName': name})
        return chains

    @staticmethod
    def _chainLetters():
        # A..Z then a..z, then AA, AB ... (very unlikely to be needed)
        import string, itertools
        singles = string.ascii_uppercase + string.ascii_lowercase
        yield from singles
        for a, b in itertools.product(singles, repeat=2):
            yield a + b

    def _applyPdbChainNames(self, chains):
        """If the user supplied an original PDB with explicit chain identifiers
        and their number matches the chains derived from the topology, rename the
        chains using those identifiers."""
        pdb = self.inputPdb.get()
        if not pdb or not os.path.exists(pdb):
            return chains
        pdbChains = parsePdbChains(pdb)
        if len(pdbChains) == len(chains) and len(pdbChains) > 0:
            self.info(f'Using chain identifiers from PDB: {[c for c, _ in pdbChains]}')
            for chain, (cid, _) in zip(chains, pdbChains):
                chain['id'] = cid
        elif pdbChains:
            self.info(f'PDB has {len(pdbChains)} chain(s) but topology yields {len(chains)}; '
                      f'keeping topology-derived chain names.')
        return chains

    def _getForcedLigands(self):
        """Set of user-specified ligand residue/molecule names (upper-cased)."""
        raw = self.ligandNames.get()
        if not raw:
            return set()
        return {n.strip().upper() for n in raw.split(',') if n.strip()}

    def _orderLigands(self, ligands, forcedLigands):
        """Order the ligand molecules so that user-named ones come first (their
        order in the parameter is preserved), keeping the rest afterwards. This
        makes the first user-named ligand the primary one."""
        if not forcedLigands:
            return ligands
        forced = [m for m in ligands if m[0].upper() in forcedLigands]
        rest = [m for m in ligands if m[0].upper() not in forcedLigands]
        return forced + rest

    def _collectResnames(self, mols, molInfo):
        """Union of the residue names of a list of (name, count) molecules,
        falling back to the molecule name when the .itp is not available."""
        resnames = set()
        for name, _ in mols:
            info = molInfo.get(name)
            if info and info['resnames']:
                resnames |= info['resnames']
            else:
                resnames.add(name)
        return sorted(resnames)

    def _pickLigandID(self, ligand, molInfo):
        """A GROMACS residue name is at most 3-4 chars; use the ligand residue
        name when it looks standard, otherwise fall back to the molecule name."""
        name = ligand[0]
        info = molInfo.get(name)
        if info and info['resnames']:
            # single-residue ligand -> use that residue name (e.g. STI)
            resnames = sorted(info['resnames'])
            if len(resnames) == 1:
                return resnames[0]
        return name

    def _getMoleculeItp(self, mol, includes):
        """Return the copied .itp file that defines the given molecule."""
        name = mol[0]
        for relInc, absPath in includes.items():
            defined = parseItpMolecules(absPath)
            if name in defined:
                return absPath
        return None

    def _logClassification(self, classified, chains, ligands, lipidResnames):
        self.info('Imported system composition:')
        self.info(f"  Protein chains : {len(chains)} "
                  f"({', '.join(c['id'] + ':' + str(c['length']) + 'res' for c in chains) or 'none'})")
        self.info(f"  Ligand(s)      : {', '.join(n for n, _ in ligands) or 'none'}")
        self.info(f"  Membrane       : {'yes -> ' + ', '.join(lipidResnames) if lipidResnames else 'no'}")
        self.info(f"  Ions           : {', '.join(n for n, _ in classified['ion']) or 'none'}")
        self.info(f"  Water          : {', '.join(n for n, _ in classified['water']) or 'none'}")

    # --------------------------- Index generation ----------------------------
    def _generateIndexFile(self, groSystem, chains, lipidResnames, ligResnames, inputIndex=None):
        """Build the .ndx file for the imported system.

        Strategy (mirrors Plugin.firstIndexCreation but for prepared systems):
          1. make_ndx with the default 'q' -> standard groups (System, Protein, Water, Ion ...).
          2. Optionally merge the groups of a user-provided .ndx into the base.
          3. In a single make_ndx pass on top of the base file, create the
             derived groups (per-chain, Membrane, Protein_<lig>).
        """
        baseIndex = self._getExtraPath('indexes_base.ndx')
        gromacsPlugin.createIndexFile(self, groSystem, outIndex=baseIndex)

        if inputIndex and os.path.exists(inputIndex):
            mergeIndexGroups(baseIndex, inputIndex)

        builder = IndexGroupBuilder(gromacsPlugin.parseIndexFile(self, baseIndex))   # {num: name}
        self._addChainGroups(builder, chains)
        protNum = self._addProteinGroup(builder, chains)
        membNum = builder.add('r ' + ' '.join(lipidResnames), 'Membrane') if lipidResnames else None
        # Ligand groups (one per ligand residue name).
        ligNums = [builder.add(f'r {resname}', resname) for resname in ligResnames]
        if protNum is not None:
            self._addComplexGroups(builder, protNum, membNum, ligResnames, ligNums)

        outIndex = self._getExtraPath('indexes.ndx')
        if builder.cmds:
            gromacsPlugin.createIndexFile(self, groSystem, inIndex=baseIndex,
                                          outIndex=outIndex, inputCommands=builder.cmds)
        else:
            shutil.copy(baseIndex, outIndex)
        return outIndex

    @staticmethod
    def _addChainGroups(builder, chains):
        """Per-chain groups -- only meaningful with more than one chain."""
        if len(chains) <= 1 or not all(c['length'] for c in chains):
            return
        start = 1
        for chain in chains:
            end = start + chain['length'] - 1
            builder.add(f"ri {start}-{end}", f"chain{chain['id']}")
            start = end + 1

    @staticmethod
    def _addProteinGroup(builder, chains):
        """Number of the Protein group, created from the chain lengths if make_ndx did not define it"""
        protNum = builder.name2num.get('Protein')
        if protNum is None and chains and all(c['length'] for c in chains):
            totalRes = sum(c['length'] for c in chains)
            protNum = builder.add(f"ri 1-{totalRes}", 'Protein')
        return protNum

    @staticmethod
    def _addComplexGroups(builder, protNum, membNum, ligResnames, ligNums):
        """Combined groups expected by the downstream protocols/viewer."""
        # Protein_<ligandID> for every ligand present.
        for resname, num in zip(ligResnames, ligNums):
            builder.add(f'{protNum} | {num}', f'Protein_{resname}')

        if membNum is not None:
            builder.add(f'{protNum} | {membNum}', 'Protein_membrane')
            # Protein_membrane_lig only when BOTH membrane and ligand exist.
            if ligNums:
                union = ' | '.join(str(n) for n in [protNum, membNum] + ligNums)
                builder.add(union, 'Protein_membrane_lig')

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        vals = []
        coords = self.inputCoords.get()
        if coords and os.path.splitext(coords)[1].lower() not in ['.gro', '.pdb']:
            vals.append('Input coordinates must be a .gro or .pdb file.')
        top = self.inputTopology.get()
        if top and not top.lower().endswith('.top'):
            vals.append('Input topology must be a .top file.')
        if top and os.path.exists(top):
            molecules = parseTopologyMolecules(top)
            if not molecules:
                vals.append('Could not find a "[ molecules ]" section in the topology file. '
                            'Is it a valid GROMACS .top file?')
            else:
                forced = self._getForcedLigands()
                present = {n.upper() for n, _ in molecules}
                unknown = [n for n in forced if n not in present]
                if unknown:
                    vals.append('Ligand name(s) not found in the topology "[ molecules ]" '
                                'section: {}. Available molecules: {}.'.format(
                                    ', '.join(unknown), ', '.join(n for n, _ in molecules)))
        return vals

    def _warnings(self):
        warns = []
        top = self.inputTopology.get()
        if top and os.path.exists(top):
            missing = []
            topDir = os.path.dirname(os.path.abspath(top))
            for relInc in parseTopologyIncludes(top):
                # Use the same resolution logic as the actual import (topology
                # dir -> optional .itp dir -> GROMACS's own force-field dir),
                # so standard force fields (amber*.ff, charmm*.ff, oplsaa.ff...)
                # shipped with GROMACS itself aren't reported as missing.
                if self._resolveInclude(relInc, topDir) is None:
                    missing.append(relInc)
            if missing:
                warns.append('The following files #included by the topology were not found and '
                             'will be skipped (the imported topology may be incomplete):\n  ' +
                             '\n  '.join(missing))
        return warns

    def _summary(self):
        summary = []
        if self.isFinished() and hasattr(self, 'outputSystem'):
            system = self.outputSystem
            summary.append('Imported GROMACS system: %s' % os.path.basename(system.getSystemFile()))
            try:
                summary.append('Chains: %s' % ', '.join(system.getChainNames()))
            except Exception:
                pass
            if system.hasLig():
                summary.append('Ligand: %s' % system.getLigandID())
            if system.getIndexFile():
                summary.append('Index groups generated: %s' %
                               ', '.join(gromacsPlugin.parseIndexFile(self, system.getIndexFile()).values()))
        else:
            summary.append('The protocol has not finished.')
        return summary

    def _methods(self):
        return ['Imports an already-prepared GROMACS system (topology, coordinates and optional '
                'trajectory), classifies its molecules from the topology, derives the chain '
                'names/lengths and generates a complete index (.ndx) file (standard GROMACS groups '
                'plus per-chain, ligand and membrane groups) so it can be used by the simulation, '
                'MMPBSA, trajectory and viewer protocols.']


class IndexGroupBuilder:
    """Queue of make_ndx selections + renames on top of an existing index, tracking the number make_ndx
    will assign to each new group. A group whose name already exists (case-insensitive) is reused
    instead of duplicated."""
    def __init__(self, groups):
        self.name2num = {v: k for k, v in groups.items()}
        self.ciName2num = {v.lower(): k for k, v in groups.items()}
        self.nextNum = max(groups.keys()) + 1
        self.cmds = []

    def add(self, selectionCmd, groupName):
        if groupName.lower() in self.ciName2num:
            return self.ciName2num[groupName.lower()]
        assigned = self.nextNum
        self.cmds += [selectionCmd, f'name {assigned} {groupName}']
        self.name2num[groupName] = self.ciName2num[groupName.lower()] = assigned
        self.nextNum += 1
        return assigned


_SECTION_RE = re.compile(r'\[([^\]]+)\]')
_INCLUDE_RE = re.compile(r'#include\s+"([^"]+)"')

def stripComment(line):
    return line.split(';', 1)[0]

def parseTopologyMolecules(topFile):
    """Return the ordered list of (moleculeName, count) from the topology's
    '[ molecules ]' section."""
    molecules, inMolecules = [], False
    with open(topFile) as f:
        for raw in f:
            line = stripComment(raw).strip()
            if not line:
                continue
            section = _SECTION_RE.match(line)
            if section:
                inMolecules = section.group(1).strip().lower() == 'molecules'
                continue
            if inMolecules:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        molecules.append((parts[0], int(parts[1])))
                    except ValueError:
                        molecules.append((parts[0], 1))
    return molecules

def parseTopologyIncludes(topFile):
    """Return the list of relative include paths ('#include "..."') found in a
    .top/.itp file (in file order, de-duplicated)."""
    includes = []
    with open(topFile) as f:
        for raw in f:
            line = raw.strip()
            if line.startswith('#include'):
                m = _INCLUDE_RE.search(line)
                if m and m.group(1) not in includes:
                    includes.append(m.group(1))
    return includes

def rewriteTopologyIncludes(topFile, includes):
    """Rewrite the '#include "rel"' directives of ``topFile`` with its absolute location. """
    with open(topFile) as f:
        lines = f.readlines()

    changed = False
    for i, line in enumerate(lines):
        if line.lstrip().startswith('#include'):
            m = _INCLUDE_RE.search(line)
            if m and m.group(1) in includes:
                absPath = includes[m.group(1)]
                lines[i] = line[:m.start(1)] + absPath + line[m.end(1):]
                changed = True

    if changed:
        with open(topFile, 'w') as f:
            f.writelines(lines)
    return topFile

def parseItpMolecules(itpFile):
    """Parse an .itp file and return {moleculeName: info} for every
    '[ moleculetype ]' defined in it. ``info`` holds: resnames (set), resseq
    (ordered residue names), nres (total residues), natoms, and chains (residue
    count of each chain segment."""
    mols, current, awaitingName = {}, None, False
    for section, line in iterTopologyLines(itpFile):
        if line is None:
            # Section header: the next line of a [ moleculetype ] is the molecule name
            awaitingName = section == 'moleculetype'
        elif section == 'moleculetype' and awaitingName:
            current, awaitingName = line.split()[0], False
            mols[current] = {'resnames': set(), 'resseq': [], 'nres': 0, 'natoms': 0, 'chains': [],
                             '_lastResnr': None, '_lastResnrInt': None, '_curChain': 0}
        elif section == 'atoms' and current in mols and len(line.split()) >= 4:
            # nr type resnr residue atom cgnr [charge] [mass]
            parts = line.split()
            addItpAtom(mols[current], parts[2], parts[3])

    # Finalise the last (or only) chain segment of every moleculetype.
    for info in mols.values():
        if info['_curChain']:
            info['chains'].append(info['_curChain'])
    return mols

def iterTopologyLines(topFile):
    """Yield (section, line) for the data lines of a .top/.itp file with comments stripped, and
    (section, None) at every section header. Yields nothing if the file cannot be read."""
    try:
        with open(topFile) as f:
            lines = f.readlines()
    except OSError:
        return
    section = None
    for raw in lines:
        line = stripComment(raw).strip()
        sectionMatch = _SECTION_RE.match(line)
        if sectionMatch:
            section = sectionMatch.group(1).strip().lower()
            yield section, None
        elif line:
            yield section, line

def addItpAtom(info, resnr, resname):
    """Account an '[ atoms ]' line of a moleculetype in its info dictionary"""
    info['natoms'] += 1
    info['resnames'].add(resname)
    if resnr == info['_lastResnr']:
        return
    resnrInt, lastResnrInt = toInt(resnr), info['_lastResnrInt']
    # A decrease in residue number => a new chain segment
    # (CHARMM-GUI can merge several chains into one moleculetype while keeping per-chain numbering).
    if lastResnrInt is not None and resnrInt is not None and resnrInt < lastResnrInt:
        info['chains'].append(info['_curChain'])
        info['_curChain'] = 0
    info['_curChain'] += 1
    info['nres'] += 1
    info['resseq'].append(resname)
    info['_lastResnr'], info['_lastResnrInt'] = resnr, resnrInt

def toInt(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None

def parsePdbChains(pdbFile):
    """Return the ordered list of (chainId, residueCount) of the protein chains
    of a PDB, reading the chain identifier column (22) of ATOM records. Blank
    chain identifiers (as CHARMM-GUI writes) yield an empty list, which makes the
    caller fall back to the topology-derived chains."""
    chains, order, lastResKey = {}, [], None
    with open(pdbFile) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            chainId = line[21].strip()
            if not chainId:
                continue
            resKey = (chainId, line[22:27])  # resSeq + insertion code
            if chainId not in chains:
                chains[chainId] = 0
                order.append(chainId)
            if resKey != lastResKey:
                chains[chainId] += 1
                lastResKey = resKey
    return [(cid, chains[cid]) for cid in order]

def isProteinResidue(resname):
    if resname in STANDARD_AA:
        return True
    # CHARMM-GUI terminal patches: N/C prefix (NGLY, CVAL ...) or P suffix (GLYP)
    if len(resname) >= 4 and resname[0] in ('N', 'C') and resname[1:] in STANDARD_AA:
        return True
    if len(resname) >= 4 and resname[-1] in ('P', 'T') and resname[:-1] in STANDARD_AA:
        return True
    return False

def classifyMolecule(name, info, forcedLigands=None):
    """Classify a single molecule into 'protein', 'lipid', 'ligand', 'ion' or
    'water' using its name and (when available) its residue composition.

    ``forcedLigands`` is an optional set of upper-cased names that the user
    explicitly declared as ligands; a match short-circuits the heuristic."""
    up = name.upper().rstrip('+-')
    resnames = {r.upper() for r in info['resnames']} if info else set()

    # User-declared ligands override every automatic rule.
    if forcedLigands and ({name.upper(), up} | resnames) & forcedLigands:
        return 'ligand'

    # Water first (most abundant, unambiguous names).
    if up in WATER_NAMES or (resnames and resnames <= WATER_NAMES):
        return 'water'

    if isProteinMolecule(info):
        return 'protein'

    # Lipids / membrane.
    if up in LIPID_NAMES or (resnames and resnames <= LIPID_NAMES):
        return 'lipid'

    # Monoatomic ions (by name or by being a single-atom molecule).
    if up in ION_NAMES or (info and info['natoms'] == 1):
        return 'ion'

    # Anything else that is not water/ion/lipid/protein is treated as a ligand.
    return 'ligand'

def isProteinMolecule(info):
    """Protein: residues are (mostly) amino acids and there is more than one"""
    if not info or info['nres'] <= 1:
        return False
    aa = sum(1 for r in info['resseq'] if isProteinResidue(r.upper()))
    return aa / info['nres'] > 0.5

def classifyMolecules(molecules, molInfo, forcedLigands=None):
    """Classify every (name, count) molecule from the topology. Returns a dict
    with keys protein/ligand/lipid/ion/water, each a list of (name, count).
    Names in ``forcedLigands`` (upper-cased set) are always classified as ligands."""
    result = {'protein': [], 'ligand': [], 'lipid': [], 'ion': [], 'water': []}
    for name, count in molecules:
        info = molInfo.get(name)
        result[classifyMolecule(name, info, forcedLigands)].append((name, count))
    return result

def detectForceFields(topFile, waterMols):
    """Best-effort detection of the (main, water) force field names of an
    imported system, used only as metadata (import does not re-run pdb2gmx)."""
    mainFF, waterFF = None, None

    header = ''
    try:
        with open(topFile) as f:
            header = ''.join(f.readline() for _ in range(40)).upper()
    except OSError:
        pass
    if 'CHARMM' in header:
        mainFF = 'charmm36-feb2021'
    elif 'AMBER' in header:
        mainFF = 'amber99sb-ildn'
    elif 'OPLS' in header:
        mainFF = 'oplsaa'

    waterMap = {'TIP3': 'tip3p', 'TIP3P': 'tip3p', 'TIP4': 'tip4p', 'TIP4P': 'tip4p',
                'TIP5': 'tip5p', 'TIP5P': 'tip5p', 'SPC': 'spc', 'SPCE': 'spce',
                'SOL': 'spc', 'OPC': 'tip4p', 'OPC3': 'tip3p'}
    for name, _ in waterMols:
        if name.upper() in waterMap:
            waterFF = waterMap[name.upper()]
            break
    return mainFF, waterFF

def mergeIndexGroups(baseNdx, extraNdx):
    """Append the groups of ``extraNdx`` that are not already present in
    ``baseNdx`` to ``baseNdx``"""
    def readGroups(path):
        groups, order, current = {}, [], None
        with open(path) as f:
            for line in f:
                m = _SECTION_RE.match(line.strip())
                if m:
                    current = m.group(1).strip()
                    groups.setdefault(current, [])
                    order.append(current)
                elif current is not None:
                    groups[current].append(line.rstrip('\n'))
        return groups, order

    baseGroups, _ = readGroups(baseNdx)
    extraGroups, extraOrder = readGroups(extraNdx)

    baseLower = {g.lower() for g in baseGroups}
    newGroups = [g for g in extraOrder if g.lower() not in baseLower]
    if not newGroups:
        return baseNdx

    with open(baseNdx, 'a') as f:
        for g in newGroups:
            f.write(f'\n[ {g} ]\n')
            for line in extraGroups[g]:
                f.write(line + '\n')
    return baseNdx
