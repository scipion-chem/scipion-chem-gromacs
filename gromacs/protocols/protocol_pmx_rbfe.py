# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Joaquin Algorta (joaquin.algorta@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
Relative binding free energy (RBFE) between two or more ligands docked/aligned
in the same pocket of a shared receptor, using pmx's non-equilibrium
fast-growth alchemical FEP, following the standard double-leg thermodynamic
cycle (confirmed against pmx's own ligand_tutorial.ipynb and Gapsys/
Perez-Benito et al. 2020, Chem. Sci. 11:1140:

  0. Build ligand A's solvated/ionized COMPLEX system (reusing
     GromacsSystemPrep's own ligand-mode steps), ligand A's solvated/ionized
     WATER-only system (ligand alone, no protein), and parametrize ligand B
     (ACPYPE).
  1. pmx atomMapping: find the morphable atoms between ligand A and B
     (ligand-pair-intrinsic - computed once, shared by both legs below).
  2. pmx ligandHybrid: build the dual-topology hybrid ligand (.itp + dummy
     atomtypes) - also shared by both legs.
  3. Splice the hybrid ligand into *both* ligand A systems' topology/
     coordinates (box, solvent and ions of each are left otherwise untouched):
     the "complex" leg (protein+ligand+solvent) and the "water" leg
     (ligand+solvent only).
  4. For *each* leg: equilibrate (EM->NVT->NPT) the hybrid system at both end
     states (lambda=0 == ligand A, lambda=1 == ligand B).
  5. For *each* leg: extract snapshots from each equilibrium trajectory and
     launch short forward (0->1) / reverse (1->0) switching ("fast growth") runs.
  6. For *each* leg: pmx analyse estimates dG_leg(A->B) from the collected
     dH/dl work values (Crooks Gaussian Intersection / Bennett Acceptance
     Ratio / Jarzynski).
  7. Combine via the thermodynamic cycle: ddG_bind(A->B) = dG_complex - dG_water.
     Both leg values are always reported (never just the combined number), so
     a wrong leg is auditable rather than silently baked into one opaque total.

When more than 2 ligands are selected, steps 0-7 above run once per edge of a
similarity-ordered A->B->C->... chain (see _computeEdgeChain and
claude/decisions/pmx_RBFE.md §15) - each edge sequentially (no two edges'
steps ever run concurrently, since several of the steps above write to fixed,
non-edge-scoped paths that the next edge's identical steps would otherwise
overwrite), though the steps *within* a single edge keep their existing
intra-edge parallelism.

Reference: pmx, Gapsys et al. (https://github.com/deGrootLab/pmx)
"""

import os
import re
import shutil
from glob import glob
from multiprocessing import cpu_count

from pyworkflow.protocol import params
from pyworkflow.utils import Message

from pwchem import Plugin as pwchemPlugin

from gromacs import Plugin as gromacsPlugin
from gromacs.objects import GromacsSystem
from gromacs.constants import PMX_SC_ALPHA, PMX_SC_SIGMA, PMX_SW_TIME, PMX_N_SNAPSHOTS
from gromacs.protocols.protocol_system_prep import (
    GromacsSystemPrep, replaceInFile, LIGAND, TOPOL_TOP, GAPS_OPTIONS,
)

try:
    from rdkit import Chem
    from rdkit import DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

MAP_DIR = 'mapping'
HYBRID_DIR = 'hybrid'
SYSTEM_DIR = 'hybrid_system'
FWD, BWD = 'fwd', 'bwd'
STATE_DIRNAME = {0: 'stateA', 1: 'stateB'}
DIRECTION_STATE = {FWD: 0, BWD: 1}

# The two legs of the standard RBFE thermodynamic cycle: ligand A->B morphed while
# bound to the protein ("complex") and morphed alone in solvent ("water"). Atom
# mapping and hybrid-ligand construction (MAP_DIR/HYBRID_DIR above) are ligand-pair-
# intrinsic and shared by both; only the surrounding system (and everything built on
# top of it - equilibration, snapshots, transitions, analysis) differs per leg.
COMPLEX, WATER = 'complex', 'water'
LEGS = (COMPLEX, WATER)
WATER_SYSTEM_DIR = 'water_system'
LEG_SYSTEM_DIR = {COMPLEX: SYSTEM_DIR, WATER: 'water_hybrid_system'}

# Whether to run every ligand in the input set, or only an explicitly-picked subset -
# see the ligandSelection param.
LIGANDS_ALL, LIGANDS_SUBSET = 0, 1


class GromacsPmxRBFE(GromacsSystemPrep):
    """
    Relative binding free energy via pmx non-equilibrium fast-growth alchemical
    FEP, for two or more ligands.

    Ligands are picked (via SelectMultiLigandsWizard) from the same
    SetOfSmallMolecules, docked/aligned in the same pocket of a shared
    receptor. With exactly 2 selected, a single A->B edge runs (this
    protocol's original behaviour, output 'outputSystem'). With 3 or more, a
    similarity-ordered A->B->C->... chain of edges runs sequentially, one
    'outputSystem_edge<i>' output per edge (see _computeEdgeChain).

    Ligand A's solvated/ionized system is built by this protocol itself,
    reusing GromacsSystemPrep's own ligand-mode preparation steps
    (PDB2GMX/editconf/solvate/genion) - subclassing it is exactly this reuse,
    the same "compose, don't duplicate" pattern GromacsMmpbsa already uses.
    """
    _label = 'pmx relative binding free energy (RBFE)'
    stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        cpus = cpu_count() // 2
        form.addParallelSection(threads=cpus, mpi=1)

        form.addHidden(params.USE_GPU, params.BooleanParam, default=False,
                       label='Use GPU for execution',
                       help='GPU may have several cores. Set it one if you don\'t know what we are '
                            'talking about but you have a GPU.')
        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       expertLevel=params.LEVEL_ADVANCED,
                       label='Choose GPU IDs',
                       help='Add a list of GPU devices that can be used (Comma separated)')
        # RBFE is always ligand-vs-ligand; the AtomStruct-only input mode GromacsSystemPrep
        # also supports doesn't apply here, so this stays fixed and hidden.
        form.addHidden('inputFrom', params.EnumParam, choices=['AtomStruct', 'SetOfSmallMolecules'],
                       default=LIGAND)

        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSetOfMols', params.PointerParam, pointerClass='SetOfSmallMolecules',
                      label='Docked molecules: ', allowsNull=False,
                      help='Set of docked/aligned molecules sharing a common receptor. Ligands are '
                           'picked from this same set below.')
        form.addParam('ligandSelection', params.EnumParam, choices=['All', 'Select subset'],
                      default=LIGANDS_SUBSET, display=params.EnumParam.DISPLAY_HLIST,
                      label='Ligands to use: ',
                      help='"All": every molecule in the set above is used (see the warning below '
                           'about the cost of this). "Select subset": pick specific ligands below.')
        form.addParam('selectedLigands', params.TextParam, width=70, allowsNull=False,
                      condition=f'ligandSelection=={LIGANDS_SUBSET}',
                      label='Ligands: ',
                      help='2 or more ligands, picked from the set above via the wizard '
                           '(Ctrl+Click or Shift+Click for multiple), one per line. With '
                           'exactly 2, a single A->B edge runs. With 3 or more, this protocol '
                           'orders them into an A->B->C->... chain by pairwise structural '
                           'similarity (RDKit Morgan/Tanimoto - falls back to the order '
                           'selected above if RDKit or a structure is unusable) and runs one '
                           'edge per consecutive pair, sequentially. This is a simplified '
                           'stand-in for real published multi-ligand RBFE networks, which use '
                           'a hub-and-spoke/bridge topology with redundant cycle-closure edges '
                           '(built with LOMAP or curated by hand) rather than a plain chain - '
                           'see claude/decisions/pmx_RBFE.md §15.')
        # Derived from selectedLigands by the wizard above (first two names) - kept as real,
        # named params (not just parsed out of selectedLigands.get() everywhere) so every
        # existing method that already reads self.inputLigand/self.ligandB (inherited from
        # GromacsSystemPrep, or this protocol's own) keeps working unchanged for the
        # single-edge (2-ligand) case; for 3+, _insertMultiEdgeSteps repoints them per edge.
        form.addHidden('inputLigand', params.StringParam, label='Ligand A: ')
        form.addHidden('ligandB', params.StringParam, label='Ligand B: ')

        group = form.addGroup('Force field')
        self._defineFFParams(group)
        self._defineACPYPEparams(form, condition=True)

        form.addSection('System preparation')
        group = form.addGroup('Boundary box')
        self._defineBoxParams(group)

        group = form.addGroup('Ions')
        self._defineIonsParams(group)

        form.addParam('addCaps', params.EnumParam, choices=GAPS_OPTIONS, default=0,
                      label='Add ACE and NME caps: ',
                      help='Add acetyl (ACE) and N-methylamide (NME) capping groups to protein '
                           'termini before building ligand A\'s syste.')

        group = form.addGroup('SS bonds')
        self._defineSSBondsParams(group)

        form.addSection(label='RBFE settings')
        grp = form.addGroup('End-state equilibration')
        grp.addParam('nStepsMin', params.IntParam, default=10000,
                     label='Max EM steps: ',
                     help='Maximum number of (minimization) steps to perform')
        grp.addParam('emTol', params.FloatParam, default=1000.0,
                     label='EM max force objective: ',
                     help='Stop minimization when the maximum force < x kJ/mol/nm.')
        grp.addParam('nvtTime', params.FloatParam, default=200.0, label='NVT time (ps): ',
                     help='Default (200 ps) is enough for thermalization on typical systems; raise '
                          'it if temperature/pressure have not visibly plateaued in the log.')
        grp.addParam('nptTime', params.FloatParam, default=2000.0, label='NPT time (ps): ',
                     help='Snapshots for the switching runs are taken evenly from the second half of '
                          'this trajectory (the first half is burn-in). Default (2000 ps, i.e. 1000 ps '
                          'of usable trajectory) is set so that, at the default 20 snapshots per '
                          'end-state, consecutive snapshots are ~50 ps apart - enough to decorrelate '
                          'for most drug-like ligands. If you raise "Snapshots per end-state" '
                          'substantially, raise this proportionally too, or the extra snapshots will '
                          'just be correlated repeats rather than independent samples.')
        grp.addParam('temperature', params.FloatParam, default=300.0, label='Temperature (K): ')
        grp.addParam('pressure', params.FloatParam, default=1.0, label='Pressure (bar): ')
        grp.addParam('timeStep', params.FloatParam, default=0.002, expertLevel=params.LEVEL_ADVANCED,
                     label='Time step (ps)[dt]: ',
                     help='MD integration time step (mdp "dt"), used for every step EM/NVT/NPT/'
                          'switching-run.')

        grp = form.addGroup('Fast-growth transitions')
        grp.addParam('nStructs', params.IntParam, default=PMX_N_SNAPSHOTS,
                     label='Snapshots per end-state: ',
                     help='Number of frames taken from each end-state equilibration to launch '
                          'independent forward/reverse switching runs from - i.e. this many forward '
                          '+ this many reverse non-equilibrium transitions get run and fed to the '
                          'CGI/BAR/JARZ estimators. Default (20, i.e. 40 total transitions) is pmx\'s '
                          'own commonly-used baseline and should give a usable estimate with a '
                          'reasonably tight error bar for a typical drug-like ligand pair; the error '
                          'bar narrows (roughly) with the square root of this number, so doubling it '
                          'to 40 will noticeably tighten it at ~2x the transition-running cost.')
        grp.addParam('swTime', params.FloatParam, default=PMX_SW_TIME,
                     label='Switching time (ps): ',
                     help='Length of the linear lambda-ramp transition run (pmx default and this '
                          'protocol\'s default: 50 ps.')
        grp.addParam('scAlpha', params.FloatParam, default=PMX_SC_ALPHA, expertLevel=params.LEVEL_ADVANCED,
                     label='Soft-core alpha: ',
                     help='Soft-core alpha (mdp "sc-alpha"): softens the vdW/Coulomb potential for '
                          'appearing/disappearing (dummy) atoms so it stays finite as they turn on/off, '
                          'avoiding the singularities a normal Lennard-Jones potential would hit.')
        grp.addParam('scSigma', params.FloatParam, default=PMX_SC_SIGMA, expertLevel=params.LEVEL_ADVANCED,
                     label='Soft-core sigma: ',
                     help='Soft-core sigma: minimum effective interaction radius used '
                          'by the soft-core potential above for appearing/disappearing atoms.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        edges = self._computeEdgeChain()
        if len(edges) <= 1:
            # Exactly 2 ligands: single-edge behaviour, output keeps its original name
            # ('outputSystem') rather than an edge suffix. Explicitly (re)set inputLigand/
            # ligandB from the resolved pair rather than trusting they're already correct:
            # true when picked via the subset wizard (which sets them as a side effect), but
            # NOT when "All" resolves to exactly 2 ligands - nothing populates them in that
            # mode otherwise.
            if edges:
                self.inputLigand.set(edges[0][0])
                self.ligandB.set(edges[0][1])
                self._store(self.inputLigand, self.ligandB)
            legAnalyseSteps = self._insertEdgeStepGraph()
            self._insertFunctionStep(self.createOutputStep, prerequisites=legAnalyseSteps)
            return

        prevStep = None
        outputSteps = []
        for edgeIdx, (ligA, ligB) in enumerate(edges):
            # Edges run strictly sequentially: several steps below (PDB2GMXStep,
            # editConfStep, solvateStep, addIonsStep, and every leg-scoped path under
            # extra/) write to fixed, non-edge-scoped paths - see archiveEdgeStep. Running
            # two edges concurrently would have them clobber each other's files (and race
            # on the shared self.inputLigand/self.ligandB set below).
            setStep = self._insertFunctionStep(self.setEdgeLigandsStep, edgeIdx, ligA, ligB,
                                               prerequisites=[prevStep] if prevStep else [])
            legAnalyseSteps = self._insertEdgeStepGraph(prerequisites=[setStep])
            archiveStep = self._insertFunctionStep(self.archiveEdgeStep, edgeIdx, ligA, ligB,
                                                   prerequisites=legAnalyseSteps)
            outputSteps.append(self._insertFunctionStep(self.createEdgeOutputStep, edgeIdx, ligA, ligB,
                                                         prerequisites=[archiveStep]))
            prevStep = archiveStep

        self._insertFunctionStep(self.createChainSummaryStep, edges, prerequisites=outputSteps)

    def _insertEdgeStepGraph(self, prerequisites=None):
        """Insert one edge's full ligand-A-vs-ligand-B step graph (system builds, atom
        mapping/hybrid ligand construction, both legs' equilibration/snapshot/transition/
        analysis steps) and return the two analyseStep ids. `prerequisites` gates only the
        first two steps (ligand A/B parametrization); everything downstream is chained off
        those exactly as in this protocol's original single-edge _insertAllSteps.

        Must default to None, NOT an empty tuple/list: pyworkflow's own __insertStep only
        special-cases `prerequisites is None` (falling back to "no prerequisites" for the very
        first step, or implicitly chaining off the previously-inserted step otherwise - the
        single-edge path below relies on exactly that implicit chaining, matching this
        protocol's original, pre-multi-edge behaviour). Passing `()` instead - confirmed by
        actually hitting this - makes pyworkflow's `isinstance(prerequisites, list)` check fail
        (a tuple isn't a list) and wrap the whole tuple as one bogus "prerequisite id", crashing
        with "TypeError: unsupported operand type(s) for -: 'tuple' and 'int'" the moment the
        step executor tries to resolve it."""
        # Ligand A's COMPLEX system build, reusing GromacsSystemPrep's own ligand-mode
        # steps unchanged (self.inputLigand/inputSetOfMols/inputFrom==LIGAND drive them).
        ligAStep = self._insertFunctionStep(self.parametrizeLigandStep, prerequisites=prerequisites)
        ligBStep = self._insertFunctionStep(self.parametrizeLigandBStep,
                                            prerequisites=prerequisites)  # independent of ligand A

        pdbStep = self._insertFunctionStep(self.PDB2GMXStep, prerequisites=[ligAStep])
        ecStep = self._insertFunctionStep(self.editConfStep, prerequisites=[pdbStep])
        solvStep = self._insertFunctionStep(self.solvateStep, prerequisites=[ecStep])
        prevStep = solvStep
        if self.placeIons.get() != 0:
            prevStep = self._insertFunctionStep(self.addIonsStep, prerequisites=[solvStep])

        # Ligand A's WATER-only system - independent of the protein/complex build above,
        # only needs ligand A's own standalone ACPYPE parametrization.
        wbStep = self._insertFunctionStep(self.buildLigandAWaterSystemStep, prerequisites=[ligAStep])

        # Atom mapping / hybrid ligand construction: ligand-pair-intrinsic, done ONCE,
        # shared by both legs below.
        mStep = self._insertFunctionStep(self.atomMappingStep, prerequisites=[prevStep, ligBStep])
        hStep = self._insertFunctionStep(self.buildHybridLigandStep, prerequisites=[mStep])

        csStep = self._insertFunctionStep(self.buildHybridSystemStep, prerequisites=[hStep])
        wsStep = self._insertFunctionStep(self.buildHybridWaterSystemStep, prerequisites=[hStep, wbStep])
        legSysStep = {COMPLEX: csStep, WATER: wsStep}

        legAnalyseSteps = []
        for leg in LEGS:
            eqSteps, extractSteps = {}, {}
            for state in (0, 1):
                eqSteps[state] = self._insertFunctionStep(self.equilibrateStateStep, leg, state,
                                                          prerequisites=[legSysStep[leg]])
                extractSteps[state] = self._insertFunctionStep(self.extractSnapshotsStep, leg, state,
                                                                prerequisites=[eqSteps[state]])

            transitionSteps = []
            for direction, state in DIRECTION_STATE.items():
                for i in range(self.nStructs.get()):
                    tStep = self._insertFunctionStep(self.runTransitionStep, leg, direction, i,
                                                     prerequisites=[extractSteps[state]])
                    transitionSteps.append(tStep)

            legAnalyseSteps.append(self._insertFunctionStep(self.analyseStep, leg,
                                                             prerequisites=transitionSteps))
        return legAnalyseSteps

    # -- Multi-ligand chain construction -----------------------------------------
    def _computeEdgeChain(self):
        """Order the selected ligand pool into an A->B->C->... chain by greedy
        nearest-neighbor structural similarity (RDKit Morgan fingerprints, Tanimoto):
        start from the first ligand in selection order, then repeatedly append whichever
        not-yet-used ligand is most similar to the last one added. Returns a list of
        (ligA, ligB) name pairs, one per edge (empty if fewer than 2 are selected).

        This is a deliberately simple stand-in, not a from-scratch reimplementation of
        how real published multi-ligand RBFE campaigns build their network: the RBFE
        methodology paper this protocol otherwise follows (Gapsys et al. 2020, Chem. Sci.
        11:1140) doesn't build one either - its "Selected datasets" section explicitly
        reuses pre-existing benchmark edges (from Wang et al. 2015 and later FEP studies)
        rather than deriving them itself. Those real networks (confirmed against the
        actual published JNK1 edge list, see claude/decisions/pmx_RBFE.md §15) are
        hub-and-spoke-plus-bridge topologies with redundant cycle-closure edges,
        typically built with LOMAP or curated by hand - out of scope here. A
        similarity-ordered chain at least keeps each individual edge's two ligands
        structurally close (the property that most directly affects a single edge's own
        convergence), though it gives up the cycle-closure consistency check a real
        network provides."""
        names = self._getSelectedLigandNames()
        if len(names) < 2:
            return []

        fps = self._computeLigandFingerprints(names)
        if fps is None:
            # RDKit unavailable, or a structure couldn't be read/fingerprinted - fall
            # back to the order the user picked them in, never crash.
            return list(zip(names[:-1], names[1:]))

        remaining = names[1:]
        chain = [names[0]]
        while remaining:
            last = fps[chain[-1]]
            best = max(remaining, key=lambda n: DataStructs.TanimotoSimilarity(last, fps[n]))
            chain.append(best)
            remaining.remove(best)
        return list(zip(chain[:-1], chain[1:]))

    def _computeLigandFingerprints(self, names):
        if not RDKIT_AVAILABLE:
            return None
        mols = {mol.__str__(): mol for mol in self.inputSetOfMols.get()}
        fps = {}
        # rdFingerprintGenerator, not the older AllChem.GetMorganFingerprintAsBitVect (deprecated
        # upstream in favor of this, confirmed by the runtime deprecation warning it prints).
        morganGen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        try:
            for name in names:
                # Reuse this protocol's own, already-battle-tested pose-file->mol2
                # conversion (same one ACPYPE parametrization uses) rather than
                # re-solving format handling (.cif/.pdb/... ) here.
                mol2File = self.addHydrogens(os.path.abspath(mols[name].getPoseFile()))
                rdmol = Chem.MolFromMol2File(mol2File, sanitize=True)
                if rdmol is None:
                    rdmol = Chem.MolFromMol2File(mol2File, sanitize=False)
                if rdmol is None:
                    raise ValueError(f'RDKit could not parse the structure for "{name}"')
                fps[name] = morganGen.GetFingerprint(rdmol)
        except Exception as e:
            self.warning(f'Could not compute structural fingerprints for the edge chain '
                        f'({e}); falling back to the order ligands were selected in.')
            return None
        return fps

    # -- Multi-edge helpers --------------------------------------------------------
    def _getEdgeDir(self, edgeIdx):
        return self._getExtraPath(f'edge_{edgeIdx}')

    def _molNameFor(self, ligandName):
        for mol in self.inputSetOfMols.get():
            if mol.__str__() == ligandName:
                return mol.getMolName()
        raise ValueError(f'Ligand "{ligandName}" not found in the input set of molecules')

    def setEdgeLigandsStep(self, edgeIdx, ligA, ligB):
        """Point every already-'ligand A'/'ligand B'-scoped method (inherited from
        GromacsSystemPrep, or this protocol's own - all of them read self.inputLigand/
        self.ligandB rather than taking a ligand as a direct argument) at this edge's
        pair. Safe only because edges run strictly sequentially (see _insertAllSteps):
        this runs once the previous edge's archiveEdgeStep has fully moved that edge's
        outputs out of the fixed paths these same methods write to."""
        self.inputLigand.set(ligA)
        self.ligandB.set(ligB)
        self._store(self.inputLigand, self.ligandB)

    @staticmethod
    def _moveIfExists(src, dst):
        if os.path.exists(src):
            shutil.move(src, dst)

    def archiveEdgeStep(self, edgeIdx, ligA, ligB):
        """Move this edge's outputs out of the fixed, non-edge-scoped paths the inherited
        system-build steps (PDB2GMXStep/editConfStep/solvateStep/addIonsStep, all under
        self._getPath()) and this protocol's own leg-scoped methods (under
        self._getExtraPath(), see the module docstring's step list) write to - otherwise
        the next edge's identical steps, run for a different ligand A/B, would silently
        overwrite them. Deliberately does NOT touch self.inputLigand/self.ligandB (ligA/
        ligB are passed in explicitly) so it stays correct regardless of exactly when it
        runs relative to the next edge's setEdgeLigandsStep."""
        edgeDir = self._getEdgeDir(edgeIdx)
        os.makedirs(edgeDir, exist_ok=True)

        # Written here (not just left for createEdgeOutputStep) because it must be
        # computed while the live per-leg analysis/ dirs still exist, before they're
        # moved below.
        self._writeCombinedResultsFile(self._getExtraPath('results_summary.txt'))

        ligAMolName = self._molNameFor(ligA)
        systemBasename = self.getSystemName()
        for name in (f'{systemBasename}_processed.gro', f'{systemBasename}_newbox.gro',
                    f'{systemBasename}_solv.gro', f'{systemBasename}_solv_ions.gro',
                    TOPOL_TOP, 'posre.itp', f'{ligAMolName}_GMX.itp'):
            self._moveIfExists(self._getPath(name), os.path.join(edgeDir, name))

        for name in (MAP_DIR, HYBRID_DIR, SYSTEM_DIR, WATER_SYSTEM_DIR, LEG_SYSTEM_DIR[WATER],
                    COMPLEX, WATER, 'ligandA.gro', 'ligandB_translated.gro', 'indexes.ndx',
                    'results_summary.txt'):
            self._moveIfExists(self._getExtraPath(name), os.path.join(edgeDir, name))

    def createEdgeOutputStep(self, edgeIdx, ligA, ligB):
        edgeDir = self._getEdgeDir(edgeIdx)
        ligAMolName = self._molNameFor(ligA)

        groBaseName = (f'{self.getSystemName()}_solv_ions.gro' if self.placeIons.get() != 0
                      else f'{self.getSystemName()}_solv.gro')
        outSystem = GromacsSystem(filename=os.path.join(edgeDir, groBaseName),
                                  topoFile=os.path.join(edgeDir, TOPOL_TOP),
                                  restrFile=os.path.join(edgeDir, 'posre.itp'),
                                  ff=self.getEnumText('mainForceField'),
                                  wff=self.getEnumText('waterForceField'))
        outSystem.setChainNames(','.join(self.getModelChains()))
        chains, lengthsDic = self.getModelChainsAndLengths()
        outSystem.setChainLengths(','.join(str(lengthsDic[c]) for c in chains))
        outSystem.setLigTopologyFile(os.path.join(edgeDir, f'{ligAMolName}_GMX.itp'))
        outSystem.setIndexFile(os.path.join(edgeDir, 'indexes.ndx'))

        _, _, _, _, ddGBind = self._combineLegs(baseDir=edgeDir)
        if ddGBind is not None:
            outSystem.setFreeEnergy(ddGBind)
        outSystem.setFreeEnergyFile(os.path.join(edgeDir, 'results_summary.txt'))

        self._defineOutputs(**{f'outputSystem_edge{edgeIdx}': outSystem})
        self._defineSourceRelation(self.inputSetOfMols, outSystem)

    def createChainSummaryStep(self, edges):
        lines = [f'RBFE chain: {" -> ".join([edges[0][0]] + [b for _, b in edges])} '
                f'({len(edges)} edge(s))']
        total, allOk = 0.0, True
        for edgeIdx, (ligA, ligB) in enumerate(edges):
            _, _, _, _, ddGBind = self._combineLegs(baseDir=self._getEdgeDir(edgeIdx))
            if ddGBind is None:
                allOk = False
                lines.append(f'  edge {edgeIdx}: {ligA} -> {ligB}: ddG_bind not available')
            else:
                total += ddGBind
                lines.append(f'  edge {edgeIdx}: {ligA} -> {ligB}: ddG_bind = {ddGBind:.2f} kJ/mol')
        if allOk:
            lines.append(f'Sum of edge ddG_bind along the chain = {total:.2f} kJ/mol - a '
                        'simple additive estimate, NOT a real network cycle-closure result '
                        '(no redundant edges were run to cross-check consistency; see the '
                        'module docstring / claude/decisions/pmx_RBFE.md §15).')
        with open(self._getExtraPath('all_edges_summary.txt'), 'w') as f:
            f.write('\n'.join(lines) + '\n')

    # -- 0. Ligand A's system (GromacsSystemPrep.parametrizeLigandStep/PDB2GMXStep/
    #       editConfStep/solvateStep/addIonsStep, all inherited unchanged) and ligand B --
    def parametrizeLigandBStep(self):
        mol = self.getLigandBMol()
        molFile = self.addHydrogens(os.path.abspath(mol.getPoseFile()))

        kwargs = self.getParameters()
        kwargs['molName'] = mol.getMolName()

        args = f'-i {molFile} -b {kwargs["molName"]} -c {kwargs["chargeMethod"]} ' \
               f'-m {kwargs["multip"]} -a {kwargs["atomType"]} -q {kwargs["qprog"]} -o gmx'
        if 'netCharge' in kwargs:
            args += f' -n {kwargs["netCharge"]}'
        pwchemPlugin.runACPYPE(self, args=args, cwd=self._getExtraPath())

    def getLigandBMol(self):
        myMol = None
        for mol in self.inputSetOfMols.get():
            if mol.__str__() == self.ligandB.get():
                myMol = mol.clone()
                break
        if myMol is None:
            raise ValueError(f'Ligand B "{self.ligandB.get()}" not found in the input set of molecules')
        return myMol

    def getLigandBPath(self, path=''):
        molName = self.getLigandBMol().getMolName()
        return self._getExtraPath(f'{molName}.acpype', path)

    def _ligBGroFile(self):
        """Ligand B's ACPYPE .gro, rigidly translated into the same frame as ligand A's
        real, final built system.

        ACPYPE's output preserves ligand B's real absolute docked-pose coordinates as-is
        (confirmed for real: its centroid matches its own Vina pose almost exactly - ACPYPE/
        openbabel don't recentre a single-molecule conversion). But ligand A's system went
        through a further whole-system recentring at editConfStep (`gmx editconf -c`, never
        called with `-princ` here, so it is a pure translation, no rotation) that ligand B's
        standalone .gro never went through. Mixing the two frames is exactly what broke this
        the first time it was tried with two genuinely different ligands (never caught by the
        earlier self-transform test, where A and B are the same molecule with zero unique
        atoms - this code path was simply never exercised): pmx's ligandHybrid placed ligand
        B's unique/dummy atoms straight from its own untranslated coordinates, landing tens of
        nm from ligand A's real position in the final box and making mdrun's domain
        decomposition fail outright ("no domain decomposition ... minimum cell size ~50 nm").

        Fix: apply the same translation ligand A underwent (computed from its own before/
        after positions) to ligand B before handing it to pmx."""
        rawGro = os.path.abspath(self.getLigandBPath(f'{self.getLigandBMol().getMolName()}_GMX.gro'))
        translatedGro = self._getExtraPath('ligandB_translated.gro')
        if not os.path.exists(translatedGro):
            translation = self._getLigandATranslation()
            self._translateGro(rawGro, translation, translatedGro)
        return os.path.abspath(translatedGro)

    def _getLigandATranslation(self):
        origGro = self.getLigandPath(f'{self.getLigandName()}_GMX.gro')
        finalGro = self._getExtraPath('ligandA.gro')
        ox, oy, oz = self._groCentroid(origGro)
        fx, fy, fz = self._groCentroid(finalGro)
        return fx - ox, fy - oy, fz - oz

    def _getWaterFrameTranslation(self, waterDir):
        """Translation from ligand A's COMPLEX-system frame (what mergedA.gro's real-atom
        coordinates are already in, since they derive from ligandA.gro, extracted from the
        complex system) into the water-only box's own, independently-recentred frame -
        see buildHybridWaterSystemStep for why this is needed."""
        complexGro = self._getExtraPath('ligandA.gro')
        waterGro = os.path.join(waterDir, 'ligA_newbox.gro')
        cx, cy, cz = self._groCentroid(complexGro)
        wx, wy, wz = self._groCentroid(waterGro)
        return wx - cx, wy - cy, wz - cz

    @staticmethod
    def _groCentroid(groPath):
        with open(groPath) as f:
            lines = f.readlines()
        n = int(lines[1])
        xs = ys = zs = 0.0
        for line in lines[2:2 + n]:
            xs += float(line[20:28])
            ys += float(line[28:36])
            zs += float(line[36:44])
        return xs / n, ys / n, zs / n

    def _fixPmxDummyAtomUnits(self, mergedAGro):
        """Correct a real unit bug in pmx's own ligandHybrid (pmx/ligand_alchemy.py's
        adjustCoords(), confirmed by reading the installed source): it fits ligand B onto
        ligand A via an RDKit round-trip through a temporary PDB file (Angstroms), then
        copies the fitted RDKit coordinates straight into pmx's Model atoms with no A->nm
        conversion. Ligand A's real atoms are never touched by this path (they keep their
        correct pmx/nm coordinates throughout) and pmx's own PDB writer + this protocol's
        subsequent `editconf` PDB->gro step apply the standard x10/1/10 conversions
        correctly for them - net zero. But the already-wrong (10x too large) dummy atoms
        pmx builds from that fit go through the SAME write+editconf chain once more,
        compounding to a net 10x-too-large final value.

        Confirmed for real: with two structurally different, independently-docked ligands
        (never exercised before - the earlier self-transform A==B test has zero dummy
        atoms, so it could never expose this), the dummy/unique atoms pmx adds for ligand
        B landed 26-48 nm from ligand A's real ~2.5-4.5 nm cluster in the same box,
        breaking mdrun's domain decomposition outright ("no domain decomposition ...
        compatible with ... minimum cell size ~50-60 nm"). Dividing those dummy atoms'
        coordinates by 10 lands them squarely back in ligand A's real cluster (checked by
        hand against the real failing run's own numbers), confirming this is exactly the
        bug and exactly the fix - not a coincidental rescaling.

        Dummy atoms are exactly the ones beyond ligand A's own real atom count: this
        protocol's hybrid ligand always lists ligand A's real atoms first (unchanged) and
        appends the new dummy atoms after (confirmed by inspecting a real mergedA.gro)."""
        with open(self._getExtraPath('ligandA.gro')) as f:
            nRealAtoms = int(f.readlines()[1])

        with open(mergedAGro) as f:
            lines = f.readlines()
        n = int(lines[1])
        outLines = lines[:2]
        for i, line in enumerate(lines[2:2 + n]):
            if i >= nRealAtoms:
                x = float(line[20:28]) / 10.0
                y = float(line[28:36]) / 10.0
                z = float(line[36:44]) / 10.0
                line = f'{line[:20]}{x:8.3f}{y:8.3f}{z:8.3f}{line[44:]}'
            outLines.append(line)
        outLines.extend(lines[2 + n:])
        with open(mergedAGro, 'w') as f:
            f.writelines(outLines)

    @staticmethod
    def _translateGro(inGro, translation, outGro):
        dx, dy, dz = translation
        with open(inGro) as f:
            lines = f.readlines()
        n = int(lines[1])
        outLines = lines[:2]
        for line in lines[2:2 + n]:
            x = float(line[20:28]) + dx
            y = float(line[28:36]) + dy
            z = float(line[36:44]) + dz
            outLines.append(f'{line[:20]}{x:8.3f}{y:8.3f}{z:8.3f}{line[44:]}')
        outLines.extend(lines[2 + n:])
        with open(outGro, 'w') as f:
            f.writelines(outLines)

    def _ligBItpFile(self):
        return os.path.abspath(self.getLigandBPath(f'{self.getLigandBMol().getMolName()}_GMX.itp'))

    def _getLigandASystem(self):
        """Fresh (non-persisted) GromacsSystem wrapper around ligand A's just-built
        solvated/ionized system. Reconstructed from disk every time it's needed rather
        than cached on self, consistent with how the rest of this plugin re-derives
        state across steps (Scipion may reload the protocol instance between steps)."""
        systemBasename = self.getSystemName()
        groBaseName = (f'{systemBasename}_solv_ions.gro' if self.placeIons.get() != 0
                       else f'{systemBasename}_solv.gro')
        system = GromacsSystem(filename=self._getPath(groBaseName), topoFile=self._getPath(TOPOL_TOP),
                               restrFile=self._getPath('posre.itp'),
                               ff=self.getEnumText('mainForceField'), wff=self.getEnumText('waterForceField'))
        system.setLigTopologyFile(self._getPath(f'{self.getLigandName()}_GMX.itp'))
        return system

    def _ensureLigandAIndexFile(self):
        indexFile = self._getExtraPath('indexes.ndx')
        if not os.path.exists(indexFile):
            chains, lengthsDic = self.getModelChainsAndLengths()
            gromacsPlugin.firstIndexCreation(self, self._getLigandASystem(), ligandName=self.getLigandName(),
                                             modelChains=chains, chainLengths=lengthsDic)
        return indexFile

    # -- 1. Atom mapping --------------------------------------------------------
    def atomMappingStep(self):
        mapDir = self._getExtraPath(MAP_DIR)
        os.makedirs(mapDir, exist_ok=True)

        ligAGro = self._extractLigandAGro()
        ligBGro = self._ligBGroFile()

        args = f'-i1 {ligAGro} -i2 {ligBGro} -o1 pairs1.dat -o2 pairs2.dat -log mapping.log'
        gromacsPlugin.runPmx(self, 'atomMapping', args, cwd=mapDir)

    def _extractLigandAGro(self):
        """Extract ligand A alone from its just-built system, via GROMACS' own
        auto-generated "Other" index group (the standard catch-all for
        non-protein/nucleic/water residues).

        Deliberately NOT a stored ligand-ID label: GROMACS' own auto-detected residue
        name for a real ligand (e.g. a genuine PDB code like 'RET') is not reliably
        the same string a naming heuristic might guess (e.g. 'LIG') - confirmed by
        hitting "Error: No such group 'LIG'" when this used to rely on one."""
        systemFile = os.path.abspath(self._getLigandASystem().getSystemFile())
        outGro = os.path.abspath(self._getExtraPath('ligandA.gro'))
        indexFile = os.path.abspath(self._ensureLigandAIndexFile())

        groups = gromacsPlugin.parseIndexFile(self, indexFile)
        invGroups = {v: k for k, v in groups.items()}
        if 'Other' not in invGroups:
            raise ValueError(f'No "Other" index group found in {indexFile}; cannot isolate '
                             f'ligand A (the system may have more than one non-protein/water '
                             f'residue type, which this protocol does not support)')

        args = f'trjconv -f {systemFile} -s {systemFile} -n {indexFile} -o {outGro}'
        gromacsPlugin.runGromacsPrintf(self, printfValues=[invGroups['Other']], args=args,
                                       cwd=self._getExtraPath())
        return outGro

    @staticmethod
    def _readResName(groFile):
        """Real residue name physically written in a single-residue .gro file's atom lines
        (ground truth) - see _extractLigandAGro for why nothing here trusts a naming guess."""
        with open(groFile) as f:
            lines = f.readlines()
        return lines[2][5:10].strip()

    # -- 2. Hybrid ligand build --------------------------------------------------
    def buildHybridLigandStep(self):
        mapDir = self._getExtraPath(MAP_DIR)
        hybDir = self._getExtraPath(HYBRID_DIR)
        os.makedirs(hybDir, exist_ok=True)

        ligAGro = os.path.abspath(self._getExtraPath('ligandA.gro'))
        ligBGro = self._ligBGroFile()
        ligAItp = os.path.abspath(self._getLigandASystem().getLigTopologyFile())
        ligBItp = self._ligBItpFile()
        pairsFile = os.path.abspath(os.path.join(mapDir, 'pairs1.dat'))

        args = (f'-i1 {ligAGro} -i2 {ligBGro} -itp1 {ligAItp} -itp2 {ligBItp} '
                f'-pairs {pairsFile} -oA mergedA.pdb -oB mergedB.pdb '
                f'-oitp merged.itp -offitp ffmerged.itp -log hybrid.log')
        gromacsPlugin.runPmx(self, 'ligandHybrid', args, cwd=hybDir)

    # -- 3a. Ligand A alone, solvated/ionized - the "water" leg's own system ------
    def buildLigandAWaterSystemStep(self):
        """Box/solvate/ionize ligand A alone (no protein) - the starting point for the
        water/solvent leg of the standard RBFE double-leg cycle
        (ddG_bind = dG_complex - dG_water, see module docstring and
        claude/decisions/pmx_RBFE.md §12/§14). Mirrors GromacsSystemPrep's own
        editConfStep/solvateStep/addIonsStep command-for-command (box size/type, water
        model, ion type/concentration all reuse the *same* form params as the complex
        leg), but can't reuse those methods directly: they're hardcoded to
        self._getPath()/getSystemName()-prefixed files, i.e. the complex leg's own
        protein+ligand system - this is a separate, ligand-only system in its own
        extra/water_system/ directory."""
        waterDir = self._getExtraPath(WATER_SYSTEM_DIR)
        os.makedirs(waterDir, exist_ok=True)
        topFile = self._buildLigandAWaterTop(waterDir)
        ligGro = os.path.abspath(self.getLigandPath(f'{self.getLigandName()}_GMX.gro'))

        boxType = self.getEnumText('boxType').lower() if self.boxType.get() != 1 else 'triclinic'
        ecArgs = (f'editconf -f {ligGro} -o ligA_newbox.gro -c -bt {boxType}'
                  + self.getDistanceArgs())
        gromacsPlugin.runGromacs(self, 'gmx', ecArgs, cwd=waterDir)

        waterModel = self.getEnumText('waterForceField')
        if waterModel in ('spc', 'spce', 'tip3p'):
            waterModel = 'spc216'
        solvArgs = (f'solvate -cp ligA_newbox.gro -cs {waterModel}.gro -o ligA_solv.gro '
                    f'-p {os.path.basename(topFile)}')
        gromacsPlugin.runGromacs(self, 'gmx', solvArgs, cwd=waterDir)

        if self.placeIons.get() != 0:
            ionsMdp = os.path.abspath(self.buildIonsMDP())
            grArgs = (f'grompp -f {ionsMdp} -c ligA_solv.gro -p {os.path.basename(topFile)} '
                      f'-o ions.tpr')
            gromacsPlugin.runGromacsPrintf(self, printfValues=['SOL'], args=grArgs, cwd=waterDir)

            cation, cc = self.parseIon(self.getEnumText('cationType'))
            anion, ac = self.parseIon(self.getEnumText('anionType'))
            genStr = (f'genion -s ions.tpr -o ligA_solv_ions.gro -p {os.path.basename(topFile)} '
                      f'-pname {cation} -nname {anion}')
            if cc == 2:
                genStr += f' -pq {cc}'
            if ac == 2:
                genStr += f' -nq {ac}'
            if self.placeIons.get() == 1:
                genStr += ' -neutral'
            elif self.placeIons.get() == 2:
                genStr += f' -np {self.cationNum.get()} -nn {self.anionNum.get()}'
            if self.addSalt:
                genStr += f' -conc {self.saltConc.get()}'
            gromacsPlugin.runGromacsPrintf(self, printfValues=['SOL'], args=genStr, cwd=waterDir)

    def _buildLigandAWaterTop(self, waterDir):
        """Hand-built topol.top for "ligand A alone in water" - mirrors exactly what
        `gmx pdb2gmx -ff <mainFF> -water <waterFF>` normally writes for a protein system
        (confirmed by inspecting a real generated topol.top from this same protocol's
        complex leg: `#include "<ff>.ff/forcefield.itp"` first, then the ligand's own
        itp, then `<ff>.ff/<waterModel>.itp` + `<ff>.ff/ions.itp` right before
        [ system ]/[ molecules ]) - just without any protein moleculetype block."""
        mainFF = self.getEnumText('mainForceField')
        waterFF = self.getEnumText('waterForceField')
        ligItpFile = self.getLigandPath(f'{self.getLigandName()}_GMX.itp')
        ligItpBase = os.path.basename(ligItpFile)
        ligMolName = self._readMoleculetypeName(ligItpFile)

        dst = os.path.join(waterDir, ligItpBase)
        if not os.path.exists(dst):
            os.link(ligItpFile, dst)

        topPath = os.path.join(waterDir, 'topol.top')
        text = (f'#include "{mainFF}.ff/forcefield.itp"\n'
                f'#include "{ligItpBase}"\n'
                f'#include "{mainFF}.ff/{waterFF}.itp"\n'
                f'#include "{mainFF}.ff/ions.itp"\n\n'
                f'[ system ]\n Ligand A alone in water\n\n'
                f'[ molecules ]\n{ligMolName}    1\n')
        with open(topPath, 'w') as f:
            f.write(text)
        return topPath

    # -- 3b. Splice the hybrid ligand into ligand A's systems (both legs) --------
    def buildHybridSystemStep(self):
        hybDir = self._getExtraPath(HYBRID_DIR)
        sysDir = self._getExtraPath(SYSTEM_DIR)
        os.makedirs(sysDir, exist_ok=True)
        system = self._getLigandASystem()
        mergedAGro = self._mergeHybridLigandGro(hybDir, sysDir)

        outTop = self._buildHybridTopology(system.getTopologyFile(), system.getLigTopologyFile(),
                                           hybDir, sysDir)
        # Real residue name (not a naming guess - see _extractLigandAGro), read
        # straight off the ligand-A.gro this same protocol already extracted for mapping.
        ligResName = self._readResName(self._getExtraPath('ligandA.gro'))
        self._spliceLigandGro(system.getSystemFile(), mergedAGro, ligResName,
                              os.path.join(sysDir, 'system.gro'))
        return outTop

    def buildHybridWaterSystemStep(self):
        hybDir = self._getExtraPath(HYBRID_DIR)
        waterDir = self._getExtraPath(WATER_SYSTEM_DIR)
        sysDir = self._getExtraPath(LEG_SYSTEM_DIR[WATER])
        os.makedirs(sysDir, exist_ok=True)
        mergedAGro = self._mergeHybridLigandGro(hybDir, sysDir)

        # The hybrid ligand's real-atom coordinates come from ligandA.gro (extracted from
        # the COMPLEX system, in that box's own post-editconf-recentring frame) - but the
        # water-only box was built and recentred completely independently, as its own
        # separate `gmx editconf -c`. Splicing the complex-frame hybrid ligand straight into
        # the water box - without correcting for this - lands it wherever the complex box's
        # recentring happened to put it, unrelated to the water box's own small volume.
        # Confirmed for real: this is exactly what caused a LINCS blowup/segfault the first
        # time the water leg actually ran (a frame mismatch the complex-only protocol could
        # never have hit, since it only ever had one box). Same class of bug, same fix
        # pattern as _ligBGroFile's ligand-B translation: move the hybrid ligand by the same
        # rigid translation ligand A itself underwent between these two independently-built
        # boxes, computed from its own before/after centroids.
        translation = self._getWaterFrameTranslation(waterDir)
        self._translateGro(mergedAGro, translation, mergedAGro)

        ligItpFile = self.getLigandPath(f'{self.getLigandName()}_GMX.itp')
        outTop = self._buildHybridTopology(os.path.join(waterDir, 'topol.top'), ligItpFile,
                                           hybDir, sysDir)
        ligResName = self._readResName(self._getExtraPath('ligandA.gro'))
        waterSystemGro = os.path.join(
            waterDir, 'ligA_solv_ions.gro' if self.placeIons.get() != 0 else 'ligA_solv.gro')
        self._spliceLigandGro(waterSystemGro, mergedAGro, ligResName,
                              os.path.join(sysDir, 'system.gro'))
        return outTop

    def _mergeHybridLigandGro(self, hybDir, sysDir):
        # -f crosses into hybDir (a sibling of the cwd below) so it must be absolute; -o stays
        # a bare filename since it lands inside cwd=sysDir itself (a joined sysDir/... string
        # here would be re-resolved against sysDir again by the subprocess and double-nest).
        mergedAPdb = os.path.abspath(os.path.join(hybDir, 'mergedA.pdb'))
        gromacsPlugin.runGromacs(self, 'gmx', f'editconf -f {mergedAPdb} -o mergedA.gro', cwd=sysDir)
        mergedAGro = os.path.join(sysDir, 'mergedA.gro')
        self._fixPmxDummyAtomUnits(mergedAGro)
        return mergedAGro

    def _buildHybridTopology(self, srcTopFile, ligItpFile, hybDir, sysDir):
        outTop = os.path.join(sysDir, 'topol.top')
        shutil.copy(srcTopFile, outTop)
        shutil.copy(os.path.join(hybDir, 'merged.itp'), os.path.join(sysDir, 'merged.itp'))

        # merged.itp only carries the *new* dummy atomtypes (ffmerged.itp): the real GAFF/AMBER
        # types it also references (e.g. ca, ha, oh...) come from each ligand's own itp, which
        # would otherwise be lost once their #include is replaced below ("Atomtype X not found"
        # in grompp, confirmed by actually running it).
        self._writeCombinedAtomTypes(
            [ligItpFile, self._ligBItpFile(), os.path.join(hybDir, 'ffmerged.itp')],
            os.path.join(sysDir, 'atomtypes.itp'))

        ligItpBase = os.path.basename(ligItpFile)
        oldMolName = ligItpBase.replace('_GMX.itp', '')
        replaceInFile(outTop, f'#include "{ligItpBase}"',
                      '#include "atomtypes.itp"\n#include "merged.itp"')

        # pmx ligandHybrid keeps mol1's original moleculetype name; only patch the
        # [ molecules ] entry in the rare case a future pmx version changes that.
        hybridMolName = self._readMoleculetypeName(os.path.join(sysDir, 'merged.itp'))
        if hybridMolName != oldMolName:
            replaceInFile(outTop, f'{oldMolName} 1', f'{hybridMolName} 1')
        return outTop

    @staticmethod
    def _extractAtomTypesLines(itpFile):
        with open(itpFile) as f:
            lines = f.readlines()
        block, inBlock = [], False
        for line in lines:
            stripped = line.strip().replace(' ', '')
            if stripped == '[atomtypes]':
                inBlock = True
                continue
            if inBlock:
                if stripped.startswith('[') and stripped.endswith(']'):
                    break
                block.append(line)
        return block

    @classmethod
    def _writeCombinedAtomTypes(cls, itpFiles, outFile):
        """Concatenate the [ atomtypes ] blocks of several itp files, deduplicated by
        atom-type name (first occurrence wins; ligands sharing a force field define
        identical params for shared types anyway)."""
        seen, combined = set(), []
        for itpFile in itpFiles:
            for line in cls._extractAtomTypesLines(itpFile):
                stripped = line.strip()
                if not stripped or stripped.startswith(';'):
                    combined.append(line)
                    continue
                name = stripped.split()[0]
                if name in seen:
                    continue
                seen.add(name)
                combined.append(line)

        with open(outFile, 'w') as f:
            f.write('[ atomtypes ]\n')
            f.writelines(combined)
        return outFile

    @staticmethod
    def _readMoleculetypeName(itpFile):
        with open(itpFile) as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if line.strip().replace(' ', '') == '[moleculetype]':
                for nxt in lines[i + 1:]:
                    nxt = nxt.strip()
                    if not nxt or nxt.startswith(';'):
                        continue
                    return nxt.split()[0]
        raise ValueError(f'No [ moleculetype ] section found in {itpFile}')

    @staticmethod
    def _spliceLigandGro(systemGro, hybridLigGro, ligResName, outGro):
        """Replace the ligand-A coordinate block (found by residue name) in systemGro
        with the hybrid ligand's coordinates, keeping box/solvent/ions untouched.
        The hybrid ligand's own residue number/name (pmx writes its own placeholder,
        e.g. "UN") are overwritten to match the original block, so index groups and
        downstream analysis keyed on the original ligand identity keep working."""
        with open(systemGro) as f:
            lines = f.readlines()
        header, boxLine = lines[0], lines[-1]
        atomLines = lines[2:-1]

        with open(hybridLigGro) as f:
            hybLines = f.readlines()
        hybAtomLines = hybLines[2:-1]

        startIdx = endIdx = None
        for i, line in enumerate(atomLines):
            if line[5:10].strip() == ligResName:
                if startIdx is None:
                    startIdx = i
                endIdx = i
        if startIdx is None:
            raise ValueError(f'Ligand residue {ligResName} not found in {systemGro}')

        resField = atomLines[startIdx][0:10]
        hybAtomLines = [resField + line[10:] for line in hybAtomLines]

        newAtomLines = atomLines[:startIdx] + hybAtomLines + atomLines[endIdx + 1:]
        with open(outGro, 'w') as f:
            f.write(header)
            f.write(f'{len(newAtomLines)}\n')
            f.writelines(newAtomLines)
            f.write(boxLine)
        return outGro

    # -- 4. End-state equilibration -----------------------------------------------
    def _getStateDir(self, leg, state):
        return self._getExtraPath(leg, STATE_DIRNAME[state])

    def _linkTopologyExtras(self, stageDir, topFile, sysDir):
        localTop = os.path.join(stageDir, os.path.basename(topFile))
        if not os.path.exists(localTop):
            os.link(topFile, localTop)
        # atomtypes.itp (not ffmerged.itp - that one's only read to build atomtypes.itp in
        # _buildHybridTopology, never copied into sysDir on its own) + merged.itp, matching
        # exactly what _buildHybridTopology writes and topol.top's #include lines reference.
        for extra in ('merged.itp', 'atomtypes.itp'):
            src = os.path.join(sysDir, extra)
            dst = os.path.join(stageDir, extra)
            if not os.path.exists(dst):
                os.link(src, dst)
        return localTop

    def _writeFepMdp(self, outFile, phase, state, **kw):
        # No couple-moltype: the A/B states are already encoded per-atom in the hybrid
        # itp (dual topology), so GROMACS' separate molecule-decoupling feature isn't
        # used here. Confirmed by actually running grompp: 'couple-moltype = none' is
        # NOT a valid "disable" sentinel, it's read as a literal molecule name to look
        # up ("Did not find any molecules of type 'none' for coupling").
        alpha, sigma = self.scAlpha.get(), self.scSigma.get()
        feBlock = (f'free-energy      = yes\n'
                   f'init-lambda      = {state}\n'
                   f'delta-lambda     = {kw.get("deltaLambda", 0)}\n'
                   f'sc-alpha         = {alpha}\n'
                   f'sc-sigma         = {sigma}\n'
                   f'sc-power         = 1\n'
                   f'sc-coul          = yes\n')
        # rcoulomb/rvdw/vdw-type/DispCorr match pmx's own validated ligand_tutorial.ipynb
        # mdp files (eq_l0.mdp/ti_l0.mdp) as closely as reasonable here - confirmed by
        # reading those files directly rather than guessing: PME with rcoulomb=1.1,
        # switched vdW between 1.0-1.1 nm, and a dispersion correction that this
        # protocol was previously missing entirely (DispCorr defaults to "no" in
        # GROMACS if unset, silently omitting a real, systematic long-range vdW
        # energy/pressure contribution).
        common = ('nstlist          = 10\n'
                  'cutoff-scheme    = Verlet\n'
                  'coulombtype      = PME\n'
                  'rcoulomb         = 1.1\n'
                  'fourierspacing   = 0.12\n'
                  'ewald-rtol       = 1e-5\n'
                  'vdw-type         = switch\n'
                  'rvdw-switch      = 1.0\n'
                  'rvdw             = 1.1\n'
                  'DispCorr         = EnerPres\n'
                  'pbc              = xyz\n'
                  'constraints      = h-bonds\n'
                  'constraint_algorithm = lincs\n')

        if phase == 'em':
            text = ('integrator       = steep\n'
                   f'emtol            = {self.emTol.get()}\n'
                   'emstep           = 0.01\n'
                   f'nsteps           = {self.nStepsMin.get()}\n'
                   + common + feBlock)
        else:
            nsteps = round(kw['simTime'] / self.timeStep.get())
            interval = kw.get('outInterval', 0)
            text = ('integrator       = md\n'
                   f'dt               = {self.timeStep.get()}\n'
                   f'nsteps           = {nsteps}\n'
                   + common +
                   f'continuation     = {kw.get("continuation", "no")}\n'
                   f'nstxout          = {interval}\n'
                   f'nstvout          = {interval}\n'
                   f'nstdhdl          = {kw.get("nstdhdl", 100)}\n'
                   'tcoupl           = V-rescale\n'
                   'tc-grps          = System\n'
                   f'ref_t            = {self.temperature.get()}\n'
                   'tau_t            = 0.1\n'
                   'nsttcouple       = 10\n')
            if kw.get('pcoupl', False):
                # Parrinello-Rahman/tau_p=5/compressibility=4.6e-5 match pmx's own
                # ligand_tutorial.ipynb mdp files exactly (this protocol previously used
                # C-rescale/tau_p=2.0/4.5e-5 - a different, also-valid barostat choice,
                # but not the one the published/benchmarked non-equilibrium RBFE workflow
                # this protocol follows actually validated its results against).
                text += ('pcoupl           = Parrinello-Rahman\n'
                        'pcoupltype       = isotropic\n'
                        f'ref_p            = {self.pressure.get()}\n'
                        'tau_p            = 5.0\n'
                        'compressibility  = 4.6e-5\n'
                        'nstpcouple       = 10\n')
            else:
                text += 'pcoupl           = no\n'
            if kw.get('genVel', False):
                text += (f'gen_vel          = yes\n'
                        f'gen_temp         = {self.temperature.get()}\n'
                        'gen_seed         = -1\n')
            else:
                text += 'gen_vel          = no\n'
            text += feBlock

        with open(outFile, 'w') as f:
            f.write(text)
        return outFile

    def _runGrompp(self, stageDir, mdpFile, groFile, topFile, tprName, sysDir, prevCpt=None, extraArgs=''):
        localTop = self._linkTopologyExtras(stageDir, topFile, sysDir)
        prevStr = f' -t {os.path.abspath(prevCpt)}' if prevCpt else ''
        outFile = f'{tprName}.tpr'
        command = (f'grompp -f {os.path.abspath(mdpFile)} -c {os.path.abspath(groFile)} '
                   f'-r {os.path.abspath(groFile)} -p {os.path.basename(localTop)}{prevStr} '
                   f'-o {outFile} -maxwarn 2 {extraArgs}')
        gromacsPlugin.runGromacs(self, 'gmx', command, cwd=stageDir, numberOfMpi=0)
        return os.path.join(stageDir, outFile)

    def _runMdrun(self, stageDir, deffnm):
        if getattr(self, params.USE_GPU):
            gpuList = getattr(self, params.GPU_LIST).get().replace(' ', '')
            gpuStr = f' -nb gpu -gpu_id {gpuList}'
        else:
            gpuStr = ' -nb cpu'
        # -ntmpi 1 (all requested threads as OpenMP, not thread-MPI ranks): a bare -nt leaves
        # gmx to auto-guess a rank count, and on a high-core-count machine it can pick more
        # thread-MPI ranks than a small box supports - confirmed for real (GromacsPmxABFE hit
        # this exact failure on its own small ligand-alone-in-solvent leg on a 64-thread remote
        # run: "no domain decomposition ... compatible with ... minimum cell size"). This
        # protocol's own water leg (ligand alone in solvent, no protein) is the same kind of
        # small box, so it's equally exposed. A single rank never needs domain decomposition,
        # sidestepping this regardless of window/box size.
        command = f'mdrun -v -deffnm {deffnm}{gpuStr} -ntmpi 1 -ntomp {self.numberOfThreads.get()}'
        gromacsPlugin.runGromacs(self, 'gmx', command, cwd=stageDir)

    def equilibrateStateStep(self, leg, state):
        stateDir = self._getStateDir(leg, state)
        os.makedirs(stateDir, exist_ok=True)
        sysDir = self._getExtraPath(LEG_SYSTEM_DIR[leg])
        topFile = os.path.join(sysDir, 'topol.top')
        groFile = os.path.join(sysDir, 'system.gro')

        emMdp = self._writeFepMdp(os.path.join(stateDir, 'em.mdp'), 'em', state)
        self._runGrompp(stateDir, emMdp, groFile, topFile, 'em', sysDir)
        self._runMdrun(stateDir, 'em')

        # Frames saved often enough to yield >= nStructs snapshots from the 2nd half of NPT
        nptSteps = round(self.nptTime.get() / self.timeStep.get())
        interval = max(1, round(nptSteps / (2 * max(self.nStructs.get(), 1))))

        nvtMdp = self._writeFepMdp(os.path.join(stateDir, 'nvt.mdp'), 'nvt', state,
                                   simTime=self.nvtTime.get(), continuation='no', genVel=True,
                                   pcoupl=False, outInterval=0)
        self._runGrompp(stateDir, nvtMdp, os.path.join(stateDir, 'em.gro'), topFile, 'nvt', sysDir)
        self._runMdrun(stateDir, 'nvt')

        nptMdp = self._writeFepMdp(os.path.join(stateDir, 'npt.mdp'), 'npt', state,
                                   simTime=self.nptTime.get(), continuation='yes', genVel=False,
                                   pcoupl=True, outInterval=interval)
        self._runGrompp(stateDir, nptMdp, os.path.join(stateDir, 'nvt.gro'), topFile, 'npt', sysDir,
                        prevCpt=os.path.join(stateDir, 'nvt.cpt'))
        self._runMdrun(stateDir, 'npt')

    # -- 5. Snapshot extraction ---------------------------------------------------
    def _getFrameDir(self, leg, state):
        return self._getExtraPath(leg, 'frames', STATE_DIRNAME[state])

    def extractSnapshotsStep(self, leg, state):
        stateDir = self._getStateDir(leg, state)
        frameDir = self._getFrameDir(leg, state)
        os.makedirs(frameDir, exist_ok=True)

        burnIn = self.nptTime.get() / 2.0
        args = (f'trjconv -f npt.trr -s npt.tpr -b {burnIn} -sep '
                f'-o {os.path.abspath(os.path.join(frameDir, "frame.gro"))}')
        gromacsPlugin.runGromacsPrintf(self, printfValues=['0'], args=args, cwd=stateDir)

        nStructs = self.nStructs.get()
        frames = sorted(glob(os.path.join(frameDir, 'frame*.gro')),
                        key=lambda p: int(re.search(r'frame(\d+)\.gro', p).group(1)))
        if len(frames) < nStructs:
            self.warning(f'Only {len(frames)} post-burn-in frames available for state {state}, '
                         f'requested {nStructs}')
        else:
            for extra in frames[:-nStructs] if nStructs else frames:
                os.remove(extra)

    def _getFrameFile(self, leg, state, i):
        frames = sorted(glob(os.path.join(self._getFrameDir(leg, state), 'frame*.gro')),
                        key=lambda p: int(re.search(r'frame(\d+)\.gro', p).group(1)))
        return frames[i] if i < len(frames) else None

    # -- 6. Fast-growth transitions ------------------------------------------------
    def _getTransitionDir(self, leg, direction, i):
        return self._getExtraPath(leg, 'transitions', direction, f'frame_{i}')

    def _failedFlagPath(self, leg, direction, i):
        return os.path.abspath(self._getExtraPath(leg, f'.failed_{direction}_{i}.flag'))

    def runTransitionStep(self, leg, direction, i):
        state = DIRECTION_STATE[direction]
        frameFile = self._getFrameFile(leg, state, i)
        if frameFile is None:
            with open(self._failedFlagPath(leg, direction, i), 'w') as fh:
                fh.write(f'no frame {i} available for state {state}\n')
            return

        tDir = self._getTransitionDir(leg, direction, i)
        os.makedirs(tDir, exist_ok=True)
        sysDir = self._getExtraPath(LEG_SYSTEM_DIR[leg])
        topFile = os.path.join(sysDir, 'topol.top')

        nsteps = round(self.swTime.get() / self.timeStep.get())
        deltaLambda = (1.0 / nsteps) if direction == FWD else (-1.0 / nsteps)
        mdp = self._writeFepMdp(os.path.join(tDir, 'ti.mdp'), 'ti', state,
                                simTime=self.swTime.get(), continuation='yes', genVel=False,
                                pcoupl=True, outInterval=0, nstdhdl=1, deltaLambda=deltaLambda)
        try:
            self._runGrompp(tDir, mdp, frameFile, topFile, 'ti', sysDir)
            self._runMdrun(tDir, 'ti')
        except Exception as e:
            self.warning(f'[{leg}] Transition {direction} frame {i} failed: {e}')
            with open(self._failedFlagPath(leg, direction, i), 'w') as fh:
                fh.write(str(e) + '\n')

    def _getDhdlFiles(self, leg, direction):
        dhdlFiles = []
        for i in range(self.nStructs.get()):
            if os.path.exists(self._failedFlagPath(leg, direction, i)):
                continue
            # mdrun names the dH/dl output <deffnm>.xvg (not <deffnm>.dhdl.xvg), confirmed
            # by actually running a free-energy mdrun on this machine. Absolute: analyseStep
            # runs with cwd=anDir, a sibling of the transitions dir these files live under.
            dhdl = os.path.abspath(os.path.join(self._getTransitionDir(leg, direction, i), 'ti.xvg'))
            if os.path.exists(dhdl):
                dhdlFiles.append(dhdl)
        return dhdlFiles

    # -- 7. pmx analyse ------------------------------------------------------------
    def analyseStep(self, leg):
        anDir = self._getExtraPath(leg, 'analysis')
        os.makedirs(anDir, exist_ok=True)

        fwdFiles = self._getDhdlFiles(leg, FWD)
        bwdFiles = self._getDhdlFiles(leg, BWD)
        if not fwdFiles or not bwdFiles:
            self.warning(f'[{leg}] Not enough successful transitions to analyse '
                         f'(forward: {len(fwdFiles)}, reverse: {len(bwdFiles)})')
            return

        args = (f'-fA {" ".join(fwdFiles)} -fB {" ".join(bwdFiles)} '
                f'-m CGI BAR JARZ -t {self.temperature.get()} -o results.txt -w wplot.png')
        try:
            gromacsPlugin.runPmx(self, 'analyse', args, cwd=anDir)
        except Exception as e:
            # pmx writes results.txt incrementally, one estimator at a time (BAR is known to
            # crash outright on numpy>=2 environments - see Plugin.addPmx's numpy<2 pin -
            # confirmed directly: a run here still produced a correct CGI dG before BAR's
            # crash aborted the whole process). Don't discard an otherwise-complete result
            # over one optional estimator failing.
            dG, _ = self._parseResultsFile(leg)
            if dG is None:
                raise
            self.warning(f'[{leg}] "pmx analyse" exited with an error ({e}), but at least one '
                         f'estimator in results.txt succeeded; continuing with that.')

    def _parseResultsFile(self, leg, baseDir=None):
        """Best-effort parse of pmx analyse's results.txt (BAR preferred, then CGI, then JARZ).
        Never raises: the raw results.txt is kept as an output regardless of whether this
        regex still matches a future pmx output-format tweak.

        `baseDir` lets this be read from an already-archived edge directory
        (extra/edge_<i>/, see archiveEdgeStep) instead of the live extra/<leg>/analysis/
        path - needed once a multi-edge run has moved a finished edge's files out of the
        way for the next edge."""
        base = baseDir if baseDir is not None else self._getExtraPath()
        resultsFile = os.path.join(base, leg, 'analysis', 'results.txt')
        if not os.path.exists(resultsFile):
            return None, None
        with open(resultsFile) as f:
            text = f.read()

        for estimator in ('BAR', 'CGI', 'JARZ'):
            dgMatch = re.search(rf'{estimator}:\s*dG\s*=\s*([-\d.]+)\s*kJ/mol', text)
            if dgMatch:
                errMatch = re.search(rf'{estimator}:\s*Std Err[^=]*=\s*([-\d.]+)\s*kJ/mol', text)
                return float(dgMatch.group(1)), float(errMatch.group(1)) if errMatch else None
        return None, None

    # -- 8. Output -------------------------------------------------------------------
    def _combineLegs(self, baseDir=None):
        """ddG_bind(A->B) = dG_complex(A->B) - dG_water(A->B) - the standard RBFE double-leg
        thermodynamic cycle (see module docstring and claude/decisions/pmx_RBFE.md §12/§14).
        Both legs' components are always returned alongside the combined value (never just
        the combined number), so a wrong leg is auditable rather than silently baked into one
        opaque total - the same transparency principle GromacsPmxABFE's 4-component reporting
        already follows for its own double-decoupling cycle.

        `baseDir` is forwarded to _parseResultsFile - see there."""
        dGComplex, errComplex = self._parseResultsFile(COMPLEX, baseDir)
        dGWater, errWater = self._parseResultsFile(WATER, baseDir)
        ddGBind = (dGComplex - dGWater) if (dGComplex is not None and dGWater is not None) else None
        return dGComplex, errComplex, dGWater, errWater, ddGBind

    def _writeCombinedResultsFile(self, outFile, baseDir=None):
        dGComplex, errComplex, dGWater, errWater, ddGBind = self._combineLegs(baseDir)

        def fmt(val, err):
            if val is None:
                return 'not available (see the leg\'s own analysis/results.txt for details)'
            errStr = f' +/- {err}' if err is not None else ''
            return f'{val}{errStr} kJ/mol'

        text = (f'dG_complex (bound leg, A->B)  = {fmt(dGComplex, errComplex)}\n'
                f'dG_water (solvent leg, A->B)  = {fmt(dGWater, errWater)}\n')
        if ddGBind is not None:
            text += f'ddG_bind = dG_complex - dG_water = {ddGBind} kJ/mol\n'
        else:
            text += ('ddG_bind: not available - both legs must succeed to combine them '
                     '(see the missing leg\'s own analysis/results.txt).\n')
        with open(outFile, 'w') as f:
            f.write(text)
        return outFile

    def createOutputStep(self):
        system = self._getLigandASystem()
        _, _, _, _, ddGBind = self._combineLegs()

        outSystem = GromacsSystem(filename=system.getSystemFile(), topoFile=system.getTopologyFile(),
                                  ff=system.getForceField(), wff=system.getWaterForceField())
        outSystem.setChainNames(','.join(self.getModelChains()))
        chains, lengthsDic = self.getModelChainsAndLengths()
        outSystem.setChainLengths(','.join(str(lengthsDic[c]) for c in chains))
        outSystem.setLigTopologyFile(system.getLigTopologyFile())
        outSystem.setIndexFile(self._getExtraPath('indexes.ndx'))

        if ddGBind is not None:
            outSystem.setFreeEnergy(ddGBind)
        summaryFile = self._writeCombinedResultsFile(self._getExtraPath('results_summary.txt'))
        outSystem.setFreeEnergyFile(summaryFile)

        self._defineOutputs(outputSystem=outSystem)
        self._defineSourceRelation(self.inputSetOfMols, outSystem)

    # --------------------------- INFO functions -----------------------------------
    def _getSelectedLigandNames(self):
        if self.ligandSelection.get() == LIGANDS_ALL:
            mols = self.inputSetOfMols.get()
            return [mol.__str__() for mol in mols] if mols is not None else []
        raw = self.selectedLigands.get()
        return [s.strip() for s in raw.split('\n') if s.strip()] if raw else []

    def _validate(self):
        errors = super()._validate()
        mols = self.inputSetOfMols.get()
        names = [mol.__str__() for mol in mols] if mols is not None else []

        selected = self._getSelectedLigandNames()
        if len(selected) < 2:
            errors.append('At least 2 ligands are needed (via "All" above, or the wizard next to '
                          f'"Ligands" for a subset) - found {len(selected)}.')
        for name in selected:
            if mols is not None and name not in names:
                errors.append(f'Selected ligand "{name}" not found in the input set of molecules.')
        return errors

    def _warnings(self):
        ws = []
        nLigands = len(self._getSelectedLigandNames())
        if nLigands > 2:
            ws.append(f'You are running a {nLigands - 1}-edge RBFE chain over {nLigands} ligands. '
                      'This calculation is computationally heavy: each edge repeats the full '
                      'pipeline (system builds, atom mapping, both legs\' equilibration and '
                      'switching runs) from scratch, and edges run sequentially, one at a time - '
                      'so total running time scales roughly linearly with the number of ligands, '
                      'not with the number of edges you might expect from a real network.')
        return ws

    def _summary(self):
        summary = []
        if self.isFinished():
            edges = self._computeEdgeChain()
            if len(edges) <= 1:
                summary.append(f'Ligand A: {self.inputLigand.get()}\nLigand B: {self.ligandB.get()}')
                dGComplex, errComplex, dGWater, errWater, ddGBind = self._combineLegs()
                if dGComplex is not None:
                    errStr = f' +/- {errComplex}' if errComplex is not None else ''
                    summary.append(f'dG_complex (bound leg, A->B) = {dGComplex}{errStr} kJ/mol')
                if dGWater is not None:
                    errStr = f' +/- {errWater}' if errWater is not None else ''
                    summary.append(f'dG_water (solvent leg, A->B) = {dGWater}{errStr} kJ/mol')
                if ddGBind is not None:
                    summary.append(f'Relative binding free energy ddG_bind(A->B) = {ddGBind:.2f} kJ/mol')
                elif dGComplex is None and dGWater is None:
                    summary.append('Finished, but neither leg\'s free energy estimate could be parsed '
                                   'from pmx analyse\'s results.txt (see the raw files in each leg\'s '
                                   'analysis output).')
            else:
                chainStr = ' -> '.join([edges[0][0]] + [b for _, b in edges])
                summary.append(f'RBFE chain ({len(edges)} edges): {chainStr}')
                for edgeIdx, (ligA, ligB) in enumerate(edges):
                    _, _, _, _, ddGBind = self._combineLegs(baseDir=self._getEdgeDir(edgeIdx))
                    if ddGBind is not None:
                        summary.append(f'  edge {edgeIdx} ({ligA} -> {ligB}): '
                                       f'ddG_bind = {ddGBind:.2f} kJ/mol')
                    else:
                        summary.append(f'  edge {edgeIdx} ({ligA} -> {ligB}): not available')
                summary.append('See extra/all_edges_summary.txt for the additive chain total.')
        else:
            summary.append('The protocol has not finished.')
        return summary

    def _methods(self):
        methods = []
        if self.isFinished():
            methods.append('Relative binding free energy was calculated with the non-equilibrium '
                           'fast-growth alchemical FEP approach implemented in pmx: atom mapping and '
                           'hybrid dual-topology ligand construction ("pmx atomMapping", "pmx '
                           'ligandHybrid"), end-state equilibration and short forward/reverse switching '
                           'runs with GROMACS ("gmx grompp"/"gmx mdrun"), and free energy estimation '
                           'with "pmx analyse" (Crooks Gaussian Intersection / BAR / Jarzynski) - '
                           'independently for the ligand-protein complex and for the ligand alone in '
                           'solvent, combined via the standard thermodynamic cycle ddG_bind = '
                           'dG_complex - dG_water.')
        return methods
