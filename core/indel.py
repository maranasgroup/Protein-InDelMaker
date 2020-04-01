from pyrosetta import *
init(extra_options = "-ex1 -ex2 -ex2aro mute all")

from pyrosetta.rosetta.protocols.relax import relax_pose
from rosetta.protocols import minimization_packing as pack_min

grafting = pyrosetta.rosetta.protocols.grafting
relax = pyrosetta.rosetta.protocols.relax.FastRelax()
relax.set_scorefxn(get_fa_scorefxn())
from core.files import InputFile, PDBFile
from core.objects import OrderedDict
from core.make_mover_from_xml import MakePyRosettaMoverFromXML
from core.kic import KICMove
from pyrosetta import Pose, Vector1, pose_from_file, create_score_function
from pyrosetta.rosetta.core.select import residue_selector as selections
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task import operation

def redock(pose,partners):
    dock_jump = 1 # jump number 1 is the inter-body jump
    pyrosetta.rosetta.protocols.docking.setup_foldtree(pose,partners,Vector1([dock_jump]))

    # Create ScoreFunctions for centroid and fullatom docking
    scorefxn = pyrosetta.create_score_function("ref2015.wts")

    # Setup the high resolution (fullatom) docking protocol using DockMCMProtocol.
    docking = pyrosetta.rosetta.protocols.docking.DockMCMProtocol()
    # Many of its options and settings can be set using the setter methods.
    docking.set_scorefxn(scorefxn)

    # # Set the native pose so that the output scorefile contains the pose rmsd metric
    # native_pose = pose 

    # # Optional: setup a PyMOLObserver
    # # pyrosetta.rosetta.protocols.moves.AddPyMOLObserver(test_pose, True)

    # # Perform protein-ligand docking
    # # counter = 0 # for pretty output to PyMOLObserver

    # test_pose = pose.clone() # Reset test pose to original structure
        
    # counter += 1 # Change the pose name, for pretty output to PyMOLObserver
    # test_pose.pdb_info().name(job_output + '_' + str(counter))
        
    # Perform docking and output to PyMOL:
    print('DOCKING'*100)
    docking.apply(pose) 

class Error(Exception):
   """Base class for other exceptions"""
   pass

class ResidueAbsent(Error):
   """Raised when the residue to delete or insert or mutate is absent"""
   pass

class NextResidueAbsent(Error):
   pass

def get_pose_number(pose,chain,res):
    info = pose.pdb_info()
    return info.pdb2pose(chain,int(res))

def get_pose(filepath): return pose_from_pdb(filepath)

def get_pose_with_ligand(filepath,LIGAND_PARAMS=[]):
    pose = Pose()
    ligand_params = Vector1(LIGAND_PARAMS)
    res_set = pose.conformation().modifiable_residue_type_set_for_conf()
    res_set.read_files_for_base_residue_types( ligand_params )
    pose.conformation().reset_residue_type_set_for_conf( res_set )
    return pose_from_file(filepath
        )

def determine_Cterm(chain,res): return False

def determine_Nterm(chain,res): return False

def get_loop_residues(flank1,flank2):
    return [flank1-2,flank1-1,flank1,flank2,flank2+1,flank2+2]

class Mover(object):
    """docstring for Move"""
    def __init__(self,indel):
        super(Mover, self).__init__()
        self.type = ''
        self.chain = ''
        self.move = None
        self.indel = indel

    def load_from_input_line(self,line):
        head = line[0]
        chain = line[1]
        tasks = line[2]
        self.type = head
        self.chain = chain

        if head=='DELETE':
            dels = []
            for each in tasks:
                if '-' in each:
                    a,b = each.split('-')
                    dels.append(list(range(int(a),int(b)+1)))
                else:
                    dels.append([int(each)])
            move = DelMove(self,self.type,self.chain,dels)
            self.move = move

        elif head=='INSERT':
            inserts = []
            # print(tasks)
            for each in tasks:
                res,aas = each.split('_')
                aas = list(aas)
                res = int(res)
                temp = []
                for i,aa in enumerate(aas):
                    # print(i,aa)
                    temp.append((res,aa))

                inserts.append(temp)
            self.move = InsMove(self,self.type,self.chain,inserts)

        elif head=='MUTATE':
            mutations = []
            for each in tasks:
                res,aa = each.split('_')
                res = int(res)
                mutations.append((res,aa))

            if self.indel.ligand:
                self.move = MutMoveLigand(self,self.type,self.chain,mutations)
            else:
                self.move = MutMove(self,self.type,self.chain,mutations)


class DelMove(object):
    """docstring for DelMove"""
    def __init__(self,mover,type,chain,dels):
        super(DelMove, self).__init__()
        self.type = type
        self.chain = chain
        self.deletions = dels
        self.mover = mover

    def apply(self,just_mutate=False):
        DEBUG = self.mover.indel.DEBUG
        pdbf = self.mover.indel.get_file('pdb')
        numbering = self.mover.indel.numbering

        for each_set in self.deletions:
            if len(each_set)>=1:
                multi = True

            for each in each_set:
                cterm = determine_Cterm(self.chain,each)
                nterm = determine_Nterm(self.chain,each)
                my_res = numbering.get_pose_number(self.chain,each)
                next_res = numbering.get_pose_number(self.chain,each+1)

                if cterm:
                    pdbf.delete_res_list([each],self.chain)

                elif nterm:
                    pdbf.delete_res_list([each],self.chain)

                elif not my_res:
                    raise ResidueAbsent

                elif not next_res:
                    raise NextResidueAbsent

                else:
                    flank1 = my_res - 1
                    flank2 = flank1+1
                    loop_residues = get_loop_residues(flank1,flank2)
                    if DEBUG: print('*'*100)
                    if DEBUG: print('Deleting pose_residue:{} pdb_residue:{}'.format(my_res,each))
                    grafting.delete_region(self.mover.indel.pose,my_res,my_res)
                    self.mover.indel.pose.dump_pdb('{}/current.pdb'.format(self.mover.indel.temp_path))
                    self.mover.indel.reload_pdbf('{}/current.pdb'.format(self.mover.indel.temp_path))
                    kic_move = KICMove(self.mover.indel,self.chain,flank1,flank2,loop_residues)
                    if DEBUG: print('Running KIC for Insertion')
                    kic_move.apply(each)
                    self.mover.indel.pose.dump_pdb('{}/current.pdb'.format(self.mover.indel.temp_path))
                    self.mover.indel.reload_pdbf('{}/current.pdb'.format(self.mover.indel.temp_path))
                    # self.mover.indel._apply_pymol()
                    if DEBUG: print(self.mover.indel.numbering.maps)
                    self.mover.indel.numbering._reload(deletion=(self.chain,my_res))
                    if DEBUG: print(self.mover.indel.numbering.maps)

        return self.mover.indel.pose.scores['total_score']

class InsMove(object):
    """docstring for InsMove"""
    def __init__(self,mover,type,chain,inserts):
        super(InsMove, self).__init__()
        self.type = type
        self.chain = chain
        self.inserts = inserts
        self.mover = mover

    def apply(self,just_mutate=False):
        DEBUG = self.mover.indel.DEBUG
        def _insert(my_res):
            if DEBUG: print('Trying to insert {} after {}'.format(aa,my_res))
            ins_pose = pose_from_sequence(aa)
            if DEBUG: print('New pose made for {}'.format(aa))
            pose = grafting.insert_pose_into_pose(self.mover.indel.pose,ins_pose,my_res,my_res+1)
            if DEBUG: print('Grafting successful!')
            
            grafting.idealize_combined_pose(pose,movemap,my_res,my_res,my_res+2,1,nres-2,True)
            # WHY nres-2 ?
            
            if DEBUG: print('Idealize successful!')
            relax.apply(pose,get_fa_scorefxn())
            return pose

        pdbf = self.mover.indel.get_file('pdb')
        numbering = self.mover.indel.numbering

        pdbf = self.mover.indel.get_file('pdb')
        nres = self.mover.indel.numbering.nres[self.chain]
        # nres = self.mover.indel.pose.pdb_info().nres()
        movemap = MoveMap()
        movemap.set_bb(True)
        movemap.set_chi(True)

        for each_set in self.inserts:
            for each in each_set:
                res,aa = each
                if DEBUG: print('Insertion: {} {}'.format(res,aa))
                # cterm = determine_Cterm(self.chain,res)
                # nterm = determine_Nterm(self.chain,res)
                cterm = False
                nterm = False
                my_res = numbering.get_pose_number(self.chain,res)
                next_res = numbering.get_pose_number(self.chain,res+1)

                if cterm:
                    self.mover.indel.pose.assign(_delete())

                elif nterm:
                    self.mover.indel.pose.assign(_delete())

                elif not my_res:
                    raise ResidueAbsent

                elif not next_res:
                    raise NextResidueAbsent

                else:
                    pose = _insert(my_res)
                    pose.dump_pdb('{}/current.pdb'.format(self.mover.indel.temp_path))
                    self.mover.indel.pose.assign(pose)

                    flank1 = my_res+1
                    flank2 = flank1+1
                    loop_residues = get_loop_residues(flank1,flank2)
                    self.mover.indel.reload_pdbf('{}/current.pdb'.format(self.mover.indel.temp_path))
                    kic_move = KICMove(self.mover.indel,self.chain,flank1,flank2,loop_residues)
                    if DEBUG: print('Running KIC for Insertion')
                    kic_move.apply(my_res)
                    self.mover.indel.pose.dump_pdb('{}/current.pdb'.format(self.mover.indel.temp_path))
                    self.mover.indel.reload_pdbf('{}/current.pdb'.format(self.mover.indel.temp_path))
                    # self.mover.indel._apply_pymol()
                    if DEBUG: print(self.mover.indel.numbering.maps)
                    self.mover.indel.numbering._reload(insertion=(self.chain,res))
                    if DEBUG: print(self.mover.indel.numbering.maps)
        return self.mover.indel.pose.scores['total_score']

class MutMoveLigand(object):
    """docstring for MutMoveLigand"""
    def __init__(self,mover,type,chain,mutations):
        super(MutMoveLigand, self).__init__()
        self.type = type
        self.chain = chain
        self.mutations = mutations
        self.mover = mover


    def apply(self,just_mutate=True):
        for num,aa in self.mutations:
            number = get_pose_number(self.mover.indel.pose,self.chain,num)
            pyrosetta.toolbox.mutate_residue(self.mover.indel.pose,number,aa)
            
            tf = TaskFactory()
            tf.push_back(operation.InitializeFromCommandline())
            tf.push_back(operation.RestrictToRepacking())
            packer = pack_min.PackRotamersMover()
            packer.task_factory(tf)

            if just_mutate: 
                packer.apply(self.mover.indel.pose)
                return self.mover.indel.pose.scores['total_score']

        packer.apply(self.mover.indel.pose)

        #relax.apply(self.mover.indel.pose)
        redock(self.mover.indel.pose,self.mover.indel.partners)
        self.mover.indel.pose.dump_pdb('{}/current.pdb'.format(self.mover.indel.temp_path))
        self.mover.indel.reload_pdbf('{}/current.pdb'.format(self.mover.indel.temp_path))
        # self.mover.indel._apply_pymol()
        return self.mover.indel.pose.scores['total_score']

class MutMove(object):
    """docstring for MutMove"""
    def __init__(self,mover,type,chain,mutations):
        super(MutMove, self).__init__()
        self.type = type
        self.chain = chain
        self.mutations = mutations
        self.mover = mover

    def apply(self,just_mutate=False):
        for num,aa in self.mutations:
            number = get_pose_number(self.mover.indel.pose,self.chain,num)
            pyrosetta.toolbox.mutate_residue(self.mover.indel.pose,number,aa)
            if just_mutate: return self.mover.indel.pose.scores['total_score']

        tf = TaskFactory()
        tf.push_back(operation.InitializeFromCommandline())
        tf.push_back(operation.RestrictToRepacking())
        packer = pack_min.PackRotamersMover()
        packer.task_factory(tf)

        packer.apply(self.mover.indel.pose)

        #relax.apply(self.mover.indel.pose)
        self.mover.indel.pose.dump_pdb('{}/current.pdb'.format(self.mover.indel.temp_path))
        self.mover.indel.reload_pdbf('{}/current.pdb'.format(self.mover.indel.temp_path))
        # self.mover.indel._apply_pymol()
        return self.mover.indel.pose.scores['total_score']

class Numbering(object):
    """docstring for Numbering."""
    def __init__(self, indel):
        super(Numbering, self).__init__()
        self.indel = indel
        self.maps = OrderedDict()
        self._load()

    def _load(self):
        pose = self.indel.pose
        pose.dump_pdb('{}/current.pdb'.format(self.indel.temp_path))
        self.indel.reload_pdbf('{}/current.pdb'.format(self.indel.temp_path))
        pdb = self.indel.get_file('pdb').pdb
        info = pose.pdb_info()
        all_residues = pdb.get_residues()

        nres = OrderedDict()
        for c in pdb.chains:
            nres_c = 0
            map = OrderedDict()
            c_residues = all_residues[c]
            for resnum,residue in c_residues.items():
                nres_c+=1
                resnum = int(resnum)
                map[resnum] = info.pdb2pose(c,resnum)
            self.maps[c] = map
            nres[c] = nres_c
        self.nres = nres

    def _reload(self,insertion=('A',0),deletion=('A',0)):
        DEBUG = self.indel.DEBUG
        
        chain_i = insertion[0]
        res_i = int(insertion[1])

        chain_d = deletion[0]
        res_d = int(deletion[1])

        if res_i:
            nres = self.indel.numbering.nres[chain_i]
            nres+=1
            self.indel.numbering.nres[chain_i] = nres
            if DEBUG: print('Changing due to insertion at {}'.format(res_i))
            new = OrderedDict()
            i = 0
            for pdb_resnum,pose_resnum in self.maps[chain_i].items():
                i+=1
                if pdb_resnum>res_i:
                    new[pdb_resnum] = 1 + pose_resnum
                    if DEBUG: print('Loop {} : {} Changed as {}'.format(i,pdb_resnum,pose_resnum+1))
                else:
                    new[pdb_resnum] = pose_resnum
                    if DEBUG: print(' {} not changed'.format(pdb_resnum))
            self.maps[chain_i] = new

        if res_d:
            nres = self.indel.numbering.nres[chain_d]
            nres-=1
            self.indel.numbering.nres[chain_d] = nres
            if DEBUG: print('Changing due to deletion at {}'.format(res_d))
            i = 0
            new = OrderedDict()
            for pdb_resnum,pose_resnum in self.maps[chain_d].items():
                i+=1
                if DEBUG: print(pdb_resnum,pdb_resnum)
                if pdb_resnum>res_d:
                    new[pdb_resnum] = -1 + pose_resnum
                    if DEBUG: print('Loop {} : {} Changed as {}'.format(i,pdb_resnum,pose_resnum-1))
                else:
                    new[pdb_resnum] = pose_resnum
                    if DEBUG: print(' {} not changed'.format(pdb_resnum))
            self.maps[chain_d] = new

    def get_pose_number(self, chain, pdb_resnum):
        return self.maps[chain][pdb_resnum]

class InDelMut(object):
    """docstring for InDel"""
    def __init__(self, temp_path,results_path,pdb_path,input_path,\
        close_cycles = 20, refine_cycles = 1,relax_cycles=3,
        DEBUG=True,watch = True,ligand=False,ligand_params=[],partners='A_X'):

        self.DEBUG = DEBUG

        super(InDelMut, self).__init__()
        if temp_path.strip()[-1]=='/':
            temp_path = temp_path.strip()[:-1]
        if results_path.strip()[-1]=='/':
            results_path = results_path.strip()[:-1]
        self.temp_path = temp_path
        self.results_path = results_path
        if DEBUG: print('temp_path:',self.temp_path)
        if DEBUG: print('results_path:',self.results_path)
        self.pdb_path = pdb_path
        self.input_path = input_path
        self.refine_cycles = refine_cycles
        self.close_cycles = close_cycles
        self.relax_cycles = relax_cycles
        self.ligand = ligand
        self.partners = partners

        self._load_files()

        # constant pose - initial pose
        if ligand:
            self.POSE = get_pose_with_ligand(self.pdb_path,ligand_params)
        else:
            self.POSE = get_pose(self.pdb_path)


        pose = Pose()
        # keeps current pose
        self.pose = pose.assign(self.POSE)

        movemap = MoveMap()
        movemap.set_bb(True)
        movemap.set_chi(True)
        minmover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
        minmover.movemap(movemap)
        minmover.score_function(get_fa_scorefxn())
        minmover.apply(self.pose)

        self._load_movers()
        self.scores[0] = pose.scores['total_score']

        self._load_numbering()

        # pymol = PyMOLMover()
        # pymol.keep_history(True)
        # self.pymol = pymol

        # if watch: self._apply_pymol()

        self.mover_index = 0

    def _load_numbering(self):
        self.numbering = Numbering(self)
        DEBUG = self.DEBUG
        if DEBUG: print('** Finished loading numbering')

    def run(self):
        DEBUG = self.DEBUG

        i = 0
        for mover in self.movers:
            move = mover.move
            i+=1
            if DEBUG: print('** Running MOVER {}'.format(i))
            score = move.apply()

            if DEBUG: print('** Finished MOVER {}'.format(i))
            self.pose.dump_pdb('{}/result_{}.pdb'.format(self.results_path,self.mover_index+1))
            self.scores[self.mover_index+1] = score
            self.status[self.mover_index+1] = True
            self.mover_index+=1

            if self.ligand:
                if DEBUG and ligand: print('** Running DOCKING {}'.format(i))
                redock(self.pose,self.partners)

            for i in range(self.relax_cycles):
                relax.apply(self.pose)

        # input_lines = open(self.input_path).readlines()

        f = open('{}/scores.txt'.format(self.results_path),'w')
        for index,score in self.scores.items():
            f.write('{}\t{}\n'.format(index,score))
        f.close()

    def make_next_move(self):
        move = self.movers[self.mover_index].move
        score = move.apply()
        self.scores[self.mover_index+1] = score
        self.status[self.mover_index+1] = True
        self.mover_index+=1

    def _apply_pymol(self): self.pymol.apply(self.pose)

    def reload_pdbf(self,newpdb):
        self.files['pdb'] = PDBFile(newpdb)

    def _load_files(self):
        DEBUG = self.DEBUG

        pdbf = PDBFile(self.pdb_path)
        inputf = InputFile(self.input_path)

        self.files = {
        'input' : inputf,
        'pdb' : pdbf
        }

        if DEBUG: print('** Finished loading pdb and input')
        return None

    def _load_movers(self):
        lines = self.get_file('input').lines

        movers_all = []
        status = OrderedDict()
        scores = OrderedDict()

        # load native score and status
        status[0] = True
        scores[0] = self.pose.scores['total_score']

        for line in lines:
            mover = Mover(self)
            mover.load_from_input_line(line)
            movers_all.append(mover)

        self.movers = movers_all
        self.status = status
        self.scores = scores

        DEBUG = self.DEBUG
        if DEBUG: print('** Finished loading movers')


    def get_file(self,name):
        if name in ['input', 'pdb']:
            return self.files[name]
        else:
            return None
