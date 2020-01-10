from core.make_mover_from_xml import MakePyRosettaMoverFromXML
import os

KIC_FULL = './xml/kic_close.xml'

def _get_flags(loop_residues,flank1,flank2,close_cycles,refine_cycles):
	FLAGS = {
	6:
	'''
	-in:file:fullatom
	-ignore_zero_occupancy false
	-overwrite
	-ignore_unrecognized_res
	-parser:script_vars start={}
	-parser:script_vars end={}
	-parser:script_vars cut={}
	-parser:script_vars flank1={}
	-parser:script_vars flank2={}
	-parser:script_vars r1={}
	-parser:script_vars r2={}
	-parser:script_vars r3={}
	-parser:script_vars r4={}
	-parser:script_vars r5={}
	-parser:script_vars r6={}
	-parser:script_vars close_cycles={}
	-parser:script_vars refine_cycles={}
	'''
	}
	start = loop_residues[0]
	end = loop_residues[-1]
	cut = flank1
	r1,r2,r3,r4,r5,r6 = loop_residues
	return FLAGS[6].format(start,end,cut,flank1,flank2,r1,r2,r3,r4,r5,r6,\
													close_cycles,refine_cycles)

def _resnum_inpose(pose,chain,res):
    info = pose.pdb_info()
    return info.pdb2pose(chain,res)

def _get_residues_to_close(pose,c,res):
    info = pose.pdb_info()
    nres = info.nres()
    curr_res = _resnum_inpose(pose,c,res)

    if curr_res==1:
        print('Encountered N-term residue! Nothing to do')
        return 1,1

    elif curr_res==nres:
        print('Encountered C-term residue! Nothing to do')
        return 1,1

    prev_res = _resnum_inpose(pose,c,res-1)
    if not prev_res:
        # res is N-terminal residue or next residue is absent
        prev_res = None
    next_res = _resnum_inpose(pose,c,res+1)
    if not next_res:
        # res is C-terminal residue or next residue is absent
        next_res = None

    return prev_res, next_res

def _get_loop_residues(pose,c,res):
    info = pose.pdb_info()
    nres = info.nres()
    flank1, flank2 = _get_residues_to_close(pose,c,res)

    if (nres - flank2)>=2 and (flank1 - 1)>=2:
        loop_residues = range(flank1-2,flank2+3)
    elif (nres - flank2)==1 and (flank1 - 1)>=2:
        loop_residues = range(flank1-2,flank2+2)
    elif (nres - flank2)==0 and (flank1 - 1)>=2:
        loop_residues = range(flank1-2,flank2+1)
    elif (nres - flank2)>=2 and (flank1 - 1)==1:
        loop_residues = range(flank1-1,flank2+3)
    elif (nres - flank2)>=2 and (flank1 - 1)==0:
        loop_residues = range(flank1,flank2+3)

    return loop_residues, flank1, flank2


class KICMove(object):
	"""docstring for KIC"""
	def __init__(self, indel, chain, flank1, flank2, loop_residues):
		super(KICMove, self).__init__()
		self.chain = chain
		self.indel = indel
		self.pose = indel.pose
		self.flank1 = flank1
		self.flank2 = flank2
		self.loop_residues = loop_residues

	def _write_flags(self):
		flags = _get_flags(self.loop_residues,self.flank1,self.flank2,
						self.indel.close_cycles,self.indel.refine_cycles)
		f = open('{}/current.flags'.format(self.indel.temp_path),'w')
		f.write(flags)
		f.close()

	def apply(self,id):
		self._write_flags()
		print(os.getcwd())
		kicmover,pose = MakePyRosettaMoverFromXML(
			KIC_FULL,'{}/current.flags'.format(self.indel.temp_path),
			'{}/current.pdb'.format(self.indel.temp_path))
		kicmover.apply(pose)
		pose.dump_pdb('{}/current.pdb'.format(self.indel.temp_path))
		# pose.dump_pdb('./tests/current_{}.pdb'.format(id))
		self.indel.pose.assign(pose)
		self.indel.reload_pdbf('{}/current.pdb'.format(self.indel.temp_path))
		# return pose.scores
