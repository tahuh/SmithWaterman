#!/usr/bin/python

"""
SmithWaterman.py

Using affine gap penalty as described in https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm

How to use

from SmithWaterman import SWAligner

#using default parameters for alignment
# Not specifying sequences to align at the start of alignment
aligner = SWAlinger()
aligner.set_reference(some sequence)
aligner.set_query(some query)
alignment = aligner.align()

#Specify sequences at initial
"""

import sys
import time

class SWAlignmentResult:
	def __init__(self,p=None,r=None,q=None,c=None):
		self.map_pos = p
		self.ref = r
		self.query = q
		self.cigar = c # Cigar string usually used in SAM file
		pass
	def get_mapping_position(self):
		return self.map_pos
	def get_ref_seq(self):
		return self.ref
	def get_query_seq(self):
		return self.query
	def get_cigar_string(self):
		return self.cigar
		
class SWAligner:
	def __init__(self, ref=None, query=None, gap_open=2, gap_ext=3, verbose=False):
		self.seq1 = ref
		self.seq2 = query
		self.gap_open = gap_open
		self.gap_ext = gap_ext
		self.verbose = verbose
	def set_reference(self, seq):
		self.seq1 = seq
	def set_query(self, seq):
		self.seq2 = seq
	def align(self):
		if (self.seq1 == None) or (self.seq2 == None):
			raise AttributeError("reference and query sequences must be specified")
		self._build_matrix()
		self._fill_matrix()
		alignment = self._backtrace()
		return alignment
	def _build_matrix(self):
		if self.verbose:
			sys.stderr.write("[SWAligner] Initialize score matrix\n")
			st = time.time()
		self.matrix = []
		length1 = len(self.seq1) # reference sequence length
		length2 = len(self.seq2) # query sequence length
		matrix_row_size = length1 + 1
		matrix_column_size = length2 + 1
		for _ in range(matrix_row_size):
			new_cols = [0] * matrix_column_size
			self.matrix.append(new_cols)
		if self.verbose:
			ed = time.time()
			sys.stderr.write("[SWAligner] Done initialization. Elapsed %.3lf sec\n"%(ed-st))
	def _fill_matrix(self):
		### Assumes that matrix is already initialized
		if self.verbose:
			st = time.time()
			sys.stderr.write("[SWAligner] Filling score matrix before alignment\n")
		colsize = len(self.matrix[0])
		rowsize = len(self.matrix)
		for i in range(1, rowsize):
			for j in range(1,colsize):
				ref_base = self.seq1[i-1]
				query_base = self.seq2[j-1]
				s = self.compute_similarity(ref_base, query_base)
				column_way_weights = []
				row_way_weights = []
				for k in range(1, i):
					wk = self.compute_gap_score(k)
					h = self.matrix[i-k][j]
					column_way_weights.append(h - wk)
				for l in range(1,j):
					wl = self.compute_gap_score(l)
					h = self.matrix[i][j-l]
					row_way_weights.append(h - wl)
				if len(column_way_weights) == 0:
					max_col = 0
				else:
					max_col = max(column_way_weights)
				if len(row_way_weights) == 0 :
					max_row = 0
				else:
					max_row = max(row_way_weights)
				L = max_row
				U = max_col
				D = self.matrix[i-1][j-1] + s
				current_maxima = self.findMax(L,U,D)
				self.matrix[i][j] = current_maxima
		if self.verbose:
			ed = time.time()
			sys.stderr.write("[SWAligner] Done filling score matrix before alignment. Elapsed %.3f sec\n"%(ed-st))
			
	def _backtrace(self):
		if self.verbose:
			st = time.time()
			sys.stderr.write("[SWAligner] Backtrace to find alignment\n")
		
		# find maximum score over alignment matrix
		max_i = -1
		max_j = -1
		max_val = -1
		rows = len(self.matrix)
		cols = len(self.matrix[0])
		
		for i in range(1,rows):
			for j in range(1,cols):
				mat_val = self.matrix[i][j]
				if mat_val > max_val:
					max_val = mat_val
					max_i = i
					max_j = j
		
		### Now we start backtrace from given value
		current_value = -1
		alignment_end_pos = max_j ### Where the query ends at this alignment
		alignment_start_pos = 0
		trace_result = []
		while current_value != 0:
			left = self.matrix[max_i][max_j-1]
			up = self.matrix[max_i-1][max_j]
			diag = self.matrix[max_i-1][max_j-1]
			max_char = self.findMaxDirection(left,up,diag)
			if max_char == -1:
				### Moving left. meaning query sequence has insertion
				max_j -= 1
				trace_result.append("<")
			elif max_char == 0:
				max_j -= 1
				max_i -= 1
				trace_result.append(".")
			else:
				### Moving up. meaning query sequence has deletion
				max_i -= 1
				trace_result.append("^")
			current_value = self.findMax(left, up, diag)
		### Now the alignment ends
		alignment_start_pos = max_i + 1
		alignment = SWAlignmentResult()
		alignment.map_pos = alignment_start_pos
		alignment.ref = self.seq1
		alignment.query = self.seq2
		
		### Generate cigar string
		numbers = []
		alphabets = []
		for letter in trace_result:
			if letter == '<':
				if len(alphabets) == 0:
					alphabets.append('I')
					numbers.append(1)
				else:
					if alphabets[-1] == 'I':
						numbers[-1] = numbers[-1] + 1
					else:
						alphabets.append('I')
						numbers.append(1)
			elif letter == '.':
				if len(alphabets) == 0 :
					alphabets.append('M')
					numbers.append(1)
				else:
					if alphabets[-1] == 'M':
						numbers[-1] = numbers[-1] + 1
					else:
						alphabets.append('M')
						numbers.append(1)
			else:
				if len(alphabets) == 0:
					alphabets.append('D')
					numbers.append(1)
				else:
					if alphabets[-1] == 'D':
						numbers[-1] = numbers[-1] + 1
					else:
						alphabets.append('D')
						numbers.append(1)
		if alignment_start_pos != 1:
			numbers.insert(0,alignment_start_pos)
			alphabets.insert(0,'H')
		if alignment_end_pos != len(self.seq2):
			diff = len(self.seq2) - alignment_end_pos
			numbers.append(diff)
			alphabets.append('H')
		
		zips = zip(numbers , alphabets)
		cig = ''
		for z in zips:
			here = str(z[0]) + z[1]
			cig += here
		alignment.cigar = cig
		
		if self.verbose:
			ed = time.time()
			sys.stderr.write("[SWAligner] Done backtrace to find alignment. Elapsed %.3f sec\n"%(ed-st))
		return alignment
	def compute_similarity(self, a,b):
		if a == b:
			return 1
		else:
			return -1
			
	def compute_gap_score(self,k):
		return self.gap_open * (k-1) + self.gap_ext
	def findMax(self, L,U, D):
		return max([L,U,D])
	def findMaxDirection(self,L,U,D):
		l = [ L,U,D]
		maxima = self.findMax(L,U,D)
		if maxima == D:
			return 0
		idx = l.index(maxima)
		if idx == 0 :
			## To the left
			return -1
		elif idx == 1:
			### To the diagonal
			return 0
		else:
			### To the upper direction
			return 1
			
			
			
#### Test case
def test():
	ref = 'TGTTACGG'
	query = 'GGTTGACTA'
	gap_e = 1
	gap_o = 1
	aligner = SWAligner(ref=ref,query=query,gap_ext=gap_e, gap_open=gap_o,verbose=True)
	alignment = aligner.align()
	sys.stdout.write("R:%s\nQ:%s\n"%(aligner.seq1,aligner.seq2))
	sys.stdout.write(alignment.cigar + "\n")
	
if __name__ == "__main__":
	print("Run test")
	test()