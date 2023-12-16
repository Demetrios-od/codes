# import sys
import random
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import bch

MAX_METRIC = 1000


def sort_llrs(data, p):
	n = len(data)
	idx = list(range(n))
	for i in range(p):
		for j in range(i+1, n):
			if abs(data[idx[j]]) < abs(data[idx[i]]):
				b = idx[i]
				idx[i] = idx[j]
				idx[j] = b
	return idx[:p].copy()


def create_flip_patterns(p):
	flip_patterns = [[],[]]
	for i in range(1 << p):
		w = 0
		ii = i
		while ii > 0:
			w += ii & 1
			ii >>= 1
		flip_patterns[w & 1] += [i]
	return flip_patterns


def chase(code, data, p, flip_patterns):
	synd = code.calc_syndrome(data)
	if all(s == 0 for s in synd):
		return [[[], 0]]
	
	idx = sort_llrs(data, p)
	# e = 1 if synd[2] == gf_pow(synd[1], 3, bch['GField']) else 0
	# pattern_parity = synd[0] ^ e

	cand_list = []
	for pattern_parity in [0,1]:
		for pattern in flip_patterns[pattern_parity]:
			tmp_synd = synd.copy()

			# apply testpattern
			i = 0
			pt = pattern
			while pt > 0:
				if pt & 1 == 1:
					code.update_syndrome(tmp_synd, idx[i])
				pt >>= 1
				i += 1
			
			# decode
			res, pos = code.decode_hard(tmp_synd)
			if not res:
				continue
			
			# merge testpattern and decoder decision
			i = 0
			pt = pattern
			while pt > 0:
				if pt & 1 == 1:
					if idx[i] in pos:
						pos.remove(idx[i])
					else:
						pos += [idx[i]]
				pt >>= 1
				i += 1
			
			# deduplication
			pos.sort()
			if any([cand[0] == pos for cand in cand_list]):
				continue

			# calculate metric and add candidate
			mtr = sum([abs(data[i]) for i in pos])
			if mtr > MAX_METRIC:
				mtr = MAX_METRIC
			cand_list += [[pos, mtr]]
	
	# order candidates by metric
	n = len(cand_list)
	idx = list(range(n))
	for i in range(n):
		for j in range(i+1,n):
			if cand_list[idx[j]][1] < cand_list[idx[i]][1]:
				b = idx[j]
				idx[j] = idx[i]
				idx[i] = b

	return [cand_list[i] for i in idx]


def map_process(data, cand_list, alpha, beta):
	map_list = []   # position differences and metrics relative to the main candidate
	for cand in cand_list:
		pos = cand_list[0][0].copy()
		for x in cand[0]:
			if x in pos:
				pos.remove(x)
			else:
				pos += [x]
		pos.sort()
		mtr = sum([abs(data[i]) for i in pos])
		map_list += [[pos, mtr]]

	n = len(data)
	ls = len(cand_list)
	v = [-1]*n
	for i in range(n):
		cand2_idx = 1
		while cand2_idx < ls and i not in map_list[cand2_idx][0]:
			cand2_idx += 1
		if cand2_idx < ls:
			v[i] = map_list[cand2_idx][1] - abs(data[i])
	idx = [i for i,vi in enumerate(v) if vi >= 0]
	v = [vi/len(idx) if i in idx else beta for i,vi in enumerate(v)]
	v = [alpha*vi*bch.SIGN(di) for di,vi in zip(data, v)]
	for i in cand_list[0][0]:
		v[i] = -v[i]

	return v


def test_hard_decoder(params):
	N = params['N']
	code = bch.ebch(N, params['T'], params['SPC'])
	print('Test hard decoder. BCH(', code.N, ',', code.K, ',', code.D, '), t=', code.T, sep='')

	info = [random.randint(0,1) for i in range(code.K)]
	# info = [0]*code.K
	par = code.encode(info)
	cw_enc = info + par
	print('CW weight:', sum(cw_enc))
	cw = [1-2*i for i in cw_enc]
	synd = code.calc_syndrome(cw)
	if not all([s == 0 for s in synd]):
		print('Clear CW has syndrome', synd)
	miscorr = 0

	for e1 in range(N):
		print('*', end='', flush=True)
		# print(e1)
		cw[e1] *= -1
		synd = code.calc_syndrome(cw)
		res, pos = code.decode_hard(synd)
		if params['T'] == 0 and res:
			miscorr += 1
		if params['T'] > 0 and not (res and pos == [e1]):
			print('Fail: error was [%d], found' % e1, pos)

		if params['T'] > 0:
			for e2 in range(e1+1, N):
				# print(e1, e2)
				cw[e2] *= -1
				synd = code.calc_syndrome(cw)
				res, pos = code.decode_hard(synd)
				if params['T'] == 1 and res:
					miscorr += 1
				if params['T'] > 1 and not (res and pos == [e1, e2]):
					print('Fail: errors were [%d %d], found' % (e1, e2), pos)

				if params['T'] > 1:
					for e3 in range(e2+1, N):
						# print(e1, e2, e3)
						cw[e3] *= -1
						synd = code.calc_syndrome(cw)
						res, pos = code.decode_hard(synd)
						if params['T'] == 2 and res:
							miscorr += 1
						if params['T'] > 2 and not (res and pos == [e1, e2, e3]):
							print('Fail: errors were [%d %d %d], found' % (e1, e2, e3), pos)
					
						if params['T'] > 2:
							for e4 in range(e3+1, N):
								# print(e1, e2, e3, e4)
								cw[e4] *= -1
								synd = code.calc_syndrome(cw)
								res, pos = code.decode_hard(synd)
								if params['T'] == 3 and res:
									miscorr += 1
								if params['T'] > 3 and not (res and pos == [e1, e2, e3, e4]):
									print('Fail: errors were [%d %d %d %d], found' % (e1, e2, e3, e4), pos)
								
								if params['T'] > 3:
									for e5 in range(e4+1, N):
										# print(e1, e2, e3, e4, e5)
										cw[e5] *= -1
										synd = code.calc_syndrome(cw)
										res, pos = code.decode_hard(synd)
										if params['T'] == 4 and res:
											miscorr += 1
										if params['T'] > 4 and not (res and pos == [e1, e2, e3, e4, e5]):
											print('Fail: errors were [%d %d %d %d %d], found' % (e1, e2, e3, e4, e5), pos)

										cw[e5] *= -1
								cw[e4] *= -1
						cw[e3] *= -1
				cw[e2] *= -1
		cw[e1] *= -1
	print('\nMiscorrections:', miscorr, '\n')


def test_cw_hard(params):
	N = params['N']
	code = bch.ebch(N, params['T'], params['SPC'])
	print('Test hard codeword. BCH(', code.N, ',', code.K, ',', code.D, '), t=', code.T, '. Errors: ', params['errpos'], sep='')

	info = [random.randint(0,1) for i in range(code.K)]
	# info = [0]*code.K
	par = code.encode(info)
	cw_enc = info + par
	cw = [1-2*i for i in cw_enc]
	for i in params['errpos']:
		cw[i] *= -1
	synd = code.calc_syndrome(cw)
	res, pos = code.decode_hard(synd)
	print('Errors were:', params['errpos'])
	print('Syndrome:', synd)
	print('Errors found:', pos, '  res =', res)


def test_chase(params, map=False):
	print('SNR = %.2f  p = %d' % (params['snr'], params['p']))
	if map:
		print('MAP decoding')

	code = bch.ebch(params['N'], params['T'], params['SPC'])
	random.seed(params['seed'])
	sigma = 10**(-params['snr']/20)

	flip_patterns = create_flip_patterns(params['p'])
	stat_distr_size = 20
	stat = {
		'frames_tx': 0,
		'bits_tx': 0,

		'in_berr': 0,
		'in_ferr': 0,
		'in_berr_distr': [0]*stat_distr_size,

		'soft_berr': 0,
		'soft_ferr': 0,
		'soft_fails': 0,
		'soft_berr_distr': [0]*stat_distr_size,
		'list_size': [0]*(1 << params['p']),

		'hard_berr': 0,
		'hard_ferr': 0,
		'hard_fails': 0,
		'hard_berr_distr': [0]*stat_distr_size,
	}

	for frame in range(params['Nframes']):
		info = [random.randint(0,1) for i in range(code.K)]
		# info = [1]*code.K
		par = code.encode(info)
		cw_enc = info + par
		cw_s = [random.normalvariate(1-2*i, sigma)  for i in cw_enc]
		cw_hd = [bch.HARD(i) for i in cw_s]
		numerr = sum([x^y for x,y in zip(cw_enc, cw_hd)])
		stat['in_berr'] += numerr
		stat['in_berr_distr'][numerr] += 1
		if numerr > 0:
			stat['in_ferr'] += 1

		# hard decoding
		synd = code.calc_syndrome(cw_s)
		res, pos = code.decode_hard(synd)
		cw_h_out = cw_hd.copy()
		if res:
			for i in pos:
				cw_h_out[i] ^= 1
		else:
			stat['hard_fails'] += 1
		numerr = sum([x^y for x,y in zip(cw_enc, cw_h_out)])
		stat['hard_berr'] += numerr
		stat['hard_berr_distr'][numerr] += 1
		if numerr > 0:
			stat['hard_ferr'] += 1

		# soft decoding
		cand_list = chase(code, cw_s, params['p'], flip_patterns)
		stat['list_size'][len(cand_list)] += 1
		cw_s_out = cw_hd.copy()
		if len(cand_list) == 0:
			stat['soft_fails'] += 1
		else:
			if map:
				cw_s_upd = map_process(cw_s, cand_list, params['alpha'], params['beta'])
				soft_out = [x+y for x,y in zip(cw_s, cw_s_upd)]
				cw_s_out = [bch.HARD(x) for x in soft_out]
			else:
				for i in cand_list[0][0]:
					cw_s_out[i] ^= 1
		
		numerr = sum([x^y for x,y in zip(cw_enc, cw_s_out)])
		stat['soft_berr'] += numerr
		stat['soft_berr_distr'][numerr] += 1
		if numerr > 0:
			stat['soft_ferr'] += 1
		
		stat['frames_tx'] += 1
		stat['bits_tx'] += params['N']

	print('\n         BER        FER      B.ers   F.ers   Fails   Err.distr')
	print('IN    %.3e  %.3e %7d %7d       -   ' % (stat['in_berr']/stat['bits_tx'], stat['in_ferr']/stat['frames_tx'], stat['in_berr'], stat['in_ferr']), end='')
	print(stat['in_berr_distr'])
	print('HARD  %.3e  %.3e %7d %7d %7d   ' % (stat['hard_berr']/stat['bits_tx'], stat['hard_ferr']/stat['frames_tx'], stat['hard_berr'], stat['hard_ferr'], stat['hard_fails']), end='')
	print(stat['hard_berr_distr'])
	print('SOFT  %.3e  %.3e %7d %7d %7d   ' % (stat['soft_berr']/stat['bits_tx'], stat['soft_ferr']/stat['frames_tx'], stat['soft_berr'], stat['soft_ferr'], stat['soft_fails']), end='')
	print(stat['soft_berr_distr'])
	print('List size distr:', stat['list_size'], '\n')


def performance_chase(params, map=False):
	random.seed(params['seed'])
	code = bch.ebch(params['N'], params['T'], params['SPC'])
	flip_patterns = create_flip_patterns(params['p'])

	stat = []
	stat_elem = {
		'snr': 0.0,
		'frames_tx': 0,
		'bits_tx': 0,

		'in_berr': 0,
		'in_ferr': 0,
		'soft_berr': 0,
		'soft_ferr': 0,
		'soft_fails': 0,
	}
	
	print('SNR     inBER       BER      inFER        FER      Frames    Be.in    B.ers    Fe.in    F.ers    Fails')
	def print_simrow(snr, stat_part):
		print('%.2f  %.3e  %.3e  %.3e  %.3e  %7d  %7d  %7d  %7d  %7d  %7d' % (
			snr,
			stat_part['in_berr']/stat_part['bits_tx'],
			stat_part['soft_berr']/stat_part['bits_tx'],
			stat_part['in_ferr']/stat_part['frames_tx'],
			stat_part['soft_ferr']/stat_part['frames_tx'],
			stat_part['frames_tx'],
			stat_part['in_berr'],
			stat_part['soft_berr'],
			stat_part['in_ferr'],
			stat_part['soft_ferr'],
			stat_part['soft_fails'],
		), end='\r')
	
	snr = params['snr_min']
	while snr <= params['snr_max']:
		sigma = 10**(-snr/20)
		stat += [stat_elem.copy()]
		stat_part = stat[-1]
		stat_part['snr'] = snr
		while stat_part['frames_tx'] < params['min_nframes'] or stat_part['soft_berr'] < params['max_berr']:
			info = [random.randint(0,1) for i in range(code.K)]
			par = code.encode(info)
			cw_enc = info + par
			cw_s = [random.normalvariate(1-2*i, sigma)  for i in cw_enc]
			cw_hd = [bch.HARD(i) for i in cw_s]
			numerr = sum([x[0]^x[1] for x in zip(cw_enc, cw_hd)])
			stat_part['in_berr'] += numerr
			if numerr > 0:
				stat_part['in_ferr'] += 1
			
			# soft decoding
			# first candidate is assumed as the hard Chase output
			cand_list = chase(code, cw_s, params['p'], flip_patterns)
			cw_s_out = cw_hd.copy()
			if len(cand_list) > 0:
				for i in cand_list[0][0]:
					cw_s_out[i] ^= 1
			else:
				stat_part['soft_fails'] += 1
			
			# error count
			numerr = sum([x[0]^x[1] for x in zip(cw_enc, cw_s_out)])
			stat_part['soft_berr'] += numerr
			if numerr > 0:
				stat_part['soft_ferr'] += 1
			
			stat_part['frames_tx'] += 1
			stat_part['bits_tx'] += params['N']
			if stat_part['frames_tx'] % params['print_interval'] == 0:
				print_simrow(snr, stat_part)
		print_simrow(snr, stat_part)
		print()
		snr += params['snr_step']
	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	snrs = [s['snr'] for s in stat]
	bers = [s['soft_berr']/s['bits_tx'] for s in stat]
	ax.plot(snrs, bers, label='p%d' % params['p'])
	ax.set(
		title = 'BCH performance',
		yscale = 'log',
		xlabel = 'SNR, dB',
		ylabel = 'BER',
		xlim = [snrs[0], snrs[-1]],
		ylim = [1e-7, 0.01],
	)
	ax.grid(
		which = 'major',
		color = [0.7, 0.7, 0.7],
		linestyle = '-',
	)
	ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
	ax.yaxis.set_major_locator(ticker.FixedLocator([10**i for i in range(-2,-16,-1)]))
	ax.legend()
	fig.set_figwidth(10)
	fig.set_figheight(9)
	plt.show()
		


if __name__ == '__main__':

	test_hard_decoder({
		'N'     : 32,
		'T'     : 4,
		'SPC'   : 1,
	})

	# test_cw_hard({
	# 	'N'     : 63,
	# 	'T'     : 4,
	# 	'SPC'   : 0,
	# 	'errpos': [0,1,2,26],
	# })


	# test_chase({
	# 	'N'  : 256,
	# 	'T'  : 2,
	# 	'SPC': 1,
	# 	'p'  : 4,

	# 	'seed'   : 100,
	# 	'snr'    : 7.5,
	# 	'Nframes': 10000,

	# 	'alpha': 3.5,
	# 	'beta' : 1,
	# }, map=True)


	# performance_chase({
	# 	'N'  : 256,
	# 	'T'  : 2,
	# 	'SPC': 1,
	# 	'p'  : 4,

	# 	'seed'    : 100,
	# 	'snr_min' : 7.0,
	# 	'snr_max' : 8.5,
	# 	'snr_step': 0.5,

	# 	'min_nframes'   : 20000,
	# 	'max_berr'      : 1000,
	# 	'print_interval': 1000,
	# }, map=False)
	
