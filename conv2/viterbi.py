import sys
import numpy as np

class convcodec():
    G = [0xb, 0xf]

    def __init__(self, *, rate_denom=len(G), dec_depth=None):
        assert 2 <= rate_denom <= len(self.G), f'Code rate 1/{rate_denom} is not supported. Must be from 1/2 to 1/{len(self.G)}'
        self.RATE = rate_denom
        self.G = self.G[:rate_denom]

        m = max(self.G)
        self.REG_L = 0   # length of the shift register
        while m > 0:
            self.REG_L += 1
            m >>= 1
        self.REG_MASK = (1 << self.REG_L) - 1
        self.DEC_DEPTH = dec_depth   # never used

        # for encoder
        self.buf = 0

        # for decoder

        # Transition table
        # Index = state after
        # States before = [index>>1, (index>>1) | (1<<(L-2))]
        # item = transmitted bits for each state_before (first bit is highest)
        self.trs = []
        for i in range(1 << (self.REG_L-1)):
            item = [0, 0]
            for j in range(2):
                state_before = (i >> 1) | (j << (self.REG_L-2))
                for Gx in self.G:
                    item[j] <<= 1
                    item[j] |= (((state_before<<1) | (i&1)) & Gx).bit_count() &1
            self.trs.append(item)
        # the table may be reduced twice

    def reset(self):
        self.buf = 0

    def encode(self, data_in):
        data_out = [0] * len(data_in) * self.RATE
        for i,b in enumerate(data_in):
            self.buf = ((self.buf << 1) | b) & self.REG_MASK   # add new bit, shift the buffer and fit it to REG_L bits
            for p in range(self.RATE):
                data_out[i*self.RATE + p] = (self.buf & self.G[p]).bit_count() & 0x1
        return data_out

    def decode(self, data_in, *, tailzeros=False):
        assert len(data_in) % self.RATE == 0, 'Number of input bits is inconsistent with code rate'
        N = len(self.trs)   # number of states (trellis nodes)
        mtr  = [0] * N      # node metrics
        path = [0] * N      # paths to each node

        inp_ind = 0
        while inp_ind < len(data_in):
            rec = 0
            for i in range(self.RATE):
                rec <<= 1
                rec |= data_in[inp_ind] # &1
                inp_ind += 1
            # now rec contains RATE bits from input

            _mtr  = [0] * N   # for new metrics and paths
            _path = [0] * N
            # one step over trellis
            for i in range(N):   # for each node
                s1 = i >> 1   # state_before[0]
                m1 = mtr[s1] + (rec ^ self.trs[i][0]).bit_count()
                s2 = s1 | (1 << (self.REG_L-2))   # state_before[1]
                m2 = mtr[s2] + (rec ^ self.trs[i][1]).bit_count()
                if m1 <= m2:
                    _mtr[i] = m1
                    _path[i] = ((path[s1]<<1) | (i&1))
                else:
                    _mtr[i] = m2
                    _path[i] = ((path[s2]<<1) | (i&1))
            mtr = _mtr
            path = _path

        ind = 0
        if not tailzeros:
            minmtr = mtr[0]
            for i in range(1, len(mtr)):
                if mtr[i] < minmtr:
                    ind = i
                    minmtr = mtr[i]
                    
        best_path = path[ind]
        data_out = [0] * (len(data_in) // self.RATE)
        for i in range(len(data_out)):
            data_out[-i-1] = best_path & 1
            best_path >>= 1
            
        return data_out


def perf_conv(*, rate=2, snr_min=-1.0, snr_max=8.0):
    snr = snr_min
    snr_step = 0.2
    cwlen = 100
    min_frames = 1000
    max_frames = 100000
    frame_errs_max = 100
    print_step = 1000

    codec = convcodec(rate_denom=rate)

    print('  SNR    InBER     frames   b.ers   f.ers       BER         FER')

    def print_report(snr, frames, in_errs, out_errs, frame_errs, *, end='\n'):
        inber = float(in_errs) / (frames * cwlen * codec.RATE)
        outber = float(out_errs) / (frames * cwlen)
        fer = float(frame_errs) / frames
        print(f'{snr:5.2f}  {inber:8.3e} {frames:8d} {out_errs:7d} {frame_errs:7d}   {outber:8.3e}   {fer:8.3e}', end=end)
        sys.stdout.flush()
    
    while snr < snr_max + snr_step / 2:
        sigma = 10**(-snr/20)
        in_errs = 0
        out_errs = 0
        frame_errs = 0
        frame = 0
        while frame < min_frames or (frame_errs < frame_errs_max and frame < max_frames):
            data = np.concatenate((np.random.randint(0, 2, cwlen), np.zeros(codec.REG_L)))
            data_enc = np.array(codec.encode(data.astype(int)))
            noise = np.random.normal(0, sigma, len(data_enc))
            data_rec = -2*data_enc + 1 + noise
            data_rec = (data_rec < 0).astype(int)
            in_errs += sum(data_rec != data_enc)
            data_dec = codec.decode(data_rec, tailzeros=True)
            errs = sum(data.astype(int) != np.array(data_dec))
            out_errs += errs
            if errs > 0: frame_errs += 1
            frame += 1
            if frame % print_step == 0:
                print_report(snr, frame, in_errs, out_errs, frame_errs, end='\r')
                
        print_report(snr, frame, in_errs, out_errs, frame_errs)
        snr += snr_step



if __name__ == '__main__':
    perf_conv()
