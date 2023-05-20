class Alignment():
    def __init__(self, ref_start, ref_end, read_start, read_end, forward, seq, read_seq, align_id):
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.read_start = read_start
        self.read_end = read_end
        self.forward = forward
        self.seq = seq
        self.read_seq = read_seq
        self.align_id = align_id

    def to_string(self):
        return str(self.align_id) + '-' + str(self.ref_start) + '-' + str(self.ref_end) + '-' + str(self.read_start) + '-' + str(self.read_end) + '-' + str(self.forward) + '-' + self.seq + '-' + self.read_seq