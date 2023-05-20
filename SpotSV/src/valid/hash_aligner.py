from valid.segment import Segment

class HashAligner():
    def __init__(self, k, windowSize, mismatchNum):
        '''
        生成哈希值
        :param k:k-mer
        :param windowSize:扩展长度
        :param mismatchNum:最大不匹配数量，避免测序错误等，默认应该为零
        '''
        self.k = k
        self.windowSize = windowSize
        self.mismatchNum = mismatchNum
        self.segments = []
        self.selfDiffSegs = []
        self.compareDiffSegs = []
        self.N = 'N'
        self.G = 'G'
        self.A = 'A'
        self.T = 'T'
        self.C = 'C'

    def getSegments(self):
        return self.segments

    def getSelfDiffSegs(self):
        return self.selfDiffSegs

    def run(self, x, y):
        '''
         x indicates read, while y ref
        :param x:
        :param y:
        :return:
        '''
        self.ref_length = len(y.get_bases())
        self.makePairwiseAlignment(x, y)

    def extendKmersForward(self, xBases, yBases, matchPositions, p, i, segId):
        '''
        正向扩展碱基
        :param xBases:
        :param yBases:
        :param matchPositions:
        :param p:
        :param i:
        :param segId:
        :return:
        '''
        matchLength = self.k  # kmer的长度
        mismatch = self.mismatchNum
        break_flag = 0
        while mismatch <= self.mismatchNum:
            # beyond the length of x
            if matchPositions[p] + matchLength >= len(xBases):
                break_flag = 1
                break
            # beyond the length of y
            if i + matchLength >= len(yBases):
                break_flag = 1
                break
            xBase = xBases[matchPositions[p] + matchLength]
            yBase = yBases[i + matchLength]
            # stop extension when meet 'N'
            if xBase == self.N or yBase == self.N:
                break_flag = 1
                break
            # stop extension when different
            if xBase != yBase:
                mismatch += 1
            matchLength += 1
        # when longer than windowSize
        if matchLength >= self.windowSize:
            if break_flag == 0:
                d = Segment(matchPositions[p], i, matchLength - 1, True, segId)
            if break_flag == 1:
                d = Segment(matchPositions[p], i, matchLength, True, segId)
            self.segments.append(d)

    def extendKmersReverse(self, reverseXbases, yBases, matchPosition, i, segId):
        matchLength = self.k
        mismatch = 0
        break_flag = 0
        while mismatch <= self.mismatchNum:
            if matchPosition + matchLength >= len(reverseXbases) - 1:
                break_flag = 1
                break
            if i + matchLength >= len(yBases) - 1:
                break_flag = 1
                break
            xBase = reverseXbases[matchPosition + matchLength]
            yBase = yBases[i + matchLength]

            if xBase == self.N or yBase == self.N:
                break_flag = 1
                break
            # stop extension when different
            if xBase != yBase:
                mismatch += 1
            matchLength += 1

        if matchLength >= self.windowSize:
            if break_flag == 0:
                d = Segment((len(reverseXbases) - 1) - matchPosition, i, matchLength - 1, False, segId)
            if break_flag == 1:
                d = Segment((len(reverseXbases) - 1) - matchPosition, i, matchLength, False, segId)
            self.segments.append(d)

    def calHash(self, kmer):
        hashValue = kmer
        return hashValue

    def makePairwiseAlignment(self, x, y):
        xBases = x.get_bases()
        # 获取碱基
        reverseXbases = x.get_reverse_complement()
        # 碱基取反
        hashedPositions = {}
        # hashPositions = {key:value} = {kmer :pos}
        # for sequence x, make hash
        for i in range(0, len(xBases) - (self.k + 1)):
            kmer = xBases[i: i + self.k]
            hashValue = self.calHash(kmer)  # hashValue = kmer
            if hashValue not in hashedPositions.keys():
                # hashPositions = {key:value} = {kmer :[pos1,pos2...]}
                hashedPositions[hashValue] = []
            hashedPositions[hashValue].append(i)
        # get reverse seq, and make hash table
        for j in range(0, len(reverseXbases) - (self.k + 1)):
            kmer = reverseXbases[j: j + self.k]
            hashValue = self.calHash(kmer)

            if hashValue not in hashedPositions.keys():
                hashedPositions[hashValue] = []
            hashedPositions[hashValue].append(-1 - j)
        yBases = y.get_bases()
        segId = 0
        self.hashvalues = []
        for l in range(0, len(yBases) - (self.k + 1)):
            kmer = yBases[l: l + self.k]
            hashValue = self.calHash(kmer)
            # hashValue = kmer
            self.hashvalues.append(hashValue)
            if hashValue in hashedPositions.keys():
                matchPositions = hashedPositions[hashValue]  # 记录一个kmer匹配位置,[pos1,pos2...]
                # for each possible hit position
                for p in range(0, len(matchPositions)):
                    # hit in the same strand
                    # Kmer extension
                    if matchPositions[p] >= 0:
                        # Already matched in previous Kmer  ???l>0
                        if matchPositions[p] > 0 and l > 0 and xBases[matchPositions[p] - 1] == yBases[l - 1]:
                            continue
                        else:
                            self.extendKmersForward(xBases, yBases, matchPositions, p, l, segId)
                            segId += 1
                    else:
                        matchPosition = -1 - matchPositions[p]
                        if matchPosition > 0 and l > 0 and reverseXbases[matchPosition - 1] == yBases[l - 1]:
                            continue
                        self.extendKmersReverse(reverseXbases, yBases, matchPosition, l, segId)
                        segId += 1

    def getMergeSegments(self):
        return self.segments

    def getHashValues(self):
        return self.hashvalues
