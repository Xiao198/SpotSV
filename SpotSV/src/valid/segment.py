class Segment():
    # x:read   y:ref
    def __init__(self, x_start, y_start, length, forward, seg_id):
        '''
        生成segment对象
        :param x_start: read start
        :param y_start: ref start
        :param length: match length
        :param forward: True or False
        :param seg_id: segment id
        '''

        self.seg_length = length
        self.seg_forward = forward
        self.seg_id = seg_id
        if seg_id != 0:
            self.x_start = x_start + 1
            self.y_start = y_start + 1
            # print(self.y_start)
        else:
            self.x_start = x_start + 1
            self.y_start = y_start

        self.seg_flag = ""
        if forward:
            # 正向，forward = true
            self.x_end = self.x_start + self.seg_length
            # self.y_end = self.y_start + (length - 1)
        else:
            self.x_end = self.x_start - self.seg_length
            # self.y_end = self.y_start - (length - 1)

        self.y_end = self.y_start + self.seg_length

    def set_flag(self, flag):
        self.seg_flag = flag

    def set_x_end(self, xEnd):
        self.x_end = xEnd
    
    def set_y_end(self, yEnd):
        self.y_end = yEnd
    
    def set_length(self, length):
        self.seg_length = length
    
    def set_x_start(self, xStart):
        self.x_start = xStart
    
    def set_y_start(self, yStart):
        self.y_start = yStart
    
    def set_forward(self, forward):
        self.seg_forward = forward
    
    def set_seg_id(self, segId):
        self.seg_id = segId

    def to_string(self):
        # id - coord_read - coord_ref - forward - lengh -
        # if self.seg_forward:
        coord_x = "[" + str(self.x_start) + "," + str(self.x_end) + "]"
        coord_y = "[" + str(self.y_start) + "," + str(self.y_end) + "]"
        # else:
        #     coord_x = "[" + str(self.x_end) + "," + str(self.x_start) + "]"
        #     coord_y = "[" + str(self.y_end) + "," + str(self.y_start) + "]"
        return str(self.seg_id) + "-" + coord_x + "-" + coord_y + "-" + str(self.seg_forward) + "-" + str(self.seg_length) + "-" + str(self.seg_flag)

    def get_x_coord(self):
        coord_x = "[" + str(self.x_start) + "," + str(self.x_end) + "]"
        return coord_x

    def get_y_coord(self):
        coord_y = "[" + str(self.y_start) + "," + str(self.y_end) + "]"
        return coord_y

    def get_xmid_coord(self):
        xmid_coord = (self.x_start + self.x_end) / 2
        return xmid_coord

    def get_ymid_coord(self):
        ymid_coord = (self.y_start + self.y_end) / 2
        return ymid_coord

    def get_xstart_coord(self):
        xstart_coord = self.x_start
        return xstart_coord

    def get_ystart_coord(self):
        ystart_coord = self.y_start
        return ystart_coord

    def get_xend_coord(self):
        xend_coord = self.x_end
        return xend_coord

    def get_yend_coord(self):
        yend_coord = self.y_end
        return yend_coord

    def get_forward(self):
        seg_forward = self.seg_forward
        return seg_forward

    def get_id(self):
        seg_id = self.seg_id
        return seg_id

    def get_length(self):
        seg_length = self.seg_length
        return seg_length