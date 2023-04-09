import torch
from Model_fir_warp import Model_fir_warp
from predict_enhance import PostProcessing

class Model_train_predict(torch.nn.Module):
    def __init__(self, is_encode, enhance, ResBlock_num=1):
        super(Model_train_predict, self).__init__()

        self.enhance_flag = enhance
        self.is_encode = is_encode
        self.warp = Model_fir_warp()
        self.enhance = PostProcessing()
        
        
    def forward(self, x_l, flo_l, mask_l, x_r, flo_r, mask_r, grid_up, grid_org, refer, alpha_level, rank):
        nearest_L, nearest_R, fir_L, fir_R = self.warp(x_l, flo_l, mask_l, x_r, flo_r, mask_r, grid_up, grid_org)
        
        predict1 = (fir_L + fir_R) / 2.0
        
        if not self.is_encode:
            H = refer / alpha_level + predict1
            shape=H.size()
            rank_map = torch.cuda.FloatTensor(shape) * 0 + rank
            if self.enhance_flag:
                H = self.enhance(H, predict1, refer/alpha_level, fir_L / 2.0, fir_R / 2.0, rank_map)

        else:
            H = refer - predict1
            H = H * alpha_level

        return H

