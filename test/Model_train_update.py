import torch
from Model_fir_warp import Model_fir_warp
from update_enhance import PostProcessing

class Model_train_update(torch.nn.Module):
    def __init__(self, is_encode, enhance, ResBlock_num=1):
        super(Model_train_update, self).__init__()

        self.enhance_flag = enhance
        self.is_encode = is_encode
        self.warp = Model_fir_warp()
        self.enhance = PostProcessing()
        
    def forward(self, x_l, flo_l, mask_l, x_r, flo_r, mask_r, grid_up, grid_org, L_frame, alpha_level, beta_level, rank):
        nearest_L, nearest_R, fir_L, fir_R = self.warp(x_l, flo_l, mask_l, x_r, flo_r, mask_r, grid_up, grid_org)
        
        predict1 = (fir_L*0.25 + fir_R*0.25)
        if not self.is_encode:
            
            L_frame_update = L_frame / beta_level - predict1/alpha_level
            shape=L_frame_update.size()
            rank_map = torch.cuda.FloatTensor(shape) * 0.0 + rank
            if self.enhance_flag:
                L_frame_update = self.enhance(L_frame_update, -predict1/alpha_level, L_frame/beta_level, fir_L*0.25/alpha_level, fir_R*0.25/alpha_level, rank_map)
            
        else:
            L_frame_update = (L_frame + predict1/alpha_level) * beta_level

        return L_frame_update

