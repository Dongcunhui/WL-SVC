import torch
from torch.autograd import Variable
import copy
import torch.nn as nn
import filter_kernel

class Model_fir_warp(torch.nn.Module):
    def __init__(self):
        super(Model_fir_warp, self).__init__()
        self.filter = filter_kernel.FIR_kernel(trainable_set=False)
        self.pixelshuffle = nn.PixelShuffle(4)
        self.upsample = torch.nn.UpsamplingNearest2d(scale_factor=4)

    def warp(self, x_l, flo_l, x_r, flo_r, grid):
        B, C, H, W = x_l.size()

        vgrid_l = grid - flo_l

        # scale grid to [-1,1]
        vgrid_l[:, 0, :, :] = 2.0 * vgrid_l[:, 0, :, :].clone() / max(W - 1, 1) - 1.0
        vgrid_l[:, 1, :, :] = 2.0 * vgrid_l[:, 1, :, :].clone() / max(H - 1, 1) - 1.0

        vgrid_l = vgrid_l.permute(0, 2, 3, 1)
        # output_l = torch.nn.functional.grid_sample(x_l, vgrid_l, mode="bilinear", padding_mode="border")
        # output_l = torch.nn.functional.grid_sample(x_l, vgrid_l, mode="nearest", padding_mode="border")
        output_l = torch.nn.functional.grid_sample(x_l, vgrid_l, mode="nearest", padding_mode="border", align_corners=True)
        # output_l = torch.nn.functional.grid_sample(x_l, vgrid_l, padding_mode="border", align_corners=True)

        vgrid_r = grid - flo_r

        # scale grid to [-1,1]
        vgrid_r[:, 0, :, :] = 2.0 * vgrid_r[:, 0, :, :].clone() / max(W - 1, 1) - 1.0
        vgrid_r[:, 1, :, :] = 2.0 * vgrid_r[:, 1, :, :].clone() / max(H - 1, 1) - 1.0

        vgrid_r = vgrid_r.permute(0, 2, 3, 1)
        # output_r = torch.nn.functional.grid_sample(x_r, vgrid_r, mode="bilinear", padding_mode="border")
        # output_r = torch.nn.functional.grid_sample(x_r, vgrid_r, mode="nearest", padding_mode="border")
        output_r = torch.nn.functional.grid_sample(x_r, vgrid_r, mode="nearest", padding_mode="border", align_corners=True)
        # output_r = torch.nn.functional.grid_sample(x_r, vgrid_r, padding_mode="border", align_corners=True)

        return output_l, output_r

    def forward(self, x_l, flo_l, mask_l, x_r, flo_r, mask_r, grid_up, grid_org):

        x_l_filter = self.filter(x_l)
        x_r_filter = self.filter(x_r)

        x_L_interpolate = self.pixelshuffle(x_l_filter)
        x_R_interpolate = self.pixelshuffle(x_r_filter)

        after_L_warp, after_R_warp = self.warp(x_L_interpolate, self.upsample(flo_l*4),
                                               x_R_interpolate, self.upsample(flo_r*4), grid_up)

        after_L_warp = after_L_warp[:,:,::4,::4] * mask_l
        after_R_warp = after_R_warp[:,:,::4,::4] * mask_r

        after_L_whole = (1 - mask_l) * after_R_warp + mask_l * after_L_warp
        after_R_whole = (1 - mask_r) * after_L_warp + mask_r * after_R_warp


        before_L, before_R = self.warp(x_l, flo_l, x_r, flo_r, grid_org)
        before_L = (1 - mask_l) * before_R + mask_l * before_L
        before_R = (1 - mask_r) * before_L + mask_r * before_R

        return before_L, before_R, after_L_whole, after_R_whole

