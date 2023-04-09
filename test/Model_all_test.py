import torch
import numpy as np
from torch.nn import functional as F

from Model_train_predict import Model_train_predict
from Model_train_update import Model_train_update
from Model_space_entropy import Model_space_entropy

alpha_level = np.sqrt(23.0/32.0)
beta_level = np.sqrt(3.0 / 2.0)

prequant_weight_high = [1.00000000000000, 0.92470128353329, 0.99029103222357, 1.12001265233249, 1.28641274725767,
  1.48343269251745, 1.71234593056217, 1.97708086092642]
prequant_weight_low = [ 1.00000000000000, 1.10554159678513, 1.26197963240006, 1.45296631451356, 1.67651411818611,
  1.93551742673391, 2.23484071726601, 2.58054224785086]
copycomp_weight_high = [ 1.00000000000000 / alpha_level,
  (prequant_weight_low[0] / prequant_weight_high[1]) / alpha_level,
  (prequant_weight_low[1] / prequant_weight_high[2]) / alpha_level,
  (prequant_weight_low[2] / prequant_weight_high[3]) / alpha_level,
  (prequant_weight_low[3] / prequant_weight_high[4]) / alpha_level,
  (prequant_weight_low[4] / prequant_weight_high[5]) / alpha_level,
  (prequant_weight_low[5] / prequant_weight_high[6]) / alpha_level,
  (prequant_weight_low[6] / prequant_weight_high[7]) / alpha_level]
copycomp_weight_low = [
  1.00000000000000 / beta_level,
  (prequant_weight_low[0] / prequant_weight_low[1]) / beta_level,
  (prequant_weight_low[1] / prequant_weight_low[2]) / beta_level,
  (prequant_weight_low[2] / prequant_weight_low[3]) / beta_level,
  (prequant_weight_low[3] / prequant_weight_low[4]) / beta_level,
  (prequant_weight_low[4] / prequant_weight_low[5]) / beta_level,
  (prequant_weight_low[5] / prequant_weight_low[6]) / beta_level,
  (prequant_weight_low[6] / prequant_weight_low[7]) / beta_level]

class Model_all(torch.nn.Module):
    def __init__(self, bitplane_number, train_step):
        super(Model_all, self).__init__()
        enhance = 1 # flag for learnable temporal wavelet inverse transform

        # the first level temporal wavelet
        self.preditct_wavelet_encode_H1 = Model_train_predict(1, enhance)
        self.space_entropy_H1 = Model_space_entropy(bitplane_number, train_step)
        self.preditct_wavelet_decode_H1 = Model_train_predict(0, enhance)
        self.update_wavelet_encode_H1 = Model_train_update(1, enhance)
        self.update_wavelet_decode_H1 = Model_train_update(0, enhance)

        # the second level temporal wavelet
        self.preditct_wavelet_encode_H2 = Model_train_predict(1, enhance)
        self.space_entropy_H2 = Model_space_entropy(bitplane_number, train_step)
        self.preditct_wavelet_decode_H2 = Model_train_predict(0, enhance)
        self.update_wavelet_encode_H2 = Model_train_update(1, enhance)
        self.update_wavelet_decode_H2 = Model_train_update(0, enhance)


        # the third level temporal wavelet
        self.preditct_wavelet_encode_H3 = Model_train_predict(1, enhance)
        self.space_entropy_H3 = Model_space_entropy(bitplane_number, train_step)
        self.preditct_wavelet_decode_H3 = Model_train_predict(0, enhance)
        self.update_wavelet_encode_H3 = Model_train_update(1, enhance)
        self.update_wavelet_decode_H3 = Model_train_update(0, enhance)

        self.space_entropy_L = Model_space_entropy(bitplane_number, train_step)

    def forward(self, frame, mask, mv, grid_up, grid_org, gop_number, bitplane_encode_number, discard_number):
        Batch_size = frame.size()[0]
        if Batch_size != 1:
            print('Batch_size must is 1!')
            exit()
        frame = frame.unsqueeze(2)
        mask = mask.unsqueeze(2)
        mv_H1 = []
        mv_H2 = []
        mv_H3 = []
        up_mv_H1 = []
        up_mv_H2 = []
        up_mv_H3 = []
        mask_H1 = []
        mask_H2 = []
        mask_H3 = []
        up_mask_H1 = []
        up_mask_H2 = []
        up_mask_H3 = []
        rank = 0

        # temporal wavelet forward transform
        # first level
        H1_list = []
        for i in range(7):
            H_tmp = self.preditct_wavelet_encode_H1(frame[:, i * 2, :, :].cuda(),
                                                                           mv[:, 4 * i:4 * i + 2, :, :].cuda(),
                                                                           mask[:, i * 2, :, :].cuda(),
                                                                           frame[:, i * 2 + 2, :, :].cuda(),
                                                                           mv[:, 4 * i + 2:4 * i + 4, :, :].cuda(),
                                                                           mask[:, i * 2 + 1, :, :].cuda(), grid_up,
                                                                           grid_org, frame[:, i * 2 + 1, :, :].cuda(), alpha_level, rank)
            H_tmp = H_tmp.cpu()
            mv_H1.append(mv[:, 4 * i:4 * i + 2, :, :])
            mv_H1.append(mv[:, 4 * i + 2:4 * i + 4, :, :])
            mask_H1.append(mask[:, i * 2, :, :])
            mask_H1.append(mask[:, i * 2 + 1, :, :])
            if torch.sum(mask[:, i * 2, :, :]) ==0 and torch.sum( mask[:, i * 2 + 1, :, :]) == 0:
                H_tmp = copycomp_weight_high[0] * frame[:, i * 2 + 1, :, :] * alpha_level
            H1_list.append(H_tmp)

        L1_list = []
        for i in range(7):
            L1 = self.update_wavelet_encode_H1(H1_list[max(i - 1, 0)].cuda(),
                                                  -mv[:, 28 + 4 * i:28 + 4 * i + 2, :, :].cuda(),
                                                  mask[:, 14 + i * 2, :, :].cuda(),
                                                  H1_list[i].cuda(),
                                                  -mv[:, 28 + 4 * i + 2:28 + 4 * i + 4, :, :].cuda(),
                                                  mask[:, 14 + i * 2 + 1, :, :].cuda(), grid_up,
                                                  grid_org, frame[:, 2 * i, :, :].cuda(), alpha_level, beta_level, rank)
            L1 = L1.cpu()
            up_mv_H1.append(mv[:, 28 + 4 * i:28 + 4 * i + 2, :, :])
            up_mv_H1.append(mv[:, 28 + 4 * i + 2:28 + 4 * i + 4, :, :])
            up_mask_H1.append(mask[:, 14 + i * 2, :, :])
            up_mask_H1.append(mask[:, 14 + i * 2 + 1, :, :])
            if torch.sum(mask[:, 14 + i * 2, :, :]) ==0 and torch.sum(mask[:, 14 + i * 2 + 1, :, :]) == 0:
                L1 = copycomp_weight_low[0]*frame[:, 2 * i, :, :] * beta_level
            L1_list.append(L1)

        # second level
        H2_list = []
        for i in range(3):
            H_tmp = self.preditct_wavelet_encode_H2(L1_list[2 * i].cuda(),
                                                   mv[:, 56 + 4 * i:56 + 4 * i + 2, :,:].cuda(),
                                                   mask[:, 28 + i * 2, :, :].cuda(),
                                                   L1_list[2 * i + 2].cuda(),
                                                   mv[:, 56 + 4 * i + 2:56 + 4 * i + 4,:, :].cuda(),
                                                   mask[:, 28 + i * 2 + 1, :, :].cuda(), grid_up,grid_org
                                                    , L1_list[2 * i + 1].cuda(), alpha_level, rank)
            H_tmp = H_tmp.cpu()
            mv_H2.append(mv[:, 56 + 4 * i:56 + 4 * i + 2, :,:])
            mv_H2.append(mv[:, 56 + 4 * i + 2:56 + 4 * i + 4,:, :])
            mask_H2.append(mask[:, 28 + i * 2, :, :])
            mask_H2.append(mask[:, 28 + i * 2 + 1, :, :])
            if torch.sum(mask[:, 28 + i * 2, :, :]) ==0 and torch.sum( mask[:, 28 + i * 2 + 1, :, :]) == 0:
                H_tmp = copycomp_weight_high[1] * L1_list[2 * i + 1] * alpha_level
            H2_list.append(H_tmp)

        L2_list = []
        for i in range(3):
            L2 = self.update_wavelet_encode_H2(H2_list[max(i - 1, 0)].cuda(),
                                              -mv[:, 68 + 4 * i:68 + 4 * i + 2, :, :].cuda(),
                                              mask[:, 34 + i * 2, :, :].cuda(),
                                              H2_list[i].cuda(),
                                              -mv[:, 68 + 4 * i + 2:68 + 4 * i + 4, :, :].cuda(),
                                              mask[:, 34 + i * 2 + 1, :, :].cuda(), grid_up,
                                              grid_org, L1_list[2 * i].cuda(), alpha_level, beta_level, rank)
            L2 = L2.cpu()
            up_mv_H2.append(mv[:, 68 + 4 * i:68 + 4 * i + 2, :, :])
            up_mv_H2.append(mv[:, 68 + 4 * i + 2:68 + 4 * i + 4, :, :])
            up_mask_H2.append(mask[:, 34 + i * 2, :, :])
            up_mask_H2.append(mask[:, 34 + i * 2 + 1, :, :])
            if torch.sum(mask[:, 34 + i * 2, :, :]) ==0 and torch.sum(mask[:, 34 + i * 2 + 1, :, :]) == 0:
                L2 = copycomp_weight_low[1]*L1_list[2 * i] * beta_level
            L2_list.append(L2)

        # third level
        H3_list = []
        i = 0
        H_tmp = self.preditct_wavelet_encode_H3(L2_list[2 * i].cuda(),
                                               mv[:, 80 + 4 * i:80 + 4 * i + 2, :, :].cuda(),
                                               mask[:, 40 + 2 * i, :, :].cuda(),
                                               L2_list[2 * i + 2].cuda(),
                                               mv[:, 80 + 4 * i + 2:80 + 4 * i + 4, :, :].cuda(),
                                               mask[:, 40 + 2 * i + 1, :, :].cuda(), grid_up,grid_org
                                                , L2_list[2 * i + 1].cuda(), alpha_level, rank)
        H_tmp = H_tmp.cpu()
        mv_H3.append(mv[:, 80 + 4 * i:80 + 4 * i + 2, :, :])
        mv_H3.append(mv[:, 80 + 4 * i + 2:80 + 4 * i + 4, :, :])
        mask_H3.append(mask[:, 40 + 2 * i, :, :])
        mask_H3.append(mask[:, 40 + 2 * i + 1, :, :])
        if torch.sum(mask[:, 40 + 2 * i, :, :]) == 0 and torch.sum(mask[:, 40 + 2 * i + 1, :, :]) == 0:
            H_tmp = copycomp_weight_high[2] * L2_list[2 * i + 1] * alpha_level
        H3_list.append(H_tmp)


        L3_list = []
        i = 0
        L3 = self.update_wavelet_encode_H3(H3_list[0].cuda(),
                                          -mv[:, 84 + 4 * i:84+ 4 * i + 2, :, :].cuda(),
                                          mask[:, 42 + i * 2, :, :].cuda(),
                                          H3_list[0].cuda(),
                                          -mv[:, 84 + 4 * i + 2:84 + 4 * i + 4, :, :].cuda(),
                                          mask[:, 42 + i * 2 + 1, :, :].cuda(), grid_up,grid_org,
                                           L2_list[0].cuda(), alpha_level, beta_level, rank)
        L3 = L3.cpu()
        up_mv_H3.append(mv[:, 84 + 4 * i:84+ 4 * i + 2, :, :])
        up_mv_H3.append(mv[:, 84 + 4 * i + 2:84 + 4 * i + 4, :, :])
        up_mask_H3.append(mask[:, 42 + i * 2, :, :])
        up_mask_H3.append(mask[:, 42 + i * 2 + 1, :, :])
        if torch.sum(mask[:, 42 + i * 2, :, :]) == 0 and torch.sum(mask[:, 42 + i * 2 + 1, :, :]) == 0:
            L3 = copycomp_weight_low[2] * L2_list[0] * beta_level
        L3_list.append(L3)

        for gop_order in range(1, gop_number):
            # print(gop_order)
            # first level
            # H1_list = []
            for i in range(4):
                H_tmp = self.preditct_wavelet_encode_H1(frame[:, 6+8*gop_order+i * 2, :, :].cuda(),
                                       mv[:,32+gop_order*56+ 4 * i:32+gop_order*56+4 * i + 2, :, :].cuda(),
                                       mask[:, 16+gop_order*28+i * 2, :, :].cuda(),
                                       frame[:, 6+8*gop_order+i * 2 + 2, :, :].cuda(),
                                       mv[:, 32+gop_order*56+4 * i + 2:32+gop_order*56+4 * i + 4, :, :].cuda(),
                                       mask[:, 16+gop_order*28+i * 2 + 1, :, :].cuda(), grid_up, grid_org
                                    , frame[:, 6+8*gop_order+i * 2 + 1, :, :].cuda(), alpha_level, rank)
                H_tmp = H_tmp.cpu()
                mv_H1.append(mv[:,32+gop_order*56+ 4 * i:32+gop_order*56+4 * i + 2, :, :])
                mv_H1.append(mv[:, 32+gop_order*56+4 * i + 2:32+gop_order*56+4 * i + 4, :, :])
                mask_H1.append(mask[:, 16+gop_order*28+i * 2, :, :])
                mask_H1.append(mask[:, 16+gop_order*28+i * 2 + 1, :, :])
                if torch.sum(mask[:, 16+gop_order*28+i * 2, :, :]) == 0 and torch.sum(mask[:, 16+gop_order*28+i * 2 + 1, :, :]) == 0:
                    H_tmp = copycomp_weight_high[0] * frame[:, 6+8*gop_order+i * 2 + 1, :, :] * alpha_level
                H1_list.append(H_tmp)

            # L1_list = []
            for i in range(4):
                L1 = self.update_wavelet_encode_H1(H1_list[3+gop_order*4+i-1].cuda(),
                              -mv[:, 32+gop_order*56+16 + 4 * i:32+gop_order*56+16 + 4 * i + 2, :, :].cuda(),
                              mask[:, 16+gop_order*28+8 + i * 2, :, :].cuda(),
                              H1_list[3+gop_order*4+i].cuda(),
                              -mv[:, 32+gop_order*56+16 + 4 * i + 2:32+gop_order*56+16 + 4 * i + 4, :, :].cuda(),
                              mask[:, 16+gop_order*28+8 + i * 2 + 1, :, :].cuda(), grid_up, grid_org,
                            frame[:, 6+8*gop_order+i * 2, :, :].cuda(), alpha_level, beta_level, rank)
                L1 = L1.cpu()
                up_mv_H1.append(mv[:, 32+gop_order*56+16 + 4 * i:32+gop_order*56+16 + 4 * i + 2, :, :])
                up_mv_H1.append(mv[:, 32+gop_order*56+16 + 4 * i + 2:32+gop_order*56+16 + 4 * i + 4, :, :])
                up_mask_H1.append(mask[:, 16+gop_order*28+8 + i * 2, :, :])
                up_mask_H1.append(mask[:, 16+gop_order*28+8 + i * 2 + 1, :, :])
                if torch.sum(mask[:, 16+gop_order*28+8 + i * 2, :, :]) == 0 \
                        and torch.sum(mask[:, 16+gop_order*28+8 + i * 2 + 1, :, :]) == 0:
                    L1 = copycomp_weight_low[0] * frame[:, 6+8*gop_order+i * 2, :, :] * beta_level
                L1_list.append(L1)

            # second level
            # H2_list = []
            for i in range(2):
                H_tmp = self.preditct_wavelet_encode_H2(L1_list[3+gop_order*4+2 * i-1].cuda(),
                               mv[:, 32+gop_order*56+32 + 4 * i:32+gop_order*56+32 + 4 * i + 2, :,:].cuda(),
                               mask[:, 16+gop_order*28+16 + i * 2, :, :].cuda(),
                               L1_list[3+gop_order*4+2 * i + 2-1].cuda(),
                               mv[:, 32+gop_order*56+32 + 4 * i + 2:32+gop_order*56+32 + 4 * i + 4, :, :].cuda(),
                               mask[:, 16+gop_order*28+16 + i * 2 + 1, :, :].cuda(), grid_up, grid_org
                                                        , L1_list[3+gop_order*4+2 * i].cuda(), alpha_level, rank)
                H_tmp = H_tmp.cpu()
                mv_H2.append(mv[:, 32+gop_order*56+32 + 4 * i:32+gop_order*56+32 + 4 * i + 2, :,:])
                mv_H2.append(mv[:, 32+gop_order*56+32 + 4 * i + 2:32+gop_order*56+32 + 4 * i + 4, :, :])
                mask_H2.append(mask[:, 16+gop_order*28+16 + i * 2, :, :])
                mask_H2.append(mask[:, 16+gop_order*28+16 + i * 2 + 1, :, :])
                if torch.sum(mask[:, 16+gop_order*28+16 + i * 2, :, :]) == 0 and torch.sum(mask[:, 16+gop_order*28+16 + i * 2 + 1, :, :]) == 0:
                    H_tmp = copycomp_weight_high[1] * L1_list[3+gop_order*4+2 * i] * alpha_level
                H2_list.append(H_tmp)

            # L2_list = []
            for i in range(2):
                L2 = self.update_wavelet_encode_H2(H2_list[1+gop_order*2+i - 1].cuda(),
                              -mv[:, 32+gop_order*56+40 + 4 * i:32+gop_order*56+40 + 4 * i + 2, :, :].cuda(),
                              mask[:, 16+gop_order*28+20 + i * 2, :, :].cuda(),
                              H2_list[1+gop_order*2+i].cuda(),
                              -mv[:, 32+gop_order*56+40 + 4 * i + 2:32+gop_order*56+40 + 4 * i + 4, :, :].cuda(),
                              mask[:, 16+gop_order*28+20 + i * 2 + 1, :, :].cuda(), grid_up,
                              grid_org, L1_list[3+gop_order*4+2 * i-1].cuda(), alpha_level, beta_level, rank)
                L2 = L2.cpu()
                up_mv_H2.append(mv[:, 32+gop_order*56+40 + 4 * i:32+gop_order*56+40 + 4 * i + 2, :, :])
                up_mv_H2.append(mv[:, 32+gop_order*56+40 + 4 * i + 2:32+gop_order*56+40 + 4 * i + 4, :, :])
                up_mask_H2.append(mask[:, 16+gop_order*28+20 + i * 2, :, :])
                up_mask_H2.append(mask[:, 16+gop_order*28+20 + i * 2 + 1, :, :])
                if torch.sum(mask[:, 16+gop_order*28+20 + i * 2, :, :]) == 0 \
                        and torch.sum(mask[:, 16+gop_order*28+20 + i * 2 + 1, :, :]) == 0:
                    L2 = copycomp_weight_low[1] * L1_list[3+gop_order*4+2 * i-1] * beta_level
                L2_list.append(L2)

            # third level
            # H3_list = []
            i = 0
            H_tmp = self.preditct_wavelet_encode_H3(L2_list[1+gop_order*2+2 * i-1].cuda(),
                                       mv[:, 32+gop_order*56+48 + 4 * i:32+gop_order*56+48 + 4 * i + 2, :, :].cuda(),
                                       mask[:, 16+gop_order*28+24 + 2 * i, :, :].cuda(),
                                       L2_list[1+gop_order*2+2 * i + 2-1].cuda(),
                                       mv[:, 32+gop_order*56+48 + 4 * i + 2:32+gop_order*56+48 + 4 * i + 4, :, :].cuda(),
                                       mask[:, 16+gop_order*28+24 + 2 * i + 1, :, :].cuda(), grid_up,
                                       grid_org, L2_list[1+gop_order*2+2 * i].cuda(), alpha_level, rank)
            H_tmp = H_tmp.cpu()
            mv_H3.append(mv[:, 32+gop_order*56+48 + 4 * i:32+gop_order*56+48 + 4 * i + 2, :, :])
            mv_H3.append(mv[:, 32+gop_order*56+48 + 4 * i + 2:32+gop_order*56+48 + 4 * i + 4, :, :])
            mask_H3.append(mask[:, 16+gop_order*28+24 + 2 * i, :, :])
            mask_H3.append(mask[:, 16+gop_order*28+24 + 2 * i + 1, :, :])
            if torch.sum(mask[:, 16+gop_order*28+24 + 2 * i, :, :]) == 0 and torch.sum(
                    mask[:, 16+gop_order*28+24 + 2 * i + 1, :, :]) == 0:
                H_tmp = copycomp_weight_high[2] * L2_list[1+gop_order*2+2 * i] * alpha_level
            H3_list.append(H_tmp)
            # L3_list = []
            i = 0
            L3 = self.update_wavelet_encode_H3(H3_list[gop_order-1].cuda(),
                                      -mv[:, 32+gop_order*56+52 + 4 * i:32+gop_order*56+52 + 4 * i + 2, :, :].cuda(),
                                      mask[:, 16+gop_order*28+26 + i * 2, :, :].cuda(),
                                      H3_list[gop_order].cuda(),
                                      -mv[:, 32+gop_order*56+52 + 4 * i + 2:32+gop_order*56+52 + 4 * i + 4, :, :].cuda(),
                                      mask[:, 16+gop_order*28+26 + i * 2 + 1, :, :].cuda(), grid_up, grid_org
                                               , L2_list[gop_order*2].cuda(), alpha_level, beta_level, rank)
            L3 = L3.cpu()
            up_mv_H3.append(mv[:, 32+gop_order*56+52 + 4 * i:32+gop_order*56+52 + 4 * i + 2, :, :])
            up_mv_H3.append(mv[:, 32+gop_order*56+52 + 4 * i + 2:32+gop_order*56+52 + 4 * i + 4, :, :])
            up_mask_H3.append(mask[:, 16+gop_order*28+26 + i * 2, :, :])
            up_mask_H3.append(mask[:, 16+gop_order*28+26 + i * 2 + 1, :, :])
            if torch.sum(mask[:, 16+gop_order*28+26 + i * 2, :, :]) == 0 \
                    and torch.sum(mask[:, 16+gop_order*28+26 + i * 2 + 1, :, :]) == 0:
                L3 = copycomp_weight_low[2] * L2_list[gop_order*2] * beta_level
            L3_list.append(L3)


        # temporal wavelet frames coding
        if bitplane_encode_number == 1:
            part_bitplane = 0
            H1 = torch.cat(H1_list[:], 0)
            rec_H1 = []
            bits_H1 = []
            for i in range(H1.size()[0]):
                tmp_rec_H1, tmp_bits_H1, wavelet_scale_H1 = self.space_entropy_H1(H1[i].unsqueeze(0).cuda(), part_bitplane, discard_number)  # 重建时域第二帧
                rec_H1.append(tmp_rec_H1.cpu())
                bits_H1.append(tmp_bits_H1.cpu())

            H2 = torch.cat(H2_list[:], 0)
            rec_H2 = []
            bits_H2 = []
            for i in range(H2.size()[0]):
                tmp_rec_H2, tmp_bits_H2, wavelet_scale_H2 = self.space_entropy_H2(H2[i].unsqueeze(0).cuda(), part_bitplane, discard_number)  # 重建时域第二帧
                rec_H2.append(tmp_rec_H2.cpu())
                bits_H2.append(tmp_bits_H2.cpu())

            H3 = torch.cat(H3_list[:], 0)
            rec_H3 = []
            bits_H3 = []
            for i in range(H3.size()[0]):
                tmp_rec_H3, tmp_bits_H3, wavelet_scale_H3 = self.space_entropy_H3(H3[i].unsqueeze(0).cuda(), part_bitplane, discard_number)  # 重建时域第二帧
                rec_H3.append(tmp_rec_H3.cpu())
                bits_H3.append(tmp_bits_H3.cpu())

            L3 = torch.cat(L3_list[:], 0)
            rec_L3 = []
            bits_L3 = []
            for i in range(L3.size()[0]):
                tmp_rec_L3, tmp_bits_L3, wavelet_scale_L3 = self.space_entropy_L(L3[i].unsqueeze(0).cuda(), part_bitplane, discard_number)  # 重建时域第二帧
                rec_L3.append(tmp_rec_L3.cpu())
                bits_L3.append(tmp_bits_L3.cpu())
        else:
            part_bitplane = 0
            H1 = torch.cat(H1_list[:], 0)
            rec_H1 = []
            bits_H1 = []
            for i in range(H1.size()[0]):
                tmp_rec_H1, tmp_bits_H1, wavelet_scale_H1 = self.space_entropy_H1(H1[i].unsqueeze(0).cuda(),
                                                                                  part_bitplane, 0)  # 重建时域第二帧
                rec_H1.append(tmp_rec_H1.cpu())
                bits_H1.append(tmp_bits_H1.cpu())

            H2 = torch.cat(H2_list[:], 0)
            rec_H2 = []
            bits_H2 = []
            for i in range(H2.size()[0]):
                tmp_rec_H2, tmp_bits_H2, wavelet_scale_H2 = self.space_entropy_H2(H2[i].unsqueeze(0).cuda(),
                                                                                  part_bitplane, 0)  # 重建时域第二帧
                rec_H2.append(tmp_rec_H2.cpu())
                bits_H2.append(tmp_bits_H2.cpu())

            H3 = torch.cat(H3_list[:], 0)
            rec_H3 = []
            bits_H3 = []
            for i in range(H3.size()[0]):
                tmp_rec_H3, tmp_bits_H3, wavelet_scale_H3 = self.space_entropy_H3(H3[i].unsqueeze(0).cuda(),
                                                                                  part_bitplane, 0)  # 重建时域第二帧
                rec_H3.append(tmp_rec_H3.cpu())
                bits_H3.append(tmp_bits_H3.cpu())

            L3 = torch.cat(L3_list[:], 0)
            rec_L3 = []
            bits_L3 = []
            for i in range(L3.size()[0]):
                tmp_rec_L3, tmp_bits_L3, wavelet_scale_L3 = self.space_entropy_L(L3[i].unsqueeze(0).cuda(),
                                                                                 part_bitplane, 0)  # 重建时域第二帧
                rec_L3.append(tmp_rec_L3.cpu())
                bits_L3.append(tmp_bits_L3.cpu())

        for part_bitplane in range(1, bitplane_encode_number-1):
            for i in range(H1.size()[0]):
                tmp_rec_H1, tmp_bits_H1, wavelet_scale_H1 = self.space_entropy_H1(H1[i].unsqueeze(0).cuda(), part_bitplane, 0, rec_H1[i].cuda())
                rec_H1[i] = tmp_rec_H1.cpu()
                bits_H1[i] += tmp_bits_H1.cpu()

            for i in range(H2.size()[0]):
                tmp_rec_H2, tmp_bits_H2, wavelet_scale_H2 = self.space_entropy_H2(H2[i].unsqueeze(0).cuda(), part_bitplane, 0, rec_H2[i].cuda())
                rec_H2[i] = tmp_rec_H2.cpu()
                bits_H2[i] += tmp_bits_H2.cpu()

            for i in range(H3.size()[0]):
                tmp_rec_H3, tmp_bits_H3, wavelet_scale_H3 = self.space_entropy_H3(H3[i].unsqueeze(0).cuda(), part_bitplane, 0, rec_H3[i].cuda())
                rec_H3[i] = tmp_rec_H3.cpu()
                bits_H3[i] += tmp_bits_H3.cpu()

            for i in range(L3.size()[0]):
                tmp_rec_L3, tmp_bits_L3, wavelet_scale_L3 = self.space_entropy_L(L3[i].unsqueeze(0).cuda(), part_bitplane, 0, rec_L3[i].cuda())
                rec_L3[i] = tmp_rec_L3.cpu()
                bits_L3[i] += tmp_bits_L3.cpu()
        if bitplane_encode_number != 1:
            part_bitplane = bitplane_encode_number-1
            for i in range(H1.size()[0]):
                tmp_rec_H1, tmp_bits_H1, wavelet_scale_H1 = self.space_entropy_H1(H1[i].unsqueeze(0).cuda(), part_bitplane,
                                                                                  discard_number, rec_H1[i].cuda())
                rec_H1[i] = tmp_rec_H1.cpu()
                bits_H1[i] += tmp_bits_H1.cpu()

            for i in range(H2.size()[0]):
                tmp_rec_H2, tmp_bits_H2, wavelet_scale_H2 = self.space_entropy_H2(H2[i].unsqueeze(0).cuda(), part_bitplane,
                                                                                  discard_number, rec_H2[i].cuda())
                rec_H2[i] = tmp_rec_H2.cpu()
                bits_H2[i] += tmp_bits_H2.cpu()

            for i in range(H3.size()[0]):
                tmp_rec_H3, tmp_bits_H3, wavelet_scale_H3 = self.space_entropy_H3(H3[i].unsqueeze(0).cuda(), part_bitplane,
                                                                                  discard_number, rec_H3[i].cuda())
                rec_H3[i] = tmp_rec_H3.cpu()
                bits_H3[i] += tmp_bits_H3.cpu()

            for i in range(L3.size()[0]):
                tmp_rec_L3, tmp_bits_L3, wavelet_scale_L3 = self.space_entropy_L(L3[i].unsqueeze(0).cuda(), part_bitplane,
                                                                                 discard_number, rec_L3[i].cuda())
                rec_L3[i] = tmp_rec_L3.cpu()
                bits_L3[i] += tmp_bits_L3.cpu()

        rec_H1 = torch.cat(rec_H1, 0)
        rec_H2 = torch.cat(rec_H2, 0)
        rec_H3 = torch.cat(rec_H3, 0)
        rec_L3 = torch.cat(rec_L3, 0)

        # temporal wavelet inverse transform
        height = rec_L3.size()[2]
        width = rec_L3.size()[3]
        pad = torch.zeros((1,1,height, width))
        rec_H1 = torch.cat((rec_H1, pad), 0)
        rec_H2 = torch.cat((rec_H2, pad), 0)
        rec_H3 = torch.cat((rec_H3, pad), 0)
        rec_L3 = torch.cat((rec_L3, pad), 0)
        rec_H1 = rec_H1.unsqueeze(1)
        rec_H2 = rec_H2.unsqueeze(1)
        rec_H3 = rec_H3.unsqueeze(1)
        rec_L3 = rec_L3.unsqueeze(1)

        L3_rec_list = []
        L3_rec_list.append(rec_L3[0])

        # third level
        L2_rec_list= []
        i = 0
        L2 = self.update_wavelet_decode_H3(rec_H3[0].cuda(),-up_mv_H3[0].cuda(),up_mask_H3[0].cuda(),
                                              rec_H3[0].cuda(),-up_mv_H3[1].cuda(),up_mask_H3[1].cuda(), grid_up, grid_org, L3_rec_list[0].cuda(), alpha_level, beta_level, rank)
        if torch.sum(up_mask_H3[0]) ==0 and torch.sum(up_mask_H3[1]) == 0:
            L2 = L3_rec_list[0] / beta_level / copycomp_weight_low[2]
        L2_rec_list.append(L2.cpu())

        # second level
        L1_rec_list = []
        i = 0
        L1 = self.update_wavelet_decode_H2(rec_H2[0].cuda(), -up_mv_H2[0].cuda(), up_mask_H2[0].cuda(),
                                              rec_H2[0].cuda(), -up_mv_H2[1].cuda(), up_mask_H2[1].cuda(), grid_up, grid_org, L2_rec_list[0].cuda(), alpha_level, beta_level, rank)
        if torch.sum(up_mask_H2[0]) ==0 and torch.sum(up_mask_H2[1]) == 0:
            L1 = L2_rec_list[0] / beta_level /copycomp_weight_low[1]
        L1_rec_list.append(L1.cpu())

        # first level
        rec_frame_list = []
        i=0
        frame = self.update_wavelet_decode_H1(rec_H1[0].cuda(), -up_mv_H1[0].cuda(),
                                              up_mask_H1[0].cuda(), rec_H1[0].cuda(),
                                              -up_mv_H1[1].cuda(), up_mask_H1[1].cuda(), grid_up, grid_org, L1_rec_list[0].cuda(), alpha_level, beta_level, rank)
        if torch.sum(up_mask_H1[0]) ==0 and torch.sum(up_mask_H1[1]) == 0:
            frame = L1_rec_list[0] / beta_level /copycomp_weight_low[0]
        rec_frame_list.append(frame.cpu())

        for gop_order in range(gop_number):
            # print(gop_order)
            L3_rec_list.append(rec_L3[gop_order+1])

            # third level
            i = 0
            if gop_order != gop_number-1:
                L2 = self.update_wavelet_decode_H3(rec_H3[gop_order].cuda(),
                          -up_mv_H3[gop_order*2+2].cuda(),
                          up_mask_H3[gop_order*2+2].cuda(),
                          rec_H3[gop_order+1].cuda(),
                          -up_mv_H3[gop_order*2+2+1].cuda(),
                          up_mask_H3[gop_order*2+2+1].cuda(), grid_up, grid_org, L3_rec_list[gop_order + 1].cuda(), alpha_level, beta_level, rank)
                if torch.sum(up_mask_H3[gop_order*2+2]) == 0 and torch.sum(up_mask_H3[gop_order*2+2+1]) == 0:
                    L2 = L3_rec_list[gop_order + 1] / beta_level / copycomp_weight_low[2]
            else:
                L2 = L3_rec_list[gop_order + 1] / beta_level / copycomp_weight_low[2]
            L2_rec_list.append(L2.cpu())
            i = 0
            H_tmp = self.preditct_wavelet_decode_H3(L2_rec_list[gop_order*2].cuda(),
                                   mv_H3[gop_order*2].cuda(),
                                   mask_H3[gop_order*2].cuda(),
                                   L2_rec_list[1 + gop_order*2].cuda(),
                                   mv_H3[gop_order*2+1].cuda(),
                                   mask_H3[gop_order*2+1].cuda(), grid_up,grid_org, rec_H3[gop_order].cuda(), alpha_level, rank)
            if torch.sum(mask_H3[gop_order*2])==0 and torch.sum(mask_H3[gop_order*2+1])==0:
                H_tmp = rec_H3[gop_order] / alpha_level / copycomp_weight_high[2]
            L2_rec_list.insert(gop_order*2+1, H_tmp.cpu())

            # second level
            for i in range(2):
                L1 = self.update_wavelet_decode_H2(rec_H2[1 + gop_order * 2 + i - 1].cuda(),
                          -up_mv_H2[gop_order*4+2+2*i].cuda(),
                          up_mask_H2[gop_order*4+2+2*i].cuda(),
                          rec_H2[1 + gop_order * 2 + i].cuda(),
                          -up_mv_H2[gop_order*4+2+2*i+1].cuda(),
                          up_mask_H2[gop_order*4+2+2*i+1].cuda(), grid_up, grid_org, L2_rec_list[gop_order*2 + 1 + i].cuda(), alpha_level, beta_level, rank)
                if torch.sum(up_mask_H2[gop_order*4+2+2*i]) == 0 and torch.sum(up_mask_H2[gop_order*4+2+2*i+1]) == 0:
                    L1 = L2_rec_list[gop_order*2 + 1 + i] / beta_level / copycomp_weight_low[1]
                L1_rec_list.append(L1.cpu())
            for i in range(2):
                H_tmp = self.preditct_wavelet_decode_H2(L1_rec_list[gop_order*4+i+i].cuda(),
                               mv_H2[gop_order*4+i*2].cuda(),
                               mask_H2[gop_order*4+2*i].cuda(),
                               L1_rec_list[gop_order*4+i + 1+i].cuda(),
                               mv_H2[gop_order*4+i*2+1].cuda(),
                               mask_H2[gop_order*4+2*i+1].cuda(), grid_up, grid_org, rec_H2[1 + gop_order * 2 + i - 1].cuda(), alpha_level, rank)
                if torch.sum(mask_H2[gop_order*4+2*i]) == 0 and torch.sum(mask_H2[gop_order*4+2*i+1]) == 0:
                    H_tmp = rec_H2[1 + gop_order * 2 + i - 1] / alpha_level / copycomp_weight_high[1]
                L1_rec_list.insert(gop_order*4+2*i+1, H_tmp.cpu())

            # first level
            for i in range(4):
                frame = self.update_wavelet_decode_H1(rec_H1[1 + gop_order * 4 + i - 1].cuda(),
                              -up_mv_H1[gop_order*8+2*i+2].cuda(),
                              up_mask_H1[gop_order*8+2+2*i].cuda(),
                              rec_H1[1 + gop_order * 4 + i].cuda(),
                              -up_mv_H1[gop_order*8+2*i+2+1].cuda(),
                              up_mask_H1[gop_order*8+2+2*i+1].cuda(), grid_up, grid_org, L1_rec_list[gop_order * 4 + 1 + i].cuda(), alpha_level, beta_level, rank)
                if torch.sum(up_mask_H1[gop_order*8+2+2*i]) == 0 and torch.sum(up_mask_H1[gop_order*8+2+2*i+1]) == 0:
                    frame = L1_rec_list[gop_order * 4 + 1 + i] / beta_level / copycomp_weight_low[0]
                rec_frame_list.append(frame.cpu())
            for i in range(4):
                H_tmp = self.preditct_wavelet_decode_H1(
                                rec_frame_list[8 * gop_order + 2*i].cuda(),
                                mv_H1[gop_order*8+2*i].cuda(),
                                mask_H1[gop_order*8+2*i].cuda(),
                                rec_frame_list[8 * gop_order + 2*i + 1].cuda(),
                                mv_H1[gop_order*8+2*i+1].cuda(),
                                mask_H1[gop_order*8+2*i+1].cuda(), grid_up, grid_org, rec_H1[1 + gop_order * 4 + i - 1].cuda(), alpha_level, rank)
                if torch.sum(mask_H1[gop_order*8+2*i]) == 0 and torch.sum(mask_H1[gop_order*8+2*i+1]) == 0:
                    H_tmp = rec_H1[1 + gop_order * 4 + i - 1] / alpha_level / copycomp_weight_high[0]
                rec_frame_list.insert(gop_order*8+1+i*2, H_tmp.cpu())
        return torch.cat(rec_frame_list, 1), [bits_H1, bits_H2, bits_H3, bits_L3]
